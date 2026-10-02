/* solids/grain_physics_loop.h — Spec definitions for the runner-template port
 * of the grain_physics neighbor loops.
 *
 * Three independent Specs:
 *   GrainBackrxSpec    — grain→gas momentum backreaction (GRAIN_BACKREACTION).
 *                        ONEWAY grain sources, j-side atomics + ghost writeback.
 *   GrainRTGasSpec     — gas→grain RT opacity coupling, direction 1
 *                        (RT_OPACITY_FROM_EXPLICIT_GRAINS). SYMMETRIC, pure i-side.
 *   GrainRTGrainSpec   — grain→gas radiation acceleration, direction 2
 *                        (RT_OPACITY_FROM_EXPLICIT_GRAINS). ONEWAY, pure i-side.
 *
 * SSOT for per-pair physics: grain_physics_functions.h (included below).
 * No physics formulas live in this header; pair_kernel bodies call into it.
 *
 * Caller TUs must `#include "../mesh/kernel.h"` BEFORE this header (kernel.h
 * has no include guards; the KOKKOS_INLINE_FUNCTION pair kernels below call
 * kernel_main / kernel_hinv). Same pattern as difffilter_loop.h et al.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef GRAIN_PHYSICS_LOOP_H
#define GRAIN_PHYSICS_LOOP_H

/* Kokkos_Core before allvars.h (its macros may conflict with stdlib names). */
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"  /* FLUID_DEFAULT for j-side fluid filter in pair_kernel */
#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"   /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

#ifdef DO_FLUID_ALTSPECIES_DRAG_CALCULATION
#include "grain_physics_functions.h"   /* SSOT: pair-kernel bodies + LocalIn types */


/* ============================================================================
 * GrainBackrxSpec — grain→gas momentum backreaction.
 * (GRAIN_BACKREACTION)
 * ============================================================================
 * Active sources: grain particles (GRAIN_PTYPES).
 * Neighbors: gas particles (Type==0).
 * Writes: atomic j-side updates to P[j].Vel, P[j].dp, P[j].Grain_AccelTimeMin,
 *         CellP[j].VelPred (+ GRAIN_EVOLUTION fields). Ghost writeback required.
 * AccumData: empty (no i-side reduction) — validation relies on cross-rank
 * ghost-writeback record count and field comparison.
 * ========================================================================== */
#if defined(GRAIN_BACKREACTION)

/* ActiveData for GrainBackrxSpec — wraps GrainBackrxLocalIn for the runner.
 * pos/h_search carry the walker-facing fields; local carries all pair inputs. */
struct GrainBackrxActiveState {
    GrainBackrxLocalIn local;
    Vec3<double>        pos;
    double              h_search;
};

/* AccumData: empty — GrainBackrxSpec does no per-source reduction. */
struct GrainBackrxAccum {
    /* Placeholder: this loop scatters to neighbours and reduces nothing on the
     * source side. Sized to a whole int so it satisfies the shuffle-width rule
     * the device lane reduction asserts (see NlrAccumReducer). */
    char _unused[4];
};

/* DeviceContext: extends base.
 * applying backreaction twice). Same pattern as ThermalFBDeviceContext. */
struct GrainBackrxDeviceContext : NeighborLoopDeviceContextBase {
};

/* NeighborData: carries the array-level P/CellP pointers + neighbor index j
 * so pair_kernel can forward to grain_backrx_pair_kernel(local, j, P, CellP,
 * dp, r2) without physics duplication. */
struct GrainBackrxNeighborData {
    struct particle_data *P_arr;
    struct gas_cell_data *CellP_arr;
    int  j;
};

struct GrainBackrxSpec {
    static constexpr const char *loop_name = "grain_backrx";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::EpsilonAtomic; /* EpsilonAtomic: drag back-reaction atomic_add(Pj.Vel/dp)+atomic_min(AccelTimeMin), stable Pj snapshot, no read-back -> ulp */

    static constexpr int                    search_mode        = MODE_B_SEARCH_ONEWAY;
    static constexpr unsigned int           neighbor_type_mask = 1u << 0;  /* gas */
    static constexpr mode_b_radius_policy_t radius_policy      = MODE_B_RADIUS_DEFAULT;

    static constexpr WritePattern  write_pattern             = WritePattern::ActiveReduceOnly;
    /* GasOnly: all j-neighbors are Type==0 gas. */
    static constexpr SidxCacheKind sidx_cache_kind           = SidxCacheKind::GasOnly;
    static constexpr bool mode_a_active_sources_in_sidx_pool = false; /* non-pool active sources (sink/star/grain) -> runner stages explicit P[].Pos */
    static constexpr bool          uses_ghost_writeback      = true;
    static constexpr bool writes_neighbour_motion = true;   /* the grain back-reaction kicks the gas */
    static constexpr bool          uses_ghost_write_detector = false;

    using CallScalars    = NlrCommonScalars;
    using ActiveData     = GrainBackrxActiveState;
    using AccumData      = GrainBackrxAccum;
    using NeighborData   = GrainBackrxNeighborData;
    using DeviceContext  = GrainBackrxDeviceContext;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

    struct Aux {};

    static_assert(uses_ghost_writeback == true,
        "GrainBackrxSpec: j-writes to gas require ghost writeback for cross-rank deltas.");

    /* Host hooks — bodies in grain_physics_loop.cc. */
    static bool        is_active(int i);
    static double      search_radius(const neighbor_loop_args& args, int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);
    static void        apply_active_writeback(const neighbor_loop_args& args,
                                              int active_slot, int i,
                                              const AccumData& accum);
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& /*local_accum*/,
                            const AccumData& /*peer_accum*/)
    {}
    static void        ghost_writeback_begin(const neighbor_loop_args& args,
                                             const NeighborLoopPlan& plan);
    static void        ghost_writeback_end  (const neighbor_loop_args& args,
                                             const NeighborLoopPlan& plan);
    /* populate_device_context: required because DeviceContext extends the base.
     * No UVM allocation needed. */
    static void populate_device_context(const neighbor_loop_args& /*args*/,
                                        DeviceContext& /*ctx*/) {}

    /* Device hooks — header-inline. */

    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& /*accum*/) {}  /* nothing to zero */

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx,
                                  int /*active_slot*/, int i,
                                  double h_search, const CallScalars& /*scalars*/) {
        ActiveData a{};
        struct particle_data& Pi = ctx.P[i];
        a.local.Pos                 = Pi.Pos;
        a.local.KernelRadius        = Pi.KernelRadius;
        a.local.Grain_DeltaMomentum = Pi.Grain_DeltaMomentum;
        a.local.Gas_Density         = Pi.Gas_Density;
        a.local.Grain_AccelTimeMin  = Pi.Grain_AccelTimeMin;
#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
        for(int kv = 0; kv < GRAIN_NUM_VOLATILE_SPECIES; kv++) {
            a.local.Grain_DeltaVolatileMass[kv] = Pi.Grain_DeltaVolatileMass[kv];
        }
        a.local.Grain_DeltaInternalEnergyHeating = Pi.Grain_DeltaInternalEnergyHeating;
#endif
        a.pos[0]   = (double)Pi.Pos[0];
        a.pos[1]   = (double)Pi.Pos[1];
        a.pos[2]   = (double)Pi.Pos[2];
        a.h_search = h_search;
        return a;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx, int j,
                                      const IdentitySidecar& /*id*/,
                                      const ActiveData& /*active*/) {
        NeighborData n{};
        n.P_arr         = ctx.P;
        n.CellP_arr     = ctx.CellP;
        n.j             = j;
        return n;
    }

    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                            const NeighborData& nb,
                            AccumData& /*accum*/,
                            NoScatter& /*scatter*/,
                            const CallScalars& /*cs*/) {
        if(active.h_search <= 0) return;
        if(nb.P_arr == nullptr || nb.CellP_arr == nullptr) return;
#ifdef HYDRO_MULTIFLUID
        /* Backreaction targets only default-fluid gas; specialty-fluid neighbors
         * (e.g. FLUID_DUST_GRAIN, FLUID_ION) are sources, not receivers. Mirrors
         * the j-side filter that lived in the legacy grain_backrx_evaluate_gpu. */
        if(nb.P_arr[nb.j].FluidType != FLUID_DEFAULT) return;
#endif
        Vec3<double> dp;
        dp[0] = active.pos[0] - (double)nb.P_arr[nb.j].Pos[0];
        dp[1] = active.pos[1] - (double)nb.P_arr[nb.j].Pos[1];
        dp[2] = active.pos[2] - (double)nb.P_arr[nb.j].Pos[2];
        nearest_xyz(dp);
        double r2 = dp.norm_sq();
        grain_backrx_pair_kernel(active.local, nb.j, nb.P_arr, nb.CellP_arr, dp, r2);
    }
};

#endif /* GRAIN_BACKREACTION */


/* ============================================================================
 * GrainRTGasSpec — gas→grain RT opacity coupling, direction 1.
 * (RT_OPACITY_FROM_EXPLICIT_GRAINS)
 * ============================================================================
 * Active sources: gas particles (Type==0).
 * Neighbors: grain particles (GRAIN_PTYPES). SYMMETRIC at P[j].KernelRadius.
 * Writes: i-side only. Accumulates InterpolatedGeometricDustCrossSection,
 *         Interpolated_Opacity[], and (MHD_BATTERY_MECHANISMS & 8) J_dust_contribution.
 * ========================================================================== */
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)

/* Shared ActiveData for both RT directions — GasGrainRTLocalIn carries all
 * pair-kernel inputs; pos/h_search are the walker-facing fields. */
struct GasGrainRTActiveState {
    GasGrainRTLocalIn local;
    Vec3<double>      pos;
    double            h_search;
};

/* AccumData for GrainRTGasSpec — gas-direction fields only. */
struct GrainRTGasAccum {
    MyDouble InterpolatedGeometricDustCrossSection;
    MyDouble Interpolated_Opacity[N_RT_FREQ_BINS];
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
    Vec3<MyDouble> J_dust_contribution;
#endif
};

/* AccumData for GrainRTGrainSpec — grain-direction field only. */
struct GrainRTGrainAccum {
    Vec3<MyDouble> Interpolated_Radiation_Acceleration;
};

/* Shared NeighborData for both RT directions. */
struct GasGrainRTNeighborData {
    struct particle_data *P_arr;
    struct gas_cell_data *CellP_arr;
    int j;
};

struct GrainRTGasSpec {
    static constexpr const char *loop_name = "grain_rt_gas";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only (out merged to accum post-pair); no j-write, no atomics -> threaded eval bit-identical */

    /* SYMMETRIC: grain kernel radius at j filters pairs (h_to_use = P[j].KernelRadius
     * inside gasgrain_rt_gas_search_pair_kernel). SYMMETRIC mode ensures pairs
     * with h_j > h_i are not missed. */
    static constexpr int                    search_mode        = MODE_B_SEARCH_SYMMETRIC;
    /* GRAIN_PTYPES: neither GasOnly (tbm=1) nor AllTypes (tbm=0x3f);
     * do not use a typed cache. */
    static constexpr unsigned int           neighbor_type_mask = GRAIN_PTYPES;
    /* grain_rt_gas pair physics: h_to_use = P[j].KernelRadius for non-gas grain
     * particles (see grain_physics_functions.h:214).  Symmetric search reach
     * needs NONGAS_KERNEL; gas excluded by mask. */
    static constexpr mode_b_radius_policy_t radius_policy      = MODE_B_RADIUS_NONGAS_KERNEL;

    static constexpr WritePattern  write_pattern             = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind sidx_cache_kind           = SidxCacheKind::None;
    static constexpr bool          uses_ghost_writeback      = false;
    static constexpr bool          uses_ghost_write_detector = false;

    using CallScalars    = NlrCommonScalars;
    using ActiveData     = GasGrainRTActiveState;
    using AccumData      = GrainRTGasAccum;
    using NeighborData   = GasGrainRTNeighborData;
    using DeviceContext  = NeighborLoopDeviceContextBase;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

    struct Aux {};

    /* Host hooks — bodies in grain_physics_loop.cc. */
    static bool        is_active(int i);
    static double      search_radius(const neighbor_loop_args& args, int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);
    static void        apply_active_writeback(const neighbor_loop_args& args,
                                              int active_slot, int i,
                                              const AccumData& accum);
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local_accum, const AccumData& peer_accum)
    {
        local_accum.InterpolatedGeometricDustCrossSection +=
            peer_accum.InterpolatedGeometricDustCrossSection;
        for(int k = 0; k < N_RT_FREQ_BINS; k++) {
            local_accum.Interpolated_Opacity[k] += peer_accum.Interpolated_Opacity[k];
        }
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
        for(int k = 0; k < 3; k++) {
            local_accum.J_dust_contribution[k] += peer_accum.J_dust_contribution[k];
        }
#endif
    }

    /* Device hooks — header-inline. */

    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        accum.InterpolatedGeometricDustCrossSection = 0;
        for(int k = 0; k < N_RT_FREQ_BINS; k++) accum.Interpolated_Opacity[k] = 0;
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
        accum.J_dust_contribution = {};
#endif
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx,
                                  int /*active_slot*/, int i,
                                  double h_search, const CallScalars& /*scalars*/) {
        ActiveData a{};
        struct particle_data& Pi = ctx.P[i];
        a.local.Type        = Pi.Type;
        a.local.Pos         = Pi.Pos;
        a.local.Vel         = Pi.Vel;
        a.local.KernelRadius = Pi.KernelRadius;
        a.local.Mass        = Pi.Mass;
        a.local.Grain_Size  = 0;
        for(int k = 0; k < N_RT_FREQ_BINS; k++) a.local.Grain_Abs_Coeff[k] = 0;
        /* Gas source: no grain-direction fields needed. */
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
        /* Gas source: pack electron fraction + density for J_dust Z_grain estimate. */
        a.local.Ne_gas      = (Pi.Type == 0 && ctx.CellP != nullptr)
                              ? ctx.CellP[i].Ne : 0.0f;
        a.local.Density_gas = (Pi.Type == 0 && ctx.CellP != nullptr)
                              ? ctx.CellP[i].Density : 0.0f;
#endif
        a.pos[0]   = (double)Pi.Pos[0];
        a.pos[1]   = (double)Pi.Pos[1];
        a.pos[2]   = (double)Pi.Pos[2];
        a.h_search = h_search;
        return a;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx, int j,
                                      const IdentitySidecar& /*id*/,
                                      const ActiveData& /*active*/) {
        GasGrainRTNeighborData n{};
        n.P_arr    = ctx.P;
        n.CellP_arr = ctx.CellP;
        n.j        = j;
        return n;
    }

    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                            const NeighborData& nb,
                            AccumData& accum,
                            NoScatter& /*scatter*/,
                            const CallScalars& /*cs*/) {
        if(active.h_search <= 0) return;
        if(nb.P_arr == nullptr) return;
        Vec3<double> dp;
        dp[0] = active.pos[0] - (double)nb.P_arr[nb.j].Pos[0];
        dp[1] = active.pos[1] - (double)nb.P_arr[nb.j].Pos[1];
        dp[2] = active.pos[2] - (double)nb.P_arr[nb.j].Pos[2];
        nearest_xyz(dp);
        double r2 = dp.norm_sq();

        /* Construct a local GasGrainRTOut, call SSOT body, then merge into accum. */
        GasGrainRTOut pair_out{};
        gasgrain_rt_gas_search_pair_kernel(active.local, nb.j, nb.P_arr, nb.CellP_arr,
                                           dp, r2, pair_out);
        accum.InterpolatedGeometricDustCrossSection += pair_out.InterpolatedGeometricDustCrossSection;
        for(int k = 0; k < N_RT_FREQ_BINS; k++) {
            accum.Interpolated_Opacity[k] += pair_out.Interpolated_Opacity[k];
        }
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
        for(int k = 0; k < 3; k++) {
            accum.J_dust_contribution[k] += pair_out.J_dust_contribution[k];
        }
#endif
    }
};


/* ============================================================================
 * GrainRTGrainSpec — grain→gas radiation acceleration, direction 2.
 * (RT_OPACITY_FROM_EXPLICIT_GRAINS)
 * ============================================================================
 * Active sources: grain particles (GRAIN_PTYPES).
 * Neighbors: gas particles (Type==0). ONEWAY at local.KernelRadius.
 * Writes: i-side only. Accumulates Interpolated_Radiation_Acceleration.
 * ========================================================================== */

struct GrainRTGrainSpec {
    static constexpr const char *loop_name = "grain_rt_grain";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only (out.Interpolated_Radiation_Acceleration merged post-pair); no j-write, no atomics -> threaded eval bit-identical */

    static constexpr int                    search_mode        = MODE_B_SEARCH_ONEWAY;
    static constexpr unsigned int           neighbor_type_mask = 1u << 0;  /* gas */
    static constexpr mode_b_radius_policy_t radius_policy      = MODE_B_RADIUS_DEFAULT;

    static constexpr WritePattern  write_pattern             = WritePattern::ActiveReduceOnly;
    /* GasOnly: all j-neighbors are Type==0 gas. */
    static constexpr SidxCacheKind sidx_cache_kind           = SidxCacheKind::GasOnly;
    static constexpr bool mode_a_active_sources_in_sidx_pool = false; /* non-pool active sources (sink/star/grain) -> runner stages explicit P[].Pos */
    static constexpr bool          uses_ghost_writeback      = false;
    static constexpr bool          uses_ghost_write_detector = false;

    using CallScalars    = NlrCommonScalars;
    using ActiveData     = GasGrainRTActiveState;
    using AccumData      = GrainRTGrainAccum;
    using NeighborData   = GasGrainRTNeighborData;
    using DeviceContext  = NeighborLoopDeviceContextBase;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

    struct Aux {};

    /* Host hooks — bodies in grain_physics_loop.cc. */
    static bool        is_active(int i);
    static double      search_radius(const neighbor_loop_args& args, int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);
    static void        apply_active_writeback(const neighbor_loop_args& args,
                                              int active_slot, int i,
                                              const AccumData& accum);
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local_accum, const AccumData& peer_accum)
    {
        for(int k = 0; k < 3; k++) {
            local_accum.Interpolated_Radiation_Acceleration[k] +=
                peer_accum.Interpolated_Radiation_Acceleration[k];
        }
    }

    /* Device hooks — header-inline. */

    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        accum.Interpolated_Radiation_Acceleration = {};
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx,
                                  int /*active_slot*/, int i,
                                  double h_search, const CallScalars& /*scalars*/) {
        ActiveData a{};
        struct particle_data& Pi = ctx.P[i];
        a.local.Type         = Pi.Type;
        a.local.Pos          = Pi.Pos;
        a.local.Vel          = Pi.Vel;
        a.local.KernelRadius = Pi.KernelRadius;
        a.local.Mass         = Pi.Mass;
        a.local.Grain_Size   = 0;
        for(int k = 0; k < N_RT_FREQ_BINS; k++) a.local.Grain_Abs_Coeff[k] = 0;
        /* Grain source: fill Grain_Size and Grain_Abs_Coeff for grain→gas direction. */
        if(((1 << Pi.Type) & GRAIN_PTYPES) && Pi.Gas_Density > 0) {
            a.local.Grain_Size = Pi.Grain_Size;
            double R_grain_code   = (double)Pi.Grain_Size / UNIT_LENGTH_IN_CGS;
            double rho_grain_code = All.Grain_Internal_Density / UNIT_DENSITY_IN_CGS;
            double rho_gas_code   = (double)Pi.Gas_Density * All.cf_a3inv;
            for(int k_freq = 0; k_freq < N_RT_FREQ_BINS; k_freq++) {
                double Q_abs = grain_extinction_Q_inline(i, k_freq, ctx.P);
                a.local.Grain_Abs_Coeff[k_freq] = (MyFloat)(Q_abs * 3.0 /
                    (4.0 * C_LIGHT_CODE_REDUCED * rho_grain_code * R_grain_code * rho_gas_code));
            }
        }
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
        a.local.Ne_gas      = 0.0f;
        a.local.Density_gas = 0.0f;
#endif
        a.pos[0]   = (double)Pi.Pos[0];
        a.pos[1]   = (double)Pi.Pos[1];
        a.pos[2]   = (double)Pi.Pos[2];
        a.h_search = h_search;
        return a;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx, int j,
                                      const IdentitySidecar& /*id*/,
                                      const ActiveData& /*active*/) {
        GasGrainRTNeighborData n{};
        n.P_arr    = ctx.P;
        n.CellP_arr = ctx.CellP;
        n.j        = j;
        return n;
    }

    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                            const NeighborData& nb,
                            AccumData& accum,
                            NoScatter& /*scatter*/,
                            const CallScalars& /*cs*/) {
        if(active.h_search <= 0) return;
        if(nb.P_arr == nullptr) return;
        Vec3<double> dp;
        dp[0] = active.pos[0] - (double)nb.P_arr[nb.j].Pos[0];
        dp[1] = active.pos[1] - (double)nb.P_arr[nb.j].Pos[1];
        dp[2] = active.pos[2] - (double)nb.P_arr[nb.j].Pos[2];
        nearest_xyz(dp);
        double r2 = dp.norm_sq();

        GasGrainRTOut pair_out{};
        gasgrain_rt_grain_search_pair_kernel(active.local, nb.j, nb.P_arr, nb.CellP_arr,
                                             dp, r2, pair_out);
        for(int k = 0; k < 3; k++) {
            accum.Interpolated_Radiation_Acceleration[k] +=
                pair_out.Interpolated_Radiation_Acceleration[k];
        }
    }
};

#endif /* RT_OPACITY_FROM_EXPLICIT_GRAINS */

#endif /* DO_FLUID_ALTSPECIES_DRAG_CALCULATION */

#endif /* GRAIN_PHYSICS_LOOP_H */
