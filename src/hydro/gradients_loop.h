/* hydro/gradients_loop.h — GradientsSpec for the neighbor-loop runner.
 *
 * Ports the legacy GPU walker `gradient_evaluate_gpu` (hydro/density_gpu.cc:
 * 93-279) to the runner Spec contract. Computes per-active gradients
 * (Density / Pressure / Velocity + many conditional fields) via symmetric
 * gas-gas neighbor topology, accumulating into the standard
 * `GasGraddata_out_` struct. Wraps the unchanged inline pair body
 * `gradient_accumulate_neighbor` (hydro/gradient_functions.h:189) — pair
 * physics is byte-for-byte preserved.
 *
 * Broad active row policy matching the legacy GPU walker (Type==0 && Mass>0);
 * narrow GasGrad_isactive gate stays at the neighbor side inside the pair
 * body.
 *
 * Written by Philip F. Hopkins (phopkins@caltech.edu) for GIZMO. */

#ifndef GRADIENTS_LOOP_H
#define GRADIENTS_LOOP_H

/* Kokkos_Core.hpp must precede allvars.h (its macros may conflict with stdlib
 * names). */
#include <Kokkos_Core.hpp>
#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"
#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"      /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
#include "gradient_functions.h"               /* Quantities_for_Gradients,
                                                * GasGraddata_in_/out_,
                                                * kernel_GasGrad,
                                                * gradient_accumulate_neighbor,
                                                * GasGrad_isactive_gpu,
                                                * SHOULD_I_USE_SPH_GRADIENTS */
/* NOTE: caller translation units must include "../mesh/kernel.h" BEFORE
 * this header. kernel.h has no include guards (defines static inline
 * kernel_main / kernel_hinv used by the pair body); double-include triggers
 * redefinition errors. Matches the cellcorrections_loop.h / sink_env1_loop.h
 * convention. gradient_functions.h is similarly include-guard-free for the
 * same reason. */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* Number of MHD-CG outer iterations (legacy: see hydro/gradients.cc top).
 * Used by hydro_gradient_calc()'s outer loop. */
#if defined(MHD_CONSTRAINED_GRADIENT)
#if (MHD_CONSTRAINED_GRADIENT > 1)
#define NUMBER_OF_GRADIENT_ITERATIONS 3
#else
#define NUMBER_OF_GRADIENT_ITERATIONS 2
#endif
#else
#define NUMBER_OF_GRADIENT_ITERATIONS 1
#endif

/* Host-side narrow predicate (defined in gradients_loop.cc; referenced by
 * cellcorrections_loop.h via extern forward-decl). */
int GasGrad_isactive(int i, struct particle_data *pp, struct gas_cell_data *cell);

/* ============================================================================
 * Per-active scratch carried across MHD-CG iterations + into finalization.
 * Migrated from the file-static `temporary_data_topass` in the legacy
 * hydro/gradients.cc (lines 209-240). Named type so GradientsAux can hold a
 * base pointer and apply_active_writeback can scatter into it.
 * ========================================================================== */
struct temporary_data_topass
{
    struct Quantities_for_Gradients Maxima;
    struct Quantities_for_Gradients Minima;
    MyFloat MaxDistance;
#if defined(KERNEL_CRK_FACES)
    MyDouble m0;
    MyDouble m1[3];
    MyDouble m2[6];
    MyDouble dm0[3];
    MyDouble dm1[3][3];
    MyDouble dm2[6][3];
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
    Vec3<MyFloat> GlassAcc;
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
    MyDouble FaceDotB;
    MyDouble FaceCrossX[3][3];
    MyDouble BGrad[3][3];
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
    Vec3<MyDouble> PhiGrad;
#endif
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
    Vec3<MyFloat> Gradients_Rad_E_gamma[N_RT_FREQ_BINS];
#endif
#ifdef TURB_DIFF_DYNAMIC
    Vec3<MyDouble> GradVelocity_bar[3];
#endif
};

/* ============================================================================
 * Per-pair physics types.
 * ========================================================================== */

/* ActiveData = GasGraddata_in_ `local` plus per-active scratch the pair body
 * needs (computed once in load_active so it isn't recomputed per neighbor:
 * Mass sign trick result + kernel_mode + V_i + hinv triplet). Legacy Mass
 * sign convention from density_gpu.cc:144-208 preserved byte-for-byte. */
struct GradientsActiveData
{
    /* Runner Mode B remote walker requires `pos` (Vec3<double>) + `h_search`
     * (double) on ActiveData — peer ranks need geometry without unpacking
     * the rest of the struct. Same shape as cellcorrections / sink_env1. */
    Vec3<double> pos;
    double       h_search;

    struct GasGraddata_in_ local;
    double hinv;
    double hinv3;
    double hinv4;
    double V_i;
    int    sph_gradients_flag_i;
    int    kernel_mode_i;
    bool   enabled;   /* row gate: false for a gas cell driven to Mass<=0 after
                       * the shared row list was frozen (feedback/swallow can
                       * kill a cell mid-corridor). Legacy rebuilt row lists
                       * post-feedback with a Mass>0 filter, so dead rows were
                       * simply ABSENT; this gate reproduces that row-absence
                       * (pair_kernel early-returns, accum stays zero). Struct
                       * member, not a computed local int, so it is safe from
                       * device constant-propagation of gate variables. */
#ifdef HYDRO_MULTIFLUID
    unsigned char FluidType;   /* packed P[i].FluidType — for same_lagrangian_fluid_id() */
#endif
};

/* NeighborData carries the integer index j plus the base P/CellP pointers.
 * gradient_accumulate_neighbor indexes ~30 P[j]/CellP[j] fields by j, so the
 * pair_kernel needs both. Different shape from sink_env1's `*particle_data`
 * because the gradient pair body accesses j by integer throughout. */
struct GradientsNeighborData
{
    int                   j;
    struct particle_data *P;
    struct gas_cell_data *CellP;
};

/* Aux carries the legacy `GasGradDataPasser` base pointer + the current
 * MHD-CG iteration index. apply_active_writeback dispatches on grad_iter to
 * replay either out2particle_GasGrad (iter==0) or out2particle_GasGrad_iter
 * (iter>0). load_active reads grad_iter through DeviceContext below (no
 * access to args from device-side load_active). */
struct GradientsAux
{
    struct temporary_data_topass *passer;    /* sized N_gas (host) */
    int                            grad_iter;
};

/* DeviceContext extension: ferries grad_iter from host (args.aux->grad_iter)
 * to device-side load_active. Trivially copyable — runner captures by value
 * into the Kokkos device lambda. populate_device_context body in
 * gradients_loop.cc. No UVM allocation -> no cleanup_device_context. */
struct GradientsDeviceContext : NeighborLoopDeviceContextBase
{
    int grad_iter;
};

/* ============================================================================
 * GradientsSpec — NeighborLoopSpec contract.
 * ========================================================================== */
struct GradientsSpec
{
    static constexpr const char *loop_name = "gradients";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only (delegates to gradient_accumulate_neighbor: out-> + local kernel scratch); no j-write, no atomics -> threaded eval bit-identical */

    /* Symmetric gas-gas topology (matches gizmo_sym_neighbor_list + the
     * legacy gradient_evaluate_gpu CSR consumer). */
    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int            neighbor_type_mask = (1u << 0);
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;

    /* Pure i-side accumulate (writes via apply_active_writeback into host
     * CellP[i].Gradients / GasGradDataPasser[i]). No ghost-writeback. */
    static constexpr WritePattern   write_pattern              = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind            = SidxCacheKind::GasOnly;
    static constexpr bool           uses_ghost_writeback       = false;
    static constexpr bool           uses_ghost_write_detector  = false;

#ifndef MHD_CONSTRAINED_GRADIENT
    /* run_mode_a may chunk this Spec's per-active staging. Audited safe: the
     * pair kernel reads only primitive GQuant fields; apply_active_writeback
     * writes only OUTPUT gradients (CellP[i].Gradients.*) + the separate
     * GasGradDataPasser -- never a pair-read primitive -- so a per-chunk
     * writeback is disjoint from later chunks' reads. Excluded under
     * MHD_CONSTRAINED_GRADIENT, where the pair reads CellP[j].Gradients.B
     * (constrained_facedotb_delta) that the writeback writes: a read-after-write
     * the chunk interleave would expose. */
    static constexpr bool           mode_a_chunked_active_staging = true;
#endif

    /* Broad active predicate matching legacy GPU walker's row source
     * `gizmo_sym_active_indices` (Type==0 && Mass>0). The narrow filter
     * (GasGrad_isactive: KernelRadius>0, Density>0, DelayTime>0 gates) is
     * applied at the neighbor side inside gradient_accumulate_neighbor via
     * GasGrad_isactive_gpu(j) — same as today. Adding the narrow gate here
     * would silently filter active-i rows the legacy GPU walker processed
     * (directive #5 violation). */
    static bool is_active(int i) {
        if(P[i].Type != 0) return false;
        if(P[i].Mass <= 0) return false;
        return true;
    }

    /* Per-pair physics types */
    using CallScalars   = NlrCommonScalars;
    using ActiveData    = GradientsActiveData;
    using AccumData     = struct GasGraddata_out_;
    using NeighborData  = GradientsNeighborData;
    using Aux           = GradientsAux;

    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;
    using DeviceContext  = GradientsDeviceContext;

    /* Per-active host hook: pre-arena search radius from P[i].KernelRadius. */
    static double search_radius(const neighbor_loop_args& /*args*/,
                                 int /*active_slot*/, int i)
    {
        return (double)P[i].KernelRadius;
    }

    /* Per-call scalars: only the common cosmology block. The pair body reads
     * All.* directly (bare `All.*` is fine under the All-mirror refactor —
     * see feedback_bare_all_in_pair_body_ok.md), so scalars is unused. */
    static CallScalars populate_call_scalars(const neighbor_loop_args& /*args*/)
    {
        return nlr_common_scalars_from_all();
    }

    /* DeviceContext extension hook: copy grad_iter from Aux into ctx.
     * Body in gradients_loop.cc. */
    static void populate_device_context(const neighbor_loop_args& args,
                                         DeviceContext& ctx);

    /* ----- Device hooks ----- */

    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum)
    {
        /* GasGraddata_out_ is POD — byte-zero matches legacy
         * `memset(&out, 0, sizeof(out))` from density_gpu.cc:212. */
        for(size_t b = 0; b < sizeof(accum); b++) ((char*)&accum)[b] = 0;
    }

    /* Build ActiveData per active i. Byte-exact reconstruction of the
     * particle2in-equivalent block at density_gpu.cc:138-227, including the
     * Mass sign trick (negate as flag for SPH gradients), MHD-CG iter>0
     * Mass=0 sentinel (read from ctx.grad_iter), and per-active
     * kernel/mode/V_i precompute. No defensive safety early-outs added —
     * preserves current GPU walker behavior exactly. */
    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx,
                                   int /*active_slot*/, int i,
                                   double h_search,
                                   const CallScalars& /*scalars*/)
    {
        ActiveData active;
        for(size_t b = 0; b < sizeof(active); b++) ((char*)&active)[b] = 0;

        struct GasGraddata_in_ &local = active.local;
        struct particle_data   *kp    = ctx.P;
        struct gas_cell_data   *kc    = ctx.CellP;

        active.pos      = kp[i].Pos;
        active.h_search = h_search;

        local.Pos          = kp[i].Pos;
        local.KernelRadius = kp[i].KernelRadius;
        local.Mass         = kp[i].Mass;
        if(local.Mass < 0) { local.Mass = 0; }

        /* Dead-row gate (see GradientsActiveData::enabled): also guards the
         * Mass/Density divisions below against a mid-corridor-killed cell. */
        active.enabled = (kp[i].Mass > 0 && kc[i].Density > 0);

        active.sph_gradients_flag_i = SHOULD_I_USE_SPH_GRADIENTS(kc[i].ConditionNumber);
        if(active.sph_gradients_flag_i) { local.Mass *= -1; }

#ifdef MHD_CONSTRAINED_GRADIENT
        if(ctx.grad_iter > 0) {
            if(kc[i].FlagForConstrainedGradients <= 0) { local.Mass = 0; }
        }
#endif

        local.GQuant.Density  = kc[i].Density;
        local.GQuant.Pressure = kc[i].Pressure;
        local.GQuant.Velocity = kc[i].VelPred;
#ifdef MAGNETIC
        local.GQuant.B = active.enabled ? (kc[i].BPred * (kc[i].Density / kp[i].Mass)) : Vec3<double>{};
#ifdef DIVBCLEANING_DEDNER
        local.GQuant.Phi = active.enabled ? (kc[i].PhiPred / kp[i].Mass) : 0;
#endif
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
        local.GQuant.InternalEnergy = kc[i].InternalEnergyPred;
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
        local.GQuant.ElectronNumberDensity = kc[i].n_e();
        local.GQuant.ElectronTemperature   = kc[i].T_e();
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
        local.GQuant.E_battery_T2 = kc[i].E_battery_T2_cell;
#endif
#ifdef DOGRAD_SOUNDSPEED
        local.GQuant.SoundSpeed = kc[i].effective_soundspeed();
#endif
#ifdef COSMIC_RAY_FLUID
        for(int k = 0; k < N_CR_PARTICLE_BINS; k++) {
            local.GQuant.CosmicRayPressure[k] = Get_Gas_CosmicRayPressure(i, k, kc);
        }
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        for(int k = 0; k < NUM_METAL_SPECIES; k++) {
            local.GQuant.Metallicity[k] = kp[i].Metallicity[k];
        }
#endif
#if defined(RT_COMPGRAD_EDDINGTON_TENSOR) && (N_RT_FREQ_BINS > 0)
        for(int k = 0; k < N_RT_FREQ_BINS; k++) {
            local.GQuant.Rad_E_gamma[k]    = kc[i].Rad_E_gamma_Pred[k];
            local.GQuant.Rad_E_gamma_ET[k] = kc[i].ET[k];
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
            for(int k2 = 0; k2 < 3; k2++) {
                local.GQuant.Rad_Flux[k][k2] = kc[i].Rad_Flux_Pred[k][k2];
            }
#endif
        }
#endif
#ifdef TURB_DIFF_DYNAMIC
        local.GQuant.Velocity_bar = kc[i].Velocity_bar;
        local.Norm_hat            = kc[i].Norm_hat;
#ifdef GALSF_SUBGRID_WINDS
        local.DelayTime = kc[i].DelayTime;
#endif
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
        local.ConditionNumber = kc[i].ConditionNumber;
        local.NV_T            = kc[i].NV_T;
        for(int k = 0; k < 3; k++) {
            for(int k2 = 0; k2 < 3; k2++) {
                local.BGrad[k][k2] = kc[i].Gradients.B[k][k2];
            }
        }
#ifdef MHD_MODIFIED_GRADIENT
        local.MG_cgcoeff = kc[i].MG_cgcoeff;
#endif
#ifdef MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV
        local.PhiGrad = kc[i].Gradients.Phi;
#endif
#endif

        if(active.sph_gradients_flag_i) { local.Mass *= -1; }    /* negate as flag */

        /* Per-active kernel triple + V_i + kernel_mode (legacy lines 215-226). */
        double h_i = local.KernelRadius;
        kernel_hinv(h_i, &active.hinv, &active.hinv3, &active.hinv4);
        if(local.Mass < 0) { local.Mass *= -1; }                 /* restore for V_i */
        active.V_i = active.enabled ? (local.Mass / local.GQuant.Density) : 0;
        if(active.sph_gradients_flag_i) { local.Mass *= -1; }    /* re-negate for kernel */

        active.kernel_mode_i = -1;
        if(active.sph_gradients_flag_i) active.kernel_mode_i = 0;
#if defined(HYDRO_SPH) || defined(KERNEL_CRK_FACES)
        active.kernel_mode_i = 0;
#endif
#ifdef HYDRO_MULTIFLUID
        active.FluidType = ctx.P[i].FluidType;
#endif
        return active;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx,
                                       int j,
                                       const IdentitySidecar& /*id*/,
                                       const ActiveData& /*active*/)
    {
        NeighborData nb{j, ctx.P, ctx.CellP};
        return nb;
    }

    /* The physics — forwards to the unchanged inline pair body. const_cast
     * on &active.local because gradient_accumulate_neighbor takes
     * GasGraddata_in_* non-const (legacy signature) but only reads. */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                             const NeighborData& neighbor,
                             AccumData& accum,
                             NoScatter& /*scatter*/,
                            const CallScalars& /*cs*/)
    {
        if(!active.enabled) return;   /* dead row (Mass<=0 mid-corridor): zero contribution */
#if defined(HYDRO_MULTIFLUID)
        if (!same_lagrangian_fluid_id(active.FluidType, neighbor.P[neighbor.j].FluidType)) return;
#endif
        struct kernel_GasGrad kernel;
        kernel.h_i = active.local.KernelRadius;
        gradient_accumulate_neighbor(
            const_cast<struct GasGraddata_in_*>(&active.local),
            &accum, &kernel, neighbor.j,
            active.sph_gradients_flag_i, active.V_i,
            active.hinv, active.hinv3, active.hinv4, active.kernel_mode_i,
            neighbor.P, neighbor.CellP);
    }

    /* Host writebacks (post-dispatch). Body in gradients_loop.cc: replays
     * out2particle_GasGrad (grad_iter==0) or out2particle_GasGrad_iter
     * (grad_iter>0), keyed on nlr_aux<GradientsSpec>(args)->grad_iter. */
    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int active_slot, int i,
                                        const AccumData& accum);

    /* Peer-rank accum merge (Mode B remote). Body in gradients_loop.cc:
     * additive for Gradients[k].*, MAX/MIN for Maxima/Minima, MAX for
     * MaxDistance, additive for everything else. */
    /* ============================================================================
     * merge_accum — peer-rank accum reduction (Mode B remote). Field-by-field
     * additive for Gradients[k].*, MAX/MIN for Maxima/Minima, MAX for
     * MaxDistance, additive for all CRK/face/SPH/AGS accumulators. Must match
     * the merge semantics that the pair body implicitly performs across passes.
     * ========================================================================== */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& dst, const AccumData& src)
    {
#define MERGE_ADD(field)  do { dst.field += src.field; } while(0)
#define MERGE_MAX(field)  do { if(src.field > dst.field) dst.field = src.field; } while(0)
#define MERGE_MIN(field)  do { if(src.field < dst.field) dst.field = src.field; } while(0)

        /* Gradients[k].* — additive for every quantity in Quantities_for_Gradients. */
        for(int k = 0; k < 3; k++) {
            MERGE_ADD(Gradients[k].Density);
            MERGE_ADD(Gradients[k].Pressure);
            for(int j = 0; j < 3; j++) { MERGE_ADD(Gradients[k].Velocity[j]); }
#ifdef MAGNETIC
            for(int j = 0; j < 3; j++) { MERGE_ADD(Gradients[k].B[j]); }
#ifdef DIVBCLEANING_DEDNER
            MERGE_ADD(Gradients[k].Phi);
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
            for(int j = 0; j < NUM_METAL_SPECIES; j++) { MERGE_ADD(Gradients[k].Metallicity[j]); }
#endif
#if defined(RT_COMPGRAD_EDDINGTON_TENSOR) && (N_RT_FREQ_BINS > 0)
            for(int j = 0; j < N_RT_FREQ_BINS; j++) {
                MERGE_ADD(Gradients[k].Rad_E_gamma[j]);
                /* SymmetricTensor2 stores 6 unique elements behind [i][j]==[j][i]
                 * aliasing — iterate raw storage to avoid double-adding the
                 * off-diagonals (xy/yz/xz appear via both index orderings). */
                for(int kd = 0; kd < 6; kd++) { MERGE_ADD(Gradients[k].Rad_E_gamma_ET[j].data[kd]); }
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
                for(int kd = 0; kd < 3; kd++) { MERGE_ADD(Gradients[k].Rad_Flux[j][kd]); }
#endif
            }
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
            MERGE_ADD(Gradients[k].InternalEnergy);
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
            MERGE_ADD(Gradients[k].ElectronNumberDensity);
            MERGE_ADD(Gradients[k].ElectronTemperature);
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
            for(int j = 0; j < 3; j++) { MERGE_ADD(Gradients[k].E_battery_T2[j]); }
#endif
#ifdef COSMIC_RAY_FLUID
            for(int j = 0; j < N_CR_PARTICLE_BINS; j++) { MERGE_ADD(Gradients[k].CosmicRayPressure[j]); }
#endif
#ifdef DOGRAD_SOUNDSPEED
            MERGE_ADD(Gradients[k].SoundSpeed);
#endif
#ifdef TURB_DIFF_DYNAMIC
            for(int j = 0; j < 3; j++) { MERGE_ADD(Gradients[k].Velocity_bar[j]); }
#endif
        }

        /* Maxima/Minima — element-wise MAX/MIN. */
        MERGE_MAX(Maxima.Density);  MERGE_MIN(Minima.Density);
        MERGE_MAX(Maxima.Pressure); MERGE_MIN(Minima.Pressure);
        for(int j = 0; j < 3; j++) { MERGE_MAX(Maxima.Velocity[j]); MERGE_MIN(Minima.Velocity[j]); }
#ifdef MAGNETIC
        for(int j = 0; j < 3; j++) { MERGE_MAX(Maxima.B[j]); MERGE_MIN(Minima.B[j]); }
#ifdef DIVBCLEANING_DEDNER
        MERGE_MAX(Maxima.Phi); MERGE_MIN(Minima.Phi);
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        for(int j = 0; j < NUM_METAL_SPECIES; j++) { MERGE_MAX(Maxima.Metallicity[j]); MERGE_MIN(Minima.Metallicity[j]); }
#endif
#if defined(RT_COMPGRAD_EDDINGTON_TENSOR) && (N_RT_FREQ_BINS > 0)
        for(int j = 0; j < N_RT_FREQ_BINS; j++) {
            MERGE_MAX(Maxima.Rad_E_gamma[j]); MERGE_MIN(Minima.Rad_E_gamma[j]);
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
            for(int kd = 0; kd < 3; kd++) { MERGE_MAX(Maxima.Rad_Flux[j][kd]); MERGE_MIN(Minima.Rad_Flux[j][kd]); }
#endif
        }
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
        MERGE_MAX(Maxima.InternalEnergy); MERGE_MIN(Minima.InternalEnergy);
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
        MERGE_MAX(Maxima.ElectronNumberDensity); MERGE_MIN(Minima.ElectronNumberDensity);
        MERGE_MAX(Maxima.ElectronTemperature);   MERGE_MIN(Minima.ElectronTemperature);
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
        for(int j = 0; j < 3; j++) { MERGE_MAX(Maxima.E_battery_T2[j]); MERGE_MIN(Minima.E_battery_T2[j]); }
#endif
#ifdef COSMIC_RAY_FLUID
        for(int j = 0; j < N_CR_PARTICLE_BINS; j++) { MERGE_MAX(Maxima.CosmicRayPressure[j]); MERGE_MIN(Minima.CosmicRayPressure[j]); }
#endif
#ifdef DOGRAD_SOUNDSPEED
        MERGE_MAX(Maxima.SoundSpeed); MERGE_MIN(Minima.SoundSpeed);
#endif
#ifdef TURB_DIFF_DYNAMIC
        for(int j = 0; j < 3; j++) { MERGE_MAX(Maxima.Velocity_bar[j]); MERGE_MIN(Minima.Velocity_bar[j]); }
#endif

        MERGE_MAX(MaxDistance);

#if defined(KERNEL_CRK_FACES)
        MERGE_ADD(m0);
        for(int k = 0; k < 3; k++) {
            MERGE_ADD(m1[k]); MERGE_ADD(dm0[k]);
            for(int kx = 0; kx < 3; kx++) { MERGE_ADD(dm1[k][kx]); }
        }
        for(int k = 0; k < 6; k++) {
            MERGE_ADD(m2[k]);
            for(int kx = 0; kx < 3; kx++) { MERGE_ADD(dm2[k][kx]); }
        }
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
        for(int j = 0; j < 3; j++) { MERGE_ADD(GlassAcc[j]); }
#endif
#ifdef HYDRO_SPH
#ifdef MAGNETIC
        for(int j = 0; j < 3; j++) { MERGE_ADD(DtB[j]); }
#ifdef DIVBCLEANING_DEDNER
        MERGE_ADD(divB);
#endif
#endif
        MERGE_ADD(alpha_limiter);
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
        for(int j = 0; j < 3; j++) {
            MERGE_ADD(Face_Area[j]);
            for(int k = 0; k < 3; k++) { MERGE_ADD(FaceCrossX[j][k]); }
        }
        MERGE_ADD(FaceDotB);
#endif
#ifdef TURB_DIFF_DYNAMIC
        for(int j = 0; j < 3; j++) { MERGE_ADD(Velocity_hat[j]); }
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
        MERGE_ADD(AGS_zeta);
#endif

#undef MERGE_ADD
#undef MERGE_MAX
#undef MERGE_MIN
    }

};


#ifdef MHD_CONSTRAINED_GRADIENT
/* Slim Spec for the constrained-gradient iterations (grad_iter>0). Identical to
 * GradientsSpec except it accumulates the slim GasGraddata_out_iter_ (FaceDotB +
 * MIDPOINT PhiGrad) instead of the full GasGraddata_out_ — the only fields the
 * iter>0 writeback (out2particle_GasGrad_iter) consumes. This drops the dead
 * ~1968 B reply payload (Mode B) / accumulator traffic (Mode A) on every CG
 * iteration. Physics-preserving: pass output is bitwise-identical (see
 * gradient_accumulate_neighbor_iter). All non-accumulator contract is aliased
 * from GradientsSpec (SSOT); only AccumData + zero/merge/compare/writeback +
 * populate_device_context + the slim pair body differ. Dispatch: grad_iter==0
 * -> GradientsSpec, grad_iter>0 -> GradientsIterSpec (gradients_loop.cc). The
 * runner learns nothing about passes. */
struct GradientsIterSpec
{
    static constexpr const char *loop_name = "gradients";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only; no j-write, no atomics */

    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int            neighbor_type_mask = (1u << 0);
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;
    static constexpr WritePattern   write_pattern              = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind            = SidxCacheKind::GasOnly;
    static constexpr bool           uses_ghost_writeback       = false;
    static constexpr bool           uses_ghost_write_detector  = false;

    using CallScalars   = NlrCommonScalars;
    using ActiveData    = GradientsActiveData;
    using AccumData     = struct GasGraddata_out_iter_;
    using NeighborData  = GradientsNeighborData;
    using Aux           = GradientsAux;

    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;
    using DeviceContext  = GradientsDeviceContext;

    /* ----- Aliased host hooks (SSOT = GradientsSpec) ----- */
    static bool is_active(int i) { return GradientsSpec::is_active(i); }
    static double search_radius(const neighbor_loop_args& args, int active_slot, int i) {
        return GradientsSpec::search_radius(args, active_slot, i);
    }
    static CallScalars populate_call_scalars(const neighbor_loop_args& args) {
        return GradientsSpec::populate_call_scalars(args);
    }

    /* DeviceContext extension hook: ferry grad_iter (own body, gradients_loop.cc). */
    static void populate_device_context(const neighbor_loop_args& args, DeviceContext& ctx);

    /* ----- Device hooks ----- */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum)
    {
        for(size_t b = 0; b < sizeof(accum); b++) ((char*)&accum)[b] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx, int active_slot, int i,
                                   double h_search, const CallScalars& scalars)
    {
        return GradientsSpec::load_active(ctx, active_slot, i, h_search, scalars);
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx, int j,
                                       const IdentitySidecar& id, const ActiveData& active)
    {
        return GradientsSpec::load_neighbor(ctx, j, id, active);
    }

    /* Slim pair body — forwards to gradient_accumulate_neighbor_iter. */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active, const NeighborData& neighbor,
                             AccumData& accum, NoScatter& /*scatter*/,
                            const CallScalars& /*cs*/)
    {
        if(!active.enabled) return;
#if defined(HYDRO_MULTIFLUID)
        if (!same_lagrangian_fluid_id(active.FluidType, neighbor.P[neighbor.j].FluidType)) return;
#endif
        struct kernel_GasGrad kernel;
        kernel.h_i = active.local.KernelRadius;
        gradient_accumulate_neighbor_iter(
            const_cast<struct GasGraddata_in_*>(&active.local),
            &accum, &kernel, neighbor.j,
            active.sph_gradients_flag_i, active.V_i,
            active.hinv, active.hinv3, active.hinv4, active.kernel_mode_i,
            neighbor.P, neighbor.CellP);
    }

    /* Slim writeback / merge (own bodies, gradients_loop.cc). */
    static void apply_active_writeback(const neighbor_loop_args& args, int active_slot, int i,
                                        const AccumData& accum);
    /* Peer-rank accum merge (Mode B remote): additive, matching GradientsSpec's
     * MERGE_ADD(FaceDotB) + MERGE_ADD(Gradients[k].Phi). */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& dst, const AccumData& src)
    {
        dst.FaceDotB += src.FaceDotB;
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
        for(int k = 0; k < 3; k++) { dst.PhiGrad[k] += src.PhiGrad[k]; }
#endif
    }
};
#endif /* MHD_CONSTRAINED_GRADIENT */

#endif /* GRADIENTS_LOOP_H */
