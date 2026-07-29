/* hydro/density_loop.h — hydro density runner Spec.
 *
 * Defines DensitySpec for the iterative neighbor-loop runner
 * (mesh/neighbor_loop_runner.h). This replaces the legacy density()
 * driver and retired density_evaluate_gpu path.
 *
 * Contract notes:
 *   - pair_kernel signature matches runner contract verbatim: 4 args
 *     (active, neighbor, accum, scalars). NO `int j` — runner doesn't
 *     pass it; if a body needs j it lives in NeighborData.
 *   - AccumData fields LITERALLY match density_functions.h:50-104
 *     density_evaluate_data_out_ (literal field-by-field port):
 *     NV_T_face_weights/ParticleVel/GradH_numer = Vec3<MyDouble>;
 *     NV_D/NV_A = MyFloat[3][3]; GradRho/Gas_B = Vec3<MyFloat>/Vec3<MyDouble>;
 *     Sink_dr_to_NearestGasNeighbor nested under BH_ACCRETE_NEARESTFIRST
 *     || SINGLE_STAR_TIMESTEPPING; GasVel under TURB_DRIVING || GRAIN_FLUID
 *     (NOT GRAIN_FLUID alone); no invented Gas_Density / SmoothedVel.
 *   - ActiveData mirrors density_functions.h:33-46 density_evaluate_data_in_:
 *     Pos=Vec3<MyDouble>, Vel=Vec3<MyFloat>, Type=int, Accel=Vec3<MyFloat>
 *     under SPHAV_CD10 only, DelayTime under GALSF_SUBGRID_WINDS. No
 *     TimeBin (j-side via NeighborData); no origin_* (runner-side via
 *     envelopes, not Spec-side).
 *   - CallScalars contains only the fields NOT already in
 *     NlrCommonScalars (cf_atime/cf_a2inv/comoving_integration_on come
 *     via common.*; pair-kernel body reads cs.common.cf_a2inv). 10 fields.
 *   - zero_accum is per-field, NOT memset: Sink_TimeBinGasNeighbor =
 *     TIMEBINS sentinel; Sink_dr_to_NearestGasNeighbor = MAX_REAL_NUMBER;
 *     all others = 0.
 *   - uses_ghost_writeback = false: density is pure i-side accumulation.
 *     Compile-time static_assert holds the invariant.
 *   - apply_active_writeback is FINAL-ONLY for iterative Specs (runner
 *     header lines 1272-1282).
 *   - mode_a_csr_buffer_factor = 1.3 == legacy DENSITY_H_BUFFER_FACTOR.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef DENSITY_LOOP_H
#define DENSITY_LOOP_H

/* Kokkos_Core must precede allvars.h (its macros may conflict with stdlib names). */
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"
#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"   /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
/* NOTE: caller TUs must include "../mesh/kernel.h" before this header.
 * kernel.h has no include guards; provides kernel_main / kernel_hinv used
 * by the inline pair body. */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

#include <vector>

/* Forward decl — defined in hydro/density_loop.cc. */
int density_isactive(int i);

/* ============================================================================
 * (1) CallScalars — deduplicated against NlrCommonScalars.
 *
 * NlrCommonScalars (mesh/neighbor_loop_runner.h:350) already carries:
 *   cf_atime, cf_a2inv, cf_a3inv, cf_hubble_a, newton_G, hubble,
 *   comoving_integration_on.
 *
 * Inside iter loop / device hooks: `cs.common.cf_a2inv` (SPHAV_CD10
 * pair-kernel reads), `cs.common.comoving_integration_on` (active
 * predicate star-type gates) come from `common`. This struct adds only
 * the fields NOT in NlrCommonScalars: time-counters, bisection
 * targets/clamps, sink overrides.
 *
 * Post-runner host-only finalization (density_finalize_post_runner) does
 * NOT read CallScalars — it reads host_all->cf_atime, host_all->cf_a2inv,
 * host_all->Time, etc. directly via nlr_host_all_ptr().
 * CallScalars is the device-facing / per-iter-hot-path channel; host_all
 * is the host-only-finalize channel. Keep the split — that's the rule
 * the runner contract uses to separate device-hot-path scalars from
 * host-only finalization.
 *
 * Field types are exact legacy types from global_data_all_struct.h. */
struct DensityCallScalars {
    NlrCommonScalars common;

    /* Time counters (not in NlrCommonScalars) */
    integertime  ti_current;         /* All.Ti_Current (sink active predicate) */
    double       time_now;           /* All.Time */
    double       time_begin;         /* All.TimeBegin */
    double       time_step;          /* All.TimeStep (extra-physics-gas gate) */

    /* Bisection target + clamps */
    double       des_num_ngb;        /* All.DesNumNgb */
    double       max_num_ngb_dev;    /* All.MaxNumNgbDeviation */
    double       min_kernel_radius;  /* All.MinKernelRadius */
    double       max_kernel_radius;  /* All.MaxKernelRadius */

    /* Sink-specific bisection overrides (only present when sinks compiled). */
#ifdef SINK_PARTICLES
    double       sink_ngb_factor;            /* All.SinkNgbFactor */
    double       sink_max_accretion_radius;  /* All.SinkMaxAccretionRadius */
#endif
};

/* ============================================================================
 * (2) AccumData — literal port of density_evaluate_data_out_
 *     (density_functions.h:50-104). EVERY `#ifdef` branch present in the
 *     legacy struct ports here verbatim. Field order matches legacy.
 *     Types match legacy (MyFloat vs MyDouble; Vec3<...> vs array).
 *
 * Merge semantics:
 *   - all fields ADD-reduced EXCEPT
 *   - Sink_TimeBinGasNeighbor: MIN, init = TIMEBINS
 *   - Sink_dr_to_NearestGasNeighbor: MIN, init = MAX_REAL_NUMBER
 * zero_accum (below) is per-field, NOT memset. ========================== */
struct DensityAccumData {
    MyDouble                   Ngb;
    MyDouble                   Rho;
    MyDouble                   DrkernNgb;
    MyDouble                   Particle_DivVel;
    SymmetricTensor2<MyDouble> NV_T;
    Vec3<MyDouble>             NV_T_face_weights;

#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
    Vec3<MyDouble>             ParticleVel;
#endif

#ifdef HYDRO_SPH
    MyDouble                   DrkernHydroSumFactor;
#endif

#ifdef RT_SOURCE_INJECTION
    MyDouble                   KernelSum_Around_RT_Source;
#endif

#ifdef HYDRO_PRESSURE_SPH
    MyDouble                   EgyRho;
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    MyFloat                    NV_D[3][3];
    MyFloat                    NV_A[3][3];
#endif

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    Vec3<MyFloat>              GradRho;
#endif

#if defined(SINK_PARTICLES)
    int                        Sink_TimeBinGasNeighbor;       /* MIN, init = TIMEBINS */
  #if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
    MyDouble                   Sink_dr_to_NearestGasNeighbor; /* MIN, init = MAX_REAL_NUMBER */
  #endif
#endif

#if defined(TURB_DRIVING) || defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION)
    Vec3<MyDouble>             GasVel;
#endif

#if defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION)
    MyDouble                   AmbientGasRho;
    MyDouble                   Gas_InternalEnergy;
  #if defined(DO_FLUID_DRAG_CALCULATION_WITHBFIELDS)
    Vec3<MyDouble>             Gas_B;
  #endif
  #if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
    MyDouble                   Gas_VolatileSpecies[GRAIN_NUM_VOLATILE_SPECIES];  /* precompiler_logic.h:1247 #define = 4 */
  #endif
#endif

#ifdef HYDRO_PARTITION_UNITY_IMPROVE_FD
    Vec3<MyDouble>             GradH_numer;
    MyDouble                   GradH_denom;
#endif
};

/* ============================================================================
 * (3) IterScratch — per-active state living across iterations.
 *
 * Legacy bisection brackets (Left[i], Right[i] at density.cc:149-150) +
 * a converged latch (replaces the legacy P[i].TimeBin negation pattern
 * at density.cc:479/591/617, for oracle compatibility) + iter-count
 * tracking for the iter>10 underflow warning at density.cc:441-452.
 *
 * Mutated host-side inside after_iter only. Not visible inside pair_kernel. */
struct DensityIterScratch {
    MyFloat  left;
    MyFloat  right;
    bool     converged;
    int      iter_count;
    bool     underflow_warned;
    /* CN tracking: legacy density.cc reads CellP[i].ConditionNumber as
     * the iter-loop's "effective" CN, and at iter==0 may update CellP[i]
     * to THIS iter's value (lines 313-322). Our design forbids CellP
     * writes inside after_iter; we mirror the iter-effective CN in
     * IterScratch instead. density_finalize_post_runner writes
     * CellP[i].ConditionNumber at the end, from the converged
     * accum.NV_T re-inversion. */
    double   condition_number_current;
    bool     cn_initialized;
};

/* ============================================================================
 * (4) ActiveData — literal port of density_evaluate_data_in_
 *     (density_functions.h:33-46).
 *
 * Types match legacy verbatim where the kernel needs legacy data. Fields
 * owned by the runner, such as origin rank/index, stay out of ActiveData.
 *
 * scalars carries the CallScalars snapshot per AGS pattern so
 * pair_kernel doesn't double-thread the args. ============================ */
struct DensityActiveState {
    /* pos is lowercase to match the runner's walker which reads
     * `active.pos[k]` (mesh/neighbor_loop_runner.cc:793). All other
     * fields use legacy density_evaluate_data_in_ capitalization. */
    Vec3<MyDouble>     pos;
    Vec3<MyFloat>      Vel;
    MyFloat            h_search;       /* legacy KernelRadius; runner manages via radii_uvm */
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    Vec3<MyFloat>      Accel;          /* P[i].GravAccel + CellP[i].HydroAccel snapshot */
#endif
#ifdef GALSF_SUBGRID_WINDS
    MyFloat            DelayTime;
#endif
    int                Type;
#ifdef HYDRO_MULTIFLUID
    unsigned char      FluidType;  /* packed P[i].FluidType — for same_lagrangian_fluid_id() */
#endif
    DensityCallScalars scalars;
};

/* NeighborData — j-side per-pair sidecar. Density has NO j-side writes
 * (verified §E inventory + design §1f static_assert), so NeighborData
 * is read-only pointers to the j-side P / CellP. neighbor_index is carried
 * only for diagnostics/assertions, not as an extra pair_kernel argument. */
struct DensityNeighborData {
    struct particle_data *neighbor_particle;
    struct gas_cell_data *neighbor_cell;     /* nullptr for non-gas / when no CellP */
    int                   neighbor_index;    /* j; diagnostic / assertions */
};

/* ============================================================================
 * (5) DensitySpec — the runner-template Spec.
 *
 * Declares the density loop traits, data structs, and hook surface. ===== */
struct DensitySpec {
    /* Identity */
    static constexpr const char *loop_name = "density";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only; no j-write, no atomics -> threaded eval bit-identical */

    /* Search policy: legacy density CSR builder uses NGB_SEARCH_ONEWAY
     * with j-bitmask = 1 (gas only) per density_gpu.cc:256. Pair predicate
     * is h_i only (no h_j check). Symmetric search would inflate NumNgb. */
    static constexpr int                     search_mode        = MODE_B_SEARCH_ONEWAY;
    static constexpr unsigned int            neighbor_type_mask = (1u << 0);  /* gas-only */
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;

    /* Write policy */
    static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
    /* SIDX cache: GAS-ONLY. An earlier draft used
     * AllTypes here, which violated the cache-mask invariant — the
     * all-types cache's compact_xyzh / pool includes DM (Type=1) which
     * a gas-only neighbor_type_mask walk should never see. Walker would
     * return DM neighbors, lazy-drift on those DM particles, drift_particle
     * abort 'no prediction into past allowed' at SP4 hydro. gpu_ngb_list_build
     * now hard-aborts on cache_tbm vs caller tbm mismatch. */
    static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::GasOnly;
    static constexpr bool mode_a_active_sources_in_sidx_pool = false; /* non-pool active sources (sink/star/grain) -> runner stages explicit P[].Pos */
    static constexpr bool           uses_ghost_writeback      = false;
    static constexpr bool           uses_ghost_write_detector = false;

    /* Convergence + iter policy */
    static constexpr double accum_tolerance  = 1e-8;
    static constexpr double radius_tolerance = 1e-9;
    using                    IterControl = Iterative;
    using                    IterScratch = DensityIterScratch;
    static constexpr int     max_iters   = MAXITER;  /* legacy density.cc bisection bound */

    /* Mode A CSR policy — exact legacy constant (density_gpu.cc:103
     * DENSITY_H_BUFFER_FACTOR 1.3). */
    static constexpr double mode_a_csr_buffer_factor      = 1.3;
    /* Fallback rule per design v0.4 §1c: flip mode_a_rebuild_csr_every_iter
     * to true ONLY if PA-A vs PB-vs-legacy validation diverges and the
     * cause is the CSR reuse. Don't tune performance before correctness. */
    static constexpr bool   mode_a_rebuild_csr_every_iter = false;

    /* Subgroup policy: single j-mask subgroup, gas-only. SupportsSubgroups
     * is OMITTED (default false_type — runner runtime-asserts single
     * subgroup at entry, TRAP 9). Per-i type branching for non-gas actives
     * lives inside the pair_kernel body. */

    /* Type aliases */
    using CallScalars   = DensityCallScalars;
    using ActiveData    = DensityActiveState;
    using AccumData     = DensityAccumData;
    using NeighborData  = DensityNeighborData;
    /* DeviceContext: density needs no extension beyond the base.
     * AGS extends with need_wakeup_uvm + oracle_dry_run flag; density's
     * pure i-side accumulation + oracle-clean after_iter contract means
     * no per-call device-side scratch is needed. */
    using DeviceContext = NeighborLoopDeviceContextBase;
    /* ScatterData: required by runner contract. Density has no j-side
     * scatter (uses_ghost_writeback=false, pure i-side accumulation),
     * so NoScatter — the runner's empty marker type. */
    using ScatterData   = NoScatter;

    /* Aux — host-owned final-accum storage.
     * apply_active_writeback (final-only for iterative Specs per runner
     * header lines 1272-1282) copies the converged AccumData into
     * per_active_final_accum[active_slot]. density_finalize_post_runner
     * consumes it for final P[]/CellP[] work. Host-only — std::vector OK. */
    struct Aux {
        /* Converged AccumData per active, copied at the post-iter-loop
         * apply_active_writeback (FINAL-only). */
        std::vector<AccumData> per_active_final_accum;
        /* Converged radius per active, copied in after_iter at the
         * Converged return paths. Needed by density_finalize_post_runner
         * to write P[i].KernelRadius; the runner does not surface
         * radii_uvm to the caller and after_iter cannot write P[i] under
         * the oracle contract. */
        std::vector<double>    per_active_final_h;
    };

    /* J-side write invariant. */
    static_assert(uses_ghost_writeback == false,
        "density() must remain pure i-side accumulation; if a j-side "
        "write is added, register a ghost_writeback op explicitly");

    /* ====================================================================
     * Host hooks. Bodies in hydro/density_loop.cc.
     * ==================================================================== */
    static bool        is_active(int particle_index) { return density_isactive(particle_index) != 0; }
    static double      search_radius(const neighbor_loop_args& args, int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);

    static void apply_active_writeback(const neighbor_loop_args& args,
                                       int active_slot, int i,
                                       const AccumData& accum);
    /* Iterative variant of apply_active_writeback. Runner calls this instead of apply_active_writeback in
     * the iterative production-only post-iter-loop writeback, passing
     * the converged radius and IterScratch. Oracle-safe: runner does NOT
     * invoke this on its oracle accum buffer, so args.aux can be written
     * from here without contamination. Density implementation copies
     * both `accum` and `final_h` into Aux for density_finalize_post_runner. */
    static void apply_active_writeback_iterative(const neighbor_loop_args& args,
                                                  int active_slot, int i,
                                                  const AccumData& accum,
                                                  double final_h,
                                                  const IterScratch& scratch);
    static void merge_accum(AccumData& local, const AccumData& peer);

    /* Oracle brute-pass toggle (NlrOracleBrutePassGuard requires the
     * symbol unconditionally — runner.h:505 calls without SFINAE guard).
     * Density's contract allows oracle (clean after_iter + clean
     * apply_active_writeback_iterative); the body is a no-op since
     * density has no j-side writes to suppress during oracle. */
    static void set_oracle_brute_pass(DeviceContext& ctx, bool on);

    /* Per-active oracle comparison (runner uses this to emit
     * ORACLE MISMATCH lines). Returns max relative diff across all
     * AccumData fields. Per-field compare (not a byte-walk) because
     * AccumData mixes MyDouble/MyFloat/int. */
    static double compare_accum(const AccumData& local, const AccumData& oracle);

    static IterResult after_iter(const AfterIterContext<DensitySpec>& ctx,
                                 const AccumData& accum);
    /* after_iter_global NOT declared — density has no per-rank
     * post-iter side-effect (bisection state is per-active, no MPI
     * inside the iter loop). Runner SFINAE-detects absence. */

    /* ====================================================================
     * Device hooks. Header-inline because the runner instantiates the
     * pair body from GPU translation units.
     * ==================================================================== */

    /* Per-field init with MIN-reduction sentinels for Sink fields. */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        accum.Ngb              = 0;
        accum.Rho              = 0;
        accum.DrkernNgb        = 0;
        accum.Particle_DivVel  = 0;
        for(int a = 0; a < 6; ++a) accum.NV_T.data[a] = 0;
        accum.NV_T_face_weights[0] = 0;
        accum.NV_T_face_weights[1] = 0;
        accum.NV_T_face_weights[2] = 0;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
        accum.ParticleVel[0] = 0; accum.ParticleVel[1] = 0; accum.ParticleVel[2] = 0;
#endif
#ifdef HYDRO_SPH
        accum.DrkernHydroSumFactor = 0;
#endif
#ifdef RT_SOURCE_INJECTION
        accum.KernelSum_Around_RT_Source = 0;
#endif
#ifdef HYDRO_PRESSURE_SPH
        accum.EgyRho = 0;
#endif
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        for(int a = 0; a < 3; ++a) for(int b = 0; b < 3; ++b) {
            accum.NV_D[a][b] = 0; accum.NV_A[a][b] = 0;
        }
#endif
#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
        accum.GradRho[0] = 0; accum.GradRho[1] = 0; accum.GradRho[2] = 0;
#endif
#if defined(SINK_PARTICLES)
        accum.Sink_TimeBinGasNeighbor       = TIMEBINS;          /* MIN-reduction sentinel */
  #if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
        accum.Sink_dr_to_NearestGasNeighbor = MAX_REAL_NUMBER;   /* MIN-reduction sentinel */
  #endif
#endif
#if defined(TURB_DRIVING) || defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION)
        accum.GasVel[0] = 0; accum.GasVel[1] = 0; accum.GasVel[2] = 0;
#endif
#if defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION)
        accum.AmbientGasRho = 0;
        accum.Gas_InternalEnergy = 0;
  #if defined(DO_FLUID_DRAG_CALCULATION_WITHBFIELDS)
        accum.Gas_B[0] = 0; accum.Gas_B[1] = 0; accum.Gas_B[2] = 0;
  #endif
  #if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
        for(int kv = 0; kv < GRAIN_NUM_VOLATILE_SPECIES; ++kv) accum.Gas_VolatileSpecies[kv] = 0;
  #endif
#endif
#ifdef HYDRO_PARTITION_UNITY_IMPROVE_FD
        accum.GradH_numer[0] = 0; accum.GradH_numer[1] = 0; accum.GradH_numer[2] = 0;
        accum.GradH_denom = 0;
#endif
    }

    /* load_active — snapshot P[i] / CellP[i] fields into ActiveData.
     * Reads through ctx.P / ctx.CellP (path-correct: Mode A arena P_gpu,
     * Mode B request-driven slab).
     *
     * Two legacy-semantics requirements are important because the pair
     * kernel inherits whatever load_active produces:
     *
     *   1. Vel: legacy (density_gpu.cc particle2in) uses CellP[i].VelPred
     *      for gas (Type==0) and P[i].Vel for non-gas. Pair kernel reads
     *      `local->Vel` as the i-side velocity in divergence/curl terms.
     *      Using P[i].Vel uniformly silently corrupts gas density physics.
     *
     *   2. Accel (SPHAV_CD10 only): legacy stores the PRE-SCALED form
     *      `All.cf_a2inv * P[i].GravAccel + CellP[i].HydroAccel`
     *      (density_gpu.cc:310). The cf_a2inv multiply is component-
     *      specific (gravity term only), not applied to the sum.
     *      density_functions.h:153-161 then subtracts the j-side analogue
     *      inline. Storing the raw sum here and "applying cf_a2inv later"
     *      breaks the component separation. Pre-scale at load_active. */
    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const NeighborLoopDeviceContextBase& ctx,
                                  int active_slot, int i,
                                  double h_search, const CallScalars& scalars) {
        ActiveData a{};
        a.pos[0]    = (MyDouble)ctx.P[i].Pos[0];
        a.pos[1]    = (MyDouble)ctx.P[i].Pos[1];
        a.pos[2]    = (MyDouble)ctx.P[i].Pos[2];
        a.h_search  = (MyFloat) h_search;
        a.Type      = (int)     ctx.P[i].Type;

        /* Vel: gas uses CellP[i].VelPred; non-gas uses P[i].Vel
         * Mirrors density_gpu.cc particle2in. */
        const bool i_is_gas = (ctx.P[i].Type == 0) && (ctx.CellP != nullptr);
        if (i_is_gas) {
            a.Vel[0] = (MyFloat)ctx.CellP[i].VelPred[0];
            a.Vel[1] = (MyFloat)ctx.CellP[i].VelPred[1];
            a.Vel[2] = (MyFloat)ctx.CellP[i].VelPred[2];
        } else {
            a.Vel[0] = (MyFloat)ctx.P[i].Vel[0];
            a.Vel[1] = (MyFloat)ctx.P[i].Vel[1];
            a.Vel[2] = (MyFloat)ctx.P[i].Vel[2];
        }

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        /* Accel: PRE-SCALED legacy form for gas; zero for non-gas
         * (SPHAV_CD10 NV_A accumulator only fires on gas pairs in the
         * legacy kernel; safe zero for non-gas i means no contribution
         * if any branch fires). */
        if (i_is_gas) {
            a.Accel[0] = (MyFloat)(scalars.common.cf_a2inv * ctx.P[i].GravAccel[0] + ctx.CellP[i].HydroAccel[0]);
            a.Accel[1] = (MyFloat)(scalars.common.cf_a2inv * ctx.P[i].GravAccel[1] + ctx.CellP[i].HydroAccel[1]);
            a.Accel[2] = (MyFloat)(scalars.common.cf_a2inv * ctx.P[i].GravAccel[2] + ctx.CellP[i].HydroAccel[2]);
        } else {
            a.Accel[0] = 0; a.Accel[1] = 0; a.Accel[2] = 0;
        }
#endif

#ifdef GALSF_SUBGRID_WINDS
        a.DelayTime = i_is_gas ? (MyFloat)ctx.CellP[i].DelayTime : (MyFloat)0;
#endif

#ifdef HYDRO_MULTIFLUID
        a.FluidType = ctx.P[i].FluidType;
#endif

        a.scalars = scalars;
        (void)active_slot;
        return a;
    }

    /* load_neighbor — route ctx.P[j] / ctx.CellP[j] pointers.
     *
     * Density requests gas-only neighbors. Imported gas ghosts carry valid
     * CellP state in the runner's active P/CellP view, so do not guard by
     * local gas capacity here; doing so silently drops cross-rank gas
     * neighbors. */
    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const NeighborLoopDeviceContextBase& ctx,
                                      int j, const IdentitySidecar& id,
                                      const ActiveData& active) {
        NeighborData n{};
        n.neighbor_particle = &ctx.P[j];
        n.neighbor_cell     = (ctx.CellP && ctx.P[j].Type == 0)
                              ? &ctx.CellP[j] : nullptr;
        n.neighbor_index    = j;
        (void)id; (void)active;
        return n;
    }

    /* pair_kernel — literal port of density_accumulate_neighbor +
     * density_evaluate_extra_physics_gas (density_functions.h:104-249).
     * Every legacy `#ifdef` branch ports verbatim per §5 + the no-#error
     * rule (feedback_port_all_branches_no_error_endstate.md).
     *
     * Signature matches runner contract: 4 args (no `int j` — see
     * mesh/neighbor_loop_runner.cc:854/1825/3147). Index j is encoded
     * in neighbor.neighbor_particle/neighbor_cell pointers.
     *
     * j is always gas (neighbor_type_mask = (1u << 0)); the runner's
     * load_neighbor populates neighbor_cell for every j the runner
     * dispatches. Keep the legacy `local->Type == 0` gate around SPHAV_CD10 /
     * TURB_DRIVING / NV_T blocks exactly as-is. */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& i_active,
                            const NeighborData& neighbor,
                            AccumData& accum, NoScatter& /*scatter*/) {

        /* CallScalars are carried by the active snapshot, matching the AGS
         * runner pattern and avoiding device-side global reads. */
        const CallScalars& cs = i_active.scalars;
        struct particle_data       &Pj    = *neighbor.neighbor_particle;
        struct gas_cell_data * const CellPj = neighbor.neighbor_cell;
        /* j is always gas under neighbor_type_mask = (1u<<0); CellPj is
         * non-null for every j the runner walks. The early-return below
         * is a defensive sanity check, valid ONLY while neighbor_type_mask
         * stays gas-only. If a future contributor widens the j-mask, this
         * early return would silently drop non-gas neighbors instead of
         * processing them — rewrite the body to handle non-gas j first,
         * THEN drop the gate. */
        if (CellPj == nullptr) return;

#ifdef GALSF_SUBGRID_WINDS
        if (CellPj->DelayTime > 0) { if (i_active.DelayTime <= 0) return; }
#endif
        if (Pj.Mass <= 0) return;

        /* Position delta with periodic wrap (legacy line 195-197). */
        Vec3<double> dp;
        dp[0] = (double)i_active.pos[0] - (double)Pj.Pos[0];
        dp[1] = (double)i_active.pos[1] - (double)Pj.Pos[1];
        dp[2] = (double)i_active.pos[2] - (double)Pj.Pos[2];
        nearest_xyz(dp);
        const double r2 = dp.norm_sq();
        const double h  = (double)i_active.h_search;
        const double h2 = h * h;
        if (r2 >= h2) return;

        /* Kernel evaluation (legacy line 200-204). */
        const double r = sqrt(r2);
        double hinv, hinv3, hinv4;
        kernel_hinv(h, &hinv, &hinv3, &hinv4);
        const double u = r * hinv;
        double wk, dwk;
        kernel_main(u, hinv3, hinv4, &wk, &dwk, 0);
        const double mass_j = (double)Pj.Mass;
        const double mj_wk  = mass_j * wk;

#ifdef DO_FLUID_ALTSPECIES_DRAG_CALCULATION
        int is_valid_fordrag = IS_PARTICLE_DRAGVALID(i_active.Type, i_active.FluidType);
#ifdef HYDRO_MULTIFLUID
        is_valid_fordrag = (is_valid_fordrag && neighbor.neighbor_particle->FluidType == FLUID_DEFAULT);
#endif
        if (is_valid_fordrag) {
            /* Cross-fluid ambient gas sampling for cross-fluid sources. */
            if (r > 0) {
                accum.AmbientGasRho += mj_wk;
                Vec3<double> dv;
                dv[0] = (double)i_active.Vel[0] - (double)CellPj->VelPred[0];
                dv[1] = (double)i_active.Vel[1] - (double)CellPj->VelPred[1];
                dv[2] = (double)i_active.Vel[2] - (double)CellPj->VelPred[2];
                {
                    Vec3<MyDouble> pos_i_md;
                    pos_i_md[0] = i_active.pos[0]; pos_i_md[1] = i_active.pos[1]; pos_i_md[2] = i_active.pos[2];
                    NGB_SHEARBOX_BOUNDARY_VELCORR_(pos_i_md, Pj.Pos, dv, 1);
                }
                accum.Gas_InternalEnergy += mj_wk * (double)CellPj->InternalEnergyPred;
                accum.GasVel[0] += mj_wk * ((double)i_active.Vel[0] - dv[0]); /* shares with TURB_DRIVING SmoothedVel; combo unsupported */
                accum.GasVel[1] += mj_wk * ((double)i_active.Vel[1] - dv[1]);
                accum.GasVel[2] += mj_wk * ((double)i_active.Vel[2] - dv[2]);
#if defined(DO_FLUID_DRAG_CALCULATION_WITHBFIELDS)
                accum.Gas_B[0] += wk * (double)CellPj->BPred[0];
                accum.Gas_B[1] += wk * (double)CellPj->BPred[1];
                accum.Gas_B[2] += wk * (double)CellPj->BPred[2];
#endif
#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
                for (int kv = 0; kv < GRAIN_NUM_VOLATILE_SPECIES; ++kv) {
                    accum.Gas_VolatileSpecies[kv] += mj_wk * (double)CellPj->VolatileSpecies[kv];
                }
#endif
            }
        }
#endif
#ifdef HYDRO_MULTIFLUID
        if (!same_lagrangian_fluid_id(i_active.FluidType, neighbor.neighbor_particle->FluidType)) return;
#endif

        /* Always-on accumulators (legacy line 206-217). */
        accum.Ngb += wk;
        accum.Rho += mj_wk;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
        if (i_active.Type == 0 && r == 0) {
            accum.ParticleVel[0] += mj_wk * (double)CellPj->VelPred[0];
            accum.ParticleVel[1] += mj_wk * (double)CellPj->VelPred[1];
            accum.ParticleVel[2] += mj_wk * (double)CellPj->VelPred[2];
        }
#endif
#if defined(RT_SOURCE_INJECTION)
  #if defined(RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION)
        if (cs.time_step == 0)
  #endif
        if ((1 << i_active.Type) & (RT_SOURCES)) {
            accum.KernelSum_Around_RT_Source += 1.0 - u * u;
        }
#endif
        accum.DrkernNgb += -(NUMDIMS * hinv * wk + u * dwk);

#ifdef HYDRO_SPH
        {
            double mass_eff = mass_j;
  #ifdef HYDRO_PRESSURE_SPH
            mass_eff *= (double)CellPj->InternalEnergyPred;
            accum.EgyRho += wk * mass_eff;
  #endif
            accum.DrkernHydroSumFactor += -mass_eff * (NUMDIMS * hinv * wk + u * dwk);
        }
#endif

        /* r > 0 branch (legacy line 228-249). */
        if (r > 0) {
            if (i_active.Type == 0) {
                /* NV_T outer-product accumulation (gas i only). */
                accum.NV_T              += wk * outer_product(dp);
                accum.NV_T_face_weights[0] += wk * dp[0];
                accum.NV_T_face_weights[1] += wk * dp[1];
                accum.NV_T_face_weights[2] += wk * dp[2];
#ifdef HYDRO_PARTITION_UNITY_IMPROVE_FD
                accum.GradH_numer[0] += (dwk / r) * dp[0];
                accum.GradH_numer[1] += (dwk / r) * dp[1];
                accum.GradH_numer[2] += (dwk / r) * dp[2];
                accum.GradH_denom    += u * dwk;
#endif
            }
            /* Velocity delta (legacy line 240). j is always gas; CellP[j].VelPred valid. */
            Vec3<double> dv;
            dv[0] = (double)i_active.Vel[0] - (double)CellPj->VelPred[0];
            dv[1] = (double)i_active.Vel[1] - (double)CellPj->VelPred[1];
            dv[2] = (double)i_active.Vel[2] - (double)CellPj->VelPred[2];
            /* SHEARBOX boundary velocity correction (no-op when not defined). */
            {
                Vec3<MyDouble> pos_i_md;
                pos_i_md[0] = i_active.pos[0]; pos_i_md[1] = i_active.pos[1]; pos_i_md[2] = i_active.pos[2];
                NGB_SHEARBOX_BOUNDARY_VELCORR_(pos_i_md, Pj.Pos, dv, 1);
            }
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
            accum.ParticleVel[0] += mj_wk * ((double)i_active.Vel[0] - dv[0]);
            accum.ParticleVel[1] += mj_wk * ((double)i_active.Vel[1] - dv[1]);
            accum.ParticleVel[2] += mj_wk * ((double)i_active.Vel[2] - dv[2]);
#endif
            accum.Particle_DivVel -= dwk * (dp[0]*dv[0] + dp[1]*dv[1] + dp[2]*dv[2]) / r;

            /* ============ density_evaluate_extra_physics_gas (inlined) ===========
             * Legacy at density_functions.h:103-179. Same `local->Type` gating
             * (Type!=0 vs Type==0) wraps the per-type physics. */
            const double mj_dwk_r = mass_j * dwk / r;

            if (i_active.Type != 0) {
#if defined(SINK_PARTICLES)
                if (i_active.Type == 5) {
                    /* Read-only MIN reductions; the j-side SwallowID/SwallowTime/Sink_Ngb_Flag
                     * writes from the CPU path are sink_swk/sink_feed responsibilities,
                     * NOT density. */
                    short int TimeBin_j = Pj.TimeBin;
                    if (TimeBin_j < 0) { TimeBin_j = -TimeBin_j - 1; }
                    if (accum.Sink_TimeBinGasNeighbor > TimeBin_j) {
                        accum.Sink_TimeBinGasNeighbor = TimeBin_j;
                    }
  #if defined(SINGLE_STAR_TIMESTEPPING)
                    double dr_eff_wtd = Pj.Get_Particle_Size();
                    dr_eff_wtd = sqrt(dr_eff_wtd * dr_eff_wtd + r * r);
                    if ((dr_eff_wtd < accum.Sink_dr_to_NearestGasNeighbor) && (Pj.Mass > 0)) {
                        accum.Sink_dr_to_NearestGasNeighbor = dr_eff_wtd;
                    }
  #endif
                }
#endif
            } else { /* i_active.Type == 0 */
#if defined(TURB_DRIVING)
                accum.GasVel[0] += mj_wk * ((double)i_active.Vel[0] - dv[0]);
                accum.GasVel[1] += mj_wk * ((double)i_active.Vel[1] - dv[1]);
                accum.GasVel[2] += mj_wk * ((double)i_active.Vel[2] - dv[2]);
#endif
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
                /* i_active.Accel is the PRE-SCALED legacy form
                 * (load_active applied cf_a2inv to GravAccel component only).
                 * Per-pair, subtract the j-side analogue inline. */
                {
                    const double gax_j = cs.common.cf_a2inv * (double)Pj.GravAccel[0] + (double)CellPj->HydroAccel[0];
                    const double gay_j = cs.common.cf_a2inv * (double)Pj.GravAccel[1] + (double)CellPj->HydroAccel[1];
                    const double gaz_j = cs.common.cf_a2inv * (double)Pj.GravAccel[2] + (double)CellPj->HydroAccel[2];
                    const double a0 = (double)i_active.Accel[0] - gax_j;
                    const double a1 = (double)i_active.Accel[1] - gay_j;
                    const double a2 = (double)i_active.Accel[2] - gaz_j;
                    accum.NV_A[0][0] += a0 * dp[0] * wk;
                    accum.NV_A[0][1] += a0 * dp[1] * wk;
                    accum.NV_A[0][2] += a0 * dp[2] * wk;
                    accum.NV_A[1][0] += a1 * dp[0] * wk;
                    accum.NV_A[1][1] += a1 * dp[1] * wk;
                    accum.NV_A[1][2] += a1 * dp[2] * wk;
                    accum.NV_A[2][0] += a2 * dp[0] * wk;
                    accum.NV_A[2][1] += a2 * dp[1] * wk;
                    accum.NV_A[2][2] += a2 * dp[2] * wk;

                    accum.NV_D[0][0] += dv[0] * dp[0] * wk;
                    accum.NV_D[0][1] += dv[0] * dp[1] * wk;
                    accum.NV_D[0][2] += dv[0] * dp[2] * wk;
                    accum.NV_D[1][0] += dv[1] * dp[0] * wk;
                    accum.NV_D[1][1] += dv[1] * dp[1] * wk;
                    accum.NV_D[1][2] += dv[1] * dp[2] * wk;
                    accum.NV_D[2][0] += dv[2] * dp[0] * wk;
                    accum.NV_D[2][1] += dv[2] * dp[1] * wk;
                    accum.NV_D[2][2] += dv[2] * dp[2] * wk;
                }
#endif
            }

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
            accum.GradRho[0] += mj_dwk_r * dp[0];
            accum.GradRho[1] += mj_dwk_r * dp[1];
            accum.GradRho[2] += mj_dwk_r * dp[2];
#endif
        } /* end r > 0 */
    }
};

#endif /* DENSITY_LOOP_H */
