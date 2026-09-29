/* sinks/sink_env1_loop.h — sink "first-pass environment" neighbor loop module.
 *
 * Defines SinkEnv1Spec for the neighbor-loop runner (see
 * mesh/neighbor_loop_runner.h). This is the canonical worked example of
 * a NeighborLoopSpec — copy this file's structure when porting a new
 * physics loop onto the runner.
 *
 * Default per-loop layout is two files: this header carries the Spec
 * contract and any inline (host+device) physics that the runner must
 * inline into the device kernel; sinks/sink_env1_loop.cc carries
 * host-only hooks and diagnostics. Unusual loops may need a reason to
 * split further, but the two-file default keeps the directory readable
 * once 5+ loops are migrated.
 *
 * Layout convention inside SinkEnv1Spec:
 *   PHYSICS BLOCK    — edit when changing the loop's physics
 *   ENGINE APPARATUS — touch only when changing the runner contract
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef SINK_ENV1_LOOP_H
#define SINK_ENV1_LOOP_H

/* Kokkos_Core.hpp must precede allvars.h (its macros may conflict with stdlib
 * names); the inline pair body below emits Kokkos::atomic_min. */
#include <Kokkos_Core.hpp>
#include "../declarations/allvars.h"

#ifdef SINK_PARTICLES

#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"      /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
#include "../mesh/ghost_exchange_functions.h" /* gx_pair_accept_wrap_and_test: the canonical pair accept */
#include "sinks_gpu_decls.h"                  /* struct sink_env_gpu_out */
#include "sink_functions.h"                   /* sink_vesc_gpu, sink_check_boundedness_gpu —
                                                 used inline in pair body under
                                                 SINK_GRAVCAPTURE_GAS. */
/* NOTE: caller translation units must include "../mesh/kernel.h" BEFORE
 * this header. kernel.h has no include guards (it defines static inline
 * kernel_main, used by the pair body below). The runner and the sink
 * loop's .cc both include kernel.h first; new callers must do the same. */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* Forward-decl from sinks/sink.h (active-particle predicate). */
int sink_isactive(int i);

/* ============================================================================
 * Per-pair physics types (file scope; PascalCase).
 *
 * SinkEnv1Spec re-exports these via `using CallScalars = SinkEnv1CallScalars;`
 * and `using ActiveData = SinkEnv1ActiveState;`. File-scope so the inline
 * pair body below can reference them without a circular dependency on
 * SinkEnv1Spec.
 * ========================================================================== */

struct SinkEnv1CallScalars {
    NlrCommonScalars common;            /* cf_atime, cf_a2inv, cf_a3inv, G, ... */
    double           sink_radius_grav;  /* SinkParticle_GravityKernelRadius */
    int              owner_task;        /* rank owning the actives this call builds; see load_active */
};

/* Active-particle state passed into the pair body. Trivially copyable for
 * byte-level MPI transfer in the Mode B remote path. Conditional fields
 * mirror the same #ifdefs as the pair body so payload size matches per
 * build config.
 *
 * h_search MUST come from the runner's per-active radii staging (driven
 * by SinkEnv1Spec::search_radius), NOT from P[i].KernelRadius — they
 * happen to be equal in the current code but the pair body is
 * contractually keyed off the caller-supplied radius. */
struct SinkEnv1ActiveState {
    Vec3<double>  pos;                  /* P[i].Pos */
    Vec3<double>  vel;                  /* P[i].Vel */
    MyIDType      id;                   /* P[i].ID — physics identity, not a self-skip */
    MyIDType      claim_token;          /* sink ownership token; see gizmo_sink_claim_token */
    double        h_search;             /* per-active radius */
    double        ags_h;                /* AGS_KernelRadius if defined, else sink_radius_grav */
#if defined(SINK_GRAVCAPTURE_GAS) || (SINK_GRAVACCRETION == 8)
    double        mass;                 /* active sink mass (boundedness/Bondi paths) */
#endif
#if defined(SINK_GRAVCAPTURE_FIXEDSINKRADIUS)
    double        sink_radius;
#endif
#if defined(SINK_RETURN_ANGMOM_TO_GAS)
    Vec3<double>  sink_angmom;          /* Sink_Specific_AngMom */
#endif
    int           origin_local_idx;     /* index into requester's per_active_accum */
    int           origin_rank;          /* requester's MPI rank (Mode B remote) */
};

/* DeviceContext extension. It adds no fields of its own today, but it stays a
 * type distinct from NeighborLoopDeviceContextBase: the runner keys off that
 * difference to decide whether the Spec must supply populate_device_context,
 * so collapsing it would change the hooks this Spec is required to define.
 * Trivially copyable: runner captures by value into device lambda. */
struct SinkEnv1DeviceContext : NeighborLoopDeviceContextBase {
};

/* ============================================================================
 * Inline pair body — single source of truth for sink_env1 per-pair physics.
 * Called from SinkEnv1Spec::pair_kernel, used by all three runner paths
 * (Mode A GPU neighbor list, Mode B local host walker, Mode B remote
 * peer-to-peer).
 *
 * Caller is responsible for the per-active early-return (Mass<=0 || h<=0).
 *
 * J-side writes:
 *   - SINGLE_STAR_SINK_DYNAMICS + SINK_GRAVCAPTURE_GAS: a per-pair atomic min
 *     on neighbor_particle.SwallowTime fires inside the boundedness-passed
 *     branch of the SINK_GRAVCAPTURE_GAS body. Reverse-comm to the home rank
 *     goes through the existing GHOST_WRITEBACK_PARTICLE_MIN(SwallowTime)
 *     bundle in sinks/sink_env1_loop.cc.
 *   - SINK_GRAVCAPTURE_GAS (without SINGLE_STAR_SINK_DYNAMICS):
 *     accum.mass_to_swallow_edd accumulates i-side; no j-write.
 * ========================================================================== */

KOKKOS_INLINE_FUNCTION
static void sink_env1_pair_kernel(const SinkEnv1ActiveState& active,
                                  const SinkEnv1CallScalars& scalars,
                                  struct particle_data& neighbor_particle,
                                  const struct gas_cell_data* neighbor_cell,
                                  struct sink_env_gpu_out& accum)
{
    /* Self-skip / mass / type-5 filter. */
    /* Identity is a separation of zero, not a matching identifier: IDs repeat, so the identifier
     * form of this skip discarded genuine distinct neighbours that happened to share one -- they
     * were dropped from the sink environment entirely. The separation is computed just below. */
    if(neighbor_particle.Mass <= 0 || neighbor_particle.Type == 5) return;

    const double h_i      = active.h_search;
    const double hinv     = 1.0 / h_i;
    const double hinv3    = hinv * hinv * hinv;
    const double ags_h_i  = active.ags_h;

    /* dP, dv with periodic wrap. Sign convention: j - i. */
    Vec3<double> dP;
    dP[0] = (double)neighbor_particle.Pos[0] - active.pos[0];
    dP[1] = (double)neighbor_particle.Pos[1] - active.pos[1];
    dP[2] = (double)neighbor_particle.Pos[2] - active.pos[2];
    Vec3<double> dv;
    dv[0] = (double)neighbor_particle.Vel[0] - active.vel[0];
    dv[1] = (double)neighbor_particle.Vel[1] - active.vel[1];
    dv[2] = (double)neighbor_particle.Vel[2] - active.vel[2];
    nearest_xyz(dP, -1);
    NGB_SHEARBOX_BOUNDARY_VELCORR_(active.pos, neighbor_particle.Pos, dv, -1);
    if(!(dP.norm_sq() > 0)) return;   /* the active sink itself, or a degenerate coincident pair */

    const double wt = (double)neighbor_particle.Mass;

#ifdef SINK_REPOSITION_ON_POTMIN
    if(neighbor_particle.Type != 0 && neighbor_particle.Type != 5) {
        double rfac  = dP.norm_sq() * (10.0 / (h_i * h_i)
                                     + 0.1 / (scalars.sink_radius_grav * scalars.sink_radius_grav));
        double wtfac = wt / (1.0 + rfac);
        if((MyFloat)neighbor_particle.Mass > accum.DF_mmax_particles) accum.DF_mmax_particles = (MyFloat)neighbor_particle.Mass;
        for(int kv = 0; kv < 3; kv++) accum.DF_mean_vel[kv] += wtfac * dv[kv];
        accum.DF_rms_vel += wtfac;
        accum.DF_rms_vel += wtfac;
        accum.DF_rms_vel += wtfac;
    }
#endif

    if(neighbor_particle.Type == 0) {
        accum.Mgas_in_Kernel += wt;
        if(neighbor_cell) { accum.Sink_SurroudingGasInternalEnergy += wt * neighbor_cell->InternalEnergy; }
        Vec3<double> J_gas = cross(dP, dv);
        for(int kv = 0; kv < 3; kv++) accum.Jgas_in_Kernel[kv] += wt * J_gas[kv];
#if defined(SINK_OUTPUT_MOREINFO)
        if(neighbor_cell) { accum.Sfr_in_Kernel += neighbor_cell->Sfr; }
#endif
#if (SINK_GRAVACCRETION >= 5) || defined(SINGLE_STAR_SINK_DYNAMICS) || defined(SINGLE_STAR_TIMESTEPPING)
        for(int kv = 0; kv < 3; kv++) accum.Sink_SurroundingGasVel[kv] += wt * dv[kv];
#endif
#ifdef JET_DIRECTION_FROM_KERNEL_AND_SINK
        for(int kv = 0; kv < 3; kv++) accum.Sink_SurroundingGasCOM[kv] += wt * dP[kv];
#endif
#if defined(SINK_RETURN_ANGMOM_TO_GAS) || defined(SINK_RETURN_BFLUX)
        {
            double u_wb = dP.norm() / DMAX(h_i, (double)neighbor_particle.KernelRadius);
            double wk_wb = 0, dwk_wb = 0;
            if(u_wb < 1) { kernel_main(u_wb, 1.0, 1.0, &wk_wb, &dwk_wb, -1); }
#if defined(SINK_RETURN_ANGMOM_TO_GAS)
            double r2j = dP.norm_sq();
            double Lrj = dot(active.sink_angmom, dP);
            Vec3<double> Ang_pass = active.sink_angmom * r2j - dP * Lrj;
            for(int kv = 0; kv < 3; kv++)
                accum.angmom_prepass_sum_for_passback[kv] += wk_wb * wt * Ang_pass[kv];
#endif
#if defined(SINK_RETURN_BFLUX)
            accum.kernel_norm_topass_in_swallowloop += wk_wb;
#endif
        }
#endif
#if (SINK_GRAVACCRETION == 8)
        if(neighbor_cell) {
            double u_h = dP.norm() / h_i;
            double wk_h = 0, dwk_h = 0;
            if(u_h < 1) { kernel_main(u_h, hinv3, hinv3 * hinv, &wk_h, &dwk_h, -1); }
            double rj = u_h * h_i * scalars.common.cf_atime;
            double csj = neighbor_cell->effective_soundspeed();
            double vdotrj = -dot(dP, dv);
            double vr_mdot = 4 * M_PI * wt * (wk_h * scalars.common.cf_a3inv) * rj * vdotrj;
            if(rj < scalars.sink_radius_grav * scalars.common.cf_atime) {
                double bondi_mdot = 4 * M_PI * scalars.common.newton_G * scalars.common.newton_G * active.mass * active.mass
                    / pow(csj * csj + dv.norm_sq() * scalars.common.cf_a2inv, 1.5)
                    * wt * (wk_h * scalars.common.cf_a3inv);
                vr_mdot = DMAX(vr_mdot, bondi_mdot);
                accum.hubber_mdot_bondi_limiter += bondi_mdot;
            }
            accum.hubber_mdot_vr_estimator    += vr_mdot;
            accum.hubber_mdot_disk_estimator  += wt * wk_h * sqrt(rj) / (neighbor_cell->Density * csj * csj);
        }
#endif
    } else if(is_galsf_stellar_candidate_type(neighbor_particle.Type, scalars.common.comoving_integration_on)) {
        accum.Mstar_in_Kernel += wt;
        Vec3<double> J_star = cross(dP, dv);
        for(int kv = 0; kv < 3; kv++) accum.Jstar_in_Kernel[kv] += wt * J_star[kv];
    } else {
        accum.Malt_in_Kernel += wt;
        Vec3<double> J_alt = cross(dP, dv);
        for(int kv = 0; kv < 3; kv++) accum.Jalt_in_Kernel[kv] += wt * J_alt[kv];
    }

    /* SINK_GRAVCAPTURE_GAS path — boundedness check + SwallowID-based
     * mass-marked-swallow accumulation. Under SINGLE_STAR_SINK_DYNAMICS,
     * the boundedness-passed branch also atomic-min's a per-pair tff into
     * neighbor_particle.SwallowTime (reverse-comm via
     * GHOST_WRITEBACK_PARTICLE_MIN(SwallowTime) bundle in sink_env1_loop.cc). */
#ifdef SINK_GRAVCAPTURE_GAS
#ifdef GRAIN_FLUID
    if(neighbor_particle.Mass > 0 && (neighbor_particle.Type == 0 || ((1<<neighbor_particle.Type) & GRAIN_PTYPES)))
#else
    if(neighbor_particle.Mass > 0 && neighbor_particle.Type == 0)
#endif
    {
        double dr_code = dP.norm();
        double vrel = dv.norm() / scalars.common.cf_atime;
#if defined(MAGNETIC) && defined(GRAIN_LORENTZFORCE)
        if((1<<neighbor_particle.Type) & GRAIN_PTYPES) {
            Vec3<double> B_vec; for(int kv = 0; kv < 3; kv++) B_vec[kv] = neighbor_particle.Gas_B[kv];
            double vrel_dot = dot(dv, B_vec), bmag2 = B_vec.norm_sq();
            vrel = (fabs(vrel_dot) / sqrt(bmag2)) / scalars.common.cf_atime;
        }
#endif
        struct gas_cell_data neighbor_cell_local;
        if(neighbor_particle.Type == 0 && neighbor_cell) { neighbor_cell_local = *neighbor_cell; }
        else { /* zero-init */
            for(size_t b = 0; b < sizeof(neighbor_cell_local); b++) {
                ((char*)&neighbor_cell_local)[b] = 0;
            }
        }
        double vbound = sink_vesc_gpu(neighbor_particle, neighbor_cell_local, active.mass, dr_code, ags_h_i);
        if(vrel < vbound) {
            double local_sink_radius = scalars.sink_radius_grav;
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
            local_sink_radius = active.sink_radius;
            double spec_mom = dot(dv, dP);
            double r2 = dP.norm_sq();
            spec_mom = r2*vrel*vrel - spec_mom*spec_mom*scalars.common.cf_a2inv;
            if(spec_mom >= scalars.common.newton_G * (active.mass + (double)neighbor_particle.Mass) * local_sink_radius) { return; }
#endif
            if(sink_check_boundedness_gpu(neighbor_particle, neighbor_cell_local, vrel, vbound, dr_code, local_sink_radius) == 1) {
#ifdef SINGLE_STAR_SINK_DYNAMICS
                const double eps = DMAX(dr_code,
                                        DMAX((double)neighbor_particle.KernelRadius, ags_h_i)
                                        * KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER);
                const double tff_pair = eps * eps * eps
                                        / ((double)active.mass + (double)neighbor_particle.Mass);
                Kokkos::atomic_min(&neighbor_particle.SwallowTime, (MyFloat)tff_pair);
#endif
                if(neighbor_particle.SwallowID < active.claim_token) { accum.mass_to_swallow_edd += (MyFloat)neighbor_particle.Mass; }
            }
        }
    }
#endif  /* SINK_GRAVCAPTURE_GAS */
}

/* ============================================================================
 * SinkEnv1Spec — the NeighborLoopSpec contract for sink_env1.
 *
 * NeighborData lifetime contract:
 *   neighbor.neighbor_particle MUST point into ctx.P (= P_gpu UVM array
 *     under Mode A; = P_host directly in Mode B walker).
 *   neighbor.neighbor_cell MUST be either nullptr OR point into ctx.CellP
 *     (UVM under Mode A).
 *   Pointers are valid only for the duration of one pair_kernel call.
 *   DO NOT store; DO NOT cross arena scope.
 * ========================================================================== */

struct SinkEnv1Spec {
    /* ====================================================================
     * PHYSICS BLOCK — edit when changing the loop's physics
     * ==================================================================== */

    /* (1) Loop identity. Supplies the PHASE0_NLR `loop=` field. */
    static constexpr const char *loop_name = "sink_env1";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::EpsilonAtomic; /* EpsilonAtomic: sole j-write gated atomic_min(SwallowTime) (order-invariant), no read-back; tiny-N usually below threshold */
    /* Legacy detector label preserved; differs from loop_name. Predates the
     * runner-template loop_name convention. Without this, the runner default
     * would label this Spec's detector "sink_env1" instead of the historical
     * "sink_environment", silently changing diagnostic output. */
    static constexpr const char *ghost_write_detector_name = "sink_environment";

    /* (2) Search policy. */
    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int            neighbor_type_mask = (unsigned int)SINK_NEIGHBOR_BITFLAG;
    /* sink_env1 pair physics: h_j reach is P[j].KernelRadius for gas AND non-gas
     * (pair body uses DMAX(h_i, P[j].KernelRadius) per sink_env1_loop.h:185,274).
     * Restores legacy KernelRadius-for-all-types symmetric reach. */
    static constexpr mode_b_radius_policy_t  radius_policy      =
        MODE_B_RADIUS_GAS_KERNEL | MODE_B_RADIUS_NONGAS_KERNEL;

    /* (3) Writeback policy. */
    static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::AllTypes;
    static constexpr bool mode_a_active_sources_in_sidx_pool = true; /* active sources (incl. non-gas) are in the AllTypes SIDX pool */

    /* (4) Active-particle predicate, passed by the caller to
     *     nlr_build_active_list. */
    static bool is_active(int particle_index) { return sink_isactive(particle_index) != 0; }

    /* (5) Per-pair physics types — file-scope structs above. */
    using CallScalars   = SinkEnv1CallScalars;
    using ActiveData    = SinkEnv1ActiveState;
    using AccumData     = struct sink_env_gpu_out;
    using DeviceContext = SinkEnv1DeviceContext;        /* file-scope struct above */

    struct NeighborData {
        struct particle_data       *neighbor_particle;  /* non-const: pair body atomic-min's SwallowTime */
        const struct gas_cell_data *neighbor_cell;      /* nullptr for non-gas / when no CellP */
    };

    /* (6) Per-active aux passed by caller through neighbor_loop_args::aux.
     *     Recover with nlr_aux<SinkEnv1Spec>(args) inside hooks. */
    struct Aux {
        struct sink_env_gpu_out *per_active_accum;   /* [num_active] — host buffer */
    };

    /* (7) Imported-ghost lifecycle traits.
     *
     * Two independent traits, each gating a hook-method pair below. The
     * runner calls each pair only when (a) the trait is true AND (b) the
     * chosen path imports ghosts (nlr_path_uses_imported_ghosts(plan.path)).
     * Mode B paths run zero lifecycle hooks (no ghost import to attend to).
     *
     * Detector vs writeback are distinct concepts:
     *   - detector  : audit/debug — catches illegal kernel writes to
     *                 imported ghosts (compiles to no-op outside
     *                 GIZMO_GPU_ARENA_DEBUG).
     *   - writeback : physics state propagation — snapshot ghost state
     *                 before the kernel and reverse-communicate any
     *                 j-side writes back to home ranks afterwards.
     *
     * uses_ghost_writeback is set true permanently for sink_env1: the actual
     * work is governed by the contents of the manifest in
     * sinks/sink_env1_loop.cc (GHOST_WRITEBACK_BUNDLE_BEGIN/END). When no
     * physics flag inside that manifest is active, the bundle has zero
     * callbacks and the scaffold short-circuits to a strict no-op.
     * Adding a new ghost-written field is a one-line edit in the manifest;
     * this trait does not change. */
    static constexpr bool uses_ghost_write_detector  = true;
    static constexpr bool uses_ghost_writeback       = true;

    /* Hook method declarations. The runner invokes each only when the
     * corresponding `uses_*` trait above is true AND the chosen path
     * imports ghosts. On the Mode A path the runner preserves this exact
     * order:
     *
     *   gizmo_request_filtered_ghost_import (in runner)
     *   ghost_write_detector_begin
     *   ghost_writeback_begin                (per-flag #ifdef union body)
     *   <run_mode_a kernel>
     *   ghost_writeback_end                  (reverse order)
     *   ghost_write_detector_end
     *   ghost_exchange_cleanup               (in runner; NTask>1 only)
     */
    /* ghost_write_detector_begin/end use the runner default
     * (::ghost_write_detector_begin(ghost_write_detector_name) /
     * ::ghost_write_detector_end()) — no Spec hook needed. */
    static void ghost_writeback_begin      (const struct neighbor_loop_args&,
                                             const struct NeighborLoopPlan&);
    static void ghost_writeback_end        (const struct neighbor_loop_args&,
                                             const struct NeighborLoopPlan&);

    /* (8) Per-active host hooks (pre-arena epoch). */

    /* Per-active search radius from P[i].KernelRadius. Pre-arena, pre-drift. */
    static double      search_radius(const neighbor_loop_args& args,
                                      int active_slot, int i);

    /* Per-call scalars captured into a POD. */
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);

    /* (9) Per-active and per-pair device hooks. KOKKOS_INLINE_FUNCTION =
     *      callable from Mode A device kernel and Mode B/Brute host walker. */

    /* Build ActiveData[slot] from device-visible state (post-neighbor-list
     * build). Byte-exact reconstruction of the legacy GPU lambda's q-packing
     * block. */
    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx,
                                   int active_slot, int i,
                                   double h_search,
                                   const CallScalars& scalars)
    {
        ActiveData active;
        active.pos      = ctx.P[i].Pos;
        active.vel      = ctx.P[i].Vel;
        active.id       = ctx.P[i].ID;
        /* The rank comes from the per-call scalar, not from the global: this body is compiled into
         * the device kernel, where a host global has no value.  populate_call_scalars runs on the
         * rank whose local particles this call loads, so the scalar IS that owner's rank. */
        active.claim_token = gizmo_sink_claim_token(scalars.owner_task, i);
        active.h_search = h_search;
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
        active.ags_h    = (double)ctx.P[i].AGS_KernelRadius;
#else
        active.ags_h    = scalars.sink_radius_grav;
#endif
#if defined(SINK_GRAVCAPTURE_GAS) || (SINK_GRAVACCRETION == 8)
        active.mass     = (double)ctx.P[i].Mass;
#endif
#if defined(SINK_GRAVCAPTURE_FIXEDSINKRADIUS)
        active.sink_radius = (double)ctx.P[i].SinkRadius;
#endif
#if defined(SINK_RETURN_ANGMOM_TO_GAS)
        for(int kv = 0; kv < 3; kv++) active.sink_angmom[kv] = ctx.P[i].Sink_Specific_AngMom[kv];
#endif
        active.origin_local_idx = active_slot;
        active.origin_rank      = -1;          /* device — rank N/A in lambda */
        (void)scalars;
        return active;
    }

    /* Zero accumulator. Byte-zero of the POD AccumData. */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum)
    {
        for(size_t b = 0; b < sizeof(accum); b++) ((char*)&accum)[b] = 0;
    }

    /* Build NeighborData for j. ctx.P/ctx.CellP are UVM under Mode A; in
     * Mode B walker they are P_host/CellP_host — same shape. */
    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx,
                                       int j,
                                       const IdentitySidecar& /*id*/,
                                       const ActiveData& /*active*/)
    {
        NeighborData neighbor;
        neighbor.neighbor_particle = &ctx.P[j];
        neighbor.neighbor_cell     = (ctx.CellP && ctx.P[j].Type == 0) ? &ctx.CellP[j] : nullptr;
        return neighbor;
    }

    /* The physics — forwards to the inline pair body above. */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                             const NeighborData& neighbor,
                             AccumData& accum,
                             NoScatter& /*scatter*/,
                            const CallScalars& cs)
    {
        /* The list holds candidates, not only neighbours: the environment sums over the pairs
           inside the symmetric kernel, r < max(h_i, h_j), with h_j this loop's reach. */
        const struct particle_data &pj = *neighbor.neighbor_particle;
        if(!gx_pair_accept_wrap_and_test(active.pos[0] - (double)pj.Pos[0], active.pos[1] - (double)pj.Pos[1],
                                         active.pos[2] - (double)pj.Pos[2], active.h_search,
                                         nlr_particle_symmetric_radius(pj, radius_policy), NGB_SEARCH_SYMMETRIC)) {return;}
        sink_env1_pair_kernel(active, cs, *neighbor.neighbor_particle,
                              neighbor.neighbor_cell, accum);
    }

    /* (10) Host writebacks (post-dispatch). */

    /* Copy AccumData into args.aux->per_active_accum[active_slot]. The
     * caller's scatter loop reads per_active_accum and applies its own
     * physics-specific reductions into SinkTempInfo. */
    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int active_slot, int i,
                                        const AccumData& accum);

    /* Per-field merge of a peer-rank reply (peer_accum) into a local
     * accumulator (local_accum). Used by run_mode_b_remote at the cross-rank
     * boundary; per-field reduction op MUST match the pair_kernel writes
     * (sum for additive fields, MAX for DF_mmax_particles). Includes all
     * #ifdef-gated optional fields. */
    /* Per-field merge of a peer rank's contribution into a local accumulator.
     * Per-field op MUST match the pair_kernel writes (sum for additive fields,
     * MAX for max-reduced fields, componentwise sum for vec3). Used by
     * run_mode_b_remote at the cross-rank boundary; within a single rank,
     * accumulation is via repeated pair_kernel calls (which already encode
     * the right per-field op).
     *
     * Nothing checks this manifest against pair_kernel at runtime: drift between
     * the two is silent, and shows up only as a cross-rank difference in the
     * affected field.
     *
     * Adding a new accumulator field for this loop = ONE LINE under the
     * appropriate physics flag's #ifdef. Operations available below; extend
     * by adding new ACCUM_* defines if a field needs different op semantics. */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local_accum, const AccumData& peer_accum)
    {
        /* Local op macros — scoped to this function via #undef below. The
         * macros expand to the same statements as the prior hand-written
         * field listing (semantically identical expansion). */
#define ACCUM_ADD(field)       local_accum.field += peer_accum.field;
#define ACCUM_ADD_VEC3(field)  for(int k = 0; k < 3; k++) local_accum.field[k] += peer_accum.field[k];
#define ACCUM_MAX(field)       if(peer_accum.field > local_accum.field) local_accum.field = peer_accum.field;

        ACCUM_ADD(Sink_SurroudingGasInternalEnergy)
        ACCUM_ADD(Mgas_in_Kernel)
        ACCUM_ADD(Mstar_in_Kernel)
        ACCUM_ADD(Malt_in_Kernel)
        ACCUM_ADD_VEC3(Jgas_in_Kernel)
        ACCUM_ADD_VEC3(Jstar_in_Kernel)
        ACCUM_ADD_VEC3(Jalt_in_Kernel)
#ifdef SINK_REPOSITION_ON_POTMIN
        ACCUM_ADD(DF_rms_vel)
        ACCUM_ADD_VEC3(DF_mean_vel)
        ACCUM_MAX(DF_mmax_particles)
#endif
#if defined(SINK_OUTPUT_MOREINFO)
        ACCUM_ADD(Sfr_in_Kernel)
#endif
#if (SINK_GRAVACCRETION >= 5) || defined(SINGLE_STAR_SINK_DYNAMICS) || defined(SINGLE_STAR_TIMESTEPPING)
        ACCUM_ADD_VEC3(Sink_SurroundingGasVel)
#endif
#ifdef JET_DIRECTION_FROM_KERNEL_AND_SINK
        ACCUM_ADD_VEC3(Sink_SurroundingGasCOM)
#endif
#if (SINK_GRAVACCRETION == 8)
        ACCUM_ADD(hubber_mdot_bondi_limiter)
        ACCUM_ADD(hubber_mdot_vr_estimator)
        ACCUM_ADD(hubber_mdot_disk_estimator)
#endif
#if defined(SINK_GRAVCAPTURE_GAS)
        ACCUM_ADD(mass_to_swallow_edd)
#endif
#if defined(SINK_RETURN_ANGMOM_TO_GAS)
        ACCUM_ADD_VEC3(angmom_prepass_sum_for_passback)
#endif
#if defined(SINK_RETURN_BFLUX)
        ACCUM_ADD(kernel_norm_topass_in_swallowloop)
#endif

#undef ACCUM_ADD
#undef ACCUM_ADD_VEC3
#undef ACCUM_MAX
    }

    /* ====================================================================
     * ENGINE APPARATUS — touch only when changing the runner contract
     * ==================================================================== */

    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;
    /* DeviceContext = SinkEnv1DeviceContext is declared in the PHYSICS
     * BLOCK above (alongside the other using aliases) so all inline
     * DeviceContext-typed method declarations below can see it. The
     * derived struct itself is file-scope above SinkEnv1Spec. */

    /* Device-context lifecycle. Both bodies are empty: this Spec stages no
     * per-active device state and has no UVM allocs to free. They are declared
     * because the extended DeviceContext above makes the runner call
     * populate_device_context, and for symmetric contract pairing. */
    static void populate_device_context(const neighbor_loop_args& args, DeviceContext& ctx);
    static void cleanup_device_context (const neighbor_loop_args& args, DeviceContext& ctx);
};

#endif /* SINK_PARTICLES */

#endif /* SINK_ENV1_LOOP_H */
