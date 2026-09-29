/* sinks/sink_env2_loop.h — sink "second environment pass" (Bulge-Disk) module.
 *
 * Defines SinkEnv2Spec for the neighbor-loop runner. SINK_GRAVACCRETION==0
 * branch only; pure i-side aggregator (no j-writes). Replaces
 * sinks/sink_environment_gpu.cc::sink_environment_second_evaluate_gpu.
 *
 * Mechanical-after-pattern port, mirroring sink_env1's
 * template (Pass A+B). The only material additions vs sink_env1 are:
 *   - per-active input vectors (Jgas, Jstar) staged via the
 *     DeviceContext extension (UVM-resident SinkEnv2PerActiveIn array).
 *     load_active reads them on device and stamps into ActiveData.
 *   - SinkEnv2-specific CallScalars is just NlrCommonScalars (common
 *     cosmology only); no sink-specific globals beyond the comoving flag,
 *     which is already in NlrCommonScalars.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef SINK_ENV2_LOOP_H
#define SINK_ENV2_LOOP_H

/* Kokkos_Core.hpp must precede allvars.h (its macros may conflict with stdlib
 * names). */
#include <Kokkos_Core.hpp>
#include "../declarations/allvars.h"

#ifdef SINK_PARTICLES
#if defined(SINK_GRAVACCRETION) && (SINK_GRAVACCRETION == 0)

#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"   /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
#include "../mesh/ghost_exchange_functions.h" /* gx_pair_accept_wrap_and_test: the canonical pair accept */
#include "sinks_gpu_decls.h"               /* struct sink_env_second_gpu_out */
/* NOTE: caller TUs must include "../mesh/kernel.h" BEFORE this header
 * (NEAREST_XYZ etc. — same convention as sink_env1_loop.h). */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* Forward decls. */
int sink_isactive(int i);

/* ============================================================================
 * Per-pair physics types (file scope; PascalCase).
 * ========================================================================== */

/* Per-active host-fill struct, populated by populate_device_context() on
 * the host (copies SinkTempInfo[t].Jgas_in_Kernel + Jstar_in_Kernel into
 * a UVM array). Read on device by Spec::load_active via
 * SinkEnv2DeviceContext::per_active_in. Trivially copyable. */
struct SinkEnv2PerActiveIn {
    Vec3<MyFloat> Jgas;
    Vec3<MyFloat> Jstar;
};

/* Active-particle state passed into the pair body. Top-level pos +
 * h_search are runner-walker-facing fields (sink_env1 precedent —
 * Mode B walker reads these directly from the active struct). The
 * Jgas / Jstar fields come from the per-active UVM array (the
 * DeviceContext extension). Trivially copyable for byte-level MPI
 * transfer in the Mode B remote path. */
struct SinkEnv2ActiveState {
    Vec3<double>      pos;
    double            h_search;
    Vec3<double>      vel;
    Vec3<double>      Jgas;
    Vec3<double>      Jstar;
    int               origin_local_idx;
    int               origin_rank;
};

/* DeviceContext extension. Carries the UVM-resident
 * per-active inputs (Jgas / Jstar) so load_active can read per-slot.
 * Allocated by populate_device_context, freed by cleanup_device_context
 * (RAII guard at runner exit on every dispatch path). Trivially
 * copyable: runner captures by value into Kokkos device lambdas. */
struct SinkEnv2DeviceContext : NeighborLoopDeviceContextBase {
    const SinkEnv2PerActiveIn *per_active_in;   /* UVM, [num_active] */
};

/* ============================================================================
 * Inline pair body — single source of truth for sink_env2 per-pair physics.
 * Pure i-side: writes only to `accum`, never to neighbor_particle. No
 * atomics, no Kokkos:: refs (so no Kokkos_Core.hpp include required from
 * this header. sink_env2_loop.cc IS in GPU_OBJS because populate_device_context
 * uses Kokkos::kokkos_malloc / kokkos_free for the per-active Jgas/Jstar
 * UVM staging — the CUDA backend requires nvcc_wrapper compilation for
 * the .cc TU even though the inline pair body itself is Kokkos-free).
 * ========================================================================== */

KOKKOS_INLINE_FUNCTION
static void sink_env2_pair_kernel(const SinkEnv2ActiveState& active,
                                  const NlrCommonScalars& scalars,
                                   const struct particle_data& neighbor_particle,
                                   struct sink_env_second_gpu_out& accum)
{
    if(neighbor_particle.Mass <= 0) return;
    if(neighbor_particle.KernelRadius <= 0) return;
    if(neighbor_particle.Type == 5) return;

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

    Vec3<double> J_tmp = cross(dP, dv);

    if(neighbor_particle.Type == 0) {
        if(dot(J_tmp, active.Jgas) < 0) {
            accum.MgasBulge_in_Kernel += (MyFloat)(2 * neighbor_particle.Mass);
        }
    }
    if(is_galsf_stellar_candidate_type(neighbor_particle.Type, scalars.comoving_integration_on)) {
        if(dot(J_tmp, active.Jstar) < 0) {
            accum.MstarBulge_in_Kernel += (MyFloat)(2 * neighbor_particle.Mass);
        }
    }
}

/* ============================================================================
 * SinkEnv2Spec — the NeighborLoopSpec contract for sink_environment_second.
 * NeighborData lifetime: same as sink_env1.
 * ========================================================================== */

struct SinkEnv2Spec {
    /* ====================================================================
     * PHYSICS BLOCK
     * ==================================================================== */

    static constexpr const char *loop_name = "sink_env2";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only (const j-pointer, MgasBulge/MstarBulge sums); no j-write -> threaded eval bit-identical */
    /* Legacy detector label preserved; differs from loop_name. Predates the
     * runner-template loop_name convention. Without this, the runner default
     * would label this Spec's detector "sink_env2" instead of the historical
     * "sink_environment_second". */
    static constexpr const char *ghost_write_detector_name = "sink_environment_second";

    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int            neighbor_type_mask = (unsigned int)SINK_NEIGHBOR_BITFLAG;
    /* sink_env2 pair physics: pair body requires P[j].KernelRadius > 0 and
     * accumulates Mstar_in_Kernel without an internal r > h_j gate, so the
     * symmetric reach must admit j with r < P[j].KernelRadius for non-gas too.
     * See sink_env2_loop.h:100. */
    static constexpr mode_b_radius_policy_t  radius_policy      =
        MODE_B_RADIUS_GAS_KERNEL | MODE_B_RADIUS_NONGAS_KERNEL;

    /* WritePattern + uses_ghost_writeback are orthogonal axes. sink_env2
     * is pure i-side (no j-writes); writeback bundle is empty. The
     * detector still fires to AUDIT that the kernel honors the pure-
     * aggregator contract. */
    static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::AllTypes;
    static constexpr bool mode_a_active_sources_in_sidx_pool = true; /* active sources (incl. non-gas) are in the AllTypes SIDX pool */

    static bool is_active(int particle_index) { return sink_isactive(particle_index) != 0; }

    using CallScalars   = NlrCommonScalars;          /* common only — no extra globals */
    using ActiveData    = SinkEnv2ActiveState;
    using AccumData     = struct sink_env_second_gpu_out;
    using DeviceContext = SinkEnv2DeviceContext;

    /* NeighborData uses const pointers — sink_env2 is read-only on j-side
     * (mirrors sink_env1; differs from sink_feed which needs non-const for
     * atomic writes). */
    struct NeighborData {
        const struct particle_data *neighbor_particle;
        const struct gas_cell_data *neighbor_cell;     /* nullptr unless gas */
    };

    /* Per-active aux. host_inputs is a host mirror of the SinkEnv2PerActiveIn
     * array; populate_device_context copies it into ctx.per_active_in. */
    struct Aux {
        struct sink_env_second_gpu_out *per_active_accum;  /* [num_active] host */
        const SinkEnv2PerActiveIn      *host_inputs;       /* [num_active] host */
    };

    /* sink_env2 has no j-writes — no ghost-writeback bundle entries needed.
     * Detector stays on for auditing (mirrors sink_env1 + sink_environment.cc
     * legacy ghost_write_detector_begin/end pair). */
    static constexpr bool uses_ghost_write_detector = true;
    static constexpr bool uses_ghost_writeback      = false;

    /* ghost_write_detector_begin/end use the runner default
     * (::ghost_write_detector_begin(ghost_write_detector_name) /
     * ::ghost_write_detector_end()) — no Spec hook needed. */

    /* Per-active host hooks. */
    static double      search_radius(const neighbor_loop_args& args,
                                      int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);

    /* Device-context lifecycle. populate copies host_inputs
     * into a UVM array; cleanup frees. cleanup runs unconditionally via
     * NlrDeviceContextCleanupGuard at runner exit on every path. */
    static void populate_device_context(const neighbor_loop_args& args, DeviceContext& ctx);
    static void cleanup_device_context (const neighbor_loop_args& args, DeviceContext& ctx);

    /* Build ActiveData[slot] on device. Reads pos/vel from ctx.P[i] and
     * Jgas/Jstar from the host-staged UVM array attached to dctx in
     * populate_device_context. */
    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& dctx,
                                   int active_slot, int i,
                                   double h_search,
                                   const CallScalars& scalars)
    {
        ActiveData active;
        active.pos[0] = (double)dctx.P[i].Pos[0];
        active.pos[1] = (double)dctx.P[i].Pos[1];
        active.pos[2] = (double)dctx.P[i].Pos[2];
        active.vel[0] = (double)dctx.P[i].Vel[0];
        active.vel[1] = (double)dctx.P[i].Vel[1];
        active.vel[2] = (double)dctx.P[i].Vel[2];
        active.h_search = h_search;
        const SinkEnv2PerActiveIn zero_in{};
        const SinkEnv2PerActiveIn& src = dctx.per_active_in ? dctx.per_active_in[active_slot]
                                                            : zero_in;
        active.Jgas[0]  = (double)src.Jgas[0];
        active.Jgas[1]  = (double)src.Jgas[1];
        active.Jgas[2]  = (double)src.Jgas[2];
        active.Jstar[0] = (double)src.Jstar[0];
        active.Jstar[1] = (double)src.Jstar[1];
        active.Jstar[2] = (double)src.Jstar[2];
        (void)scalars;
        active.origin_local_idx = active_slot;
        active.origin_rank      = -1;
        return active;
    }

    /* Zero accumulator. POD AccumData; byte-zero suffices. */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum)
    {
        for(size_t b = 0; b < sizeof(accum); b++) ((char*)&accum)[b] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& dctx,
                                       int j,
                                       const IdentitySidecar& /*id*/,
                                       const ActiveData& /*active*/)
    {
        NeighborData neighbor;
        neighbor.neighbor_particle = &dctx.P[j];
        neighbor.neighbor_cell     = (dctx.CellP && dctx.P[j].Type == 0) ? &dctx.CellP[j] : nullptr;
        return neighbor;
    }

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
        sink_env2_pair_kernel(active, cs, *neighbor.neighbor_particle, accum);
    }

    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int active_slot, int i,
                                        const AccumData& accum);

    /* Per-field merge (Mode B remote). sink_env2 has only two additive fields. */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local_accum, const AccumData& peer_accum)
    {
#define ACCUM_ADD(field)  local_accum.field += peer_accum.field;

        ACCUM_ADD(MgasBulge_in_Kernel)
        ACCUM_ADD(MstarBulge_in_Kernel)

#undef ACCUM_ADD
    }

    /* ====================================================================
     * ENGINE APPARATUS
     * ==================================================================== */
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

};

#endif /* SINK_GRAVACCRETION == 0 */
#endif /* SINK_PARTICLES */

#endif /* SINK_ENV2_LOOP_H */
