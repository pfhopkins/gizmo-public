/* gravity/ags_force_loop.cc — host hooks + ghost-writeback bundle +
 * toplevel AGSForce_calc() for AgsForceSpec.
 *
 * KOKKOS_INLINE_FUNCTION hooks (load_active, load_neighbor, pair_kernel,
 * zero_accum) and the inline pair body (ags_force_pair_kernel_body) live
 * in ags_force_loop.h so they inline from device kernels (Mode A) and
 * host walkers (Mode B). This TU carries the host-only
 * hooks, the generic ghost-writeback bundle (PARTICLE_ADD_VEC3 on
 * Vel/dp, PARTICLE_ADD on NInteractions, PARTICLE_MAX on wakeup), the
 * host-side wakeup pre-zero + arena invalidate that gives the generic
 * MAX op the same event semantics as the retired snapshot-based
 * ghost_writeback_agsforce, and the toplevel caller that builds the
 * MPI-collective global-union subgroup list (mirrors ags_density_loop.cc:828-865).
 *
 * Replaces gravity/ags_force_gpu.cc and the agsforce block in
 * mesh/ghost_writeback.cc (the latter retired in the cleanup commit).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"     /* MUST precede allvars.h */
#include "../declarations/allvars.h"
#include "../declarations/gpu_numeric_macros.h"
#include "../declarations/gpu_error_check.h"
#include "../core/proto.h"
#include "../core/timestep_functions.h"
#include "../mesh/kernel.h"                     /* MUST precede ags_force_loop.h */
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_writeback_ops.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "../system/gpu_particles_arena.h"
#include "../sidm/sidm_gpu_decls.h"

/* Flux helper templates are now included by ags_force_loop.h itself, so the
 * Spec header is self-sufficient (required because neighbor_loop_runner.cc
 * explicitly instantiates run_neighbor_loop_iterative<AgsForceSpec> and the
 * runner template body transitively references the flux helpers). */
#include "ags_force_loop.h"

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE


/* ============================================================================
 * HOST HOOKS
 * ========================================================================== */

double AgsForceSpec::search_radius(const neighbor_loop_args& /*args*/,
                                    int /*active_slot*/, int i)
{
    /* The pair terms all vanish at r >= max(h_i,h_j), so the search radius is
     * the particle's own reach: the symmetric walk supplies the h_j > h_i side. */
    return (double)P[i].AGS_KernelRadius;
}

AgsForceSpec::CallScalars
AgsForceSpec::populate_call_scalars(const neighbor_loop_args& /*args*/)
{
    CallScalars scalars;
    scalars.common = nlr_common_scalars_from_all();
    for(int k = 0; k < TIMEBINS; k++) scalars.TimeBinActive[k] = TimeBinActive[k];
    scalars.rng_salt = AGS_FORCE_RNG_SALT;
    return scalars;
}

void AgsForceSpec::populate_device_context(const neighbor_loop_args& args,
                                            DeviceContext& ctx)
{

    /* Sticky single-int wakeup flag, lives across all subgroups of one
     * toplevel call. Lifecycle matches ags_density's need_wakeup_uvm. */
    ctx.need_wakeup_uvm = NULL;
#if defined(DM_SIDM)
#endif
    ctx.wakeup_dirty_base = WakeupDirty;   /* global UVM sidecar base; kernel marks WakeupDirty[j] on wakeup */
    ctx.need_wakeup_uvm = (int *) gizmo_gpu_alloc_shared(sizeof(int), NULL);
    /* Zeroed the moment it exists: cleanup_device_context reads it on every exit
     * path, including the failure return below. */
    if(ctx.need_wakeup_uvm) { *ctx.need_wakeup_uvm = 0; }
    if(!ctx.need_wakeup_uvm || ctx.populate_failed) {
        ctx.populate_failed = 1;
        gizmo_request_controlled_stop(7723,
            "adaptive-softening forces: could not stage the wakeup flag; "
            "the forces are not computed",
            __FILE__, __LINE__, __FUNCTION__);
        return;
    }

    /* CBE gradients are persistent on P[i].Gradients_CBE_basis_moments and
     * ghost-transported naturally by standard P[] import — no scratch
     * pointer to plumb through the DeviceContext. */
    (void)args;
}

void AgsForceSpec::cleanup_device_context(const neighbor_loop_args& /*args*/,
                                           DeviceContext& ctx)
{
    if(ctx.need_wakeup_uvm) {
        if(*ctx.need_wakeup_uvm) NeedToWakeupParticles_local = 1;
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(ctx.need_wakeup_uvm);
        ctx.need_wakeup_uvm = nullptr;
    }
#if defined(DM_SIDM)
#endif
}

/* apply_active_writeback — verbatim translation of the i-side scatter at
 * gravity/ags_rkern.cc:278-336 (legacy AGSForce_calc body). Each block is
 * exactly the legacy block reformatted into Spec-hook form. */
void AgsForceSpec::apply_active_writeback(const neighbor_loop_args& /*args*/,
                                           int /*active_slot*/, int i,
                                           const AccumData& accum)
{
#if defined(DM_SIDM)
    for(int k = 0; k < 3; k++) {
        P[i].Vel[k] += accum.sidm_kick[k];
        P[i].dp[k]  += accum.sidm_kick[k] * P[i].Mass;
    }
    /* MIN-merge against the AGSForce_calc preamble seed
     * (P[i].dtime_sidm = 10*get_particle_timestep_in_physical(i)). */
    if(accum.dtime_sidm < P[i].dtime_sidm) P[i].dtime_sidm = accum.dtime_sidm;
    P[i].NInteractions += accum.si_count;
#endif
#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & 7)
    if(P[i].Mass > 0) {
        if(accum.Grain_DeltaCoagMass > 0) {
            double M_old = (double)P[i].Mass;
            double M_new = M_old + accum.Grain_DeltaCoagMass;
            for(int s = 0; s < GRAIN_NUM_SPECIES; s++) {
                double M_species_s = M_old * (double)P[i].Composition[s]
                                     + accum.Grain_DeltaCoag_CompositionMass[s];
                if(M_species_s < 0) M_species_s = 0;
                P[i].Composition[s] = (MyFloat)(M_species_s / M_new);
            }
            P[i].Grain_Size = (MyFloat)((double)P[i].Grain_Size * pow(M_new / M_old, 1.0 / 3.0));
            P[i].Mass       = (MyDouble)M_new;
        }
        if(accum.Grain_DeltaErosionFrac != 1.0 && accum.Grain_DeltaErosionFrac > 0.0) {
            P[i].Grain_Size = (MyFloat)((double)P[i].Grain_Size * accum.Grain_DeltaErosionFrac);
        }
    }
#endif
#ifdef DM_FUZZY
    for(int k = 0; k < 3; k++) P[i].GravAccel[k] += accum.acc[k];
    P[i].AGS_Dt_Numerical_QuantumPotential += accum.AGS_Dt_Numerical_QuantumPotential;
#if (DM_FUZZY > 0)
    P[i].AGS_Dt_Psi_Re   += accum.AGS_Dt_Psi_Re;
    P[i].AGS_Dt_Psi_Im   += accum.AGS_Dt_Psi_Im;
    P[i].AGS_Dt_Psi_Mass += accum.AGS_Dt_Psi_Mass;
#endif
#endif
#if defined(CBE_INTEGRATOR)
    if(accum.AGS_vsig > P[i].AGS_vsig) P[i].AGS_vsig = accum.AGS_vsig;
    for(int k1 = 0; k1 < CBE_INTEGRATOR_NBASIS; k1++) {
        for(int k2 = 0; k2 < CBE_INTEGRATOR_NMOMENTS; k2++) {
            P[i].CBE_basis_moments_dt[k1][k2]  += accum.CBE_basis_moments_dt[k1][k2];
            P[i].CBE_basis_out_rate_dt[k1][k2] += accum.CBE_basis_out_rate_dt[k1][k2];
        }
    }
#if defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(CBE_INTEGRATOR_OUTPUT_MOREINFO)
    cbe_step_diagnostics_observe_face(accum.cbe_face_residual_max,
                                       accum.cbe_face_residual_sum,
                                       accum.cbe_bracket_fail_count);
    cbe_step_diagnostics_observe_recon(accum.cbe_recon_rho_clamp_count,
                                        accum.cbe_recon_S_clamp_count);
    cbe_step_diagnostics_observe_pairing_free_slot(accum.cbe_pairing_free_slot_count);
#if defined(CBE_INTEGRATOR_WITHGRADIENTS)
    cbe_step_diagnostics_observe_grad_nonfinite(accum.cbe_grad_nonfinite_count);
#endif
#endif
#endif
    (void)accum; (void)i;
}


/* ============================================================================
 * GHOST-WRITEBACK BUNDLE (generic ops; replaces the 164-line bespoke
 * ghost_writeback_{zero_,}agsforce in mesh/ghost_writeback.cc).
 *
 * Coverage:
 *   P[j].Vel             — DM_SIDM scatter kicks       (PARTICLE_ADD_VEC3)
 *   P[j].dp              — DM_SIDM scatter kicks       (PARTICLE_ADD_VEC3)
 *   P[j].NInteractions   — DM_SIDM scatter counter     (PARTICLE_ADD)
 *   P[j].wakeup          — DM_SIDM/CBE wakeup writes   (PARTICLE_MAX)
 *
 * The wakeup-event semantics needs a host-side pre-zero of ghost
 * P[j].wakeup BEFORE the bundle snapshot — see ghost_writeback_begin
 * below. Without that, an imported ghost wakeup ≥ what the pair body
 * writes would let the snapshot-diff op (post>snap) silently drop the
 * event, where legacy ghost_writeback_zero_agsforce + the post-kernel
 * `wakeup != 0` predicate would have sent it.
 * ========================================================================== */

GHOST_WRITEBACK_BUNDLE_BEGIN(ags_force)
#if defined(DM_SIDM)
    GHOST_WRITEBACK_PARTICLE_ADD_VEC3(Vel)
    GHOST_WRITEBACK_PARTICLE_ADD_VEC3(dp)
    GHOST_WRITEBACK_PARTICLE_ADD(NInteractions)
#endif
#if defined(DM_SIDM) || defined(CBE_INTEGRATOR)
    GHOST_WRITEBACK_PARTICLE_MAX(wakeup)
#endif
GHOST_WRITEBACK_BUNDLE_END(ags_force)


/* ghost_write_detector_begin/end: runner default (loop_name = "ags_force"). */

void AgsForceSpec::ghost_writeback_begin(const neighbor_loop_args& args,
                                           const NeighborLoopPlan& /*plan*/)
{
    (void)args;   /* may be unused depending on physics flags below */
#if defined(DM_SIDM) || defined(CBE_INTEGRATOR)
    /* Pre-zero ghost P[j].wakeup BEFORE the bundle snapshot. Restores the
     * event semantics of the retired ghost_writeback_zero_agsforce: any
     * nonzero wakeup written by the pair body to an imported ghost ships
     * home regardless of whether the home rank's wakeup is already higher
     * (home then MAX-merges, so a no-op home update is fine; the missing
     * case was a kernel write that didn't trip post>snap because the
     * ghost was imported with an already-larger wakeup). */
    int num_ghosts = ghost_get_num_ghosts();
    int num_local  = ghost_get_num_local();
    for(int g = 0; g < num_ghosts; g++) {
        P[num_local + g].wakeup = 0;
    }
    /* Host mutated ghost particles; force the next arena acquire to re-
     * stage P. Mirrors the legacy invalidate at the end of
     * the retired ghost_writeback_zero_agsforce path (now deleted). */
    gpu_particles_arena_invalidate();
#endif
    /* CBE gradients are persistent on P[i].Gradients_CBE_basis_moments,
     * so the standard P[] ghost import (driven by the runner's ghost
     * lifecycle, not this hook) carries them onto ghost slots
     * automatically. No custom CBE import call lives here under the v5
     * corrective architecture. */
    ghost_writeback_begin_bundle(ags_force_ghost_writeback_bundle_ptr());
}

void AgsForceSpec::ghost_writeback_end(const neighbor_loop_args& /*args*/,
                                         const NeighborLoopPlan& /*plan*/)
{
    ghost_writeback_end_bundle(ags_force_ghost_writeback_bundle_ptr());
}


/* ============================================================================
 * TOPLEVEL CALLER — AGSForce_calc()
 *
 * Replaces the body that used to live in gravity/ags_rkern.cc:204-350.
 * Builds the MPI-collective global-union subgroup partition by Allreduce-BOR
 * of local bitmask presence, mirroring gravity/ags_density_loop.cc:828-865.
 * All ranks enter run_neighbor_loop_iterative<AgsForceSpec> with the SAME
 * subgroup list ordering; ranks with no actives for a bm contribute an
 * empty subgroup entry (num_active_local=0, active_indices=nullptr) so
 * the collective is symmetric.
 *
 * Iterative-shaped wrapper (max_iters=1, after_iter always Converged) is
 * a fit-the-runner pattern; physics is single-pass per subgroup. Same
 * idiom mechfb / radfb_rp / dm_dispersion use.
 * ========================================================================== */

void AGSForce_calc(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    double t00_truestart = my_second(); double child0_span = CPU_ChildCharged;
    PRINT_STATUS(" ..entering AGS-Force calculation [as hydro loop for non-gas elements]\n");

    double timecomp = 0, timecomm = 0, timewait = 0;
    double t0 = my_second(); const double child0_timecomp = CPU_ChildCharged;

    /* SIDM preamble (legacy gravity/ags_rkern.cc:209-211): seed
     * P[i].dtime_sidm with 10*dt so the in-loop MIN reduce pulls it down
     * only on real scatter events. zero_accum uses MAX_REAL_NUMBER as the
     * MIN-merge sentinel; apply_active_writeback MIN-applies accum into
     * the seeded P[i].dtime_sidm. */
#if defined(DM_SIDM)
    for(int i : ActiveParticleList) {
        if(AGSForce_isactive(i)) {
            P[i].dtime_sidm = 10. * get_particle_timestep_in_physical(i, P);
        }
    }
#endif

    /* Build per-bm subgroup partition. Same partition rule as legacy
     * gravity/ags_rkern.cc:224-234. */
    std::map<int, std::vector<int>> bm_groups_host;
    uint64_t local_bm_presence = 0;
    for(int ii : ActiveParticleList) {
        if(AGSForce_isactive(ii)) {
            int bm = ags_gravity_kernel_shared_BITFLAG(P[ii].Type);
            if(bm > 0 && bm < 64) {
                bm_groups_host[bm].push_back(ii);
                local_bm_presence |= (1ULL << bm);
            }
        }
    }
    uint64_t global_bm_presence = local_bm_presence;
    if(NTask > 1) {
        MPI_Allreduce(&local_bm_presence, &global_bm_presence, 1,
                      MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
    }

    /* Zero per-call i-side ASSIGN-mode accumulators on local actives
     * (legacy gravity/ags_rkern.cc:236-256). These match OUTPUT fields the
     * pair body writes via ASSIGN, not via ACCUM_ADD — they MUST be zeroed
     * before the kernel so the post-kernel scatter is well-defined. */
    for(auto& kv : bm_groups_host) {
        for(int ii : kv.second) {
#ifdef DM_FUZZY
            P[ii].AGS_Dt_Numerical_QuantumPotential = 0;
#if (DM_FUZZY > 0)
            P[ii].AGS_Dt_Psi_Re = P[ii].AGS_Dt_Psi_Im = P[ii].AGS_Dt_Psi_Mass = 0;
#endif
#endif
#if defined(CBE_INTEGRATOR)
            P[ii].AGS_vsig = 0;
            for(int k1 = 0; k1 < CBE_INTEGRATOR_NBASIS; k1++) {
                for(int k2 = 0; k2 < CBE_INTEGRATOR_NMOMENTS; k2++) {
                    P[ii].CBE_basis_moments_dt[k1][k2]  = 0;
                    P[ii].CBE_basis_out_rate_dt[k1][k2] = 0;
                }
            }
#endif
            (void)ii;
        }
    }

    /* Build NlrSubgroup array in global-union ordering. ALL ranks enter
     * with the same subgroup ordering; ranks with no actives for a bm
     * contribute an empty entry (collective-symmetry filler). */
    std::vector<NlrSubgroup> subgroups;
    std::vector<std::vector<int>> subgroup_actives_storage;
    int total_active_local = 0;
    for(int bm = 1; bm < 64; bm++) {
        if(!(global_bm_presence & (1ULL << bm))) continue;
        std::vector<int> actives_for_bm;
        auto it = bm_groups_host.find(bm);
        if(it != bm_groups_host.end()) actives_for_bm = std::move(it->second);
        subgroup_actives_storage.push_back(std::move(actives_for_bm));
        NlrSubgroup sg{};
        sg.j_type_bitmask   = (unsigned int)bm;
        sg.num_active_local = (int)subgroup_actives_storage.back().size();
        sg.active_indices   = (sg.num_active_local > 0)
                              ? subgroup_actives_storage.back().data()
                              : nullptr;
        subgroups.push_back(sg);
        total_active_local += sg.num_active_local;
    }

    if(subgroups.empty()) {
        /* No globally-active AGS particles anywhere. Skip the runner call. */
        double t1 = my_second(); cpu_chain_sync(t1);
        CPU_Step[CPU_AGSDENSMISC] += cpu_minus_children(timediff(t00_truestart, t1), child0_span);
        return;
    }

    /* Concatenated active list across subgroups (runner's union view). */
    std::vector<int> active_list_concat;
    active_list_concat.reserve(total_active_local);
    for(const auto& sg_actives : subgroup_actives_storage) {
        for(int i : sg_actives) active_list_concat.push_back(i);
    }

    /* The runner owns the ghost lifecycle (see ags_density). */
    double ags_ghost_safety = gizmo_ghost_safety_factor();

    /* CBE gradients (when CBE_INTEGRATOR_WITHGRADIENTS is on) are refreshed
     * pre-force by CBEGrad_gradient_calc() (called from core/accel.cc
     * before this toplevel), with persistent storage on
     * P[i].Gradients_CBE_basis_moments. AGSForce_calc itself owns no
     * gradient scratch; the flux body reads the persistent field via
     * load_active's by-value snapshot (i side) and direct P[j] access
     * (j side, naturally ghost-imported). */
    AgsForceSpec::Aux aux{};

    neighbor_loop_args_iterative args{};
    static_cast<neighbor_loop_args&>(args) = nlr_default_args();
    args.P                   = P;
    args.CellP               = (All.TotN_gas > 0) ? CellP : nullptr;
    args.num_total           = NumPart;
    args.active_list         = (total_active_local > 0)
                               ? active_list_concat.data() : nullptr;
    args.num_active          = total_active_local;
    args.aux                 = &aux;
    args.num_subgroups       = (int)subgroups.size();
    args.subgroups           = subgroups.data();
    args.ghost_safety_factor = ags_ghost_safety;

    run_neighbor_loop_iterative<AgsForceSpec>(args);

    if(NTask > 1) ghost_exchange_cleanup();
    timecomp += cpu_minus_children(timediff(t0, my_second()), child0_timecomp);

    /* CBE post-tree-walk per-active finalization (legacy gravity/ags_rkern.cc:343-345).
       Active set is exactly ActiveParticleList — no extra gating; matches legacy. */
#ifdef CBE_INTEGRATOR
    cbe_postgravity_evaluate_gpu(P, ActiveParticleList.data(), (int)ActiveParticleList.size());
#endif

    double t1 = my_second(); cpu_chain_sync(t1);
    double timeall = cpu_minus_children(timediff(t00_truestart, t1), child0_span);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp;
    CPU_Step[CPU_AGSDENSWAIT]    += timewait;
    CPU_Step[CPU_AGSDENSCOMM]    += timecomm;
    CPU_Step[CPU_AGSDENSMISC]    += timeall - (timecomp + timewait + timecomm);
}


/* Explicit instantiation lives in mesh/neighbor_loop_runner.cc alongside
 * the template definition (per project convention — see the explicit-
 * instantiation block near the end of that file). */

#endif /* AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE */
