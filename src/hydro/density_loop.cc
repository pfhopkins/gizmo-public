/* hydro/density_loop.cc — DensitySpec host hooks and runner-driven density().
 *
 * This file owns the modern hydro density loop: active selection,
 * iterative radius convergence, final P/CellP scatter, and the downstream
 * hydro ghost handoff. The device pair body lives in density_loop.h because
 * the runner instantiates it from GPU translation units.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>            /* MPI_Datatype used in declarations/typedefs.h via gpu_all_mirror.h's transitive include */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"  /* MUST precede allvars.h: installs device-pass `#define All AllDeviceMirror` redirect before cell_data.h is parsed */
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"               /* MUST precede density_loop.h */
#include "../mesh/ghost_symlist_lifecycle.h" /* gizmo_ghost_safety_factor + hydro_density_redo (post-finalize downstream refresh) */
#include "../system/gpu_particles_arena.h"   /* gpu_particles_arena_invalidate — pre-runner prepass mirrors legacy density.cc:183 */
#include "density_loop.h"
#include "../mesh/gpu_neighbor_list.h"      /* gizmo_mark_kernel_radius_dirty_indices */

/* ====================================================================
 * density_isactive — hydro density active-particle predicate.
 *
 * Host-only. Bare All.* is avoided because this file is compiled as a GPU
 * translation unit in some builds; nlr_host_all_ptr() is the canonical
 * accessor for the real host global.
 * ==================================================================== */
int density_isactive(int n)
{
    const struct global_data_all_processes *host_all = nlr_host_all_ptr();

    /* Marker for particles done iterating — same negation pattern legacy
     * density() uses. The runner path tracks convergence in IterScratch,
     * but keeping this guard preserves compatibility with callers that
     * check the shared predicate. */
    if(P[n].TimeBin < 0) {return 0;}
    if(P[n].Type == 0) {if(CellP[n].recent_refinement_flag == 1) return 1;}

#if defined(GRAIN_FLUID)
    if((1 << P[n].Type) & (GRAIN_PTYPES)) {return 1;}
#endif

#if defined(SINK_INTERACT_ON_GAS_TIMESTEP)
    if(P[n].Type == 5) {
        if(!P[n].do_gas_search_this_timestep && host_all->Ti_Current > 0) return 0;
    }
#endif

#if defined(RT_SOURCE_INJECTION)
    if((1 << P[n].Type) & (RT_SOURCES))
    {
#if defined(GALSF)
       if(is_galsf_stellar_candidate_type(P[n].Type, host_all->ComovingIntegrationOn) && (P[n].Mass>0))
        {
            double star_age = evaluate_stellar_age_Gyr(n);
            if((star_age < 0.1)&&(star_age > 0)&&(!isnan(star_age))) return 1;
        }
#else
        if(Flag_FullStep) {return 1;}
#endif
    }
#endif

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    if(is_galsf_stellar_candidate_type(P[n].Type, host_all->ComovingIntegrationOn) && (P[n].Mass>0))
    {
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
        if(P[n].SNe_ThisTimeStep>0) return 1;
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION)
        if(P[n].MassReturn_ThisTimeStep>0) return 1;
#ifdef GALSF_FB_FIRE_RPROCESS
        if(P[n].RProcessEvent_ThisTimeStep>0) return 1;
#endif
#if defined(GALSF_FB_FIRE_AGE_TRACERS)
        if(P[n].AgeDeposition_ThisTimeStep>0) return 1;
#endif
#endif
#endif

#if defined(GALSF)
        if(P[n].DensityAroundParticle <= 0) return 1;
        if(host_all->ComovingIntegrationOn == 0)
        {
            double star_age = evaluate_stellar_age_Gyr(n);
            if(star_age < 0.035) return 1;
        }
#endif
#if (defined(GRAIN_FLUID) || defined(RADTRANSFER)) && (!defined(GALSF) && !(defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)))
        return 1;
#endif
    }
#endif

#ifdef SINK_PARTICLES
    if(P[n].Type == 5) return 1;
#endif

    if(P[n].Type == 0 && P[n].Mass > 0) return 1;
    return 0;
}

/* ====================================================================
 * search_radius — per-active radius hook.
 *
 * Returns P[i].KernelRadius (legacy density's iteration variable). The
 * runner stages this into radii_uvm at iter 0; after_iter mutates
 * per-iter via IterResult::new_h_search.
 * ==================================================================== */
double DensitySpec::search_radius(const neighbor_loop_args& args,
                                  int /*active_slot*/, int i) {
    return (double)args.P[i].KernelRadius;
}

/* ====================================================================
 * populate_call_scalars — per-call CallScalars snapshot.
 *
 * NLR host-All accessor rule: read host_all via nlr_host_all_ptr()
 * — NEVER bare All.*.
 * NlrCommonScalars provides cf_atime/cf_a2inv/comoving_integration_on
 * (etc.) via the shared helper.
 * ==================================================================== */
DensitySpec::CallScalars DensitySpec::populate_call_scalars(
        const neighbor_loop_args& /*args*/) {
    const struct global_data_all_processes *host_all = nlr_host_all_ptr();
    CallScalars scalars;
    scalars.common = nlr_common_scalars_from_all();

    scalars.ti_current               = host_all->Ti_Current;
    scalars.time_now                 = host_all->Time;
    scalars.time_begin               = host_all->TimeBegin;
    scalars.time_step                = host_all->TimeStep;

    scalars.des_num_ngb              = host_all->DesNumNgb;
    scalars.max_num_ngb_dev          = host_all->MaxNumNgbDeviation;
    scalars.min_kernel_radius        = host_all->MinKernelRadius;
    scalars.max_kernel_radius        = host_all->MaxKernelRadius;

#ifdef SINK_PARTICLES
    scalars.sink_ngb_factor          = host_all->SinkNgbFactor;
    scalars.sink_max_accretion_radius = host_all->SinkMaxAccretionRadius;
#endif

    return scalars;
}

/* ====================================================================
 * density_build_active_list — host-side active particle list.
 *
 * Walks ActiveParticleList, filters with density_isactive, and returns
 * a contiguous std::vector<int> for the runner. ======================== */
std::vector<int> density_build_active_list(void)
{
    std::vector<int> active_list_concat;
    active_list_concat.reserve(NumPart);
    for (int i : ActiveParticleList) {
        if (density_isactive(i)) active_list_concat.push_back(i);
    }
    return active_list_concat;
}

/* ====================================================================
 * density_build_subgroups — single gas-only subgroup.
 *
 * Density partitions actives into ONE subgroup with j_type_bitmask=1
 * (gas). Per-i type branching (for non-gas active particles via
 * DO_DENSITY_AROUND_NONGAS_PARTICLES etc.) lives inside the pair-kernel
 * body; the runner subgroup is about the neighbor (j) mask, not the
 * active (i) type.
 *
 * Backing storage for active_indices is the caller-owned vector;
 * subgroup.active_indices points into it. Caller
 * must keep both alive for the duration of the runner call. ============ */
std::vector<NlrSubgroup> density_build_subgroups(int *active_indices, int num_active)
{
    std::vector<NlrSubgroup> subgroups;
    NlrSubgroup sg{};
    sg.j_type_bitmask    = (1u << 0);   /* gas only */
    sg.active_indices    = active_indices;
    sg.num_active_local  = num_active;
    subgroups.push_back(sg);
    return subgroups;
}

/* ====================================================================
 * Hook bodies.
 * ==================================================================== */

/* ====================================================================
 * apply_active_writeback — required by the Spec contract, but for
 * DensitySpec the runner uses the iterative variant below
 * (apply_active_writeback_iterative) which carries the converged radius
 * + IterScratch. This hook stays declared/defined so existing runner
 * dispatch paths that fall through to the non-iterative hook (e.g.,
 * future audit/diagnostic paths) link cleanly. The body copies just
 * the accum; density_finalize_post_runner needs the radius too, which
 * only the iterative variant supplies. Keeping the non-iterative form
 * accum-only is consistent with the runner's documented contract.
 * ==================================================================== */
void DensitySpec::apply_active_writeback(const neighbor_loop_args& args,
                                          int active_slot, int /*i*/,
                                          const AccumData& accum) {
    auto* aux = static_cast<Aux*>(args.aux);
    aux->per_active_final_accum[active_slot] = accum;
}

/* ====================================================================
 * apply_active_writeback_iterative.
 *
 * Production-only writeback channel. The runner's iterative post-loop
 * writeback (mesh/neighbor_loop_runner.cc:~4121) calls this INSTEAD of
 * apply_active_writeback when the Spec opts in (SFINAE-detected via
 * nlr_spec_has_apply_active_writeback_iterative_v). That loop runs
 * exactly once per active, which is what makes writing args.aux from
 * HERE single-valued.
 *
 * ⚠ Do NOT move this into after_iter. An earlier draft wrote
 * per_active_final_h there and it was reverted: after_iter can be
 * invoked more than once for the same active, so a converged quantity
 * written from it is not well defined.
 *
 * IterScratch passed for forward-compat: future Specs might persist
 * (e.g.) ConditionNumber from the final iter through scratch into
 * downstream state. Density currently does not need it — finalize
 * re-computes ConditionNumber from accum.NV_T inversion.
 * ==================================================================== */
void DensitySpec::apply_active_writeback_iterative(const neighbor_loop_args& args,
                                                    int active_slot, int /*i*/,
                                                    const AccumData& accum,
                                                    double final_h,
                                                    const IterScratch& /*scratch*/) {
    auto* aux = static_cast<Aux*>(args.aux);
    aux->per_active_final_accum[active_slot] = accum;
    aux->per_active_final_h[active_slot]     = final_h;
}


/* ====================================================================
 * after_iter — bisection convergence + radius update.
 *
 * Literal port of hydro/density.cc:244-590 with the strict design
 * contract:
 *   - READS: accum (this iter's accumulated values), ctx.scratch
 *     (bisection state across iters), ctx.h_search_current (this iter's
 *     radius), ctx.scalars (CallScalars), ctx.iter_index, args.P[ctx.i]
 *     and args.CellP[ctx.i] (read-only field reads).
 *   - WRITES: ONLY ctx.scratch (Left/Right brackets, converged latch,
 *     iter counter, condition_number_current). NEVER writes P/CellP.
 *   - RETURNS: IterResult{ status, new_h_search }. status =
 *     Converged | NeedsMore | AdjustRadius; new_h_search honored only
 *     on AdjustRadius (runner header line 632-642).
 *
 * The legacy in-loop NV_T inversion + ConditionNumber computation
 * happens LOCALLY here (host-side matrix inversion, doesn't write the
 * inverted matrix anywhere); the eventual write of CellP[i].NV_T and
 * CellP[i].ConditionNumber happens in density_finalize_post_runner
 * from Aux's converged accum.
 * ==================================================================== */
IterResult DensitySpec::after_iter(const AfterIterContext<DensitySpec>& ctx,
                                    const AccumData& accum)
{
    const neighbor_loop_args_iterative& args = ctx.args;
    const int          i        = ctx.i;
    const CallScalars& cs       = ctx.scalars;
    DensityIterScratch& scratch = ctx.scratch;
    const int          iter     = ctx.iter_index;

    /* ---- (a) local normalization (mirrors density.cc:247-268). Computes
     * effective NumNgb / DrkernNgbFactor / Particle_DivVel from accum
     * without touching P[i]. ---- */
    double desnumngb    = cs.des_num_ngb;
    double desnumngbdev = cs.max_num_ngb_dev;
    if (cs.time_now == cs.time_begin) {
        if (cs.max_num_ngb_dev > 0.05) desnumngbdev = 0.05;
    }
    const double desnumngbdev_0 = desnumngbdev;

    double NumNgb_eff = (double)accum.Ngb;
    double DrkernNgbFactor_eff = (double)accum.DrkernNgb;
    double Particle_DivVel_eff = (double)accum.Particle_DivVel;
    const double h = (double)ctx.h_search_current;

    if (NumNgb_eff > 0) {
        DrkernNgbFactor_eff *= h / (NUMDIMS * NumNgb_eff);
        Particle_DivVel_eff /= NumNgb_eff;
        NumNgb_eff *= VOLUME_NORM_COEFF_FOR_NDIMS * pow(h, (double)NUMDIMS);
    } else {
        NumNgb_eff = 0;
        DrkernNgbFactor_eff = 0;
        Particle_DivVel_eff = 0;
    }
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
    if (ags_density_isactive(i) && (args.P[i].Type > 0)) {
        Particle_DivVel_eff = 0;
    }
#endif
    if (DrkernNgbFactor_eff > -0.9) {
        DrkernNgbFactor_eff = 1.0 / (1.0 + DrkernNgbFactor_eff);
    } else {
        DrkernNgbFactor_eff = 1.0;
    }
    Particle_DivVel_eff *= DrkernNgbFactor_eff;

    /* ---- (b) NV_T inversion + FaceClosureError computation (gas only,
     * host-side LOCAL — no write back to CellP). Mirrors density.cc:271-303.
     * Produces ConditionNumber AND face_closure_error_local — both used
     * below for ncorr_ngb. The inverted matrix and FaceClosureError are
     * discarded here; finalize re-inverts from Aux for the persistent
     * CellP writes.
     *
     * FaceClosureError nudge at legacy line 326 modifies ncorr_ngb, which
     * changes desnumngb and therefore the convergence trajectory. It
     * cannot be deferred to finalize; the nudge must fire mid-iter. Tinv
     * is kept in scope specifically so the face-leak computation can run. */
    double ConditionNumber = 0.0;
    double face_closure_error_local = 0.0;  /* gas-only; nudge gated below */
    if (args.P[i].Type == 0) {
        const double V_i = VOLUME_NORM_COEFF_FOR_NDIMS * pow(h, (double)NUMDIMS) / (NumNgb_eff > 0 ? NumNgb_eff : 1.0);
        const double dimensional_NV_T_normalizer = pow(h, (double)(2 - NUMDIMS));
        double NV_T_local[3][3];
        for (int k1 = 0; k1 < 3; ++k1) {
            for (int k2 = 0; k2 < 3; ++k2) {
                NV_T_local[k1][k2] = accum.NV_T[k1][k2] / dimensional_NV_T_normalizer;
            }
        }
        const double ConditionNumber_threshold = 10.0 * CONDITION_NUMBER_DANGER;
        const double trace_initial = NV_T_local[0][0] + NV_T_local[1][1] + NV_T_local[2][2];
        double conditioning_term_to_add = 1.05 * (trace_initial / NUMDIMS) / ConditionNumber_threshold;
        double Tinv[3][3];
        while (true) {
            ConditionNumber = matrix_invert_ndims(NV_T_local, Tinv);
            if (ConditionNumber < ConditionNumber_threshold) break;
            for (int k1 = 0; k1 < NUMDIMS; ++k1) {
                NV_T_local[k1][k1] += conditioning_term_to_add;
            }
            conditioning_term_to_add *= 1.2;
        }

        /* Face-leak / FaceClosureError, mirrors density.cc:273-303 verbatim
         * but reading accum.NV_T_face_weights (this iter's accumulator)
         * instead of CellP[i].NV_T_face_weights (which legacy had already
         * scattered into CellP). The legacy NV_T used in the face_out
         * sum is the INVERTED matrix; we use Tinv-with-normalizer-undone
         * to match.
         *
         * Note: legacy line 275 overwrites dx_i with sqrt(V_i * trace).
         * Trace here is the trace of the RAW (un-inverted) NV_T per
         * legacy semantics (CellP[i].NV_T at that program point is the
         * raw accumulated NV_T scattered back from accum, not yet
         * inverted). */
        const double trace_raw_NV_T = (double)accum.NV_T[0][0] + (double)accum.NV_T[1][1] + (double)accum.NV_T[2][2];
        const double dx_i = sqrt(V_i * trace_raw_NV_T);

        double face_in[3];
        face_in[0] = (double)accum.NV_T_face_weights[0];
        face_in[1] = (double)accum.NV_T_face_weights[1];
        face_in[2] = (double)accum.NV_T_face_weights[2];

        /* Full 3x3 inverted (Tinv stored upper triangle in legacy line
         * 292; legacy then reads CellP[i].NV_T as a full symmetric
         * matrix in the face_out sum, with implicit lower-triangle
         * symmetry — matrix_invert_ndims returns full Tinv). */
        double NV_T_inverted[3][3];
        for (int k1 = 0; k1 < 3; ++k1) {
            for (int k2 = 0; k2 < 3; ++k2) {
                NV_T_inverted[k1][k2] = Tinv[k1][k2] / dimensional_NV_T_normalizer;
            }
        }

        double face_out[3] = {0, 0, 0};
        for (int k1 = 0; k1 < 3; ++k1) {
            for (int k2 = 0; k2 < 3; ++k2) {
                face_out[k1] += 2.0 * V_i * NV_T_inverted[k1][k2] * face_in[k2];
            }
        }
        double dimless_face_leak = 0;
        for (int k1 = 0; k1 < 3; ++k1) dimless_face_leak += fabs(face_out[k1]) / (double)NUMDIMS;

#ifdef HYDRO_KERNEL_SURFACE_VOLCORR
        double closure_asymm = 0;
        for (int k1 = 0; k1 < 3; ++k1) closure_asymm += face_in[k1] * face_in[k1];
        const double particle_inverse_volume = NumNgb_eff / (VOLUME_NORM_COEFF_FOR_NDIMS * pow(h, (double)NUMDIMS));
        closure_asymm = sqrt(closure_asymm) / (h * particle_inverse_volume);
        face_closure_error_local = DMIN(DMAX(1.0259 - 2.52444 * closure_asymm, 0.344301), 1.0);
#else
        face_closure_error_local = dimless_face_leak / (2.0 * (double)NUMDIMS * pow(dx_i, (double)(NUMDIMS - 1)));
#endif
    }

    /* ---- (c) iter==0 ConditionNumber update logic (mirrors
     * density.cc:313-322). Reads previous-call stored CellP[i].ConditionNumber;
     * decides whether to "promote" this iter's CN. Result lives in
     * scratch.condition_number_current; legacy's CellP[i] write is
     * deferred to finalize. ---- */
    if (args.P[i].Type == 0) {
        const double c0 = 0.1 * (double)CONDITION_NUMBER_DANGER;
        if ((iter == 0) && (ConditionNumber > args.CellP[i].ConditionNumber)
                        && (args.CellP[i].ConditionNumber > 0)) {
            double ncorr_a = 1.0, cn_a = args.CellP[i].ConditionNumber;
            if (cn_a > c0) ncorr_a = sqrt(1.0 + (cn_a - c0) / (double)CONDITION_NUMBER_DANGER);
            if (ncorr_a > 2.0) ncorr_a = 2.0;
            const double dn_ngb     = fabs(NumNgb_eff - cs.des_num_ngb * ncorr_a) / (desnumngbdev_0 * ncorr_a);

            double ncorr_b = 1.0, cn_b = ConditionNumber;
            if (cn_b > c0) ncorr_b = sqrt(1.0 + (cn_b - c0) / (double)CONDITION_NUMBER_DANGER);
            if (ncorr_b > 2.0) ncorr_b = 2.0;
            const double dn_ngb_alt = fabs(NumNgb_eff - cs.des_num_ngb * ncorr_b) / (desnumngbdev_0 * ncorr_b);

            const double dn_ngb_min = DMIN(dn_ngb, dn_ngb_alt);
            if (dn_ngb_min < 10.0) {
                scratch.condition_number_current = ConditionNumber;  /* legacy line 322's write — but to scratch, NOT CellP */
                scratch.cn_initialized = true;
            }
        }
        if (!scratch.cn_initialized) {
            /* No iter==0 promotion fired: use stored prev-call value. */
            scratch.condition_number_current = args.CellP[i].ConditionNumber;
            scratch.cn_initialized = true;
        }
    }

    /* ---- (d) ncorr_ngb from effective CN + FaceClosureError nudge
     * (mirrors density.cc:324-327). FaceClosureError nudge fires INSIDE
     * the iter loop, modifies ncorr_ngb, changes convergence trajectory.
     * Cannot be deferred. ---- */
    double ncorr_ngb = 1.0;
    if (args.P[i].Type == 0) {
        const double c0 = 0.1 * (double)CONDITION_NUMBER_DANGER;
        const double cn = scratch.condition_number_current;
        if (cn > c0) ncorr_ngb = sqrt(1.0 + (cn - c0) / (double)CONDITION_NUMBER_DANGER);
        if (ncorr_ngb > 2.0) ncorr_ngb = 2.0;
#if !defined(HYDRO_KERNEL_SURFACE_VOLCORR)
        const double d00 = 0.35;
        if (face_closure_error_local > d00) {
            ncorr_ngb = DMAX(ncorr_ngb, DMIN(face_closure_error_local / d00, 2.0));
        }
#endif
    }
    desnumngb    = cs.des_num_ngb * ncorr_ngb;
    desnumngbdev = desnumngbdev_0 * ncorr_ngb;
#if !defined(EOS_ELASTIC)
    if (iter > 1) {
        desnumngbdev = DMIN(0.25 * desnumngb,
                            desnumngbdev * exp(0.1 * log(desnumngb / (16.0 * desnumngbdev)) * (double)iter));
    }
#endif

#ifdef SINK_PARTICLES
    if (args.P[i].Type == 5) {
        desnumngb = cs.des_num_ngb * cs.sink_ngb_factor;
  #ifdef SINGLE_STAR_SINK_DYNAMICS
        desnumngbdev = (cs.sink_ngb_factor + 1.0);
  #else
        desnumngbdev = 4.0 * (cs.sink_ngb_factor + 1.0);
  #endif
    }
#endif

#ifdef GRAIN_FLUID
    if ((1 << args.P[i].Type) & (GRAIN_PTYPES)) {
        desnumngb    = cs.des_num_ngb;
        desnumngbdev = cs.des_num_ngb / 4.0;
  #if defined(GRAIN_BACKREACTION)
        desnumngbdev = desnumngbdev_0;
  #endif
    }
#endif

    double minsoft = cs.min_kernel_radius;
    double maxsoft = cs.max_kernel_radius;

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    /* legacy density.cc:153-159 constants — Config-resolved at compile time */
    int valid_stellar_types = 2+4+8+16, invalid_stellar_types = 1+32;
  #if (defined(GRAIN_FLUID) || defined(RADTRANSFER)) && (!defined(GALSF) && !(defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)))
    valid_stellar_types = 16; invalid_stellar_types = 1+2+4+8+32;
    #ifdef RADTRANSFER
    invalid_stellar_types = 64; valid_stellar_types = RT_SOURCES;
    #endif
  #endif
    if (((1 << args.P[i].Type) & (valid_stellar_types)) && !((1 << args.P[i].Type) & (invalid_stellar_types))) {
        desnumngb = cs.des_num_ngb;
  #if defined(RT_SOURCE_INJECTION)
        if (desnumngb < 64.0) desnumngb = 64.0;
  #endif
  #ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
        desnumngb = DMAX(desnumngb, 128.0);
        if (KERNEL_FUNCTION > 3) desnumngb = DMAX(desnumngb, 256.0);
  #endif
  #ifdef GALSF
        if (desnumngb < 64.0) desnumngb = 64.0;
        const double unitlength_in_kpc = UNIT_LENGTH_IN_KPC * cs.common.cf_atime;
        maxsoft = 2.0 / unitlength_in_kpc;
    #if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && defined(SINK_PARTICLES) && (defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL))
        if ((args.P[i].SNe_ThisTimeStep > 0) || (args.P[i].MassReturn_ThisTimeStep > 0)
                || (cs.time_now == cs.time_begin)) {
            maxsoft = 2.0 / unitlength_in_kpc;
        } else {
            maxsoft = 0.1 / unitlength_in_kpc;
        }
    #endif
  #endif
        desnumngbdev = desnumngb / 2.0;
    }
#endif

#ifdef SINK_PARTICLES
    if (args.P[i].Type == 5) {
        maxsoft = cs.sink_max_accretion_radius / cs.common.cf_atime;
  #ifdef SINGLE_STAR_SINK_DYNAMICS
        minsoft = SinkParticle_GravityKernelRadius;
    #ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
        minsoft = DMAX(minsoft, DMIN(args.P[i].SinkRadius, 0.1 * SinkParticle_GravityKernelRadius));
    #endif
  #endif
    }
#endif

    /* ---- (e) Convergence check + bracket / radius update (mirrors
     * density.cc:398-589). The legacy mutates P[i].KernelRadius
     * directly; we compute new_h_search locally and return via
     * IterResult. min/max kernel flags are per-iter scratch (not
     * IterScratch fields). ---- */
    int redo_particle = 0;
    int particle_set_to_minrkern_flag = 0;
    int particle_set_to_maxrkern_flag = 0;
    double new_h = h;  /* candidate next radius; legacy mutates P[i].KernelRadius */

    /* normal range check (line 401-403) */
    if ((NumNgb_eff < (desnumngb - desnumngbdev) && new_h < 0.999 * maxsoft) ||
        (NumNgb_eff > (desnumngb + desnumngbdev) && new_h > 1.001 * minsoft)) {
        redo_particle = 1;
    }

    /* max kernel size (line 406-420) */
    if ((new_h >= 0.999 * maxsoft) && (NumNgb_eff < (desnumngb - desnumngbdev))) {
        redo_particle = 0;
        if (new_h == maxsoft) {
            particle_set_to_maxrkern_flag = 0;
        } else {
            redo_particle = 1;
            new_h = maxsoft;
            particle_set_to_maxrkern_flag = 1;
        }
    }

    /* min kernel size (line 423-438) */
    if ((new_h <= 1.001 * minsoft) && (NumNgb_eff > (desnumngb + desnumngbdev))) {
        redo_particle = 0;
        if (new_h == minsoft) {
            particle_set_to_minrkern_flag = 0;
        } else {
            redo_particle = 1;
            new_h = minsoft;
            particle_set_to_minrkern_flag = 1;
        }
    }

#ifdef GALSF
    if (cs.common.comoving_integration_on && cs.time_now > cs.time_begin) {
        if ((args.P[i].Type == 4) && (iter > 1) && (NumNgb_eff > 4) && (NumNgb_eff < 100) && (redo_particle == 1)) {
            redo_particle = 0;
        }
    }
#endif

    if ((redo_particle == 0) && (args.P[i].Type == 0)) {
        /* legacy line 447-459: condition-number diagnostic + persist CN
         * to CellP[i].ConditionNumber. Persist is deferred to finalize;
         * here we just stash ConditionNumber into scratch for finalize. */
        scratch.condition_number_current = ConditionNumber;
        /* legacy PRINT_WARNING for ConditionNumber > 1e6 * threshold is a
         * diagnostic only — preserved here for fidelity. */
        if (ConditionNumber > 1e6 * (double)CONDITION_NUMBER_DANGER) {
            PRINT_WARNING("Condition number=%g Cnum_prev=%g threshold=%g iter=%d NumNgb=%g desnumngb=%g KernelRadius=%g i=%d task=%d ID=%llu Type=%d Left=%g Right=%g",
                ConditionNumber, args.CellP[i].ConditionNumber, (double)CONDITION_NUMBER_DANGER,
                iter, NumNgb_eff, desnumngb, new_h, i, ThisTask,
                (unsigned long long) args.P[i].ID, args.P[i].Type,
                (double)scratch.left, (double)scratch.right);
        }
    }

    if (redo_particle) {
        /* legacy line 463-468: iter>=MAXITER-10 warning */
        if (iter >= MAXITER - 10) {
            PRINT_WARNING("i=%d task=%d ID=%llu iter=%d Type=%d KernelRadius=%g Drkern=%g Left=%g Right=%g NumNgb=%g Right-Left=%g maxh=%d minh=%d minsoft=%g maxsoft=%g desnum=%g desnumdev=%g",
                i, ThisTask, (unsigned long long) args.P[i].ID, iter, args.P[i].Type,
                new_h, DrkernNgbFactor_eff,
                (double)scratch.left, (double)scratch.right,
                NumNgb_eff, (double)scratch.right - (double)scratch.left,
                particle_set_to_maxrkern_flag, particle_set_to_minrkern_flag,
                minsoft, maxsoft, desnumngb, desnumngbdev);
        }

        /* tight-bracket early exit (legacy line 474-482).
         * Returns Converged — caller in this case in legacy did `npleft--;
         * P[i].TimeBin = -P[i].TimeBin - 1; continue;` which marks the
         * particle done. */
        if (scratch.left > 0 && scratch.right > 0) {
            if ((double)(scratch.right - scratch.left) < 1.0e-3 * (double)scratch.left) {
                scratch.converged = true;
                scratch.condition_number_current = ConditionNumber;
                /* Final radius travels to finalize through the single-valued
                 * apply_active_writeback_iterative channel. */
                return IterResult{IterStatus::Converged, (double)new_h};
            }
        }

        if ((particle_set_to_maxrkern_flag == 0) && (particle_set_to_minrkern_flag == 0)) {
            /* bracket update (line 486-490) */
            if (NumNgb_eff < (desnumngb - desnumngbdev)) {
                scratch.left = (MyFloat)DMAX(new_h, (double)scratch.left);
            } else {
                if (scratch.right != 0) {
                    if (new_h < (double)scratch.right) scratch.right = (MyFloat)new_h;
                } else {
                    scratch.right = (MyFloat)new_h;
                }
            }

            /* new radius computation (line 493-583) */
            if (scratch.right > 0 && scratch.left > 0) {
                /* both brackets set: geometric interpolation */
                double maxjump = 0;
                if (iter > 1) maxjump = 0.2 * log((double)scratch.right / (double)scratch.left);
                if (NumNgb_eff > 1) {
                    double jumpvar = DrkernNgbFactor_eff * log(desnumngb / NumNgb_eff) / NUMDIMS;
                    if (iter > 1) {
                        if (fabs(jumpvar) < maxjump) {
                            if (jumpvar < 0) jumpvar = -maxjump; else jumpvar = maxjump;
                        }
                    }
                    new_h *= exp(jumpvar);
                } else {
                    new_h *= 2.0;
                }
                if ((new_h < (double)scratch.right) && (new_h > (double)scratch.left)) {
                    if (iter > 1) {
                        const double hfac = exp(maxjump);
                        if (new_h > (double)scratch.right / hfac) new_h = (double)scratch.right / hfac;
                        if (new_h < (double)scratch.left  * hfac) new_h = (double)scratch.left  * hfac;
                    }
                } else {
                    if (new_h > (double)scratch.right) new_h = (double)scratch.right;
                    if (new_h < (double)scratch.left)  new_h = (double)scratch.left;
                    new_h = pow(new_h * (double)scratch.left * (double)scratch.right, 1.0 / 3.0);
                }
            } else {
                if (scratch.right == 0 && scratch.left == 0) {
                    printf("DensitySpec::after_iter: Right==0 && Left==0 && KernelRadius=%g (task=%d i=%d)\n", new_h, ThisTask, i); fflush(stdout);
                    endrun(90001009);
                    return IterResult{IterStatus::Converged, (double)new_h};   /* graceful: bad-stop set; stop iterating this active with the last finite radius; drains at runner completion -> phase poll (per-active: NO immediate collective) */
                }

                if (scratch.right == 0 && scratch.left > 0) {
                    double fac_lim;
                    if (NumNgb_eff > 1) {
                        fac_lim = log(desnumngb / NumNgb_eff) / NUMDIMS;
                    } else {
                        fac_lim = 1.4;
                    }
                    if ((NumNgb_eff < 2.0 * desnumngb) && (NumNgb_eff > 0.1 * desnumngb)) {
                        double slope = DrkernNgbFactor_eff;
                        if (iter > 2 && slope < 1.0) slope = 0.5 * (slope + 1.0);
                        double fac = fac_lim * slope;
                        if (iter >= 4) { if (DrkernNgbFactor_eff == 1.0) fac *= 10.0; }
                        if (fac < fac_lim + 0.231) {
                            new_h *= exp(fac);
                        } else {
                            new_h *= exp(fac_lim + 0.231);
                        }
                    } else {
                        new_h *= exp(fac_lim);
                    }
                }

                if (scratch.right > 0 && scratch.left == 0) {
                    double fac_lim;
                    if (NumNgb_eff > 1) {
                        fac_lim = log(desnumngb / NumNgb_eff) / NUMDIMS;
                    } else {
                        fac_lim = 1.4;
                    }
                    if (fac_lim < -1.535) fac_lim = -1.535;
                    if ((NumNgb_eff < 2.0 * desnumngb) && (NumNgb_eff > 0.1 * desnumngb)) {
                        double slope = DrkernNgbFactor_eff;
                        if (iter > 2 && slope < 1.0) slope = 0.5 * (slope + 1.0);
                        double fac = fac_lim * slope;
                        if (iter >= 4) { if (DrkernNgbFactor_eff == 1.0) fac *= 10.0; }
                        if (fac > fac_lim - 0.231) {
                            new_h *= exp(fac);
                        } else {
                            new_h *= exp(fac_lim - 0.231);
                        }
                    } else {
                        new_h *= exp(fac_lim);
                    }
                }
            }
        }
        /* min/max clamps (line 586-589) */
        if (new_h < minsoft) new_h = minsoft;
        if (particle_set_to_minrkern_flag == 1) new_h = minsoft;
        if (new_h > maxsoft) new_h = maxsoft;
        if (particle_set_to_maxrkern_flag == 1) new_h = maxsoft;

        scratch.iter_count = iter + 1;
        return IterResult{IterStatus::AdjustRadius, new_h};
    }

    /* redo_particle == 0: converged */
    scratch.converged = true;
    if (args.P[i].Type == 0) {
        scratch.condition_number_current = ConditionNumber;
    }
    /* Final radius travels to finalize through the single-valued
     * apply_active_writeback_iterative channel. */
    return IterResult{IterStatus::Converged, (double)new_h};
}

/* ====================================================================
 * density_finalize_post_runner — legacy density.cc:253-304 (in-iter
 * normalization that we deferred) + density.cc:620-795 (post-loop final
 * physics) + density_gpu.cc:336-401 (per-active accum scatter to
 * P[]/CellP[]). Runs ONCE after the runner's iter loop converges.
 *
 * Reads Aux.per_active_final_accum (converged AccumData per active)
 * and Aux.per_active_final_h (converged radius per active). Writes
 * P[i] / CellP[i] freely — this is the only path that does so for
 * density's port. host_all is the canonical accessor for All.* reads.
 *
 * Degenerate empty-active-gas cases are guarded locally where normalization
 * would otherwise divide by zero.
 * ==================================================================== */
void density_finalize_post_runner(const std::vector<int>& active_list_concat,
                                  DensitySpec::Aux& aux,
                                  const struct global_data_all_processes* host_all)
{
    const int N = (int)active_list_concat.size();

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    /* Legacy density.cc:153-159 type-validity constants. Config-resolved
     * at compile time; needed for valid_stellar_types gate around the
     * DensityAroundParticle / GradRho scatter. */
    int valid_stellar_types = 2+4+8+16, invalid_stellar_types = 1+32;
  #if (defined(GRAIN_FLUID) || defined(RADTRANSFER)) && (!defined(GALSF) && !(defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)))
    valid_stellar_types = 16; invalid_stellar_types = 1+2+4+8+32;
    #ifdef RADTRANSFER
    invalid_stellar_types = 64; valid_stellar_types = RT_SOURCES;
    #endif
  #endif
#endif

    /* Parallel-for: legacy uses #pragma omp parallel for schedule(dynamic)
     * around density.cc:629-795. Preserved here. */
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int slot = 0; slot < N; ++slot) {
        const int i = active_list_concat[slot];

        const DensitySpec::AccumData& accum = aux.per_active_final_accum[slot];
        const double h = aux.per_active_final_h[slot];

        /* ---- (1) Per-active normalization (legacy density.cc:253-268;
         * lifted out of the iter loop to honor the no-P/CellP-write
         * contract in after_iter). Writes P[i] freely here. ---- */
        P[i].KernelRadius = (MyFloat)h;
        double NumNgb_eff = (double)accum.Ngb;
        double DrkernNgbFactor_eff = (double)accum.DrkernNgb;
        double Particle_DivVel_eff = (double)accum.Particle_DivVel;
        if (NumNgb_eff > 0) {
            DrkernNgbFactor_eff *= h / (NUMDIMS * NumNgb_eff);
            Particle_DivVel_eff /= NumNgb_eff;
            NumNgb_eff *= VOLUME_NORM_COEFF_FOR_NDIMS * pow(h, (double)NUMDIMS);
        } else {
            NumNgb_eff = DrkernNgbFactor_eff = Particle_DivVel_eff = 0;
        }
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
        if (ags_density_isactive(i) && (P[i].Type > 0)) {
            Particle_DivVel_eff = 0;
        }
#endif
        if (DrkernNgbFactor_eff > -0.9) {
            DrkernNgbFactor_eff = 1.0 / (1.0 + DrkernNgbFactor_eff);
        } else {
            DrkernNgbFactor_eff = 1.0;
        }
        Particle_DivVel_eff *= DrkernNgbFactor_eff;

        P[i].NumNgb           = (MyFloat)NumNgb_eff;
        P[i].DrkernNgbFactor  = (MyFloat)DrkernNgbFactor_eff;
        P[i].Particle_DivVel  = (MyFloat)Particle_DivVel_eff;

        /* ---- (2) Scatter accum → CellP[i] / P[i]: legacy density_gpu.cc:
         * 336-401 per-active fields. Done unconditionally per legacy
         * (writes are guarded by #ifdef where the field is declared). ---- */
        if (P[i].Type == 0) {
            /* Type==0 scatters: CellP[i].* fields. */
            CellP[i].Density = (MyFloat)accum.Rho;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
            CellP[i].ParticleVel[0] = (MyFloat)accum.ParticleVel[0];
            CellP[i].ParticleVel[1] = (MyFloat)accum.ParticleVel[1];
            CellP[i].ParticleVel[2] = (MyFloat)accum.ParticleVel[2];
#endif
            for (int k = 0; k < 6; ++k) {
                CellP[i].NV_T.data[k] = (MyFloat)accum.NV_T.data[k];
            }
            CellP[i].NV_T_face_weights[0] = (MyFloat)accum.NV_T_face_weights[0];
            CellP[i].NV_T_face_weights[1] = (MyFloat)accum.NV_T_face_weights[1];
            CellP[i].NV_T_face_weights[2] = (MyFloat)accum.NV_T_face_weights[2];
#ifdef HYDRO_PARTITION_UNITY_IMPROVE_FD
            CellP[i].GradH_numer[0] = (MyFloat)accum.GradH_numer[0];
            CellP[i].GradH_numer[1] = (MyFloat)accum.GradH_numer[1];
            CellP[i].GradH_numer[2] = (MyFloat)accum.GradH_numer[2];
            CellP[i].GradH_denom    = (MyFloat)accum.GradH_denom;
#endif
#ifdef HYDRO_SPH
            CellP[i].DrkernHydroSumFactor = (MyFloat)accum.DrkernHydroSumFactor;
#endif
#ifdef HYDRO_PRESSURE_SPH
            CellP[i].EgyWtDensity = (MyFloat)accum.EgyRho;
#endif
#if defined(TURB_DRIVING)
            CellP[i].SmoothedVel[0] = (MyFloat)accum.GasVel[0];
            CellP[i].SmoothedVel[1] = (MyFloat)accum.GasVel[1];
            CellP[i].SmoothedVel[2] = (MyFloat)accum.GasVel[2];
#endif
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
            for (int k1 = 0; k1 < 3; ++k1) {
                for (int k2 = 0; k2 < 3; ++k2) {
                    CellP[i].NV_D[k1][k2] = (MyFloat)accum.NV_D[k1][k2];
                    CellP[i].NV_A[k1][k2] = (MyFloat)accum.NV_A[k1][k2];
                }
            }
#endif
        }

#if defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION)
        /* multi-fluid collisional interaction ambient gas scatter. */
        if (IS_PARTICLE_DRAGVALID(P[i].Type, P[i].FluidType)) {
            P[i].Gas_Density = (MyFloat)accum.AmbientGasRho;
            P[i].Gas_Temperature = (MyFloat)accum.Gas_Temperature;
            P[i].Gas_fion = (MyFloat)accum.Gas_fion;
            P[i].Gas_Velocity[0] = (MyFloat)accum.GasVel[0];
            P[i].Gas_Velocity[1] = (MyFloat)accum.GasVel[1];
            P[i].Gas_Velocity[2] = (MyFloat)accum.GasVel[2];
  #if defined(DO_FLUID_DRAG_CALCULATION_WITHBFIELDS)
            P[i].Gas_B[0] = (MyFloat)accum.Gas_B[0];
            P[i].Gas_B[1] = (MyFloat)accum.Gas_B[1];
            P[i].Gas_B[2] = (MyFloat)accum.Gas_B[2];
  #endif
  #if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
            for (int kv = 0; kv < GRAIN_NUM_VOLATILE_SPECIES; ++kv) {
                P[i].Gas_VolatileSpecies[kv] = (MyFloat)accum.Gas_VolatileSpecies[kv];
            }
  #endif
        }
#endif

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
        P[i].DensityAroundParticle = (MyFloat)accum.Rho;
        P[i].GradRho[0] = (MyFloat)accum.GradRho[0];
        P[i].GradRho[1] = (MyFloat)accum.GradRho[1];
        P[i].GradRho[2] = (MyFloat)accum.GradRho[2];
#endif

#if defined(RT_SOURCE_INJECTION)
  #if defined(RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION)
        if (host_all->TimeStep == 0)
  #endif
        if ((1 << P[i].Type) & (RT_SOURCES)) {
            P[i].KernelSum_Around_RT_Source = (MyFloat)accum.KernelSum_Around_RT_Source;
        }
#endif

#if defined(SINK_PARTICLES)
        if (P[i].Type == 5) {
            P[i].Sink_TimeBinGasNeighbor = (short int)accum.Sink_TimeBinGasNeighbor;
  #if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
            P[i].Sink_dr_to_NearestGasNeighbor = (MyFloat)accum.Sink_dr_to_NearestGasNeighbor;
  #endif
        }
#endif

        /* ---- (3) Persistent NV_T inversion + FaceClosureError + CN
         * (gas-i only; mirrors legacy density.cc:271-304 final state).
         * Computes and writes the persistent CellP[i].NV_T, FaceClosureError,
         * ConditionNumber from the converged accum. Recompute from final
         * Aux rather than relying on after_iter's scratch values. ---- */
        if (P[i].Type == 0) {
            const double V_i = VOLUME_NORM_COEFF_FOR_NDIMS * pow(h, (double)NUMDIMS) / (NumNgb_eff > 0 ? NumNgb_eff : 1.0);
            const double dimensional_NV_T_normalizer = pow(h, (double)(2 - NUMDIMS));
            double NV_T_local[3][3];
            for (int k1 = 0; k1 < 3; ++k1) {
                for (int k2 = 0; k2 < 3; ++k2) {
                    NV_T_local[k1][k2] = (double)accum.NV_T[k1][k2] / dimensional_NV_T_normalizer;
                }
            }
            const double trace_raw_NV_T = (double)accum.NV_T[0][0] + (double)accum.NV_T[1][1] + (double)accum.NV_T[2][2];
            const double dx_i = sqrt(V_i * trace_raw_NV_T);

            const double ConditionNumber_threshold = 10.0 * CONDITION_NUMBER_DANGER;
            const double trace_initial = NV_T_local[0][0] + NV_T_local[1][1] + NV_T_local[2][2];
            double conditioning_term_to_add = 1.05 * (trace_initial / NUMDIMS) / ConditionNumber_threshold;
            double Tinv[3][3];
            double ConditionNumber = 0;
            while (true) {
                ConditionNumber = matrix_invert_ndims(NV_T_local, Tinv);
                if (ConditionNumber < ConditionNumber_threshold) break;
                for (int k1 = 0; k1 < NUMDIMS; ++k1) {
                    NV_T_local[k1][k1] += conditioning_term_to_add;
                }
                conditioning_term_to_add *= 1.2;
            }
            /* Write back inverted matrix (legacy line 292): upper triangle only. */
            for (int k1 = 0; k1 < 3; ++k1) {
                for (int k2 = k1; k2 < 3; ++k2) {
                    CellP[i].NV_T[k1][k2] = (MyFloat)(Tinv[k1][k2] / dimensional_NV_T_normalizer);
                }
            }
            CellP[i].ConditionNumber = (MyFloat)ConditionNumber;

            /* FaceClosureError final compute (legacy line 296-303). */
            double face_in[3] = { (double)accum.NV_T_face_weights[0],
                                  (double)accum.NV_T_face_weights[1],
                                  (double)accum.NV_T_face_weights[2] };
            double NV_T_inverted_full[3][3];
            for (int k1 = 0; k1 < 3; ++k1) {
                for (int k2 = 0; k2 < 3; ++k2) {
                    NV_T_inverted_full[k1][k2] = Tinv[k1][k2] / dimensional_NV_T_normalizer;
                }
            }
            double face_out[3] = {0, 0, 0};
            for (int k1 = 0; k1 < 3; ++k1) {
                for (int k2 = 0; k2 < 3; ++k2) {
                    face_out[k1] += 2.0 * V_i * NV_T_inverted_full[k1][k2] * face_in[k2];
                }
            }
            double dimless_face_leak = 0;
            for (int k1 = 0; k1 < 3; ++k1) dimless_face_leak += fabs(face_out[k1]) / (double)NUMDIMS;
#ifdef HYDRO_KERNEL_SURFACE_VOLCORR
            double closure_asymm = 0;
            for (int k1 = 0; k1 < 3; ++k1) closure_asymm += face_in[k1] * face_in[k1];
            const double particle_inverse_volume = NumNgb_eff / (VOLUME_NORM_COEFF_FOR_NDIMS * pow(h, (double)NUMDIMS));
            closure_asymm = sqrt(closure_asymm) / (h * particle_inverse_volume);
            CellP[i].FaceClosureError = (MyFloat)DMIN(DMAX(1.0259 - 2.52444 * closure_asymm, 0.344301), 1.0);
#else
            CellP[i].FaceClosureError = (MyFloat)(dimless_face_leak / (2.0 * (double)NUMDIMS * pow(dx_i, (double)(NUMDIMS - 1))));
#endif
        }

        /* ---- (4) Legacy post-loop "real final operations" block
         * (density.cc:632-748). ---- */
        if (P[i].Type == 0 && P[i].Mass > 0) {
            if (CellP[i].Density > 0) {
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
  #if (HYDRO_FIX_MESH_MOTION==4)
                set_mesh_motion(i);
  #elif ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                {
                    const double eps_pvel = 0.3;
                    /* legacy: CellP[i].ParticleVel = CellP[i].VelPred*(1-eps) + CellP[i].ParticleVel*(eps/Density) */
                    for (int k = 0; k < 3; ++k) {
                        CellP[i].ParticleVel[k] = (MyFloat)(
                            (double)CellP[i].VelPred[k] * (1.0 - eps_pvel)
                            + (double)CellP[i].ParticleVel[k] * (eps_pvel / (double)CellP[i].Density));
                    }
                }
  #elif (HYDRO_FIX_MESH_MOTION==7)
                CellP[i].ParticleVel[0] = CellP[i].VelPred[0];
                CellP[i].ParticleVel[1] = CellP[i].VelPred[1];
                CellP[i].ParticleVel[2] = CellP[i].VelPred[2];
  #endif
#endif

#ifdef HYDRO_SPH
  #ifdef HYDRO_PRESSURE_SPH
                if (CellP[i].InternalEnergyPred > 0) {
                    CellP[i].EgyWtDensity /= CellP[i].InternalEnergyPred;
                } else {
                    CellP[i].EgyWtDensity = 0;
                }
  #endif
                if ((P[i].KernelRadius > 0) && (P[i].NumNgb > 0)) {
                    const double numden_ngb = (double)P[i].NumNgb / (VOLUME_NORM_COEFF_FOR_NDIMS * pow((double)P[i].KernelRadius, (double)NUMDIMS));
                    double DrkernHydroSumFactor_eff = (double)CellP[i].DrkernHydroSumFactor;
                    DrkernHydroSumFactor_eff *= (double)P[i].KernelRadius / (NUMDIMS * numden_ngb);
                    DrkernHydroSumFactor_eff *= -(double)P[i].DrkernNgbFactor;
                    CellP[i].DrkernHydroSumFactor = (MyFloat)DrkernHydroSumFactor_eff;
                } else {
                    CellP[i].DrkernHydroSumFactor = 0;
                }
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
                {
                    /* legacy lines 669-697: physical-units conversion + dt_DivVel/NV_trSSt computation */
                    for (int k1 = 0; k1 < 3; ++k1) {
                        for (int k2 = 0; k2 < 3; ++k2) {
                            CellP[i].NV_D[k2][k1] *= (MyFloat)host_all->cf_a2inv;
                            CellP[i].NV_A[k2][k1] /= (MyFloat)host_all->cf_atime;
                        }
                    }
                    double dtDV[3][3], A[3][3], V[3][3], S[3][3];
                    for (int k1 = 0; k1 < 3; ++k1) {
                        for (int k2 = 0; k2 < 3; ++k2) {
                            V[k1][k2] = (double)CellP[i].NV_D[k1][0] * (double)CellP[i].NV_T[0][k2]
                                      + (double)CellP[i].NV_D[k1][1] * (double)CellP[i].NV_T[1][k2]
                                      + (double)CellP[i].NV_D[k1][2] * (double)CellP[i].NV_T[2][k2];
                            A[k1][k2] = (double)CellP[i].NV_A[k1][0] * (double)CellP[i].NV_T[0][k2]
                                      + (double)CellP[i].NV_A[k1][1] * (double)CellP[i].NV_T[1][k2]
                                      + (double)CellP[i].NV_A[k1][2] * (double)CellP[i].NV_T[2][k2];
                        }
                    }
                    CellP[i].NV_DivVel = (MyFloat)(V[0][0] + V[1][1] + V[2][2]);
                    double NV_trSSt = 0;
                    for (int k1 = 0; k1 < 3; ++k1) {
                        for (int k2 = 0; k2 < 3; ++k2) {
                            dtDV[k1][k2] = A[k1][k2] - (V[k1][0]*V[0][k2] + V[k1][1]*V[1][k2] + V[k1][2]*V[2][k2]);
                            S[k1][k2] = 0.5 * (V[k1][k2] + V[k2][k1]);
                            if (k2 == k1) S[k1][k2] -= (double)CellP[i].NV_DivVel / NUMDIMS;
                            NV_trSSt += S[k1][k2] * S[k1][k2];
                        }
                    }
                    CellP[i].NV_trSSt     = (MyFloat)NV_trSSt;
                    CellP[i].NV_dt_DivVel = (MyFloat)(dtDV[0][0] + dtDV[1][1] + dtDV[2][2]);
                }
#endif

#if defined(TURB_DRIVING)
                if (CellP[i].Density > 0) {
                    CellP[i].SmoothedVel[0] /= CellP[i].Density;
                    CellP[i].SmoothedVel[1] /= CellP[i].Density;
                    CellP[i].SmoothedVel[2] /= CellP[i].Density;
                } else {
                    CellP[i].SmoothedVel[0] = 0;
                    CellP[i].SmoothedVel[1] = 0;
                    CellP[i].SmoothedVel[2] = 0;
                }
#endif
            } // CellP[i].Density > 0

#ifndef HYDRO_SPH
            /* legacy line 711-723: non-SPH density finalize */
            if ((P[i].KernelRadius > 0) && (P[i].NumNgb > 0)) {
                CellP[i].Density = (MyFloat)((double)P[i].Mass * (double)P[i].NumNgb
                    / (VOLUME_NORM_COEFF_FOR_NDIMS * pow((double)P[i].KernelRadius, (double)NUMDIMS)));
            } else {
                if (P[i].KernelRadius <= 0) {
                    CellP[i].Density = 0;
                } else {
                    CellP[i].Density = (MyFloat)((double)P[i].Mass / (VOLUME_NORM_COEFF_FOR_NDIMS * pow((double)P[i].KernelRadius, (double)NUMDIMS)));
                }
            }
#endif
            double Volume_0 = (double)P[i].Mass / (double)CellP[i].Density;

#ifdef HYDRO_PARTITION_UNITY_IMPROVE_FD
            if (CellP[i].GradH_denom != 0) {
                const double inv_denom = -1.0 / (double)CellP[i].GradH_denom;
                Vec3<double> gradH;
                gradH[0] = (double)CellP[i].GradH_numer[0] * inv_denom;
                gradH[1] = (double)CellP[i].GradH_numer[1] * inv_denom;
                gradH[2] = (double)CellP[i].GradH_numer[2] * inv_denom;
                const double fd_correction = 1.0 + 4.0 * KERNEL_AWPMHD_FD_ALPHA * gradH.norm_sq();
                CellP[i].Density /= fd_correction;
                Volume_0 = (double)P[i].Mass / (double)CellP[i].Density;
            }
#endif

#if defined(HYDRO_KERNEL_SURFACE_VOLCORR)
            CellP[i].Density /= CellP[i].FaceClosureError;
            CellP[i].FaceClosureError = (MyFloat)Volume_0;
#endif

#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
            Volume_0 = (double)P[i].Mass / (double)CellP[i].Density;
            if (host_all->Time == host_all->TimeBegin) {
                CellP[i].Density_ExplicitInt = CellP[i].Density;
            } else {
                CellP[i].Density = CellP[i].Density_ExplicitInt;
            }
            CellP[i].FaceClosureError = (MyFloat)Volume_0;
#endif

#ifdef HYDRO_VOLUME_CORRECTIONS
            CellP[i].Volume_1 = CellP[i].Volume_0 = (MyFloat)Volume_0;
#endif

            set_eos_pressure(i, P, CellP);
        } // P[i].Type == 0 && P[i].Mass > 0

#if defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION)
        if (IS_PARTICLE_DRAGVALID(P[i].Type, P[i].FluidType)) {
            if (P[i].Gas_Density > 0) {
                P[i].Gas_Temperature /= P[i].Gas_Density;
#if (defined(COOLING) && !defined(CHIMES)) || (defined(RADTRANSFER) && defined(RT_CHEM_PHOTOION))
                P[i].Gas_fion /= P[i].Gas_Density;
#else
                P[i].Gas_fion = -1; /* no chemistry to supply one; the estimator falls back on this */
#endif
                P[i].Gas_Velocity[0] /= P[i].Gas_Density;
                P[i].Gas_Velocity[1] /= P[i].Gas_Density;
                P[i].Gas_Velocity[2] /= P[i].Gas_Density;
  #if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
                for (int kv = 0; kv < GRAIN_NUM_VOLATILE_SPECIES; ++kv) {
                    P[i].Gas_VolatileSpecies[kv] /= P[i].Gas_Density;
                }
  #endif
            } else {
                P[i].Gas_Temperature = 0; P[i].Gas_fion = -1;
                P[i].Gas_Velocity[0] = P[i].Gas_Velocity[1] = P[i].Gas_Velocity[2] = 0;
  #if defined(DO_FLUID_DRAG_CALCULATION_WITHBFIELDS)
                P[i].Gas_B[0] = P[i].Gas_B[1] = P[i].Gas_B[2] = 0;
  #endif
  #if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
                for (int kv = 0; kv < GRAIN_NUM_VOLATILE_SPECIES; ++kv) {
                    P[i].Gas_VolatileSpecies[kv] = 0;
                }
  #endif
            }
        }
#endif

        /* NumNgb final transform: ^(1/NDIMS). Legacy line 778. */
        if (P[i].NumNgb > 0) {
            P[i].NumNgb = (MyFloat)pow((double)P[i].NumNgb, 1.0 / NUMDIMS);
        } else {
            P[i].NumNgb = 0;
        }

#if defined(MAGNETIC)
        if (P[i].Type == 0) {
            if (CellP[i].recent_refinement_flag == 1) {
                const MyFloat scale = (MyFloat)((double)P[i].Mass / (double)CellP[i].Density);
                CellP[i].BPred[0] = CellP[i].B[0] = CellP[i].BField_prerefinement[0] * scale;
                CellP[i].BPred[1] = CellP[i].B[1] = CellP[i].BField_prerefinement[1] * scale;
                CellP[i].BPred[2] = CellP[i].B[2] = CellP[i].BField_prerefinement[2] * scale;
                CellP[i].BField_prerefinement[0] = 0;
                CellP[i].BField_prerefinement[1] = 0;
                CellP[i].BField_prerefinement[2] = 0;
            }
        }
#endif
        if (P[i].Type == 0) {
            CellP[i].recent_refinement_flag = 0;
        }

#if defined(SINK_WIND_SPAWN_SET_BFIELD_POLTOR)
        if (P[i].Type == 0) {
            if (P[i].ID == host_all->SpawnedWindCellID && CellP[i].IniDen < 0) {
                CellP[i].IniDen = CellP[i].Density;
                const MyFloat scale = (MyFloat)(
                    (host_all->UnitMagneticField_in_gauss / UNIT_B_IN_GAUSS)
                    * ((double)P[i].Mass / (host_all->cf_a2inv * (double)CellP[i].Density)));
                CellP[i].BPred[0] = CellP[i].B[0] = CellP[i].IniB[0] * scale;
                CellP[i].BPred[1] = CellP[i].B[1] = CellP[i].IniB[1] * scale;
                CellP[i].BPred[2] = CellP[i].B[2] = CellP[i].IniB[2] * scale;
            }
        }
#endif
    } // for slot
}

/* ====================================================================
 * density() — runner-driven host driver.
 *
 * Mirrors the AGS pattern at gravity/ags_density_loop.cc::ags_density()
 * with density-specific changes:
 *   - Single gas-only subgroup (no bm partitioning).
 *   - hydro-typed ghost prep + redo (vs ags-typed).
 *   - Aux carries both per_active_final_accum AND per_active_final_h
 *     (the latter populated via apply_active_writeback_iterative).
 *
 * Driver responsibilities:
 *   (1) host_all early exit on TotN_gas <= 0 (legacy density.cc:136).
 *   (2) Hydro-typed ghost prep BEFORE active-list build (legacy 146).
 *   (3) Build active list via density_isactive.
 *   (4) Allocate both Aux vectors to active_list.size().
 *   (5) Build single gas-only subgroup; build args; drive runner.
 *   (6) density_finalize_post_runner reads Aux + writes P/CellP.
 *   (7) gizmo_hydro_density_import_ghosts_fresh_no_drift: build the broad
 *       downstream hydro-oneway ghost pool consumed by cellcorrections_calc,
 *       hydro_gradient_calc, hydro_force. the iterative
 *       runner's exact-query Mode A pool is too narrow / wrong shape to grow
 *       into the downstream one, so we build the broad pool fresh from
 *       converged radii after the runner has cleaned its internal pool.
 *       IMPORT-ONLY (no drift): post-density move_particles at All.Ti_Current
 *       breaks downstream hydro's drift contract (SP4 abort
 *       'no prediction into past allowed').
 *   (8) Timing: CPU_Step[CPU_MISC] += measure_time() at entry; PRINT_STATUS
 *       at end mirroring legacy line 800-802.
 * ==================================================================== */
void density(void)
{
    const struct global_data_all_processes *host_all = nlr_host_all_ptr();
    if (host_all->TotN_gas <= 0) return;

    CPU_Step[CPU_MISC] += measure_time();
    const double t00_truestart = my_second();

    /* (1) Ghost lifecycle: caller does NOT call gizmo_hydro_density_prep_ghosts
     *     before the runner. A pre-runner global ghost import would mutate
     *     state outside the Mode B corridor snapshot — violating the
     *     "Mode B must not touch globals on tiny-N path" rule and
     *     Instead:
     *       - Mode A iterative runner owns Mode A's ghost import,
     *         using args.ghost_safety_factor + Spec::search_mode +
     *         Spec::neighbor_type_mask + the per-active radii_uvm.
     *       - Mode B (local/remote) performs NO ghost import; lazy-drift
     *         on touched candidates only.
     *       - Downstream consumers (cellcorrections_calc, hydro_gradient_calc,
     *         hydro_force) consume a broad hydro-oneway pool built fresh
     *         post-finalize via gizmo_hydro_density_import_ghosts_fresh_no_drift
     *         at step (7) below. The runner's narrow exact-query Mode A pool
     *         is wrong shape for downstream — fresh broad rebuild from
     *         converged radii is the correct handoff. (The older idiom grew the
     *         pre-density pool instead, which only works when that pool was
     *         itself broad.) */
    const double gsl_safety = gizmo_ghost_safety_factor();

    /* (2) Build active list. */
    std::vector<int> active_list_concat = density_build_active_list();

    /* (2.5) Legacy pre-iteration active prepass — port of density.cc:183-205.
     *       Mirrors legacy verbatim: gpu_particles_arena_invalidate() so the
     *       runner re-seeds with freshly written host values, then for each
     *       active particle reset sink fields and clamp KernelRadius to
     *       0.99*maxsoft. Without this, runner's Spec::search_radius can seed
     *       gpu_ngb_list_build with unchecked / oversize / non-finite
     *       P[i].KernelRadius values; iter-1 CSR rebuild then asks for
     *       neighbor lists at radii up to All.MaxKernelRadius (~211) for
     *       ~2.3M actives → device pool overflow → cudaErrorIllegalAddress.
     *       IterScratch.{left,right} are zeroed by the runner at allocation
     *       so the legacy `Left[i] = Right[i] = 0` is already covered. */
    gpu_particles_arena_invalidate();
    {
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES) && defined(GALSF)
        int valid_stellar_types   = 2+4+8+16;
        int invalid_stellar_types = 1+32;
#ifdef BLACK_HOLES
        valid_stellar_types   = 16;
        invalid_stellar_types = 1+2+4+8+32;
#endif
#ifdef RT_SOURCES
        invalid_stellar_types = 64;
        valid_stellar_types   = RT_SOURCES;
#endif
#endif
        for (int i : active_list_concat) {
#ifdef SINK_PARTICLES
            P[i].SwallowID = 0;
#ifdef SINGLE_STAR_SINK_DYNAMICS
            P[i].SwallowTime = MAX_REAL_NUMBER;
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
            P[i].Sink_Ngb_Flag = 0;
#endif
#endif
            double maxsoft = host_all->MaxKernelRadius;
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES) && defined(GALSF)
            if (((1 << P[i].Type) & (valid_stellar_types)) &&
                !((1 << P[i].Type) & (invalid_stellar_types))) {
                maxsoft = 2.0 / (UNIT_LENGTH_IN_KPC * host_all->cf_atime);
            }
#endif
#ifdef SINK_PARTICLES
            /* Enforce the recommended-policy invariant "nothing exceeds
             * MaxKernelRadius" via DMIN (mirrors the drift-path AGS cap in
             * ags_return_maxsoft): keeps this sink's setup KernelRadius <=
             * MaxKernelRadius so the Mode-B per-type node band (capped at
             * MaxKernelRadius) stays a valid upper bound on the sink leaf reach.
             * No effect where SinkMaxAccretionRadius < MaxKernelRadius (every
             * tested config sets MaxKernelRadius huge). */
            if (P[i].Type == 5) {
                maxsoft = DMIN(host_all->MaxKernelRadius,
                               host_all->SinkMaxAccretionRadius / host_all->cf_atime);
            }
#endif
            if ((P[i].KernelRadius < 0) || !isfinite(P[i].KernelRadius) ||
                (P[i].KernelRadius > 0.99 * maxsoft)) {
                P[i].KernelRadius = 0.99 * maxsoft;
            }
        }
    }

    /* (3) Allocate converged output vectors to num_active. */
    DensitySpec::Aux aux;
    aux.per_active_final_accum.resize(active_list_concat.size());
    aux.per_active_final_h    .resize(active_list_concat.size());

    /* (4) Single gas-only subgroup (design §1b). */
    std::vector<NlrSubgroup> subgroups = density_build_subgroups(
        active_list_concat.data(), (int)active_list_concat.size());

    /* (5) Build args + drive runner. */
    neighbor_loop_args_iterative args{};
    static_cast<neighbor_loop_args&>(args) = nlr_default_args();
    args.P                   = P;
    args.CellP               = (host_all->TotN_gas > 0) ? CellP : nullptr;
    args.num_total           = NumPart;
    args.active_list         = active_list_concat.data();
    args.num_active          = (int)active_list_concat.size();
    args.aux                 = &aux;
    args.num_subgroups       = (int)subgroups.size();
    args.subgroups           = subgroups.data();
    args.ghost_safety_factor = gsl_safety;
    run_neighbor_loop_iterative<DensitySpec>(args);

    /* (6) Post-runner finalize: per-active normalization + NV_T inversion
     *     + scatter accum -> P/CellP + legacy 620-795 post-loop physics. */
    density_finalize_post_runner(active_list_concat, aux, host_all);
    /* The finalize wrote each active's converged KernelRadius: a kept neighbour index takes the new
     * reaches before its next walk. */
    gizmo_mark_kernel_radius_dirty_indices(active_list_concat.data(), (int)active_list_concat.size());

    /* (7) NO downstream ghost-pool handoff here. Every downstream consumer owns
     *     its own ghost lifecycle and cleans up any prior pool before importing,
     *     so a broad pool built here is dead work — it is destroyed before any
     *     consumer reads it: ags_density() cleans + re-imports its own pool
     *     (gravity/ags_density_loop.cc); the hydro corridor cleans + imports its
     *     own SYMMETRIC gas pool in Mode A (or nothing in Mode B) at
     *     gizmo_hydro_corridor_begin() -> gizmo_gradients_prep_symlist(); and
     *     force_update_hmax() between them reads no ghost pool. The prior
     *     unconditional gizmo_hydro_density_import_ghosts_fresh_no_drift() call
     *     was a legacy leftover from before the corridor owned the gas ghost
     *     lifecycle; under forced universal Mode-B it fired a large ONEWAY
     *     import (up to ~7s/all-active step) that nothing consumed, contaminating
     *     the density wall. (Drift contract note, still load-bearing for any
     *     future post-density ghost work: NEVER call gizmo_hydro_density_prep_ghosts
     *     here — its move_particles(All.Ti_Current) breaks the next phase's
     *     lazy-drift contract with 'no prediction into past allowed'. Post-density
     *     ghost work, if ever reintroduced, must be import-only / no-drift.) */

    /* (8) Timing. Mirrors legacy density.cc:798-802 (PRINT_STATUS at end). */
    const double t1 = my_second();
    cpu_chain_sync(t1);
    const double timeall = timediff(t00_truestart, t1);
    if (ThisTask == 0) {
        PRINT_STATUS("  density computation done (%.4f s) via runner", timeall);
    }
    /* NOTE: cellcorrections_calc is a SEPARATE step phase invoked from
     * core/accel.cc:90 after density() returns. NOT called from here. */
}
/* GPU All_dev sync stub. density_loop.cc is in GPU_OBJS (compiled with
 * nvcc_wrapper because pair_kernel uses Kokkos atomics). Density-side
 * device kernel runs inside the runner's TU (mesh/neighbor_loop_runner.cc),
 * which carries its own All_dev sync; this is a no-op stub keeping
 * cooling.cc's per-TU sync sweep happy if a future build adds an extern
 * declaration in cooling.cc. */
