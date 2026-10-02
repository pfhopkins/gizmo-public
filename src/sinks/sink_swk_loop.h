/* sinks/sink_swk_loop.h — sink "swallow-and-kick" neighbor loop module.
 *
 * Defines SinkSwkSpec for the neighbor-loop runner (mesh/neighbor_loop_runner.h).
 * Per the standard 2-file convention, this header carries the Spec contract +
 * inline (host+device) physics that the runner must inline into the device
 * kernel; sinks/sink_swk_loop.cc carries the host-only hooks.
 *
 * Replaces sinks/sink_swallow_and_kick_gpu.cc. The legacy
 * sink_swallow_pair_kernel in sinks/sink_swallow_and_kick_functions.h is
 * mirrored below as sink_swk_pair_kernel with the runner's pair contract:
 *   (const SinkSwkActiveState&, particle_data&, gas_cell_data*,
 *    SinkSwallowOut&).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef SINK_SWK_LOOP_H
#define SINK_SWK_LOOP_H

/* Kokkos_Core must precede allvars.h (its macros may conflict with stdlib
 * names). */
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"

#ifdef SINK_PARTICLES

#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"   /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
#include "sink_functions.h"                /* sink_vesc_gpu, sink_check_boundedness_gpu, sink_fb_angleweight_localcoupling_gpu */
#if defined(COSMIC_RAY_FLUID) && defined(SINK_COSMIC_RAYS)
#include "../eos/cosmic_ray_fluid/cosmic_ray_functions.h"
#endif
/* NOTE: caller TUs must include "../mesh/kernel.h" BEFORE this header
 * (kernel_main, NEAREST_XYZ; same convention as sink_feed_loop.h). */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* Forward decls. */
int sink_isactive(int i);

/* ============================================================================
 * Per-source (active sink) input + output structs. Single source of truth
 * for the sink_swk physics types; previously lived in
 * sink_swallow_and_kick_functions.h alongside a duplicate pair-kernel body
 * (now retired).
 *
 * SinkSwallowLocalIn carries per-active host-staged inputs; SinkSwallowOut
 * is the per-active accumulator (also SinkSwkSpec::AccumData). The runner's
 * Mode B remote path transports these byte-for-byte across MPI, so both must
 * be trivially copyable.
 * ========================================================================== */

struct SinkSwallowLocalIn {
    Vec3<MyDouble> Pos;
    Vec3<MyFloat>  Vel;
    MyFloat KernelRadius;
    MyFloat Mass;
    MyFloat Sink_Mass;
    MyFloat Dt;
    MyFloat Mdot;
    MyIDType ID;
    MyIDType claim_token;   /* must match what sink_feed wrote into SwallowID */
    MyIDType ID_child_number;
    MyIDType ID_generation;
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS) || defined(SINK_WIND_KICK)
    Vec3<MyFloat> Jgas_in_Kernel;
#endif
#ifdef SINK_ALPHADISK_ACCRETION
    MyFloat Sink_Mass_Reservoir;
#endif
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS)
    MyFloat Sink_angle_weighted_kernel_sum;
#endif
#if defined(SINK_RETURN_ANGMOM_TO_GAS)
    Vec3<MyFloat> Sink_Specific_AngMom;
    MyFloat angmom_norm_topass_in_swallowloop;
#endif
#if defined(SINK_RETURN_BFLUX)
    Vec3<MyFloat> B;
    MyFloat kernel_norm_topass_in_swallowloop;
#endif
#ifdef SINGLE_STAR_FB_LOCAL_RP
    MyFloat Luminosity;
#endif
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS)
    /* mom_budget precomputed on the host via sink_lum_bol (host-only).
     * Read on the device by SinkSwkSpec::load_active. */
    MyFloat mom_budget;
#endif
};

struct SinkSwallowOut {
    MyDouble accreted_Mass;
    MyDouble accreted_Sink_Mass;
    MyDouble accreted_Sink_Mass_reservoir;
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
    MyDouble Sink_AccretionDeficit;
#endif
#ifdef GRAIN_FLUID
    MyDouble accreted_dust_Mass;
#endif
#ifdef RT_REINJECT_ACCRETED_PHOTONS
    MyDouble accreted_photon_energy;
#endif
#if defined(SINK_FOLLOW_ACCRETED_MOMENTUM)
    Vec3<MyDouble> accreted_momentum;
#endif
#if defined(SINK_FOLLOW_ACCRETED_COM)
    Vec3<MyDouble> accreted_centerofmass;
#endif
#if defined(SINK_RETURN_BFLUX)
    Vec3<MyDouble> accreted_B;
#endif
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
    Vec3<MyDouble> accreted_J;
#endif
#ifdef SINK_COUNTPROGS
    int Sink_CountProgs;
#endif
#ifdef GALSF
    MyFloat Accreted_Age;
#endif
    int n_gas_swallowed;
    int n_sink_swallowed;
    int n_star_swallowed;
    int n_dm_swallowed;
    /* Per-bin deltas for TimeBin_Sink_* (sink-sink merger Type-5 swallow);
     * scattered into globals on the host since those globals are not
     * device-accessible. */
    MyDouble delta_TimeBin_Sink_mass[TIMEBINS];
    MyDouble delta_TimeBin_Sink_dynamicalmass[TIMEBINS];
    MyDouble delta_TimeBin_Sink_Mdot[TIMEBINS];
    MyDouble delta_TimeBin_Sink_Medd[TIMEBINS];
};

/* ============================================================================
 * Per-pair physics types (file scope; PascalCase). SinkSwkSpec re-exports
 * via `using CallScalars = ...; using ActiveData = ...; using AccumData = ...`.
 * ========================================================================== */

/* Per-call cosmology + sink globals captured once from All.* on the host
 * by populate_call_scalars. Passed to the inline pair body as its own
 * argument, once per call, so the body never reads All.* directly (TRAP 1).
 *
 * No rng_step — sink_swk has no RNG sites in the pair body (audit
 * confirmed via grep: no gizmo_gpu_rand* in sink_swallow_and_kick_functions.h
 * or this loop's pair body).
 *
 * The 6 sink-specific scalars below correspond to the audit table of
 * All.* reads (10 distinct reads; 3 covered by NlrCommonScalars; 1
 * reused from sink_feed; 6 NEW here). */
struct SinkSwkCallScalars {
    NlrCommonScalars common;                       /* cf_atime, cf_a2inv, cf_hubble_a, comoving_integration_on, ... */
    double sink_radius_grav;                       /* SinkParticle_GravityKernelRadius */
    /* The next four fields are populated only under their matching #ifdef gates
     * in sink_swk_loop.cc::populate_call_scalars (zero otherwise). Pair-body
     * readers are themselves under matching #ifdefs, so the fallback zero is
     * never observed. Declared unconditionally so the struct shape doesn't
     * shift with Config and populate code can do unconditional zero-init. */
    double sink_accreted_fraction;                 /* All.Sink_accreted_fraction (SINK_WIND_KICK||SINK_WIND_SPAWN) */
    double sink_feedback_factor;                   /* All.SinkFeedbackFactor */
    double sink_radiative_efficiency;              /* All.SinkRadiativeEfficiency */
    double sink_outflow_velocity;                  /* All.Sink_outflow_velocity (SINK_WIND_KICK||SINK_WIND_SPAWN) */
    double sink_rad_momentum_factor;               /* All.Sink_Rad_MomentumFactor (SINK_PHOTONMOMENTUM) */
#if defined(COSMIC_RAY_FLUID) && defined(SINK_COSMIC_RAYS)
    double sink_cr_injection_efficiency;           /* All.Sink_CosmicRay_Injection_Efficiency */
#endif
#ifdef GALSF_SUBGRID_WINDS
    double wind_free_travel_max_time_factor;       /* All.WindFreeTravelMaxTimeFactor */
#endif
};

/* Active-particle state passed into the pair body. Same shape as
 * SinkFeedActiveState — `pos` and `h_search` top-level (Mode B walker
 * reads them); `local: SinkSwallowLocalIn` carries the host-fill struct;
 * origin_* carries Mode-B-remote requester identity. Per-call scalars are
 * NOT here -- they are one value for the whole call and are passed alongside
 * this record.
 *
 * mom_budget / J_dir / sink_mass_withdisk are precomputed in load_active
 * from the host-staged local and the call scalars, and reused by the pair body
 * across all j (mirrors the per-active precomputation block in legacy
 * sink_swallow_and_kick_gpu.cc:265-285). */
struct SinkSwkActiveState {
    Vec3<double>           pos;
    double                 h_search;
    SinkSwallowLocalIn     local;
    int                    origin_local_idx;
    int                    origin_rank;
    /* Per-active precomputations. */
    double                 mom_budget;
    Vec3<double>           J_dir;
    double                 sink_mass_withdisk;
};

/* DeviceContext extension. Carries:
 *   per_active_local — UVM array of SinkSwallowLocalIn[num_active], host-staged
 *                      by populate_device_context, read on device by load_active.
 *
 * Trivially copyable. Mirrors SinkFeedDeviceContext exactly. */
struct SinkSwkDeviceContext : NeighborLoopDeviceContextBase {
    const SinkSwallowLocalIn *per_active_local;   /* UVM, [num_active] */
};

/* ============================================================================
 * Inline pair body — single source of truth for sink_swk per-pair physics.
 * Mirrors legacy sink_swallow_pair_kernel
 * (sinks/sink_swallow_and_kick_functions.h:108-440) with the new signature:
 *   (const SinkSwkActiveState& active, particle_data& neighbor_particle,
 *    gas_cell_data* neighbor_cell, SinkSwallowOut& out)
 *
 * The pair body writes j-side atomics (P[j].Mass / Vel / dp / Sink_Mass /
 * Sink_Mdot / Sink_Mass_Reservoir; CellP[j].MassTrue / VelPred /
 * InternalEnergy / InternalEnergyPred / B / BPred / DelayTime / CR fields)
 * alongside the i-side accumulator writes (out.*).
 *
 * CR fix vs legacy: the runner port unifies CR-field handling across both
 * modes. Legacy Mode A imported-ghost path silently dropped CR fields
 * (ghost_writeback_sinkswallow had no CR delta record); Mode B local
 * already applied them in place. The compound ghost-writeback callback
 * in sink_swk_loop.cc now ships CR deltas through reverse-comm in Mode A
 * too, so both modes correctly apply CR injections to home-rank gas.
 * Phil-confirmed 2026-05-09 — fix, not a divergence.
 * ========================================================================== */

KOKKOS_INLINE_FUNCTION
static void sink_swk_pair_kernel(const SinkSwkActiveState& active,
                                  const SinkSwkCallScalars& scalars,
                                  struct particle_data& neighbor_particle,
                                  struct gas_cell_data* neighbor_cell,
                                  SinkSwallowOut& out)
{
    const SinkSwallowLocalIn& local      = active.local;
    const double              h_i        = (double)local.KernelRadius;

    if(neighbor_particle.Mass <= 0) return;
#ifdef HYDRO_MULTIFLUID_DM
    /* skip dark-fluid gas: not coupled to baryonic sinks via accretion or feedback. */
    if(neighbor_particle.Type == 0 && neighbor_particle.FluidType == FLUID_DM) return;
#endif

    /* dpos / r2 with periodic wrap (mirrors sink_swallow_and_kick_gpu.cc:291-293). */
    Vec3<double> dpos;
    dpos[0] = (double)neighbor_particle.Pos[0] - (double)local.Pos[0];
    dpos[1] = (double)neighbor_particle.Pos[1] - (double)local.Pos[1];
    dpos[2] = (double)neighbor_particle.Pos[2] - (double)local.Pos[2];
    nearest_xyz(dpos, -1);
    double r2 = dpos.norm_sq();

    /* Re-bind locals for clarity (legacy used these names). The mass, velocity
     * and internal energy are read atomically: another sink's lane may be
     * kicking this same unswallowed cell through the atomic adds below at this
     * moment, and a plain read beside them is a torn value, not an old one. */
    MyIDType OriginallyMarkedSwallowID = neighbor_particle.SwallowID;
    double Mass_j = (double)Kokkos::atomic_load(&neighbor_particle.Mass);
    double Mass_j_0 = Mass_j;
    double InternalEnergy_j = 0, InternalEnergy_j_0 = 0;
    if(neighbor_particle.Type == 0 && neighbor_cell) {
        InternalEnergy_j   = (double)Kokkos::atomic_load(&neighbor_cell->InternalEnergy);
        InternalEnergy_j_0 = InternalEnergy_j;
    }
    double Vel_j[3], Vel_j_0[3];
    for(int k = 0; k < 3; k++) { Vel_j[k] = (double)Kokkos::atomic_load(&neighbor_particle.Vel[k]); Vel_j_0[k] = Vel_j[k]; }

    Vec3<double> dvel{Vel_j[0] - (double)local.Vel[0],
                      Vel_j[1] - (double)local.Vel[1],
                      Vel_j[2] - (double)local.Vel[2]};
    NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos, neighbor_particle.Pos, dvel, -1);

#if defined(SINK_RETURN_ANGMOM_TO_GAS) || defined(SINK_RETURN_BFLUX)
    double wk = 0, dwk_tmp = 0;
    if(neighbor_particle.Type == 0) {
        double hj = (double)neighbor_particle.KernelRadius;
        double h_denom = (h_i > hj) ? h_i : hj;
        double u = sqrt(r2) / h_denom;
        if(u < 1) { kernel_main(u, 1.0, 1.0, &wk, &dwk_tmp, -1); } else { wk = 0; }
    }
#endif
#if defined(SINK_RETURN_ANGMOM_TO_GAS)
    if(neighbor_particle.Type == 0) {
        Vec3<double> sam{(double)local.Sink_Specific_AngMom[0],
                         (double)local.Sink_Specific_AngMom[1],
                         (double)local.Sink_Specific_AngMom[2]};
        Vec3<double> dlv = cross(sam, dpos) * (wk * (double)local.angmom_norm_topass_in_swallowloop);
        for(int k = 0; k < 3; k++) Vel_j[k] += dlv[k];
        out.accreted_momentum[0] -= Mass_j * dlv[0];
        out.accreted_momentum[1] -= Mass_j * dlv[1];
        out.accreted_momentum[2] -= Mass_j * dlv[2];
        Vec3<double> cross_dl = cross(dpos, dlv);
        out.accreted_J[0] -= Mass_j * cross_dl[0];
        out.accreted_J[1] -= Mass_j * cross_dl[1];
        out.accreted_J[2] -= Mass_j * cross_dl[2];
    }
#endif
#if defined(SINK_RETURN_BFLUX)
    if(neighbor_particle.Type == 0 && neighbor_cell &&
       (double)local.kernel_norm_topass_in_swallowloop > 0) {
        double b_frac = DMIN(0.1, (double)local.Dt / ((double)local.Sink_Mass_Reservoir / (double)local.Mdot))
                        * wk / (double)local.kernel_norm_topass_in_swallowloop;
        for(int k = 0; k < 3; k++) {
            double dB = b_frac * (double)local.B[k];
            Kokkos::atomic_add(&neighbor_cell->B[k],     (MyFloat)dB);
            Kokkos::atomic_add(&neighbor_cell->BPred[k], (MyFloat)dB);
            out.accreted_B[k] -= dB;
        }
    }
#endif

    double f_accreted = 0;

    if(neighbor_particle.SwallowID == local.claim_token && Mass_j > 0 && r2 > 0) {
        f_accreted = 1;
#ifdef SINK_WIND_KICK
        if(neighbor_particle.Type == 0) {
            f_accreted = scalars.sink_accreted_fraction;
#ifndef SINK_GRAVCAPTURE_GAS
            if((scalars.sink_feedback_factor > 0) && (scalars.sink_feedback_factor != 1.0)) {
                f_accreted /= scalars.sink_feedback_factor;
            } else {
                if(scalars.sink_outflow_velocity > 0) {
                    f_accreted = 1.0 / (1.0 + fabs(1.0 * SINK_WIND_KICK)
                                  * scalars.sink_radiative_efficiency
                                  * C_LIGHT_CODE / scalars.sink_outflow_velocity);
                }
            }
            if((active.sink_mass_withdisk - (double)local.Mass) <= 0) { f_accreted = 0; }
#endif
        }
#endif

        double mcount_for_conserve = f_accreted * Mass_j;
#if (SINK_FOLLOW_ACCRETED_ANGMOM == 1)
        if(neighbor_particle.Type != 5) { mcount_for_conserve = 0; }
        else                            { mcount_for_conserve = (double)neighbor_particle.Sink_Mass; }
#ifdef SINK_ALPHADISK_ACCRETION
        if(neighbor_particle.Type == 5) { mcount_for_conserve += (double)neighbor_particle.Sink_Mass_Reservoir; }
#endif
#endif
#ifdef GRAIN_FLUID
        if((1 << neighbor_particle.Type) & GRAIN_PTYPES) { out.accreted_dust_Mass += Mass_j; }
#endif
#ifdef RT_REINJECT_ACCRETED_PHOTONS
        if(neighbor_particle.Type == 0 && neighbor_cell) {
            double photon_energy = 0;
            for(int kf = 0; kf < N_RT_FREQ_BINS; kf++) photon_energy += (double)neighbor_cell->Rad_E_gamma[kf];
            out.accreted_photon_energy += photon_energy;
        }
#endif
#if defined(SINK_FOLLOW_ACCRETED_MOMENTUM)
        out.accreted_momentum[0] += mcount_for_conserve * dvel[0];
        out.accreted_momentum[1] += mcount_for_conserve * dvel[1];
        out.accreted_momentum[2] += mcount_for_conserve * dvel[2];
#endif
#if defined(SINK_FOLLOW_ACCRETED_COM)
        out.accreted_centerofmass[0] += mcount_for_conserve * dpos[0];
        out.accreted_centerofmass[1] += mcount_for_conserve * dpos[1];
        out.accreted_centerofmass[2] += mcount_for_conserve * dpos[2];
#endif
#ifdef SINK_RETURN_BFLUX
        if(neighbor_cell) {
            out.accreted_B[0] += (double)neighbor_cell->BPred[0];
            out.accreted_B[1] += (double)neighbor_cell->BPred[1];
            out.accreted_B[2] += (double)neighbor_cell->BPred[2];
        }
#endif
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
        {
            Vec3<double> cdv = cross(dpos, dvel);
            out.accreted_J[0] += mcount_for_conserve * cdv[0];
            out.accreted_J[1] += mcount_for_conserve * cdv[1];
            out.accreted_J[2] += mcount_for_conserve * cdv[2];
            if(neighbor_particle.Type == 5) {
                Vec3<double> psam{(double)neighbor_particle.Sink_Specific_AngMom[0],
                                  (double)neighbor_particle.Sink_Specific_AngMom[1],
                                  (double)neighbor_particle.Sink_Specific_AngMom[2]};
                out.accreted_J[0] += mcount_for_conserve * psam[0];
                out.accreted_J[1] += mcount_for_conserve * psam[1];
                out.accreted_J[2] += mcount_for_conserve * psam[2];
            }
        }
#endif

        if(neighbor_particle.Type == 5) {
            /* Sink-sink merger */
#ifdef SINK_INCREASE_DYNAMIC_MASS
            double acc_mass = ((double)neighbor_particle.Sink_Mass > Mass_j / SINK_INCREASE_DYNAMIC_MASS)
                                ? (double)neighbor_particle.Sink_Mass : (Mass_j / SINK_INCREASE_DYNAMIC_MASS);
            out.accreted_Mass += acc_mass;
#else
            out.accreted_Mass += Mass_j;
#endif
            out.accreted_Sink_Mass += (double)neighbor_particle.Sink_Mass;
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
            out.Sink_AccretionDeficit += (double)neighbor_particle.Sink_AccretionDeficit;
#endif
#ifdef SINK_ALPHADISK_ACCRETION
            out.accreted_Sink_Mass_reservoir += (double)neighbor_particle.Sink_Mass_Reservoir;
#endif
#ifdef SINK_WIND_SPAWN
#ifdef SINGLE_STAR_FB_JETS
            const double swallowed_unspawned = (double)neighbor_particle.unspawned_wind_mass + (double)neighbor_particle.unspawned_jet_mass; // fold the swallowed sink's jet reservoir in too, so a merger cannot strand or lose it
#else
            const double swallowed_unspawned = (double)neighbor_particle.unspawned_wind_mass;
#endif
#ifdef SINK_ALPHADISK_ACCRETION
            out.accreted_Sink_Mass_reservoir += swallowed_unspawned;
#else
            out.accreted_Sink_Mass          += swallowed_unspawned;
#endif
#endif
#ifdef SINK_COUNTPROGS
            out.Sink_CountProgs += neighbor_particle.Sink_CountProgs;
#endif
            int bin = neighbor_particle.TimeBin;
            Kokkos::atomic_add(&out.delta_TimeBin_Sink_mass[bin],          -(double)neighbor_particle.Sink_Mass);
            Kokkos::atomic_add(&out.delta_TimeBin_Sink_dynamicalmass[bin], -Mass_j);
            Kokkos::atomic_add(&out.delta_TimeBin_Sink_Mdot[bin],          -(double)neighbor_particle.Sink_Mdot);
            if((double)neighbor_particle.Sink_Mass > 0) {
                Kokkos::atomic_add(&out.delta_TimeBin_Sink_Medd[bin],
                                   -(double)neighbor_particle.Sink_Mdot / (double)neighbor_particle.Sink_Mass);
            }
            Mass_j = 0;
#ifdef SINK_ALPHADISK_ACCRETION
            Kokkos::atomic_store(&neighbor_particle.Sink_Mass_Reservoir, (MyFloat)0);
#endif
            Kokkos::atomic_store(&neighbor_particle.Sink_Mdot, (MyFloat)0);
            Kokkos::atomic_store(&neighbor_particle.Sink_Mass, (MyFloat)0);
#ifdef GALSF
            /* MIN-merge with sentinel (zero_accum sets Accreted_Age = MAX_REAL_NUMBER). */
            if((MyFloat)neighbor_particle.StellarAge < out.Accreted_Age) {
                out.Accreted_Age = (MyFloat)neighbor_particle.StellarAge;
            }
#endif
            out.n_sink_swallowed++;
        }

#if defined(SINK_GRAVCAPTURE_NONGAS)
        if(neighbor_particle.Type > 0 && neighbor_particle.Type < 5) {
            out.accreted_Mass += Mass_j;
            if(!is_galsf_stellar_candidate_type(neighbor_particle.Type, scalars.common.comoving_integration_on)) { // non-star (DM/collisionless) over types 1-4
                out.accreted_Sink_Mass += Mass_j;
                out.n_dm_swallowed++;
            } else {
#ifdef SINK_ALPHADISK_ACCRETION
                out.accreted_Sink_Mass_reservoir += Mass_j;
#else
                out.accreted_Sink_Mass           += Mass_j;
#endif
                out.n_star_swallowed++;
            }
            Mass_j = 0;
        }
#endif

        if(neighbor_particle.Type == 0) {
            out.accreted_Mass += f_accreted * Mass_j;
#ifdef SINK_GRAVCAPTURE_GAS
#ifdef SINK_ALPHADISK_ACCRETION
            out.accreted_Sink_Mass_reservoir += f_accreted * Mass_j;
#else
            out.accreted_Sink_Mass           += f_accreted * Mass_j;
#endif
#endif
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
            out.Sink_AccretionDeficit -= f_accreted * Mass_j;
#endif
            Mass_j *= (1.0 - f_accreted);
#ifdef SINK_WIND_KICK
            {
                double v_kick = scalars.sink_outflow_velocity;
                Vec3<double> dir{dpos[0], dpos[1], dpos[2]};
#if (SINK_WIND_KICK < 0)
                if(dot(dir, active.J_dir) > 0) { dir = active.J_dir; } else { dir = -active.J_dir; }
#endif
#if defined(COSMIC_RAY_FLUID) && defined(SINK_COSMIC_RAYS)
                /* CR injection alongside the wind kick. The injected CR
                 * deltas land in CellP[j].CosmicRay* via atomic_add; for
                 * Mode A imported ghosts, the compound ghost-writeback
                 * callback in sink_swk_loop.cc reverse-communicates these
                 * deltas to the home rank (the CR-fields fix vs legacy).
                 *
                 * Bug fix vs legacy sink_swallow_pair_kernel
                 * (sinks/sink_swallow_and_kick_functions.h, now retired):
                 * the legacy CR helper call passed
                 *   evaluate_cr_transport_reductionfactor(j, kc, 0, &CellP[j])
                 * which inside the helper indexed cell[target] = cell[j]
                 * = (&CellP[j])[j] = CellP[2j] -- out-of-bounds for
                 * j > num_total/2, wrong-cell read otherwise. Canonical
                 * call shape (matches galaxy_sf/mechanical_fb_functions.h)
                 * is (target=0, cell=&CellP[j]), which we use here.
                 * Physics correctness wins; both runner Modes A and B now
                 * use this corrected indexing (no Mode A/Mode B asymmetry
                 * introduced).
                 *
                 * CR_energy_spectrum_injection_fraction is unaffected:
                 * the helper passes the cell arg only to
                 * return_CRbin_CR_rigidity_in_GV(target=-1, ...) whose
                 * target<0 branch never indexes cell (returns the global
                 * All.CR_global_rigidity_at_bin_center[] value). */
                if(neighbor_cell) {
                    double dEcr = scalars.sink_cr_injection_efficiency * Mass_j *
                                  (scalars.sink_accreted_fraction
                                   / (1.0 - scalars.sink_accreted_fraction)) *
                                  C_LIGHT_CODE * C_LIGHT_CODE;
                    if(dEcr > 0) {
                        double f_inj[N_CR_PARTICLE_BINS]; f_inj[0] = 1.0;
#if (N_CR_PARTICLE_BINS > 1)
                        double sum_in = 0.0;
                        for(int kc = 0; kc < N_CR_PARTICLE_BINS; kc++) {
                            f_inj[kc] = CR_energy_spectrum_injection_fraction(
                                kc, 5, scalars.sink_outflow_velocity, 0,
                                0, &neighbor_particle, neighbor_cell);
                            sum_in += f_inj[kc];
                        }
                        if(sum_in > 0) { for(int kc = 0; kc < N_CR_PARTICLE_BINS; kc++) f_inj[kc] /= sum_in; }
                        else           { for(int kc = 0; kc < N_CR_PARTICLE_BINS; kc++) f_inj[kc] = 1.0 / N_CR_PARTICLE_BINS; }
#endif
                        Vec3<double> cr_dir = {dir[0], dir[1], dir[2]};
                        double cr_dir_norm2 = cr_dir.norm_sq();
                        if(cr_dir_norm2 <= 0) { cr_dir[0] = 0; cr_dir[1] = 0; cr_dir[2] = 1; cr_dir_norm2 = 1; }
                        double cr_dir_norm = sqrt(cr_dir_norm2);
                        for(int kc = 0; kc < N_CR_PARTICLE_BINS; kc++) {
                            double dEcr_k = evaluate_cr_transport_reductionfactor(0, kc, 0, neighbor_cell)
                                            * dEcr * f_inj[kc];
                            if(dEcr_k <= 0) continue;
                            Kokkos::atomic_add(&neighbor_cell->CosmicRayEnergy[kc],     dEcr_k);
                            Kokkos::atomic_add(&neighbor_cell->CosmicRayEnergyPred[kc], dEcr_k);
                            double flux_mag = dEcr_k * CRFLUID_REDUCED_C_CODE(kc);
                            for(int kk = 0; kk < 3; kk++) {
                                double dflux = flux_mag * cr_dir[kk] / cr_dir_norm;
                                Kokkos::atomic_add(&neighbor_cell->CosmicRayFlux[kc][kk],     dflux);
                                Kokkos::atomic_add(&neighbor_cell->CosmicRayFluxPred[kc][kk], dflux);
                            }
                        }
                    }
                }
#endif
                double nrm = dir.norm_sq();
                if(nrm <= 0) { dir[0] = 0; dir[1] = 0; dir[2] = 1; }
                else         { nrm = sqrt(nrm); dir /= nrm; }
                for(int k = 0; k < 3; k++) { Vel_j[k] += v_kick * scalars.common.cf_atime * dir[k]; }
#ifdef GALSF_SUBGRID_WINDS
                if(neighbor_cell) {
                    Kokkos::atomic_store(&neighbor_cell->DelayTime,
                        (MyFloat)(scalars.wind_free_travel_max_time_factor / scalars.common.cf_hubble_a));
                }
#endif
            }
#endif
            out.n_gas_swallowed++;
        }
    } /* end of (SwallowID matches) branch */

#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS)
    if(active.mom_budget > 0 && (double)local.Dt > 0 && OriginallyMarkedSwallowID == 0 &&
       neighbor_particle.SwallowID == 0 && Mass_j > 0 && neighbor_particle.Type == 0 && neighbor_cell) {
        double r = sqrt(r2);
        if(r > 0) {
            Vec3<double> dir{dpos[0], dpos[1], dpos[2]}; dir /= r;
            double cos_t = dot(dir, active.J_dir);
            double w = sink_fb_angleweight_localcoupling_gpu(neighbor_particle, *neighbor_cell, cos_t, r, h_i);
            double mom_wt = ((double)local.Sink_angle_weighted_kernel_sum > 0)
                            ? (w / (double)local.Sink_angle_weighted_kernel_sum) : 0;
#ifdef SINK_PHOTONMOMENTUM
            double v_kick = scalars.sink_rad_momentum_factor * mom_wt * active.mom_budget / Mass_j;
            Vel_j[0] += v_kick * scalars.common.cf_atime * dir[0];
            Vel_j[1] += v_kick * scalars.common.cf_atime * dir[1];
            Vel_j[2] += v_kick * scalars.common.cf_atime * dir[2];
#endif
            (void)mom_wt;
        }
    }
#endif

    /* delta j-writes atomically. Mirrors legacy
     * sink_swallow_and_kick_functions.h:399-439. */
    if(Mass_j != Mass_j_0 ||
       Vel_j[0] != Vel_j_0[0] || Vel_j[1] != Vel_j_0[1] || Vel_j[2] != Vel_j_0[2] ||
       InternalEnergy_j != InternalEnergy_j_0) {
        if(Mass_j > 0) {
            double dmass = Mass_j - Mass_j_0;
            double current_mass = (double)Kokkos::atomic_load(&neighbor_particle.Mass);
            if(current_mass > 0) {
                Kokkos::atomic_add(&neighbor_particle.Mass, (MyDouble)dmass);
            }
        } else {
            Kokkos::atomic_store(&neighbor_particle.Mass, (MyDouble)0);
        }
        for(int k = 0; k < 3; k++) {
            double dVel = Vel_j[k] - Vel_j_0[k];
            double dpk  = Mass_j * Vel_j[k] - Mass_j_0 * Vel_j_0[k];
            Kokkos::atomic_add(&neighbor_particle.Vel[k], (MyDouble)dVel);
            Kokkos::atomic_add(&neighbor_particle.dp[k],  (MyFloat)dpk);
        }
        if(neighbor_particle.Type == 0 && neighbor_cell) {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            if(Mass_j > 0) {
                Kokkos::atomic_add(&neighbor_cell->MassTrue, (MyDouble)(Mass_j - Mass_j_0));
            } else {
                Kokkos::atomic_store(&neighbor_cell->MassTrue, (MyDouble)0);
            }
#endif
            for(int k = 0; k < 3; k++) {
                Kokkos::atomic_add(&neighbor_cell->VelPred[k], (MyDouble)(Vel_j[k] - Vel_j_0[k]));
            }
            double dIE = InternalEnergy_j - InternalEnergy_j_0;
            Kokkos::atomic_add(&neighbor_cell->InternalEnergy,     (MyDouble)dIE);
            Kokkos::atomic_add(&neighbor_cell->InternalEnergyPred, (MyDouble)dIE);
        }
    }
}

/* ============================================================================
 * SinkSwkSpec — the NeighborLoopSpec contract for sink_swallow_and_kick.
 * ========================================================================== */

struct SinkSwkSpec {
    /* ====================================================================
     * PHYSICS BLOCK
     * ==================================================================== */

    static constexpr const char *loop_name = "sink_swk";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::EpsilonAtomic; /* a cell is drained only by the one sink whose ID it carries in SwallowID, so the read-then-add on its mass has a single writer; the kicks to unswallowed neighbours are order-independent atomics */

    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int            neighbor_type_mask = (unsigned int)SINK_NEIGHBOR_BITFLAG;
    /* sink_swk pair physics: pair body uses hj = P[j].KernelRadius for SPH
     * weighting on gas-j (sink_swk_loop.h:282) and accumulates onto active i for
     * any admitted j without an internal r > h_j gate.  Symmetric reach must
     * include j-side KernelRadius for both gas and non-gas. */
    static constexpr mode_b_radius_policy_t  radius_policy      =
        MODE_B_RADIUS_GAS_KERNEL | MODE_B_RADIUS_NONGAS_KERNEL;

    static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::AllTypes;
    static constexpr bool mode_a_active_sources_in_sidx_pool = true; /* active sources (incl. non-gas) are in the AllTypes SIDX pool */

    /* Active predicate matches the legacy host caller block:
     * sink_isactive(i) && P[i].SwallowID == 0. */
    static bool is_active(int particle_index) {
        return (sink_isactive(particle_index) != 0) && (P[particle_index].SwallowID == 0);
    }

    using CallScalars   = SinkSwkCallScalars;
    using ActiveData    = SinkSwkActiveState;
    using AccumData     = SinkSwallowOut;
    using DeviceContext = SinkSwkDeviceContext;     /* device-context extension */

    /* NeighborData carries non-const pointers — sink_swk pair body does
     * many j-side atomic writes. */
    struct NeighborData {
        struct particle_data *neighbor_particle;
        struct gas_cell_data *neighbor_cell;
    };

    struct Aux {
        SinkSwallowOut          *per_active_accum;   /* [num_active] host */
        const SinkSwallowLocalIn *host_locals;       /* [num_active] host */
    };

    static constexpr bool uses_ghost_write_detector = true;
    static constexpr bool uses_ghost_writeback      = true;
    static constexpr bool writes_neighbour_motion = true;   /* a swallowing sink kicks the cells it drains */

    /* ghost_write_detector_begin/end: runner default (loop_name = "sink_swk"). */
    static void ghost_writeback_begin      (const struct neighbor_loop_args&,
                                             const struct NeighborLoopPlan&);
    static void ghost_writeback_end        (const struct neighbor_loop_args&,
                                             const struct NeighborLoopPlan&);

    /* Per-active host hooks. */
    static double      search_radius(const neighbor_loop_args& args,
                                      int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);

    /* Device-context lifecycle. */
    static void populate_device_context(const neighbor_loop_args& args, DeviceContext& ctx);
    static void cleanup_device_context (const neighbor_loop_args& args, DeviceContext& ctx);

    /* Build ActiveData[slot] on device. Reads host-staged per-active local
     * from ctx.per_active_local; computes per-active mom_budget / J_dir /
     * sink_mass_withdisk inline. */
    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& dctx,
                                   int active_slot, int /*i*/,
                                   double h_search,
                                   const CallScalars& scalars)
    {
        ActiveData active{};
        if(dctx.per_active_local) { active.local = dctx.per_active_local[active_slot]; }
        active.pos[0] = (double)active.local.Pos[0];
        active.pos[1] = (double)active.local.Pos[1];
        active.pos[2] = (double)active.local.Pos[2];
        active.h_search         = h_search;
        (void)scalars;
        active.origin_local_idx = active_slot;
        active.origin_rank      = -1;

        /* Per-active precomputation (mirrors sink_swallow_and_kick_gpu.cc:265-285). */
        active.mom_budget = 0;
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS)
        /* mom_budget is precomputed on the host and stashed in
         * SinkSwallowLocalIn.mom_budget — see sink_swk_loop.cc::sink_swk_fill_local.
         * Avoids host-only sink_lum_bol() on the device. */
        active.mom_budget = (double)active.local.mom_budget;
#endif
        active.J_dir[0] = 0; active.J_dir[1] = 0; active.J_dir[2] = 1;
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS) || defined(SINK_WIND_KICK)
        {
            Vec3<double> tmp{(double)active.local.Jgas_in_Kernel[0],
                             (double)active.local.Jgas_in_Kernel[1],
                             (double)active.local.Jgas_in_Kernel[2]};
            double n2 = tmp.norm_sq();
            if(n2 > 0) { double n = sqrt(n2); active.J_dir = tmp / n; }
        }
#endif
        active.sink_mass_withdisk = (double)active.local.Sink_Mass;
#if defined(SINK_WIND_KICK) && defined(SINK_ALPHADISK_ACCRETION)
        active.sink_mass_withdisk += (double)active.local.Sink_Mass_Reservoir;
#endif
        return active;
    }

    /* Zero accumulator — POD byte-zero, plus sentinel init for Accreted_Age
     * (ACCUM_MIN merge requires a ceiling; merge with peer takes the
     * smaller, so MAX_REAL_NUMBER means "not yet seen"). */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        for(size_t b = 0; b < sizeof(accum); b++) ((char*)&accum)[b] = 0;
#ifdef GALSF
        accum.Accreted_Age = (MyFloat)MAX_REAL_NUMBER;
#endif
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
        sink_swk_pair_kernel(active, cs, *neighbor.neighbor_particle,
                              neighbor.neighbor_cell, accum);
    }

    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int active_slot, int i,
                                        const AccumData& accum);

    /* Per-field merge of a peer rank's contribution. Per-field op MUST match
     * the pair_kernel writes; Accreted_Age uses MIN-merge with the
     * MAX_REAL_NUMBER sentinel from zero_accum. Nothing checks the two against
     * each other at runtime, so drift between them is silent.
     *
     * Adding a new accumulator field for this loop = ONE LINE under the
     * appropriate physics flag's #ifdef. */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local_accum, const AccumData& peer_accum)
    {
#define ACCUM_ADD(field)        local_accum.field += peer_accum.field;
#define ACCUM_ADD_VEC3(field)   for(int k = 0; k < 3; k++) local_accum.field[k] += peer_accum.field[k];
#define ACCUM_MIN(field)        if(peer_accum.field < local_accum.field) local_accum.field = peer_accum.field;
#define ACCUM_ADD_ARRAY(field, len)                                       \
        for(int k = 0; k < (len); k++) local_accum.field[k] += peer_accum.field[k];

        ACCUM_ADD(accreted_Mass)
        ACCUM_ADD(accreted_Sink_Mass)
        ACCUM_ADD(accreted_Sink_Mass_reservoir)
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
        ACCUM_ADD(Sink_AccretionDeficit)
#endif
#ifdef GRAIN_FLUID
        ACCUM_ADD(accreted_dust_Mass)
#endif
#ifdef RT_REINJECT_ACCRETED_PHOTONS
        ACCUM_ADD(accreted_photon_energy)
#endif
#if defined(SINK_FOLLOW_ACCRETED_MOMENTUM)
        ACCUM_ADD_VEC3(accreted_momentum)
#endif
#if defined(SINK_FOLLOW_ACCRETED_COM)
        ACCUM_ADD_VEC3(accreted_centerofmass)
#endif
#if defined(SINK_RETURN_BFLUX)
        ACCUM_ADD_VEC3(accreted_B)
#endif
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
        ACCUM_ADD_VEC3(accreted_J)
#endif
#ifdef SINK_COUNTPROGS
        ACCUM_ADD(Sink_CountProgs)
#endif
#ifdef GALSF
        ACCUM_MIN(Accreted_Age)
#endif
        ACCUM_ADD(n_gas_swallowed)
        ACCUM_ADD(n_sink_swallowed)
        ACCUM_ADD(n_star_swallowed)
        ACCUM_ADD(n_dm_swallowed)
        ACCUM_ADD_ARRAY(delta_TimeBin_Sink_mass,          TIMEBINS)
        ACCUM_ADD_ARRAY(delta_TimeBin_Sink_dynamicalmass, TIMEBINS)
        ACCUM_ADD_ARRAY(delta_TimeBin_Sink_Mdot,          TIMEBINS)
        ACCUM_ADD_ARRAY(delta_TimeBin_Sink_Medd,          TIMEBINS)

#undef ACCUM_ADD
#undef ACCUM_ADD_VEC3
#undef ACCUM_MIN
#undef ACCUM_ADD_ARRAY
    }

    /* ====================================================================
     * ENGINE APPARATUS
     * ==================================================================== */
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

};

#endif /* SINK_PARTICLES */

#endif /* SINK_SWK_LOOP_H */
