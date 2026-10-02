#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../core/wakeup_sidecar.h"
#include "../mesh/gpu_neighbor_list.h"   /* gizmo_motion_bound_raise */
#include <vector>
#include "../mesh/kernel.h"
#include "mechanical_fb_types.h"  /* provides struct MechFBGasDelta */
#include "galsf_gpu_decls.h"      /* mechfb_{alloc,free,run_iterative,fill_call_scalars} prototypes */

/* Routines for mechanical feedback/enrichment models: stellar winds, supernovae, etc
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef GALSF_FB_MECHANICAL

int addFB_evaluate_active_check(int i, int fb_loop_iteration);
int addFB_evaluate_active_check(int i, int fb_loop_iteration)
{
    if(P[i].Type <= 1) {return 0;} // note quantities used here must -not- change in the loop [hence not using mass here], b/c can change offsets for return from different processors, giving a negative mass and undefined behaviors
    if(P[i].KernelRadius <= 0) {return 0;}
    if(P[i].NumNgb <= 0) {return 0;}
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
    if(P[i].Type == 5 && !P[i].do_gas_search_this_timestep) {return 0;}
#endif
    if(P[i].SNe_ThisTimeStep>0) {if(fb_loop_iteration<0 || fb_loop_iteration==0) {return 1;}}
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
    if(P[i].MassReturn_ThisTimeStep>0) {if(fb_loop_iteration<0 || fb_loop_iteration==1) {return 1;}}
#ifdef GALSF_FB_FIRE_RPROCESS
    if(P[i].RProcessEvent_ThisTimeStep>0) {if(fb_loop_iteration<0 || fb_loop_iteration==2) {return 1;}}
#endif
#ifdef GALSF_FB_FIRE_AGE_TRACERS
    if(P[i].AgeDeposition_ThisTimeStep>0) {if(fb_loop_iteration<0 || fb_loop_iteration==3) {return 1;}}
#endif
#endif
#if defined(SINGLE_STAR_FB_WINDS) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    if(P[i].wind_mode != 2 || P[i].ProtoStellarStage != 5) {return 0;}
#endif
    return 0;
}


void determine_where_SNe_occur(void)
{
    if(All.Time<=0) return;
    int i; double dt,star_age,npossible,nhosttotal,ntotal,ptotal,dtmean,rmean;
    npossible=nhosttotal=ntotal=ptotal=dtmean=rmean=0;
    double mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean;
    mpi_npossible=mpi_nhosttotal=mpi_ntotal=mpi_ptotal=mpi_dtmean=mpi_rmean=0;
    // loop over particles //
    for (int i : ActiveParticleList)
    {
        P[i].SNe_ThisTimeStep=0;
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
        P[i].MassReturn_ThisTimeStep=0;
#ifdef GALSF_FB_FIRE_RPROCESS
        P[i].RProcessEvent_ThisTimeStep=0;
#endif
#ifdef GALSF_FB_FIRE_AGE_TRACERS
        P[i].AgeDeposition_ThisTimeStep=0;
#endif
#endif
    

#if defined(SINGLE_STAR_SINK_DYNAMICS)
        if(P[i].Type == 0) {continue;} // any non-gas type is eligible to be a 'star' here
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
        if(P[i].Type == 5)
        {   // single sink type, for type 5
            if((P[i].ProtoStellarStage != 5) && (P[i].ProtoStellarStage != 6)) {continue;} // only MS or SN-tagged particles eligible to have winds or SN in this routine
        } else { // stellar population type allowed here
            if(!is_galsf_stellar_candidate_type(P[i].Type, All.ComovingIntegrationOn)) {continue;} // stellar-candidate type (see proto.h)
        }
#endif
#else
        if(!is_galsf_stellar_candidate_type(P[i].Type, All.ComovingIntegrationOn)) {continue;} // stellar-candidate type (see proto.h)
#endif
        if(P[i].Mass<=0) {continue;}
        dt = get_particle_feedback_timestep_in_physical(i, P);
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
        if(P[i].Type == 5) {dt = P[i].dt_since_last_gas_search;}
#endif
        if(dt<=0) {continue;} // no time, no events
        star_age = evaluate_stellar_age_Gyr(i);
        if(star_age<=0) {continue;} // unphysical age, no events
        // now use a calculation of mechanical event rates to determine where/when the events actually occur //
        npossible++;
        double RSNe = mechanical_fb_calculate_eventrates(i,dt);
        rmean += RSNe; ptotal += RSNe * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR);
#ifdef GALSF_SFR_IMF_SAMPLING
        if(P[i].Type<5) {if(P[i].IMF_NumMassiveStars>0) {P[i].IMF_NumMassiveStars=DMAX(0,P[i].IMF_NumMassiveStars-P[i].SNe_ThisTimeStep);}} // lose an O-star for every SNe //
#endif
        if(P[i].SNe_ThisTimeStep>0) {ntotal+=P[i].SNe_ThisTimeStep; nhosttotal++;}
        dtmean += dt;
    } // for (int i : ActiveParticleList) //

    /* The six log-file diagnostics are independent sums to rank 0, so they ride
     * one six-element reduction rather than six single-value ones: same values,
     * one collective latency instead of six. This routine runs on every step
     * that has any feedback-eligible particle, so the five saved collectives
     * are five fewer places for rank skew to be absorbed. */
    const double loc_fb_diag[6] = {dtmean, rmean, ptotal, nhosttotal, ntotal, npossible};
    double glob_fb_diag[6] = {0,0,0,0,0,0};
    MPI_Reduce(loc_fb_diag, glob_fb_diag, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mpi_dtmean = glob_fb_diag[0]; mpi_rmean = glob_fb_diag[1]; mpi_ptotal = glob_fb_diag[2];
    mpi_nhosttotal = glob_fb_diag[3]; mpi_ntotal = glob_fb_diag[4]; mpi_npossible = glob_fb_diag[5];

    if(ThisTask == 0)
    {
        if(mpi_ntotal > 0 && mpi_nhosttotal > 0 && mpi_dtmean > 0)
        if(mpi_npossible>0)
        {
            mpi_dtmean /= mpi_npossible; mpi_rmean /= mpi_npossible;
            fprintf(FdSNeFBLogFile, "%.16g %g %g %g %g %g %g \n", All.Time,mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean); fflush(FdSNeFBLogFile);
        }
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) {fflush(FdSNeFBLogFile);}
    } // if(ThisTask == 0) //

} // void determine_where_SNe_occur() //



/* this is a temporary structure for quantities used ONLY in the loops below,
 for example for ensuring conservation if there are many overlapping events.
 Struct type MechFBGasDelta is defined in mechanical_fb_functions.h so the
 GPU port can share the exact layout. */
static struct MechFBGasDelta *LocalGasMechFBInfoTemp;






/* subroutine to check for total kinetic energy change and thermal energy change after integrating the effects of all SNe over all cells, and coupling this to particles,  */
void verify_and_assign_local_mechfb_integrals(void)
{
    /* The cells this pass kicks: their motion bounds are raised once at the
     * end, so the neighbour indexes stay valid for the velocities they now have. */
    std::vector<int> kicked;
    int j,k; for(j=0;j<N_gas;j++)
    {
        if(LocalGasMechFBInfoTemp[j].N_injected <= 0) {continue;} /* all mechanisms deposit non-zero mass, so skip if this is not >0*/
        if(P[j].Type==0 && P[j].Mass>0)
        {
            kicked.push_back(j);
            double m0=P[j].Mass, dm=LocalGasMechFBInfoTemp[j].m_injected; P[j].Mass += dm; /* update mass */ if(P[j].Type==0) {CellP[j].Mass = P[j].Mass;}
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            m0=CellP[j].MassTrue; CellP[j].MassTrue += dm; /* update conserved mass */
#endif
            double mf=m0+dm; /* save for below */
            for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k] = (m0/mf)*P[j].Metallicity[k] + (1./mf)*LocalGasMechFBInfoTemp[j].Z_injected[k];} /* update metallicity */
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            update_ISMDustChem_after_mechanical_injection(j, LocalGasMechFBInfoTemp[j].Mass_Where_Dust_Shocked, m0, mf, LocalGasMechFBInfoTemp[j].Z_injected); /* update dust chemistry quantities */
#endif
            CellP[j].Density *= mf/m0; /* update density [semi-drift] */
            double dTE=LocalGasMechFBInfoTemp[j].TE_injected;
            if(dTE > 0)
            {
                double TE_0=m0*CellP[j].InternalEnergy; dTE=DMAX(-TE_0,dTE); /* ensure against non-negative values */
                double dU = (-dm/mf)*CellP[j].InternalEnergy + (1./mf)*dTE; /* using new mass get updated internal energy */
                double dt = get_particle_timestep_in_physical(j, P), implied_heating_cgs=(dU*UNIT_SPECEGY_IN_CGS*PROTONMASS_CGS)/(dt*UNIT_TIME_IN_CGS), typical_cooling_cgs=1.e-23*(CellP[j].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS);
                if((implied_heating_cgs < 0.3*typical_cooling_cgs) && (dt > MIN_REAL_NUMBER) && ((dU < 4.*CellP[j].InternalEnergy) || ((dU < 1000.*CellP[j].InternalEnergy) && ((dU+CellP[j].InternalEnergy)*U_TO_TEMP_UNITS*2./3.*1.28 < 5.e5)))) {CellP[j].DtInternalEnergy += dU/dt;} else {
                    CellP[j].InternalEnergy += dU; CellP[j].InternalEnergyPred += dU;
                    /* Direct-dU branch: mark this receiver for positive wakeup
                     * from the most-aggressive source TimeBin that contributed.
                     * Without this, the receiver's stale post-injection state
                     * (huge u, low MaxSignalVel) is only "discovered" later by a
                     * pair-body wakeup-check, which then propagates shell-by-shell
                     * with waker_bin - offset per generation, cascading bins down.
                     * Coupled to the process_wake_ups floor at lowest_occupied_
                     * active_bin, this gets the receiver active at the right pace
                     * immediately, killing the cascade at its root. Standard
                     * hydro-convention encoding via max_source_wakeup accumulator
                     * -- no special encoding, no EOS/MaxSignalVel refresh (which
                     * would be silent suppression of legitimate wakeup). Deferred-dU
                     * branch intentionally NOT marked: u unchanged there, no
                     * perturbation. */
                    int wakeup_val = LocalGasMechFBInfoTemp[j].max_source_wakeup;
                    if(wakeup_val > 0 && wakeup_val > P[j].wakeup) {
                        P[j].wakeup = (short int)wakeup_val;
                        wakeup_sidecar_mark(j);
                        NeedToWakeupParticles_local = 1;
                    }
                }
                //CellP[j].InternalEnergy += dU; CellP[j].InternalEnergyPred += dU; /* update internal energy; simpler (old) way to do it - less accurate phase diagrams at high density, however */
            }
            double dKE=LocalGasMechFBInfoTemp[j].KE_injected, dp[3];
            if(dKE != 0 || LocalGasMechFBInfoTemp[j].p_injected[0] != 0 || LocalGasMechFBInfoTemp[j].p_injected[1] != 0 || LocalGasMechFBInfoTemp[j].p_injected[2] != 0 )
            {
                double KE_0=0, p0[3], p2=0, dp2=0, pdp=0; /* define variables for momentum update (the non-trivial update term) */
                for(k=0;k<3;k++) {p0[k]=m0*P[j].Vel[k]/All.cf_atime; dp[k]=LocalGasMechFBInfoTemp[j].p_injected[k]; p2+=p0[k]*p0[k]; dp2+=dp[k]*dp[k]; pdp+=p0[k]*dp[k];} /* define variables */
                KE_0=p2/(2.*m0); dKE=DMAX(-KE_0,dKE); /* limit to retain positive definite kinetic energy */
                double a0=dp2, b0=2.*pdp, c0=p2*dm/m0+2.*mf*dKE, sfac=b0*b0+4.*a0*c0, f0=0; /* assume coupled delta_p = f0*delta_p, and solve for the multiplier f0 that gives the desired total kinetic energy change */
                if(sfac<=0 || !isfinite(sfac)) {if((fabs(a0)<1.e-40) || !isfinite(b0/a0)) {f0=0;} else {f0=-b0/(2.*a0);}} // catches for floating-point error
                    else {if((b0>0) && fabs(4.*a0*c0)<1.e-3*b0*b0) {f0=c0/b0;} // catches for floating-point error
                    else {if(fabs(a0)<1.e-40 || !isfinite(a0)) {f0=0;} else {f0=(-b0+sqrt(b0*b0+4.*a0*c0))/(2.*a0);}}} // catches for floating-point error, if pass all of them, use exact solution for desired KE
                if(f0<0 || !isfinite(f0)) {f0=0;} else {if(f0>1.) {f0=1.;}} /* limit to physical values (should never be an issue but again because of float error it could be) */
                for(k=0;k<3;k++) {
                    double dv = (-dm/mf)*P[j].Vel[k] + f0*(1./mf)*dp[k]*All.cf_atime; /* calculate total momentum change and mass change and therefore final velocity (in code units) */
                    P[j].Vel[k] += dv; CellP[j].VelPred[k] += dv; P[j].dp[k] += f0*dp[k]; /* update velocities */
                }
            }
#if defined(COSMIC_RAY_FLUID)
            /* Apply CR deltas accumulated by the batched GPU kernel. On CPU-only
             * builds, the existing inject_cosmic_rays in the pair loop writes
             * directly to CellP, so these delta fields remain zero and this block
             * is a no-op. */
            {
                double cr_tot = 0;
                for(int kb = 0; kb < N_CR_PARTICLE_BINS; kb++) cr_tot += LocalGasMechFBInfoTemp[j].CR_energy_injected[kb];
                if(cr_tot > 0) {
                    /* compute unit injection direction from weighted sum */
                    double dmag_sq = 0;
                    for(int c = 0; c < 3; c++) dmag_sq += LocalGasMechFBInfoTemp[j].CR_dir_weighted[c] * LocalGasMechFBInfoTemp[j].CR_dir_weighted[c];
                    double dmag = sqrt(dmag_sq);
                    double dir_hat[3];
                    if(dmag > 0) { for(int c = 0; c < 3; c++) dir_hat[c] = LocalGasMechFBInfoTemp[j].CR_dir_weighted[c] / dmag; }
                    else         { dir_hat[0] = 0; dir_hat[1] = 0; dir_hat[2] = 1; }
#ifdef MAGNETIC
                    /* project onto B, keeping sign from requested direction */
                    {
                        double B_dot_d = 0;
                        for(int c = 0; c < 3; c++) B_dot_d += dir_hat[c] * CellP[j].BPred[c];
                        for(int c = 0; c < 3; c++) dir_hat[c] = B_dot_d * CellP[j].BPred[c];
                        double bmag_sq = 0;
                        for(int c = 0; c < 3; c++) bmag_sq += dir_hat[c] * dir_hat[c];
                        double bmag = sqrt(bmag_sq);
                        if(bmag > 0) { for(int c = 0; c < 3; c++) dir_hat[c] /= bmag; }
                        else         { dir_hat[0] = 0; dir_hat[1] = 0; dir_hat[2] = 1; }
                    }
#endif
                    /* apply per-bin energy / number / flux */
                    for(int kb = 0; kb < N_CR_PARTICLE_BINS; kb++) {
                        double dEcr = LocalGasMechFBInfoTemp[j].CR_energy_injected[kb];
                        if(dEcr <= 0) continue;
                        CellP[j].CosmicRayEnergy[kb]     += dEcr;
                        CellP[j].CosmicRayEnergyPred[kb] += dEcr;
#if defined(CRFLUID_EVOLVE_SPECTRUM)
                        CellP[j].CosmicRay_Number_in_Bin[kb] += LocalGasMechFBInfoTemp[j].CR_number_injected[kb];
#endif
                        double flux_mag = dEcr * CRFLUID_REDUCED_C_CODE(kb);
                        for(int c = 0; c < 3; c++) {
                            double dflux = flux_mag * dir_hat[c];
                            CellP[j].CosmicRayFlux[kb][c]     += dflux;
                            CellP[j].CosmicRayFluxPred[kb][c] += dflux;
                        }
                    }
                }
            }
#endif /* COSMIC_RAY_FLUID */
            /* apply to every home gas cell with pending deltas (N_injected>0):
               this buffer is the single source of truth for received feedback,
               populated by all paths (single/multi-rank, both dispatch modes)
               before this routine runs, so no source-side coupling count is
               consulted -- that count was rank-attributed to sources, not
               receivers, and truncated cross-rank deposits. */
        }
        /* drain-zero: re-zero this touched cell (SSOT) so the persistent buffer
         * stays all-zero between steps. Runs for every N_injected>0 cell (incl.
         * Type!=0/Mass<=0, drained-but-not-applied) so no dirty cell survives. */
        mechfb_reset_one_gas_delta(LocalGasMechFBInfoTemp, j);
    }
    if(!kicked.empty()) {gizmo_motion_bound_raise(kicked.data(), (int)kicked.size());}
    return;
}



/* parent routine which calls all of the mechanical fb loops and verifies coupling, if needed */
void mechanical_fb_calc_toplevel(void)
{
    PRINT_STATUS("Start mechanical feedback computation...");

    /* GPU port: dispatch all 6 modes at once on the GPU, sharing a single
       neighbor list across modes. See galaxy_sf/mechfb_loop.{h,cc} (runner-template
       port); the legacy galaxy_sf/mechanical_fb_gpu.cc evaluator has been deleted. */
    {
        /* Build SUPERSET active list: any star active in ANY mode is included.
           addFB_evaluate_active_check(ii, -2) with a negative iteration acts as
           "check all modes" — the fb_loop_iteration<0 branch fires for every
           event type. Canonical single-predicate form; avoids hard-coding the
           {-2..3} range in multiple call sites (MechFBSpec::is_active is the
           runner-side mirror of this same predicate). */
        int num_active = 0;
        for(int ii : ActiveParticleList) {
            if(addFB_evaluate_active_check(ii, -2)) num_active++;
        }
        /* Multi-rank correctness: ghost_prep is collective so all ranks must agree
         * on whether to call it. MPI_Allreduce first; if no rank has any source,
         * every rank skips the mech_fb path entirely (no allocation, no kernel). */
        int global_num_active = num_active;
        MPI_Allreduce(&num_active, &global_num_active, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if(global_num_active == 0) { return; }

        /* persistent grow-only buffer: no per-step alloc or O(N_gas) zero. The
         * buffer is all-zero at entry by the drain-zero invariant below.
         * Kept AFTER the active-count reduce above so a globally quiet step
         * never allocates it. If any rank could not obtain it, every rank must
         * skip together: peers deposit into this buffer through the feedback
         * exchange, so one rank cannot bow out alone. The failing rank has
         * already asked for the stop, naming the buffer. */
        LocalGasMechFBInfoTemp = mechfb_get_persistent_gas_delta(N_gas);
        int deposit_buffer_missing = (LocalGasMechFBInfoTemp == NULL) ? 1 : 0, any_missing = 0;
        MPI_Allreduce(&deposit_buffer_missing, &any_missing, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(any_missing) { return; }
        int *nl_active = (int *) mymalloc("mechfb_nl_active",
            (num_active > 0 ? num_active : 1) * sizeof(int));
        {int aa = 0; for(int ii : ActiveParticleList) {
            if(addFB_evaluate_active_check(ii, -2)) nl_active[aa++] = ii;
        }}
        mechfb_run_iterative(nl_active, num_active,
                              LocalGasMechFBInfoTemp, N_gas);
        myfree(nl_active);
    }



    verify_and_assign_local_mechfb_integrals();
    /* persistent buffer: NOT freed -- reused across steps, kept all-zero by the
     * drain-zero in verify_and_assign_local_mechfb_integrals. */
}




#endif /* GALSF_FB_MECHANICAL */
