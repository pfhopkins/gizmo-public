/*! \file sink_swallow_and_kick.c
*  \brief routines for gas accretion onto sink particles, and sink particle mergers
*/
/* Stdlib + Kokkos must precede any project header (allvars.h macros may
 * conflict with stdlib names). */
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"  /* MUST precede allvars.h: installs device-pass `#define All AllDeviceMirror` redirect before cell_data.h is parsed */
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/gpu_neighbor_list.h" /* gizmo_mark_kernel_radius_dirty_* */
#include "../core/wakeup_sidecar.h"
#include "sink_functions.h"
#include "sink_swk_loop.h"
/*
* This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
* see notes in sink.c for details on code history.
*/


#ifdef SINK_PARTICLES // top-level flag [needs to be here to prevent compiler breaking when this is not active] //


static int N_gas_swallowed, N_star_swallowed, N_dm_swallowed, N_sink_swallowed;

#ifdef SINK_ALPHADISK_ACCRETION
#define out_accreted_Sink_Mass_alphaornot out.accreted_Sink_Mass_reservoir
#else
#define out_accreted_Sink_Mass_alphaornot out.accreted_Sink_Mass
#endif





/* Host-side fill of SinkSwallowLocalIn for source particle i. Mirrors
 * legacy sink_swallow_local_fill from sinks/sink_swallow_and_kick_gpu.cc:47-87,
 * with one addition: precomputed mom_budget under SINK_CALC_LOCAL_ANGLEWEIGHTS
 * (sink_lum_bol is host-only). */
static void sink_swk_fill_local(int i, struct SinkSwallowLocalIn *loc)
{
    int j_tempinfo = P[i].IndexMapToTempStruc;
    loc->Pos             = P[i].Pos;
    loc->Vel             = P[i].Vel;
    loc->KernelRadius    = P[i].KernelRadius;
    loc->Mass            = P[i].Mass;
    loc->Sink_Mass       = P[i].Sink_Mass;
    loc->ID              = P[i].ID;
    loc->claim_token     = gizmo_sink_claim_token(ThisTask, i);
    loc->ID_child_number = P[i].ID_child_number;
    loc->ID_generation   = P[i].ID_generation;
    loc->Mdot            = P[i].Sink_Mdot;
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS) || defined(SINK_WIND_KICK)
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
    loc->Jgas_in_Kernel = P[i].Sink_Specific_AngMom;
#else
    loc->Jgas_in_Kernel = SinkTempInfo[j_tempinfo].Jgas_in_Kernel;
#endif
#endif
#ifdef SINK_ALPHADISK_ACCRETION
    loc->Sink_Mass_Reservoir = P[i].Sink_Mass_Reservoir;
#endif
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS)
    loc->Sink_angle_weighted_kernel_sum = SinkTempInfo[j_tempinfo].Sink_angle_weighted_kernel_sum;
#endif
    loc->Dt = (MyFloat)get_particle_feedback_timestep_in_physical(i, P);
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
    loc->Dt = P[i].dt_since_last_gas_search;
#endif
#if defined(SINK_RETURN_ANGMOM_TO_GAS)
    loc->Sink_Specific_AngMom = P[i].Sink_Specific_AngMom;
    loc->angmom_norm_topass_in_swallowloop = SinkTempInfo[j_tempinfo].angmom_norm_topass_in_swallowloop;
#endif
#if defined(SINK_RETURN_BFLUX)
    loc->B = P[i].B;
    loc->kernel_norm_topass_in_swallowloop = SinkTempInfo[j_tempinfo].kernel_norm_topass_in_swallowloop;
#endif
#ifdef SINGLE_STAR_FB_LOCAL_RP
    loc->Luminosity = (MyFloat)sink_lum_bol(loc->Mdot, loc->Sink_Mass, i);
#endif
#if defined(SINK_CALC_LOCAL_ANGLEWEIGHTS)
#if defined(SINGLE_STAR_FB_LOCAL_RP)
    loc->mom_budget = (MyFloat)((double)loc->Luminosity * (double)loc->Dt / C_LIGHT_CODE);
#else
    /* The particle index is required, not optional: with a protostellar-evolution
     * model the luminosity is a cached per-star quantity read from the particle,
     * and the accretion rate and mass arguments are not used at all. */
    loc->mom_budget = (MyFloat)(sink_lum_bol((double)loc->Mdot, (double)loc->Sink_Mass, i)
                                 * (double)loc->Dt / C_LIGHT_CODE);
#endif
#endif
}


/* Caller-side scatter from per-active accum into SinkTempInfo + globals.
 * Mirrors legacy sink_swallow_apply_out (sink_swallow_and_kick_gpu.cc:92-132)
 * including the StellarAge MIN-replace from out.Accreted_Age. */
static void sink_swk_scatter(const int *active_list, int num_active,
                              const struct SinkSwallowOut *per_active_accum,
                              int *N_gas_sw, int *N_sink_sw,
                              int *N_star_sw, int *N_dm_sw)
{
    *N_gas_sw = *N_sink_sw = *N_star_sw = *N_dm_sw = 0;
    for(int a = 0; a < num_active; a++) {
        int i = active_list[a];
        const struct SinkSwallowOut& out = per_active_accum[a];
        int t = P[i].IndexMapToTempStruc;

#define SCATTER_ADD(dst, src)        dst += out.src;
#define SCATTER_ADD_VEC3(dst, src)   for(int k = 0; k < 3; k++) dst[k] += out.src[k];
#define SCATTER_MIN(dst, src)        if(out.src < (dst)) (dst) = out.src;

        SCATTER_ADD(SinkTempInfo[t].accreted_Mass,                   accreted_Mass)
        SCATTER_ADD(SinkTempInfo[t].accreted_Sink_Mass,              accreted_Sink_Mass)
        SCATTER_ADD(SinkTempInfo[t].accreted_Sink_Mass_reservoir,    accreted_Sink_Mass_reservoir)
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
        SCATTER_ADD(SinkTempInfo[t].Sink_AccretionDeficit,           Sink_AccretionDeficit)
#endif
#ifdef GRAIN_FLUID
        SCATTER_ADD(SinkTempInfo[t].accreted_dust_Mass,              accreted_dust_Mass)
#endif
#ifdef RT_REINJECT_ACCRETED_PHOTONS
        SCATTER_ADD(SinkTempInfo[t].accreted_photon_energy,          accreted_photon_energy)
#endif
#if defined(SINK_FOLLOW_ACCRETED_MOMENTUM)
        SCATTER_ADD_VEC3(SinkTempInfo[t].accreted_momentum,          accreted_momentum)
#endif
#if defined(SINK_FOLLOW_ACCRETED_COM)
        SCATTER_ADD_VEC3(SinkTempInfo[t].accreted_centerofmass,      accreted_centerofmass)
#endif
#if defined(SINK_RETURN_BFLUX)
        SCATTER_ADD_VEC3(SinkTempInfo[t].accreted_B,                 accreted_B)
#endif
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
        SCATTER_ADD_VEC3(SinkTempInfo[t].accreted_J,                 accreted_J)
#endif
#ifdef SINK_COUNTPROGS
        P[i].Sink_CountProgs += out.Sink_CountProgs;
#endif
#ifdef GALSF
        SCATTER_MIN(P[i].StellarAge, Accreted_Age)
#endif
        for(int b = 0; b < TIMEBINS; b++) {
            TimeBin_Sink_mass[b]          += out.delta_TimeBin_Sink_mass[b];
            TimeBin_Sink_dynamicalmass[b] += out.delta_TimeBin_Sink_dynamicalmass[b];
            TimeBin_Sink_Mdot[b]          += out.delta_TimeBin_Sink_Mdot[b];
            TimeBin_Sink_Medd[b]          += out.delta_TimeBin_Sink_Medd[b];
        }

        *N_gas_sw  += out.n_gas_swallowed;
        *N_sink_sw += out.n_sink_swallowed;
        *N_star_sw += out.n_star_swallowed;
        *N_dm_sw   += out.n_dm_swallowed;

#undef SCATTER_ADD
#undef SCATTER_ADD_VEC3
#undef SCATTER_MIN
        (void)t;
    }
}




void sink_swallow_and_kick_loop(void)
{
    CPU_Step[CPU_SINK_FEEDSWK] += measure_time();
    (void)N_gas_swallowed; (void)N_star_swallowed; (void)N_dm_swallowed; (void)N_sink_swallowed;
    /* D1 GPU path — atomic j-writes + ghost writeback; per-source output scatter
     * and MPI_Reduce of swallow counters handled inside the launcher. */
    /* engine: build active list */
    int *active_list = nullptr;
    int  num_active = 0, num_global_active = 0;
    if(!nlr_build_active_list(SinkSwkSpec::is_active,
                               &active_list, &num_active, &num_global_active,
                               "sinkswk_active_list")) {
        CPU_Step[CPU_SINK_FEEDSWK] += measure_time();
        return;
    }

    /* physics: per-active accumulator + per-active host-fill */
    int alloc_n = (num_active > 0) ? num_active : 1;
    struct SinkSwallowOut *per_active_accum = (struct SinkSwallowOut *)
        mymalloc("sinkswk_per_active_accum", alloc_n * sizeof(struct SinkSwallowOut));
    struct SinkSwallowLocalIn *host_locals = (struct SinkSwallowLocalIn *)
        mymalloc("sinkswk_host_locals", alloc_n * sizeof(struct SinkSwallowLocalIn));

    for(int a = 0; a < num_active; a++) {
        sink_swk_fill_local(active_list[a], &host_locals[a]);
    }

    SinkSwkSpec::Aux aux;
    aux.per_active_accum = per_active_accum;
    aux.host_locals      = host_locals;

    /* engine: hand off to runner */
    neighbor_loop_args args = nlr_default_args();
    args.active_list = active_list;
    args.num_active  = num_active;
    args.aux         = &aux;
    run_neighbor_loop<SinkSwkSpec>(args);

    /* physics: caller scatter + MPI_Reduce of swallow counters */
    int N_gas_sw = 0, N_sink_sw = 0, N_star_sw = 0, N_dm_sw = 0;
    sink_swk_scatter(active_list, num_active, per_active_accum,
                     &N_gas_sw, &N_sink_sw, &N_star_sw, &N_dm_sw);

    int Ntot_gas = 0, Ntot_sink = 0, Ntot_star = 0, Ntot_dm = 0;
    MPI_Reduce(&N_gas_sw,  &Ntot_gas,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_sink_sw, &Ntot_sink, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_star_sw, &Ntot_star, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_dm_sw,   &Ntot_dm,   1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0 && (Ntot_gas || Ntot_sink || Ntot_star || Ntot_dm)) {
        printf("Accretion done: swallowed %d gas, %d star, %d dm, and %d sink particles\n",
               Ntot_gas, Ntot_star, Ntot_dm, Ntot_sink);
        fflush(stdout);
    }

    /* engine: free + return */
    myfree(host_locals);
    myfree(per_active_accum);
    nlr_free_active_list(active_list);

    CPU_Step[CPU_SINK_FEEDSWK] += measure_time();
}





#ifdef SINK_WIND_SPAWN
void spawn_sink_wind_feedback(void)
{
    if(ghost_require_no_live_pool_for_layout_change("spawn_sink_wind_feedback")) {return;}
    int i, n_particles_split = 0, MPI_n_particles_split, dummy_gas_tag=0;
    for(i = 0; i < NumPart; i++)
        if(P[i].Type==0)
        {
            dummy_gas_tag=i;
            break;
        }

    /* don't loop or go forward if there are no gas particles in the domain, or the code will crash */
    for (int i : ActiveParticleList)
    {
        long nmax = (int)(0.99*All.MaxPart); if(All.MaxPart-20 < nmax) nmax=All.MaxPart-20; int ptype_can_spawn = 0; if(P[i].Type == 5) {ptype_can_spawn = 1;}
#ifdef SNE_NONSINK_SPAWN
        if(P[i].Type == 4) {ptype_can_spawn = 1;}
#endif
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4)
        if(is_particle_a_special_zoom_target(i)) {ptype_can_spawn = 1;}
#endif
        if((NumPart+n_particles_split+(int)(2.*(SINK_WIND_SPAWN+0.1)) < nmax) && (ptype_can_spawn==1)) // basic condition: particle is a 'spawner' (sink), and code can handle the event safely without crashing.
        {
            int sink_eligible_to_spawn = 0; // flag to check eligibility for spawning
            if(*active_unspawned_mass_ptr(i) >= (SINK_WIND_SPAWN)*target_mass_for_wind_spawning(i)) {sink_eligible_to_spawn=1;} // have 'enough' mass to spawn, in whichever reservoir currently holds the discrete-spawn channel
#if defined(SINGLE_STAR_SINK_DYNAMICS)
            if(P[i].Type==5) {if((P[i].Mass <= 3.5*P[i].Sink_Formation_Mass) || (P[i].Sink_Mass*UNIT_MASS_IN_SOLAR < 0.01)) {sink_eligible_to_spawn=0;}}  // spawning causes problems in these modules for low-mass sinks, so arbitrarily restrict to this, since it's roughly a criterion on the minimum particle mass. and for <0.01 Msun, in pre-collapse phase, no jets
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
            if(P[i].Type==5) {if(P[i].ProtoStellarStage == 6) {sink_eligible_to_spawn=1;}} // spawn the SNe ejecta no matter what the sink or 'unspawned' mass flag actually is
#endif
#endif
            if(sink_eligible_to_spawn)
            {
                int j; dummy_gas_tag=-1; double r2=MAX_REAL_NUMBER;
                for(j=0; j<N_gas; j++) /* find the closest gas particle on the domain to act as the dummy */
                {
                    if(P[j].Type==0)
                    {
                        if((P[j].Mass>0) && (CellP[j].Density>0) && (CellP[j].recent_refinement_flag==0))
                        {
                            double dx2=(P[j].Pos[0]-P[i].Pos[0])*(P[j].Pos[0]-P[i].Pos[0]) + (P[j].Pos[1]-P[i].Pos[1])*(P[j].Pos[1]-P[i].Pos[1]) + (P[j].Pos[2]-P[i].Pos[2])*(P[j].Pos[2]-P[i].Pos[2]);
                            if(dx2 < r2) {r2=dx2; dummy_gas_tag=j;}
                        }
                    }
                }
                if(dummy_gas_tag >= 0)
                {
                    n_particles_split += sink_spawn_particle_wind_shell( i , dummy_gas_tag, n_particles_split);
                }
            }
        }
    }
    MPI_Allreduce(&n_particles_split, &MPI_n_particles_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(MPI_n_particles_split>0) {
#ifdef MAINTAIN_TREE_IN_REARRANGE
        All.NumForcesSinceLastTreeBuild +=  0.0001 * All.TreeRebuild_ActiveFraction * All.TotNumPart; // we can insert spawned particles in the tree, but still a good idea to rebuild the tree every now and then, so we bring the next treebuild a bit sooner; additional cost should be small
#else
        TreeReconstructFlag = 1; // otherwise just wipe and rebuild the tree next chance you get - more expensive but more accurate
#endif
        if(ThisTask == 0) {printf(" ..Sink-Spawn Event: %d particles spawned \n", MPI_n_particles_split);}
    }

    /* rearrange_particle_sequence -must- be called immediately after this routine! */
    All.TotNumPart += (long long)MPI_n_particles_split;
    All.TotN_gas   += (long long)MPI_n_particles_split;
    Gas_split       = n_particles_split;                    // specific to the local processor //
}


void get_random_orthonormal_basis(int seed, Vec3<double>& nx, Vec3<double>& ny, Vec3<double>& nz)
{
    double phi, cos_theta, sin_theta, sin_phi, cos_phi;
    phi=2.*M_PI*get_random_number(seed+1+ThisTask), cos_theta=2.*(get_random_number(seed+3+2*ThisTask)-0.5); sin_theta=sqrt(1-cos_theta*cos_theta), sin_phi=sin(phi), cos_phi=cos(phi);
    /* velocities (determined by wind velocity direction) */
    nz = {sin_theta*cos_phi, sin_theta*sin_phi, cos_theta}; // random z axis

    double norm=0;
    while(norm==0){ // necessary in case ny is parallel to nz - believe it or not this happened once!
        phi=2.*M_PI*get_random_number(seed+4+ThisTask), cos_theta=2.*(get_random_number(seed+5+2*ThisTask)-0.5); sin_theta=sqrt(1-cos_theta*cos_theta), sin_phi=sin(phi), cos_phi=cos(phi);
        ny = {sin_theta*cos_phi, sin_theta*sin_phi, cos_theta}; // random y axis, needs to have its z component deprojected
        // do Gram-Schmidt to get an orthonormal basis
        ny -= nz * dot(ny, nz); // deproject component along z
        norm = ny.norm_sq();
        if(norm==0) continue;
        ny *= 1./sqrt(norm);
    }
    nx = cross(ny, nz);
    return;
}

/* Convenience function to compute the direction to launch a wind particle                                      */
/*                                                                                                              */
/* i - index of particle doing the spawning                                                                     */
/* num_spawned_this_call - how many we have already spawned in this call of sink_spawn_particle_wind_shell */
/* mode - 0 for random, 1 for collimated, 2 for isotropized random, 3 for angular grid                          */
/* ny, nz - shape (3,) arrays containing 2 vectors in the fixed orthonormal basis - for collimated winds, nz    */
/*          is the axis                                                                                         */
/* dir - shape (3,) array containing the direction - pass as an input to remember the previous direction        */

void get_wind_spawn_direction(int i, int num_spawned_this_call, int mode, Vec3<double>& ny, Vec3<double>& nz, Vec3<double>& veldir, Vec3<double>& dpdir)
{
    int k;
    if((mode != 3) && (num_spawned_this_call % 2)) { // every second particle is spawned in the opposite direction to the last, conserving momentum and COM
        veldir = -veldir; dpdir = -dpdir;
        return; // we're done
    }
    Vec3<double> nx = cross(ny, nz);
    // now do the actual direction based on the mode we're in
    double phi, cos_theta, sin_theta, sin_phi, cos_phi;
    if(mode==0){ // fully random
        phi=2.*M_PI*get_random_number(num_spawned_this_call+1+ThisTask), cos_theta=2.*(get_random_number(num_spawned_this_call+3+2*ThisTask)-0.5); sin_theta=sqrt(1-cos_theta*cos_theta), sin_phi=sin(phi), cos_phi=cos(phi);
        veldir = {sin_theta*cos_phi, sin_theta*sin_phi, cos_theta}; dpdir = veldir;
    } else if (mode==1){ // collimated according to a conical velocity field
        double theta0=0.01, thetamax=30.*(M_PI/180.); // "flattening parameter" and max opening angle of jet velocity distribution from Matzner & McKee 1999, sets the collimation of the jets
#if !defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
        theta0=1.e-4; thetamax=10.*(M_PI/180.); // narrower opening angle distribution for agn jets
#endif
        double theta=atan(theta0*tan(get_random_number(num_spawned_this_call+7+5*ThisTask)*atan(sqrt(1+theta0*theta0)*tan(thetamax)/theta0))/sqrt(1+theta0*theta0)); // biased sampling to get collimation
        phi=2.*M_PI*get_random_number(num_spawned_this_call+1+ThisTask);
        cos_theta = cos(theta), sin_theta=sin(theta), sin_phi=sin(phi), cos_phi=cos(phi);
        veldir = nx*(sin_theta*cos_phi) + ny*(sin_theta*sin_phi) + nz*cos_theta; dpdir = veldir; //converted from angular momentum relative to into standard coordinates
    }
#if defined(SINGLE_STAR_FB_WINDS) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    else if (mode==2){ //random 3-axis isotropized - spawn along z axis, then y, then x
        if(((P[i].ID_generation-1) % 6) == 0) { // need to generate a brand new coordinate frame
            get_random_orthonormal_basis(P[i].ID_generation, nx, ny, nz);
            veldir = nz; for(k=0; k<3; k++) {P[i].Wind_direction[k]=nx[k]; P[i].Wind_direction[k+3]=ny[k];} dpdir = veldir;
        }
        else if(((P[i].ID_generation-1) % 6) == 2) {for(k=0; k<3; k++) {veldir[k] = P[i].Wind_direction[k];} dpdir = veldir;}
        else {for(k=0; k<3; k++) {veldir[k] = P[i].Wind_direction[k+3];} dpdir = veldir;}
    }
#endif
#if (defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)) || defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
    else if (mode==3) { // spawn on a specific angular grid
        int dir_ind = num_spawned_this_call % SINGLE_STAR_FB_SNE_N_EJECTA;
        veldir = nx*All.SN_Ejecta_Direction[dir_ind][0] + ny*All.SN_Ejecta_Direction[dir_ind][1] + nz*All.SN_Ejecta_Direction[dir_ind][2]; //use directions pre-computed to isotropically cover a sphere with SINGLE_STAR_FB_SNE_N_EJECTA particles
        dpdir = veldir;
    }
#endif
    return;
}


/* return desired cell launch speed for spawned cells, in physical (not comoving) units */
double get_spawned_cell_launch_speed(int i, struct particle_data *pp)
{
    double v_magnitude = All.Sink_outflow_velocity; // velocity of the jet: default mode is to set this manually to a specific value in physical units

#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES == 3)
    if(is_particle_a_special_zoom_target(i)) {return 1.e5/UNIT_VEL_IN_KMS;} // need an initial velocity for launch here //
#endif
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4)
    if(is_particle_a_special_zoom_target(i)) {return 3.e4/UNIT_VEL_IN_KMS;} // need an initial velocity for launch here //
#endif
    
#ifdef SNE_NONSINK_SPAWN
    if(pp[i].Type == 4) {
        double t_gyr = evaluate_stellar_age_Gyr(i); int SNeIaFlag=0; if(t_gyr > 0.03753) {SNeIaFlag=1;}; /* assume SNe before critical time are core-collapse, later are Ia */
        double Msne=10.5/UNIT_MASS_IN_SOLAR; if(SNeIaFlag) {Msne=1.4/UNIT_MASS_IN_SOLAR;} // average ejecta mass for single event (normalized to give total mass loss correctly)
        double SNeEgy = (1.0e51/UNIT_ENERGY_IN_CGS);
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
        if(SNeIaFlag==0) {double z_eff = pp[i].Metallicity[10]/All.SolarAbundances[10]; if(z_eff < 1) {SNeEgy *= pow(z_eff + 1.e-5 , -0.12);}} // updated to use same metallicity used for stellar evolution, rather than total metallicity, if this derives from pre-explosion winds, etc, for consistency
#if (FIRE_SNE_ENERGY_METAL_DEPENDENCE_EXPERIMENT > 1)
        if(i>0) {double z0 = pp[i].Metallicity[0]/All.SolarAbundances[0];
#if (FIRE_SNE_ENERGY_METAL_DEPENDENCE_EXPERIMENT > 2)
            SNeEgy *= pow(z0/0.1 + 1.e-3 , -0.2);
#else
            SNeEgy *= pow(z0/0.1 + 1.e-3 , -0.1);
#endif
        }
#endif
#endif
        return sqrt(2.0*SNeEgy/Msne); // v_ej in code units: assume all SNe = 1e51 erg //
    }
#endif

#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
    double Mdot_wind = pp[i].Sink_Mdot_ROI - pp[i].Sink_Mdot;
    if(Mdot_wind < 0) {return MIN_REAL_NUMBER;} // should be invalid
    double mdot = pp[i].Sink_Mdot / (pp[i].Sink_Mass / (4.e7 / UNIT_TIME_IN_YR));
    double L_over_c = evaluate_sink_radiative_efficiency(pp[i].Sink_Mdot,pp[i].Sink_Mass,i) * pp[i].Sink_Mdot * C_LIGHT_CODE;
    double Pdot_rad = 0.;
    if(mdot > 0.01) {Pdot_rad = L_over_c * DMIN(DMAX(mdot,1.),10.);}
    double sigma_ROI = sqrt(All.G * pp[i].Sink_Mass / pp[i].Sink_ROI);
    double Pdot_turb = 3. * Mdot_wind * sigma_ROI;
    double Pdot_wind = Pdot_rad + Pdot_turb;
    v_magnitude = Pdot_wind / Mdot_wind; 
    /* // (older deprecated model here)
    double MSINK_4 = pp[i].Sink_Mass * UNIT_MASS_IN_SOLAR / 1.e4; // sink mass in 1e4 Msun to scale
    double lambda_edd_eff = DMAX( pp[i].Sink_Mdot / sink_eddington_mdot(pp[i].Sink_Mass) , 1.e-10 ); // eddington ratio, with floor just to prevent unphysical behaviors
    if(lambda_edd_eff > (SINK_RIAF_SUBEDDINGTON_MODEL))
    {
        double v_eff_esc_BLR = 270. * sqrt(sqrt(MSINK_4 / lambda_edd_eff)) / UNIT_VEL_IN_KMS; // escape velocity from BLR in km/s, using canonical RBLR ~ 20 light-days * (L_bol/1e45)^(1/2)-ish scaling
        v_magnitude = DMIN(v_magnitude , v_eff_esc_BLR); // the input Sink_outflow_velocity parameter now sets the maximum efficiency/velocity this is allowed to reach, but it can be arbitrarily lower
    } else {
        v_magnitude = DMAX(v_magnitude , 1.e5 / UNIT_VEL_IN_KMS); // fast jet speed
    }
    */
#endif
    
#ifdef SINGLE_STAR_FB_JETS
    v_magnitude = single_star_jet_velocity(i); // get velocity from our more detailed function
#endif
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION) && defined(SINGLE_STAR_FB_WINDS)
    if((pp[i].ProtoStellarStage == 5) && (pp[i].wind_mode==1)) {v_magnitude = single_star_wind_velocity(i);} // only MS stars launch winds: get velocity from fancy model
#endif
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION) && defined(SINGLE_STAR_FB_SNE)
    if(pp[i].ProtoStellarStage == 6) {v_magnitude = single_star_SN_velocity(i);} // this star is about to go SNe: get velocity from fancy model
#endif
    return v_magnitude;
}


#ifdef MAGNETIC
void get_wind_spawn_magnetic_field(int j, int mode, Vec3<double>& ny, Vec3<double>& nz, Vec3<double>& dpdir, double d_r)
{
    int k; CellP[j].divB = 0; CellP[j].DtB = {};
#ifdef DIVBCLEANING_DEDNER
    CellP[j].DtPhi = CellP[j].PhiPred = CellP[j].Phi = 0; CellP[j].DtB_PhiCorr = {};
#endif

    double volume_for_BtoVB = P[j].Mass / CellP[j].Density;
#ifdef SINK_WIND_SPAWN_SET_BFIELD_POLTOR /* user manually sets the poloidal and toroidal components here */
    double inj_scale=All.Sink_spawn_injectionradius/All.cf_atime; Vec3<double> Bfield={}, nx=cross(ny,nz);
    double cos_theta=dot(nz,dpdir), sin_theta=sqrt(1-cos_theta*cos_theta), cos_phi=dot(nx,dpdir)/sin_theta, sin_phi=dot(ny,dpdir)/sin_theta;
    /* initialize poloidal component, in the nx,ny,nz coordinate frame */
    Bfield[0]+= All.B_spawn_pol*d_r*cos_theta*d_r*sin_theta/inj_scale/inj_scale*cos_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    Bfield[1]+= All.B_spawn_pol*d_r*cos_theta*d_r*sin_theta/inj_scale/inj_scale*sin_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    Bfield[2]+= All.B_spawn_pol*(1-d_r*cos_theta*d_r*sin_theta/inj_scale/inj_scale)    *exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    /* initialize toroidal component, in the nx,ny,nz coordinate frame */
    Bfield[0]+= -1*All.B_spawn_tor*(d_r/inj_scale)*sin_theta*sin_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    Bfield[1]+=    All.B_spawn_tor*(d_r/inj_scale)*sin_theta*cos_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    /* assign it back to the actual evolved B in the lab/simulation coordinate frame */
    CellP[j].IniB = nx*Bfield[0] + ny*Bfield[1] + nz*Bfield[2]; CellP[j].DtB = {};
    CellP[j].BPred = CellP[j].B = CellP[j].IniB * ((All.UnitMagneticField_in_gauss/UNIT_B_IN_GAUSS)*(volume_for_BtoVB/All.cf_a2inv));
    
#else /* set B-fields to be weak relative to local ISM values */

    double Bmag=0, Bmag_0=0;
    {auto Bphys = CellP[j].B * (All.cf_a2inv/volume_for_BtoVB); Bmag = Bphys.norm_sq(); Bmag_0 = CellP[j].B.norm_sq();} // get actual Bfield
    double Bmag_low_rel_to_progenitor = 1.e-10 * sqrt(Bmag); // set to some extremely low value relative to cloned element
    double u_internal_new_cell = All.Sink_outflow_temperature / (  MEAN_MOLECULAR_WEIGHT_DEFAULT * (GAMMA_DEFAULT-1.) * U_TO_TEMP_UNITS ); // internal energy of new wind cell
    double Bmag_low_rel_to_pressure = 1.e-3 * sqrt(2.*CellP[j].Density*All.cf_a3inv * u_internal_new_cell); // set to beta = 1e6
    Bmag = DMAX(Bmag_low_rel_to_progenitor , Bmag_low_rel_to_pressure); // pick the larger of these (still small) B-field values
#ifdef MHD_B_SET_IN_PARAMS
    double Bmag_IC = sqrt(All.BiniX*All.BiniX + All.BiniY*All.BiniY + All.BiniZ*All.BiniZ) * All.UnitMagneticField_in_gauss / UNIT_B_IN_GAUSS; // IC B-field sets floor as well
    Bmag = DMAX(Bmag , 0.1 * Bmag_IC);
#endif
#if defined(SINGLE_STAR_FB_SNE)
    if(P[j].Type==5) {if(P[j].ProtoStellarStage == 6) {Bmag *= 1.e-3;}} // No need to have flux in SN ejecta - note that this assumes we inherited this attribute from the spawning sink before calling this routine
#endif
    Bmag = DMAX(Bmag, MIN_REAL_NUMBER); // floor to prevent underflow errors
    /* add magnetic flux here to 'Bmag' if desired */
    Bmag *= volume_for_BtoVB / All.cf_a2inv; // convert back to code units
    for(k=0;k<3;k++) {if(Bmag_0>0) {CellP[j].B[k]*=Bmag/sqrt(Bmag_0);} else {CellP[j].B[k]=Bmag;}} // assign if valid values
    CellP[j].BPred=CellP[j].B; CellP[j].DtB={}; // set predicted = actual, derivative to null
#endif
    CellP[j].BField_prerefinement = CellP[j].B * (1.0 / volume_for_BtoVB); /* record the real value of B pre-split to know what we need to correctly re-initialize to once the volume partition can be recomputed */
    CellP[j].BPred = CellP[j].B; /* set predicted/drifted equal to the value above */
    return;
}
#endif


/*! this code copies what was used in merge_split.c for the gas particle split case */
int sink_spawn_particle_wind_shell( int i, int dummy_cell_i_to_clone, int num_already_spawned )
{
    if(ghost_require_no_live_pool_for_layout_change("sink_spawn_particle_wind_shell")) {return 0;}
    double *unspawned_mass_ptr = active_unspawned_mass_ptr(i); // whichever reservoir currently holds the discrete-spawn channel for this particle
    double total_mass_in_winds = *unspawned_mass_ptr;

    int n_particles_split   = (int) floor( total_mass_in_winds / target_mass_for_wind_spawning(i) ); /* if we set SINK_WIND_SPAWN we presumably wanted to do this in an exactly-conservative manner, which means we want to have an even number here. */
    int k=0, j;   /* j is a particle index bounded by NumPart, so int, matching gizmo_mark_kernel_radius_dirty_indices(const int*) */

#if defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    if(P[i].Type==5) {
        if(P[i].ProtoStellarStage == 6) {
            n_particles_split = (int) floor( total_mass_in_winds / (P[i].Sink_Formation_Mass) );
            double m_relic = single_star_relic_SN_mass(i); // get the intended relic mass //
            if(P[i].Sink_Mass <= m_relic) { // last batch to be spawned
                n_particles_split = SINGLE_STAR_FB_SNE_N_EJECTA; // we are going to spawn a bunch of low mass particles to take the last bit of mass away
                printf("Spawning last SN ejecta of star %llu with %g mass and %d particles \n",(unsigned long long) P[i].ID,total_mass_in_winds,n_particles_split);
                P[i].Mass = DMAX(0, m_relic); // set mass to zero so that this sink will get cleaned up (TreeReconstructFlag = 1 should be already set in sink.c) if(P[i].Type==0) {CellP[i].Mass = P[i].Mass;}
#ifdef SINK_ALPHADISK_ACCRETION
                P[i].Sink_Mass_Reservoir = 0; // just to be safe
#endif
                if(P[i].Sink_Mass > 0 && P[i].Sink_Mass > P[i].Sink_Formation_Mass) {P[i].ProtoStellarStage = 7;} // this is a relic now, move it to the next stage
	    }
	}
    }
#if (defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)) || defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
    if(P[i].Type==5) {if(P[i].ProtoStellarStage == 6)
    {
        if (n_particles_split<SINGLE_STAR_FB_SNE_N_EJECTA) {return 0;} // we have to wait until we get a full shell
        else {n_particles_split = n_particles_split - (n_particles_split % SINGLE_STAR_FB_SNE_N_EJECTA);} // we only eject full shells, in practice this will be one shell at a time
    }}
#endif
#endif
    if((((int)SINK_WIND_SPAWN) % 2) == 0) {if(( n_particles_split % 2 ) != 0) {n_particles_split -= 1;}} /* n_particles_split was not even. we'll wait to spawn this last particle, to keep an even number, rather than do it right now and break momentum conservation */	
    if( (n_particles_split == 0) || (n_particles_split < 1) ) {return 0;}
    int n0max = DMAX(20 , (int)(3.*(SINK_WIND_SPAWN)+0.1));
#if defined(SNE_NONSINK_SPAWN)
    n0max = DMAX(6 , (int)(3.*(SINK_WIND_SPAWN)+0.1)); // more conservative to spread over more timesteps to avoid nasty overlaps //
#if (defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)) || defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
    n0max = DMAX(SINGLE_STAR_FB_SNE_N_EJECTA , (int)(3.*(SINK_WIND_SPAWN)+0.1)); if((n0max % 2) != 0) {n0max += 1;} // should ensure n0max is always an even number //
#endif
#endif
    if((n0max % 2) != 0) {n0max += 1;} // should ensure n0max is always an even number //
#if defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    if(P[i].Type==5) {if(P[i].ProtoStellarStage == 6) {n0max = DMAX(n0max, SINGLE_STAR_FB_SNE_N_EJECTA);}} // so that we can spawn the number of wind particles we want, by setting SINK_WIND_SPAWN high it ispossible to spawn multitudes of SINGLE_STAR_FB_SNE_N_EJECTA, but in practice we usually spawn just one
#endif
    if(n_particles_split > n0max) {n_particles_split = n0max;}


    /* here is where the details of the split are coded, the rest is bookkeeping */
    //double mass_of_new_particle = total_mass_in_winds / n_particles_split; /* don't do this, as can produce particles with extremely large masses; instead wait to spawn */
    double mass_of_new_particle = target_mass_for_wind_spawning(i); double mass_of_new_particle_default,mass_of_new_particle_prev; mass_of_new_particle_default=mass_of_new_particle; mass_of_new_particle_prev=mass_of_new_particle;
#if defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    if(P[i].Type==5) {if(P[i].ProtoStellarStage == 6) {mass_of_new_particle = total_mass_in_winds/(double) n_particles_split;}} // ejecta will have the gas mass resolution except the last batch which will lower masses
#endif
    printf("Task %d wants to create %g mass in wind with %d new particles each of mass %g \n .. splitting sink %d using hydro element %d\n", ThisTask,total_mass_in_winds, n_particles_split, mass_of_new_particle, i, dummy_cell_i_to_clone);

    /* the spawned wind elements are gas, taking a slot in both P[] and CellP[] at the same index, so
       both capacities have to hold them */
    if(!gizmo_particle_index_fits_live_tree(NumPart + num_already_spawned + n_particles_split))
    {
        /* Storage has room; the tree standing right now does not reach these slots, and spawned wind
         * elements are inserted into it.  Decline this spawn the way the routine already declines
         * when it is out of room -- the wind is spawned after the next rebuild instead. */
        printf("On Task=%d we skip spawning %d wind elements from sink %d: the standing tree has only %d particle slots. They can spawn after the next tree rebuild.\n", ThisTask, n_particles_split, i, All.TreeParticleSlots);
        fflush(stdout);
        return 0;
    }
    if(NumPart + num_already_spawned + n_particles_split >= All.MaxPart ||
       NumPart + num_already_spawned + n_particles_split >= All.MaxPartGas)
    {
        printf("On Task=%d there is no room to spawn a wind element: %d particles are here (+%d already spawned this step) and this rank holds %d particle and %d gas slots. Raise ParticleNumberMemoryImbalance_Limit (PartAllocFactor) so a rank may be assigned more of them, or spread the run over more ranks.\n", ThisTask, NumPart, num_already_spawned, All.MaxPart, All.MaxPartGas);
        fflush(stdout); endrun(8888);
        return 0;   /* no space left: honor the no-split contract (return 0) instead of spawning past MaxPart */
    }
    double d_r = 0.25 * KERNEL_CORE_SIZE*P[i].KernelRadius; // needs to be epsilon*KernelRadius where epsilon<<1, to maintain stability //
    double r2 = (P[dummy_cell_i_to_clone].Pos - P[i].Pos).norm_sq();
    d_r = DMIN(d_r, 0.5*sqrt(r2));
#ifndef SELFGRAVITY_OFF
    d_r = DMAX(d_r , 2.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0]);
#endif
#ifdef SINK_WIND_SPAWN_SET_BFIELD_POLTOR
    d_r = DMIN(d_r , All.Sink_spawn_injectionradius/All.cf_atime); /* KYSu: sets spawn scale manually */
#endif
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
    d_r = DMIN(P[i].SinkRadius, d_r); //launch close to the sink
#endif
#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM) && defined(PARTICLE_EXCISION)
    if(is_particle_a_special_zoom_target(i)) {double rmin=All.ForceSoftening[3], r=sqrt(r2), r0=0.5*(rmin+r)*(1.+0.1*get_random_number(i+j)); d_r=r0;} // make sure to spawn OUTSIDE of the excision radius!
#endif
#if defined(SNE_NONSINK_SPAWN)
    if(P[i].Type == 4) {double rmin=All.ForceSoftening[4], r=sqrt(r2), r0=0.5*(rmin+r)*(0.5+1.5*get_random_number(i+j)); d_r=r0;} // need a generous padding to ensure no overlaps
#endif
    long bin, bin_0; for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinCount[bin] > 0) break;} /* gives minimum active timebin of any particle */
    bin_0 = bin; int i0 = i; /* save minimum timebin, also save ID of sink particle for use below */
    bin = P[i0].TimeBin; /* make this particle active on the BH/star timestep */
    Vec3<double> veldir, dpdir; // velocity direction to spawn in - declare outside the loop so we remember it from the last iteration
    int mode = 0; // 0 if doing totally random directions, 1 if collimated, 2 for 3-axis isotropized, and 3 if using an angular grid,  4 old collimatation script, position isotropic velicity coliminated within certain open angle (might be useful to still keep this option owing to the free open angle choice and better sampling the magnetic field geometry)
#if defined(SINGLE_STAR_FB_JETS) || defined(JET_DIRECTION_FROM_KERNEL_AND_SINK) || defined(SINK_FB_COLLIMATED)
    mode = 1; // collimated mode
#endif
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
    mode=0; // broad-angle default, but will modify below
#endif
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    mode = 0;
#if defined(SINGLE_STAR_FB_JETS)
    mode = 1; // launch polar jets
#endif
#if defined(SINGLE_STAR_FB_WINDS)
    if(P[i].Type==5) {if(P[i].ProtoStellarStage == 5) {
        mode = 2; // winds use 3-axis isotropized directions
#if defined(SINGLE_STAR_FB_JETS)
        if(P[i].wind_mode == 2) {mode = 1;} // we inject winds with the FIRE module, we only spawn polar jets here
#endif
    }}
#endif
#if defined(SINGLE_STAR_FB_SNE)
    if(P[i].Type==5) {if(P[i].ProtoStellarStage == 6) {mode = 3;}} // SNe use an angular grid
#endif
#endif // single-star evolution clause
    if(P[i].Type==3) {mode = 1;} // special particle spawn is collimated
    if(P[i].Type==4) {mode = 0;} // star particle spawn is isotropic
#ifdef SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT
    if(P[i].Type==4) {mode = 3;} // star particle spawn is isotropic but regularized
#endif

    // based on the mode we're in, let's pick a fixed orthonormal basis that all spawned elements are aware of
    Vec3<double> jz={0,0,1},jy={0,1,0},jx={1,0,0};  /* set up a coordinate system [xyz if we don't have any other information */
#ifndef JET_DIRECTION_FIXED_Z  /* TESTING ONLY: pin the launch axis to +/-z by keeping the default
                                  basis, skipping the angular-momentum reorientation and the
                                  precession below, so a launch geometry can be checked against a
                                  known axis. Real jets follow the disk/spin axis, so this is not
                                  physical for production. */
#ifdef SINK_FOLLOW_ACCRETED_ANGMOM  /* use local angular momentum to estimate preferred directions/coordinates for spawning */
    if(mode==1){ // set up so that the z axis is the angular momentum vector
#ifdef JET_DIRECTION_FROM_KERNEL_AND_SINK // Jgas stores total angmom in COM frame of sink-gas system; use this for direction
        double Jtot=P[i].Jgas_in_Kernel.norm_sq();
        if(Jtot>0) {Jtot=1/sqrt(Jtot); jz = P[i].Jgas_in_Kernel * Jtot;}
#else
        double Jtot=P[i].Sink_Specific_AngMom.norm_sq();
        if(Jtot>0) {Jtot=1/sqrt(Jtot); jz = P[i].Sink_Specific_AngMom * Jtot;}
#endif
        Jtot=jz[1]*jz[1]+jz[2]*jz[2]; if(Jtot>0) {Jtot=1/sqrt(Jtot); jy={0, jz[2]*Jtot, -jz[1]*Jtot};} else {jy={0, 1, 0};}
        jx = cross(jz, jy);
    }
#endif
    if(mode == 3){ // if doing an angular grid, need some fixed coordinates to orient it, but want to switch em up each time to avoid artifacts
        get_random_orthonormal_basis(P[i].ID_generation, jx, jy, jz);
    }
#ifdef SINK_WIND_SPAWN_SET_JET_PRECESSION /* rotate the jet angle according to the explicitly-included precession parameters */
    double degree = All.Sink_jet_precess_degree, period = All.Sink_jet_precess_period/UNIT_TIME_IN_GYR; Vec3<double> new_dir;
    new_dir[0]= jx[0]*cos(degree/180.*M_PI)-jx[2]*sin(degree/180.*M_PI); new_dir[1]= 1.0*jx[1]; new_dir[2]= jx[0]*sin(degree/180.*M_PI)+jx[2]*cos(degree/180.*M_PI);
    jx[0]= new_dir[0]*cos(2.*M_PI/period*All.Time)-new_dir[1]*sin(2.*M_PI/period*All.Time); jx[1]= new_dir[0]*sin(2.*M_PI/period*All.Time)+new_dir[1]*cos(2.*M_PI/period*All.Time); jx[2]= new_dir[2];

    new_dir[0]= jy[0]*cos(degree/180.*M_PI)-jy[2]*sin(degree/180.*M_PI); new_dir[1]= 1.0*jy[1]; new_dir[2]= jy[0]*sin(degree/180.*M_PI)+jy[2]*cos(degree/180.*M_PI);
    jy[0]= new_dir[0]*cos(2.*M_PI/period*All.Time)-new_dir[1]*sin(2.*M_PI/period*All.Time); jy[1]= new_dir[0]*sin(2.*M_PI/period*All.Time)+new_dir[1]*cos(2.*M_PI/period*All.Time); jy[2]= new_dir[2];

    new_dir[0]= jz[0]*cos(degree/180.*M_PI)-jz[2]*sin(degree/180.*M_PI); new_dir[1]= 1.0*jz[1]; new_dir[2]= jz[0]*sin(degree/180.*M_PI)+jz[2]*cos(degree/180.*M_PI);
    jz[0]= new_dir[0]*cos(2.*M_PI/period*All.Time)-new_dir[1]*sin(2.*M_PI/period*All.Time); jz[1]= new_dir[0]*sin(2.*M_PI/period*All.Time)+new_dir[1]*cos(2.*M_PI/period*All.Time); jz[2]= new_dir[2];
#endif    
#endif  /* closes !JET_DIRECTION_FIXED_Z */

    /* create the  new particles to be added to the end of the particle list :
        i is the sink particle tag, j is the new "spawed" particle's location, dummy_cell_i_to_clone is a dummy gas cell's tag to be used to init the wind particle */
    int mode_default = mode, mode_prev = mode;
    double v_magnitude_physical_default = get_spawned_cell_launch_speed(i, P), v_magnitude_physical=v_magnitude_physical_default, v_magnitude_physical_prev=v_magnitude_physical; /* call subroutine for this velocity */
    
    for(j = NumPart + num_already_spawned; j < NumPart + num_already_spawned + n_particles_split; j++)
    {   /* first, clone the 'dummy' particle so various fields are set appropriately */
        P[j] = P[dummy_cell_i_to_clone]; CellP[j] = CellP[dummy_cell_i_to_clone]; /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */
        wakeup_sidecar_mark(j);   /* whole-struct clone inherits the template particle's wakeup into the spawned slot */

#if defined(SINK_TEST_WIND_MIXED_FASTSLOW) || defined(SINK_RIAF_SUBEDDINGTON_MODEL)
        if(P[i].Type==5) {
            double masscorrfac_fast = 100.; /* ratio of spawned jet cell mass to non-jet cell mass */
            double fraction_to_spawn_in_jet = 0.1; /* fraction of spawned cells by number in jet */
            double frac_clight_jet = 0.1; /* default fraction of C for jet speed */
#ifdef SINK_TEST_WIND_MIXED_FASTSLOW
            frac_clight_jet = (SINK_TEST_WIND_MIXED_FASTSLOW/UNIT_VEL_IN_KMS) / C_LIGHT_CODE;
#endif
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
            double frac_clight_jet_max = 0.10; //1.0; // don't let the jet be too fast, for physical (superluminal) or numerical (timestep) reasons
            double frac_clight_jet_min = 0.03; //0.1; // don't let the jet be too slow, or it won't behave like a jet; lower jet mass to compensate
            double Mdot_wind = P[i].Sink_Mdot_ROI - P[i].Sink_Mdot;
            if(Mdot_wind > 0)
            {
                double a_spin = 0.33;
                double HR = 0.33;
                double mdot = P[i].Sink_Mdot / (P[i].Sink_Mass / (5.e7/UNIT_TIME_IN_YR));
                if(mdot < 0.01) {HR = 1;}
                double eff_jet = 0.1 * pow(a_spin*HR,2);
                double Mdot_jet = (fraction_to_spawn_in_jet/masscorrfac_fast) * Mdot_wind;
                double eta_jet = Mdot_jet / P[i].Sink_Mdot;
                frac_clight_jet = sqrt(2.*eff_jet/eta_jet); /* scaling so that KE of jet = desired */
                if(frac_clight_jet > frac_clight_jet_max) { /* superluminal - need more mass in jet to make this make sense */
                    frac_clight_jet = frac_clight_jet_max; // cap this at luminal
                    masscorrfac_fast = frac_clight_jet*frac_clight_jet * (fraction_to_spawn_in_jet/(2.*eff_jet)) * (Mdot_wind / P[i].Sink_Mdot); // boost this term to make up the difference
                }
                if(frac_clight_jet < frac_clight_jet_min) {
                    frac_clight_jet = frac_clight_jet_min; // cap this at minimum
                    masscorrfac_fast = frac_clight_jet*frac_clight_jet * (fraction_to_spawn_in_jet/(2.*eff_jet)) * (Mdot_wind / P[i].Sink_Mdot); // boost this term to make up the difference
                }
            }
#endif
            double m0_newparticlemass_for_target_spawnedmass = mass_of_new_particle_default / (1. - fraction_to_spawn_in_jet * (1.-1./masscorrfac_fast));
            mode = mode_default; mass_of_new_particle=m0_newparticlemass_for_target_spawnedmass; v_magnitude_physical=v_magnitude_physical_default;
            if((j - (NumPart + num_already_spawned)) % 2) {mode=mode_prev; v_magnitude_physical=v_magnitude_physical_prev; mass_of_new_particle=mass_of_new_particle_prev;  /* for every-other particle, need to match previous for conservation */
            } else { /* collimated jet */
                if(get_random_number(j) < fraction_to_spawn_in_jet) {
                    mode=1; mass_of_new_particle=m0_newparticlemass_for_target_spawnedmass/masscorrfac_fast; v_magnitude_physical = frac_clight_jet * C_LIGHT_CODE;
                } else { /* isotropic/broad-angle wind */
                    mode=0; mass_of_new_particle=m0_newparticlemass_for_target_spawnedmass; v_magnitude_physical=v_magnitude_physical_default; /* isotropic slow wind */
                }
            }
        }
#endif
        v_magnitude_physical_prev = v_magnitude_physical; mode_prev = mode; mass_of_new_particle_prev=mass_of_new_particle;

        /* now we need to make sure everything is correctly placed in timebins for the tree */
        P[j].TimeBin = bin; // get the timebin, and put this particle into the appropriate timebin
        ActiveParticleList.push_back(j);
        NumForceUpdate++;
        TimeBinCount[bin]++; TimeBinCountGas[bin]++; PrevInTimeBin[j] = i0; /* likewise add it to the counters that register how many particles are in each timebin */
        NextInTimeBin[j] = NextInTimeBin[i0]; if(NextInTimeBin[i0] >= 0) {PrevInTimeBin[NextInTimeBin[i0]] = j;} NextInTimeBin[i0] = j; if(LastInTimeBin[bin] == i0) {LastInTimeBin[bin] = j;}
        P[j].Ti_begstep = All.Ti_Current; P[j].Ti_current = All.Ti_Current;
        P[j].dt_step = GET_INTEGERTIME_FROM_TIMEBIN(bin);
        P[j].wakeup = -1; wakeup_sidecar_mark(j);
        NeedToWakeupParticles_local = 1;
        /* this is a giant pile of variables to zero out. dont need everything here because we cloned a valid particle, but handy anyways */
        P[j].Particle_DivVel = 0; CellP[j].DtInternalEnergy = 0; CellP[j].HydroAccel = {}; P[j].GravAccel = {};
        P[j].NumNgb=cbrt(All.DesNumNgb); // this gets cube rooted at the end of the density loop, so take cbrt here
#ifdef PMGRID
        P[j].GravPM = {};
#endif
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        CellP[j].MaxKineticEnergyNgb = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        CellP[j].dMass = 0; CellP[j].DtMass = 0; CellP[j].MassTrue = P[j].Mass; CellP[j].GravWorkTerm = {};
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        P[j].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        P[j].AGS_KernelRadius = P[j].KernelRadius;
#endif
#endif
#ifdef ADAPTIVE_TREEFORCE_UPDATE
        P[j].tdyn_step_for_treeforce = 0; P[j].time_since_last_treeforce = MAX_REAL_NUMBER; // make sure we get a new tree force right off the bat
#endif
#ifdef CONDUCTION
        CellP[j].Kappa_Conduction = 0;
#endif
#ifdef MHD_NON_IDEAL
        CellP[j].Eta_MHD_OhmicResistivity_Coeff = 0; CellP[j].Eta_MHD_HallEffect_Coeff = 0; CellP[j].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
#ifdef VISCOSITY
        CellP[j].Eta_ShearViscosity = 0; CellP[j].Zeta_BulkViscosity = 0;
#endif
#ifdef TURB_DIFFUSION
        CellP[j].TD_DiffCoeff = 0;
#endif
#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==1)
        CellP[j].HostHaloMass = 0;
#endif
#endif
#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
        CellP[j].DelayTimeHII = 0;
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        CellP[j].DelayTimeCoolingSNe = 0;
#endif
#ifdef GALSF
        CellP[j].Sfr = 0;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        CellP[j].alpha = 0.0;
#endif
#if defined(SINK_THERMALFEEDBACK)
        CellP[j].Injected_Sink_Energy = 0;
#endif
#ifdef RADTRANSFER
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            CellP[j].Rad_E_gamma[k] = 0;
#if defined(RT_EVOLVE_ENERGY)
            CellP[j].Rad_E_gamma_Pred[k] = 0; CellP[j].Dt_Rad_E_gamma[k] = 0;
#endif
#if defined(RT_EVOLVE_FLUX)
            int kdir; for(kdir=0;kdir<3;kdir++){CellP[j].Rad_Flux[k][kdir] = 0;}
#endif
        }
#endif
#if defined(GALSF)
        P[j].StellarAge = All.Time; // use this attibute to save the gas cell's formation time for possible subsequent checks for special behavior on its first timestep
#endif

        /* now set the real hydro variables. */
        /* set the particle ID */ // unsigned int bits; int SPLIT_GENERATIONS = 4; for(bits = 0; SPLIT_GENERATIONS > (1 << bits); bits++); /* the particle needs an ID: we give it a bit-flip from the original particle to signify the split */
        P[j].ID = All.SpawnedWindCellID; /* update:  We are using a fixed wind ID, to allow for trivial wind particle identification */
#if defined(SINGLE_STAR_SINK_DYNAMICS)
        if(mass_of_new_particle >= 0.5*P[i].Sink_Formation_Mass) {P[j].ID = All.SpawnedWindCellID + 1;} // this just has the nominal mass resolution, so no special treatment - this avoids the P[i].ID == All.SpawnedWindCellID checks throughout the code
#endif
        P[j].ID_child_number = P[i].ID_child_number + P[i].ID_generation; P[i].ID_generation++; P[j].ID_generation = P[i].ID; // this allows us to track spawned particles by giving them unique sub-IDs. Remember we MUST NEVER alter an existing particle ID OR ID_child_number!
        P[j].Mass = mass_of_new_particle; /* assign masses to both particles (so they sum correctly) */ if(P[j].Type==0) {CellP[j].Mass = P[j].Mass;}
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        CellP[j].MassTrue = P[j].Mass;
#endif
        P[i].dp -= P[j].Mass * P[i].Vel; /* track momentum change from mass loss for tree node update */
        P[i].Mass -= P[j].Mass; /* make sure the operation is mass conserving! */ if(P[i].Type==0) {CellP[i].Mass = P[i].Mass;}
        *unspawned_mass_ptr -= P[j].Mass; /* remove the mass successfully spawned, to update the remaining unspawned mass in whichever reservoir is active */

#if defined(METALS) && (defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_SNE) || defined(SNE_NONSINK_SPAWN) || (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4))
        double yields[NUM_METAL_SPECIES]={0}; get_jet_yields(yields,i); // default to jet-type
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION) && defined(SINGLE_STAR_FB_WINDS)
        if(P[i].Type==5) {if((P[i].ProtoStellarStage == 5) && (P[i].wind_mode==1)) {get_wind_yields(yields,i);}} // get abundances in wind
#endif
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION) && defined(SINGLE_STAR_FB_SNE)
        if(P[i].Type==5) {if(P[i].ProtoStellarStage == 6) {double Msne; get_SNe_yields(yields,i,stellar_lifetime_in_Gyr(i),0,&Msne);}} // get sne yields
#endif
        for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k]=yields[k];}  // update metallicity of spawned cell modules
#endif
        
        // actually lay down position and velocities using coordinate basis
        get_wind_spawn_direction(i, j - (NumPart + num_already_spawned), mode, jy, jz, veldir, dpdir);
        P[j].Pos = P[i].Pos + dpdir*d_r; P[j].Vel = P[i].Vel + veldir*(v_magnitude_physical*All.cf_atime); CellP[j].VelPred = P[j].Vel; // convert to code (comoving) velocity units
#if defined(USE_TIMESTEP_DILATION_FOR_ZOOMS)
        /* the clone inherited the template cell's factor, which belongs to a different position; this
           cell was assigned a timebin and Ti_begstep of its own above, so freeze its factor here at
           the position just laid down, as any other timestep assignment would */
        P[j].TimestepDilationFactor = return_timestep_dilation_factor(j, P);
#endif

        /* condition number, smoothing length, and density */
        CellP[j].ConditionNumber *= 100.0; /* boost the condition number to be conservative, so we don't trigger madness in the kernel */
        CellP[j].recent_refinement_flag = 1; /* tag the newly-created cell as recently-refined for all purposes */
#if defined(SINGLE_STAR_SINK_DYNAMICS)
        CellP[j].MaxSignalVel = 2.*DMAX(v_magnitude_physical, CellP[j].MaxSignalVel); // need this to satisfy the Courant condition in the first timestep after spawn; note here MaxSignalVel is now defined in physical code units
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
	    P[j].ProtoStellarStage = P[i].ProtoStellarStage; // inherit this from the spawning sink particle so we can use it in subroutines
        // need to initialize the gas density and search radius so that we get sensible CFL timesteps (which happens before density() is called and we recalculate these self-consistently)
        if(n_particles_split > All.DesNumNgb) { // we are spawning a whole "shell" together, so initialize search radii/densities assuming kernels are confined to the region of spawned material.
            CellP[j].Density = mass_of_new_particle / (4 * M_PI * d_r*d_r*d_r);
            P[j].KernelRadius = P[j].NumNgb * 2.32489404843 * d_r;
        } else { // we are spawning in the jet/wind piecemeal, so use the local density estimator around the star
            CellP[j].Density = P[i].DensityAroundParticle;
            P[j].KernelRadius = P[i].KernelRadius;
        }
        /* New particle's KernelRadius was just initialized. Mark h-dirty for
         * all caches that survive across this spawn (caches whose range covers
         * j). The new particle is appended, so a kept index built before it
         * fails its particle-count check and is rebuilt; this h-dirty mark keeps
         * the compact_xyzh.h slot fresh-on-next-build. */
        gizmo_mark_kernel_radius_dirty_indices(&j, 1);
#endif
#endif
        /* The density above is only a guess for the single timestep before this cell's first density
           pass, but it has to be positive and finite: a zero seed -- which is what the spawning sink
           hands over whenever no gas lies inside its own kernel, so its density estimator is empty --
           gives zero pressure, a 0/0 sound speed, zero conduction and viscosity timesteps, and an
           infinite cell volume in the magnetic field seed below. Fall back to, and keep a loose floor
           well below, the density the density loop itself assigns an element with no neighbors in its
           kernel. The floor is deliberately far under that value because the seeded kernel radius can
           badly under-estimate the converged one, which would over-estimate the density. */
        double rho_alone_in_kernel = 0;
        if((P[j].KernelRadius > 0) && isfinite(P[j].KernelRadius)) {rho_alone_in_kernel = P[j].Mass / (VOLUME_NORM_COEFF_FOR_NDIMS * pow((double)P[j].KernelRadius, (double)NUMDIMS));}
        if(!(rho_alone_in_kernel > 0) || !isfinite(rho_alone_in_kernel)) {rho_alone_in_kernel = CellP[dummy_cell_i_to_clone].Density;} /* the cloned element always has positive density */
        if(!(CellP[j].Density > 0) || !isfinite(CellP[j].Density)) {CellP[j].Density = rho_alone_in_kernel;}
        CellP[j].Density = DMAX(CellP[j].Density , 0.01 * rho_alone_in_kernel);

        /* note, if you want to use this routine to inject magnetic flux or cosmic rays, do this below */
#if defined(SINK_WIND_SPAWN_SET_BFIELD_POLTOR)
        CellP[j].IniDen = -1. * CellP[j].Density; /* this is essentially acting like a bitflag, to signal to the code that the density needs to be recalculated because a spawn event just occurred */
#endif
#ifdef MAGNETIC
        get_wind_spawn_magnetic_field(j, mode, jy, jz, dpdir, d_r);
#endif
#ifdef COSMIC_RAY_FLUID
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
        CellP[j].DtCREgyNewInjectionFromShocks=0;
#endif
        int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) /* initialize CR energy and other related terms to nil */
        {
            CellP[j].CosmicRayEnergyPred[k_CRegy]=CellP[j].CosmicRayEnergy[k_CRegy]=CellP[j].DtCosmicRayEnergy[k_CRegy]=0;
#ifdef CRFLUID_EVOLVE_SPECTRUM
            CellP[j].CosmicRay_Number_in_Bin[k_CRegy]=CellP[j].DtCosmicRay_Number_in_Bin[k_CRegy]=0;
#endif
            CellP[j].CosmicRayFlux[k_CRegy] = {}; CellP[j].CosmicRayFluxPred[k_CRegy] = {};
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            CellP[j].CosmicRayAlfvenEnergy[k_CRegy] = {}; CellP[j].CosmicRayAlfvenEnergyPred[k_CRegy] = {}; CellP[j].DtCosmicRayAlfvenEnergy[k_CRegy] = {};
#endif
        } /* complete CR initialization to null */
#endif
        CellP[j].InternalEnergy = All.Sink_outflow_temperature / (  MEAN_MOLECULAR_WEIGHT_DEFAULT * (GAMMA_DEFAULT-1.) * U_TO_TEMP_UNITS ); /* internal energy, determined by desired wind temperature (assume fully ionized primordial gas with gamma=5/3) */
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
        CellP[j].InternalEnergy = 0.01 * (0.5*v_magnitude_physical*v_magnitude_physical); /* set to be 1% of the kinetic energy of the ejecta, here */
#endif
#if defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
        double sne_energy_fraction_in_thermal = 1.e-3;
        if(P[i].Type==5) {if(P[i].ProtoStellarStage == 6) {CellP[j].InternalEnergy = All.MinGasTemp / (  MEAN_MOLECULAR_WEIGHT_MOLECULAR * (1.4-1.) * U_TO_TEMP_UNITS ) + sne_energy_fraction_in_thermal/(1.-sne_energy_fraction_in_thermal) * pow(single_star_SN_velocity(i),2.0);}}
#endif
        CellP[j].InternalEnergyPred = CellP[j].InternalEnergy;
        CellP[j].MeanMolecularWeight = MEAN_MOLECULAR_WEIGHT_DEFAULT; CellP[j].Gamma = GAMMA_DEFAULT; /* the ionized ideal gas assumed when the energy above was set from a desired temperature */
        CellP[j].Temperature = CellP[j].gas_temperature_from_u(CellP[j].InternalEnergy);
#ifdef EOS_ANCHOR_INTERNALENERGY_IN_DRIFTS
        CellP[j].u_anchor = CellP[j].InternalEnergy;
#endif

#if defined(COSMIC_RAY_FLUID) && defined(SINK_COSMIC_RAYS) /* inject cosmic rays alongside wind injection */
        double eps_cr = evaluate_sink_cosmicray_efficiency(P[i].Sink_Mdot,P[i].Sink_Mass,i);
        double fac_wind_corr = All.Sink_accreted_fraction / (1.-All.Sink_accreted_fraction);
        double dEcr = eps_cr * P[j].Mass * fac_wind_corr * C_LIGHT_CODE*C_LIGHT_CODE;
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
        dEcr = eps_cr * 0.5 * v_magnitude_physical*v_magnitude_physical * P[j].Mass; /* in this case, eps_cr refers to fraction relative to KE going into jet spawns */
        if(mass_of_new_particle > 2.*mass_of_new_particle_default/masscorrfac_fast) {dEcr=0;} /* only the jet cells carry CR energy */
#endif
#if defined(SINK_TEST_WIND_MIXED_FASTSLOW)
        if(mass_of_new_particle < 2.*mass_of_new_particle_default/masscorrfac_fast) {dEcr*=masscorrfac_fast;} else {dEcr=0;}
#endif
#if defined(SINK_CR_INJECTION_AT_TERMINATION)
        CellP[j].Sink_CR_Energy_Available_For_Injection = dEcr;     /* store energy for later injection */
#else
        inject_cosmic_rays(dEcr, v_magnitude_physical, 5, j, veldir.data, P, CellP); /* inject directly */
#endif
#endif
        /* Note: New tree construction can be avoided because of  `force_add_element_to_tree()' */
        force_add_element_to_tree(i0, j);// (buggy) /* we solve this by only calling the merge/split algorithm when we're doing the new domain decomposition */
    }
    if(*unspawned_mass_ptr < 0) {*unspawned_mass_ptr=0;}
    return n_particles_split;
}
#endif



#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4)
/* routine for injection from sink boundary around 'special' particle types */
void special_rt_feedback_injection(void)
{
    double L0_cgs = 7.e45, MdotJetMsunYr=1., mspecial_tot=0; int iBH0=-1, k;
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES == 3)
    L0_cgs = 1.e43; MdotJetMsunYr = 1.e-3;
#endif
    for(k=0;k<SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM;k++) {mspecial_tot += All.Mass_of_SpecialParticle[k];}
    if(mspecial_tot <= 0) {return;}
    double delta_wt_sum = 0, delta_wt_sumsum, r_min = All.ForceSoftening[3] * All.cf_atime, r_max = 5. * r_min, dt = All.TimeStep, subgrid_lum = L0_cgs / (UNIT_ENERGY_IN_CGS/UNIT_TIME_IN_CGS), de_00 = subgrid_lum * dt; if(dt <= 0) {return;}
    /* Owner-only scans: bound by the LOCAL particle count, not NumPart. A live
       ghost pool (hydro-corridor span / TRANSPORT_SUBCYCLE) extends NumPart with
       ghost COPIES of remote particles: including them here double-counts the
       MPI-summed weight normalization below (owner rank + ghost-holding rank),
       and a ghost copy of the target BH would absorb the reservoir update into
       a discarded slot. ghost_get_num_local()==NumPart when no pool is live. */
    int n_local_real = ghost_get_num_local();
    int n_wt = 0, i; for(i=0;i<n_local_real;i++) {
        if(is_particle_a_special_zoom_target(i)) {iBH0=i;}
        if(P[i].Type != 0) {continue;}
        Vec3<double> dp{All.cf_atime*(double)P[i].Pos[0], All.cf_atime*(double)P[i].Pos[1], All.cf_atime*(double)P[i].Pos[2]};
        double r2 = dp.norm_sq(), wt, wt_new=0, r;
        r = sqrt(r2); if(r < r_min || r >= r_max) {continue;}
        double vol = P[i].Mass / (CellP[i].Density*All.cf_a3inv), cos_t = dp[0] / r;
        wt = 1.e-5 * pow(fabs(cos_t),8) * vol * (r_max*r_max/(r*r)-1.); delta_wt_sum += wt;
    }
    MPI_Allreduce(&delta_wt_sum, &delta_wt_sumsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // collect the information on weight sums
    if(All.Time <= All.TimeBegin) {return;}
    if(delta_wt_sumsum <= 0) {return;}
    for(i=0;i<n_local_real;i++) {   /* owner-only: see bound rationale above */
        if(P[i].Type != 0) {continue;}
        Vec3<double> dp{All.cf_atime*(double)P[i].Pos[0], All.cf_atime*(double)P[i].Pos[1], All.cf_atime*(double)P[i].Pos[2]};
        double r2 = dp.norm_sq(), wt, wt_new=0, r, de;
        r = sqrt(r2); if(r < r_min || r >= r_max) {continue;}
        double vol = P[i].Mass / (CellP[i].Density*All.cf_a3inv), cos_t = dp[0] / r;
        wt = 1.e-5 * pow(fabs(cos_t),8) * vol * (r_max*r_max/(r*r)-1.); de = de_00 * wt / delta_wt_sumsum;
        if(de <= 0 || !isfinite(de)) {continue;}
        k=RT_FREQ_BIN_INFRARED; double T00 = 1.e5, f0 = de * C_LIGHT_CODE / r;
        if(CellP[i].Radiation_Temperature > 0 && CellP[i].Rad_E_gamma[k] > 0) {CellP[i].Radiation_Temperature = (CellP[i].Rad_E_gamma[k] + de)/(CellP[i].Rad_E_gamma[k]/CellP[i].Radiation_Temperature + de/T00);} else {CellP[i].Radiation_Temperature = T00;}
        CellP[i].Rad_E_gamma[k] += de; CellP[i].Rad_E_gamma_Pred[k] += de;
        int j; for(j=0;j<3;j++) {CellP[i].Rad_Flux[k][j] += f0 * dp[j]; CellP[i].Rad_Flux_Pred[k][j] += f0 * dp[j];}
    }
    if(iBH0 >= 0) {
        P[iBH0].unspawned_wind_mass += MdotJetMsunYr * dt * (6.304e25 * UNIT_TIME_IN_CGS/UNIT_MASS_IN_CGS); // will sent to jets subroutine, for spawning, alongside radiation injection //
        double n_unspawned = P[iBH0].unspawned_wind_mass / ((SINK_WIND_SPAWN)*target_mass_for_wind_spawning(iBH0)); // number of spawned gas cells that can be made from the mass in the reservoir
        if(n_unspawned> Max_Unspawned_MassUnits_fromSink) {Max_Unspawned_MassUnits_fromSink = n_unspawned;} // track the maximum integer number of elements this sink could spawn
        P[iBH0].Sink_Specific_AngMom[0]=1; P[iBH0].Sink_Specific_AngMom[1]=0; P[iBH0].Sink_Specific_AngMom[2]=0; // placeholder direction; not yet set to the desired physical direction
    }
    return;
}
#endif


#ifdef SINK_WIND_SPAWN
/* Which reservoir currently feeds discrete spawning. This must agree with the channel that
   target_mass_for_wind_spawning() below prices cells for, or mass banks at one resolution and spawns at
   another. Only the STARFORGE accretion-jet channel banks into unspawned_jet_mass: that is a type-5 sink
   in a build that compiles the jets, and it is the same eligibility single_star_jet_mdot() applies. SNe
   ejecta always uses unspawned_wind_mass, and at main sequence the winds use it while they hold the
   channel. The black-hole outflow producer also writes unspawned_wind_mass; the jets ARE compiled in that
   build, but it deposits onto a special zoom target, which is type 3, so it never addresses a particle this
   function routes to the jet reservoir. Keep that in mind if the special-target type ever changes. */
double* active_unspawned_mass_ptr(int i)
{
#ifdef SINGLE_STAR_FB_JETS
    if(P[i].Type == 5) {
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
        if(P[i].ProtoStellarStage >= 6) {return &P[i].unspawned_wind_mass;} // supernova ejecta, and anything folded in on the way past it. This must not require SINGLE_STAR_FB_WINDS: a jets+SNe build still reaches stage 6, and routing it to the jet reservoir would bank the ejecta in one reservoir and drain the other
#endif
#if defined(SINGLE_STAR_FB_WINDS) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
        if((P[i].ProtoStellarStage == 5) && (P[i].wind_mode == 1)) {return &P[i].unspawned_wind_mass;} // main-sequence winds hold the channel
#endif
        return &P[i].unspawned_jet_mass; // jets are the only, or the currently dominant, discrete-spawn channel
    }
#endif
    return &P[i].unspawned_wind_mass;
}
#endif


/* simple routine that evaluates the target cell mass for the spawning subroutine */
double target_mass_for_wind_spawning(int i)
{
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES == 3) // replace later as needed //
    if(is_particle_a_special_zoom_target(i)) {return 1.e-9/UNIT_MASS_IN_SOLAR;} //
#endif
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4) // replace later as needed //
    if(is_particle_a_special_zoom_target(i)) {return 1.e-6/UNIT_MASS_IN_SOLAR;} //
#endif

#if defined(SNE_NONSINK_SPAWN)
    if(P[i].Type==4) {return 0.5 / UNIT_MASS_IN_SOLAR;} // replace later as needed //
#endif
    
#ifdef SINK_WIND_SPAWN

#if defined(SINGLE_STAR_FB_WINDS) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    if(P[i].Type==5) {
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL) || defined(SINK_SCALE_SPAWNINGMASS_WITH_INITIALMASS) // we specify the value relative to Sink_Formation_Mass
        if((All.Cell_Spawn_Mass_ratio_MS>0.0)&&(P[i].ProtoStellarStage == 5)&&(P[i].wind_mode==1)) {return All.Cell_Spawn_Mass_ratio_MS * P[i].Sink_Formation_Mass;} //use different (probably lower) mass for winds than for jets (will also reduce it for MS jets, but that should be fine)
        else {return All.Sink_outflow_particlemass * P[i].Sink_Formation_Mass;}
#else // we specify the absolute value
        if((P[i].ProtoStellarStage == 5) && (P[i].wind_mode == 1)) {return (All.Cell_Spawn_Mass_ratio_MS > 0) ? All.Cell_Spawn_Mass_ratio_MS : All.Sink_outflow_particlemass;} // the winds hold the discrete-spawn channel: use their specified absolute mass resolution, falling back to the outflow particle mass if it is not set. If the jets hold the channel instead, even at main sequence, fall through below to the jet mass
        else if(P[i].ProtoStellarStage == 6) {return P[i].Sink_Formation_Mass;} // If supernova, use the nominal "average" mass resolution
#endif
    }
#endif // single-star if above 

#if defined(SINK_SCALE_SPAWNINGMASS_WITH_INITIALMASS)
    return All.Sink_outflow_particlemass * P[i].Sink_Formation_Mass;
#else
    return All.Sink_outflow_particlemass;
#endif

#endif // SINK_WIND_SPAWN clause
    return 0; // no well-defined answer, this shouldn't be called in this instance
}


#endif // top-level flag
