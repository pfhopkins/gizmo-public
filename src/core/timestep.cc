#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"
#include "../core/proto.h"
#include "../core/wakeup_sidecar.h"
#include "../mesh/kernel.h"
#if defined(DM_SIDM)
#include "../sidm/sidm_helper_functions.h"
#ifdef GRAIN_COLLISIONS
#include "../solids/grain_helper_functions.h"
#endif
#endif

/*! \file timestep.c
 *  routines for assigning new timesteps
 */
/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO; these
 * modifications include the addition of various timestep criteria, the WAKEUP
 * additions, and various changes of units and variable naming conventions throughout,
 * as well as timestep conditions for all physics and alternative solver options
 * and different timestep schemes entirely.
 */

static double dt_displacement = 0;


/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void find_timesteps(void)
{
    CPU_Step[CPU_MISC] += measure_time();

    int i, bin, binold, prev, next;
    integertime ti_step, ti_step_old, ti_min, ti_stepmax, ti_max;
    double aphys;
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    int special_particle_active_with_this_index[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM], j_specialpartical_counter=0;
    double xyz_local[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM][3], xyz_global[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM][3], special_particle_mass_local[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM]={0}, special_particle_mass_global[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM]={0};
    for(i=0;i<SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM;i++) {special_particle_active_with_this_index[i] = -1; xyz_local[i][0]=xyz_local[i][1]=xyz_local[i][2] = -MAX_REAL_NUMBER;}
#endif

    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin || dt_displacement == 0)
        find_dt_displacement_constraint(All.cf_hubble_a * All.cf_atime * All.cf_atime);

#ifdef DIVBCLEANING_DEDNER
    /* need to calculate the global fastest wave speed to manage the damping terms stably */
    if((All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)||(All.FastestWaveSpeed == 0))
    {
        double fastwavespeed = 0.0;
        double fastwavedecay = 0.0;
        double fac_magnetic_pressure = 1. / All.cf_atime;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(max:fastwavespeed,fastwavedecay)
#endif
        for(i=0;i<NumPart;i++)
        {
            if(P[i].Type==0)
            {
                double vsig2 = 0.5  * fabs(CellP[i].MaxSignalVel); // in v_phys units //
                double vA2_for_vsig = 0; // the magnetic term needs a density to divide by; a cell without one contributes only its sound speed, rather than an infinity that this max-reduction would then spread to every rank
                if(CellP[i].Density > 0) {vA2_for_vsig = fac_magnetic_pressure * CellP[i].Bfield().norm_sq() / CellP[i].Density;}
                double vsig1 = sqrt( CellP[i].effective_soundspeed()*CellP[i].effective_soundspeed() + vA2_for_vsig );
                double vsig0 = DMAX(vsig1,vsig2);

                if(vsig0 > fastwavespeed) fastwavespeed = vsig0; // physical unit
                double hsig0 = P[i].Get_Particle_Size() * All.cf_atime; // physical unit
                if(vsig0/hsig0 > fastwavedecay) fastwavedecay = vsig0 / hsig0; // physical unit
            }
        }
        /* if desired, can just do this by domain; otherwise we use an MPI call over all domains to collect */
        double fastwavespeed_max_glob=fastwavespeed;
        MPI_Allreduce(&fastwavespeed, &fastwavespeed_max_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double fastwavedecay_max_glob=fastwavedecay;
        MPI_Allreduce(&fastwavedecay, &fastwavedecay_max_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        /* now set the variables */
        All.FastestWaveSpeed = fastwavespeed_max_glob;
        All.FastestWaveDecay = fastwavedecay_max_glob;
    }
#endif

#if defined(FORCE_EQUAL_TIMESTEPS) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    ti_max = 0;
    ti_min = TIMEBASE;
    for (int i : ActiveParticleList)
    {
#if defined(FORCE_EQUAL_TIMESTEPS)
        /* get the timestep for this particle, and apply any dilation factor -- we're trying to find the minimum active BIN, not the minimum active timestep in float.
           get_timestep sets the particle's dilation factor, so it must be called in its own statement before the factor is read. */
        integertime ti_step_undilated = get_timestep(i, &aphys, 0);
        ti_step = (integertime)(((double)ti_step_undilated) / timestep_dilation_factor(i, P));
#elif defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
        if(is_particle_a_special_zoom_target(i)==0) {ti_step = P[i].dt_step;} else {ti_step = TIMEBASE;} // set the source particle to have a timestep no more than 4 bins larger than the previous smallest active particle/cell bin timestep
#endif
        if(ti_step < ti_min) {ti_min = ti_step;}
        if(ti_step > ti_max) {ti_max = ti_step;}
    }
    if(ti_min > (dt_displacement / All.Timebase_interval)) {ti_min = (dt_displacement / All.Timebase_interval);}

    ti_step = TIMEBASE;
    while(ti_step > ti_min) {ti_step >>= 1;}
    ti_stepmax = TIMEBASE;
    while(ti_stepmax > ti_max) {ti_stepmax >>= 1;}
    integertime ti_min_glob, ti_max_glob;
    MPI_Allreduce(&ti_step, &ti_min_glob, 1, MPI_TYPE_TIME, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&ti_stepmax, &ti_max_glob, 1, MPI_TYPE_TIME, MPI_MAX, MPI_COMM_WORLD);
#if !defined(FORCE_EQUAL_TIMESTEPS) && defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
#if defined(USE_TIMESTEP_DILATION_FOR_ZOOMS)
    ti_min_glob <<= 2; // 2^N times min timestep - shift to N bins higher
#else
    ti_min_glob <<= 4; // 2^N times min timestep - shift to N bins higher
#endif
    if(ti_min_glob > ti_max_glob) {ti_min_glob = ti_max_glob;}
#endif
#endif


    /* Now assign new timesteps  */
    for (int i : ActiveParticleList)
    {
#ifdef FORCE_EQUAL_TIMESTEPS
        ti_step = ti_min_glob;  /* note that the dilation factor is already applied to ti_min_glob above - re-applying here would double-count it */
#else
        ti_step = get_timestep(i, &aphys, 0);
        ti_step = (integertime)(((double)ti_step) / timestep_dilation_factor(i, P)); /* get_timestep above froze the factor for this particle */
#endif
        
#if !defined(FORCE_EQUAL_TIMESTEPS) && defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
        if(ti_step < 0) {ti_step = ti_min_glob;}
        if(is_particle_a_special_zoom_target(i)) {
            if(ti_step > ti_min_glob) {ti_step = ti_min_glob;}
            if(ti_step > ti_max_glob) {ti_step = ti_max_glob;}
        }
        //if(ti_min_glob > 0) {if(is_particle_a_special_zoom_target(i)) {while(ti_step > ti_min_glob) {ti_step >>= 1;}}} // set this per the above loop to minimum threshold relative to previous steps
#endif
        /* make it a power 2 subdivision */
        ti_min = TIMEBASE;
        while(ti_min > ti_step) {ti_min >>= 1;}
        ti_step = ti_min;
        bin = get_timestep_bin(ti_step);
        binold = P[i].TimeBin;
        if(bin > binold)		/* timestep wants to increase */
        {
            while(TimeBinActive[bin] == 0 && bin > binold) {bin--;}	/* make sure the new step is synchronized */
            ti_step = GET_INTEGERTIME_FROM_TIMEBIN(bin);
        }
        if(All.Ti_Current >= TIMEBASE) {ti_step = 0; bin = 0;} /* we here finish the last timestep. */

        if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
        {
            printf("we are beyond the end of the timeline (task=%d) -- clamping ti_step and requesting controlled stop\n", ThisTask); fflush(stdout);	/* should not happen */
            endrun(90001005);	/* graceful bad-stop; the clamp below keeps ti_step finite; drains at the find_timesteps poll (no per-particle collective) */
            ti_step = TIMEBASE - All.Ti_Current;
            ti_min = TIMEBASE;
            while(ti_min > ti_step) {ti_min >>= 1;}
            ti_step = ti_min;
        }

        if(bin != binold)
        {
            TimeBinCount[binold]--;
            if(P[i].Type == 0)
            {
                TimeBinCountGas[binold]--;
#ifdef GALSF
                TimeBinSfr[binold] -= CellP[i].Sfr;
                TimeBinSfr[bin] += CellP[i].Sfr;
#endif
            }

#ifdef SINK_PARTICLES
            if(P[i].Type == 5)
            {
                TimeBin_Sink_mass[binold] -= P[i].Sink_Mass;
                TimeBin_Sink_dynamicalmass[binold] -= P[i].Mass;
                TimeBin_Sink_Mdot[binold] -= P[i].Sink_Mdot;
                if(P[i].Sink_Mass > 0) {TimeBin_Sink_Medd[binold] -= P[i].Sink_Mdot / P[i].Sink_Mass;}
                TimeBin_Sink_mass[bin] += P[i].Sink_Mass;
                TimeBin_Sink_dynamicalmass[bin] += P[i].Mass;
                TimeBin_Sink_Mdot[bin] += P[i].Sink_Mdot;
                if(P[i].Sink_Mass > 0) {TimeBin_Sink_Medd[bin] += P[i].Sink_Mdot / P[i].Sink_Mass;}
            }
#endif
            prev = PrevInTimeBin[i];
            next = NextInTimeBin[i];

            if(FirstInTimeBin[binold] == i) {FirstInTimeBin[binold] = next;}
            if(LastInTimeBin[binold] == i) {LastInTimeBin[binold] = prev;}
            if(prev >= 0) {NextInTimeBin[prev] = next;}
            if(next >= 0) {PrevInTimeBin[next] = prev;}

            if(TimeBinCount[bin] > 0)
            {
                PrevInTimeBin[i] = LastInTimeBin[bin];
                NextInTimeBin[LastInTimeBin[bin]] = i;
                NextInTimeBin[i] = -1;
                LastInTimeBin[bin] = i;
            }
            else
            {
                FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
                PrevInTimeBin[i] = NextInTimeBin[i] = -1;
            }
            TimeBinCount[bin]++;
            if(P[i].Type == 0) {TimeBinCountGas[bin]++;}
            P[i].TimeBin = bin;
        }

        ti_step_old = P[i].dt_step;
        P[i].Ti_begstep += ti_step_old;
        P[i].dt_step = ti_step;
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
        if(P[i].Type == 5){
            if(All.Ti_Current == 0) { // first timestep
                P[i].dt_since_last_gas_search = get_physical_timestep_from_timebin(P[i].TimeBin, i, P);
                P[i].do_gas_search_this_timestep = 1;
            } else {
                P[i].dt_since_last_gas_search += get_physical_timestep_from_timebin(P[i].TimeBin, i, P);
                if(P[i].dt_since_last_gas_search > 0.49 * get_physical_timestep_from_timebin(P[i].Sink_TimeBinGasNeighbor, i, P)){
                    P[i].do_gas_search_this_timestep = 1;
                } else {P[i].do_gas_search_this_timestep = 0;}
            }
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
	    if(P[i].ProtoStellarStage == 6) {P[i].do_gas_search_this_timestep = 1;} // always do gas search if we're rapidly spawning in new gas shells
#endif
        }
#endif
        
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
        if(is_particle_a_special_zoom_target(i) && P[i].Mass > 0) {xyz_local[j_specialpartical_counter][0]=P[i].Pos[0]; xyz_local[j_specialpartical_counter][1]=P[i].Pos[1]; xyz_local[j_specialpartical_counter][2]=P[i].Pos[2]; special_particle_active_with_this_index[j_specialpartical_counter]=i; special_particle_mass_local[j_specialpartical_counter]=P[i].Mass; j_specialpartical_counter++;} // active on this processor, set
#endif
        
    }

#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    MPI_Allreduce(xyz_local, xyz_global, 3*SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // broadcast the new position of the special particle
    double mass_to_sum_local[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM],  mass_to_sum_global[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM]={0}; // define mass variables for passing
    int k; for(k=0;k<SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM;k++) {mass_to_sum_local[k] = All.Mass_Accreted_By_SpecialParticle[k];}
    MPI_Allreduce(mass_to_sum_local, mass_to_sum_global, SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // broadcast the mass update of the special particle
    for(k=0;k<SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM;k++)
    {
        if(xyz_global[k][0] > -1.e10) { // this indicates that the special particle was active on one task
            All.SpecialParticle_Position_ForRefinement[k] = {xyz_global[k][0], xyz_global[k][1], xyz_global[k][2]}; // variable was updated, update global variable as needed
            if(special_particle_active_with_this_index[k]>=0) {P[special_particle_active_with_this_index[k]].Mass += mass_to_sum_global[k]; special_particle_mass_local[k] += mass_to_sum_global[k];} // the special particle lives here with this id, so we can update it with this mass
            All.Mass_Accreted_By_SpecialParticle[k] = 0; // reset this variable on all processors because we have added it now to the special particle, to conserve mass properly
        }
    }
    MPI_Allreduce(special_particle_mass_local, special_particle_mass_global, SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // broadcast the mass of the special particle
    for(k=0;k<SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM;k++) {if(special_particle_mass_global[k] > 0) {All.Mass_of_SpecialParticle[k] = special_particle_mass_global[k];}} // update the mass of the special particle for everyone to use
    // ???
#endif

#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_TAG_ANCHOR
    update_tag_anchor_refinement_center(); // override slot 0 with the tag-derived anchor (COM or densest tagged cell); the Type-3 scan above leaves slot 0 unset in this mode since there are no Type-3 special particles
#endif


#ifdef PMGRID
    if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = TIMEBASE;
        while(ti_step > (dt_displacement / All.Timebase_interval)) {ti_step >>= 1;}
        if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
        {
            bin = get_timestep_bin(ti_step);
            binold = get_timestep_bin(All.PM_Ti_endstep - All.PM_Ti_begstep);
            while(TimeBinActive[bin] == 0 && bin > binold) {bin--;}	/* make sure the new step is synchronized */
            ti_step = GET_INTEGERTIME_FROM_TIMEBIN(bin);
        }
        if(All.Ti_Current == TIMEBASE) {ti_step = 0;} /* we here finish the last timestep. */
        All.PM_Ti_begstep = All.PM_Ti_endstep;
        All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;
    }
#endif

    process_wake_ups();

    /* Controlled-stop collective. Belongs at the end of find_timesteps
     * because every rank reaches this point exactly once per sync-point,
     * after its full active-particle timestep loop has drained. No MPI
     * inside get_timestep() — ranks complete that loop at different
     * cumulative counts so any per-particle collective would deadlock. */
    gizmo_collect_controlled_stop();

    CPU_Step[CPU_FIND_TIMESTEPS] += measure_time();
}



/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
integertime get_timestep(int p,		/*!< particle index */
                         double *aphys,	/*!< acceleration (physical units) */
                         int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding aphys */ )
{
    double ax, ay, az, ac, csnd = 0, dt = All.MaxSizeTimestep, dt_courant = 0, dt_divv = 0;
    integertime ti_step; int k; k=0;
#ifdef TRANSPORT_SUBCYCLE
    if(P[p].Type == 0) {CellP[p].Transport_Dt_Subcycle = MAX_REAL_NUMBER;}
#endif

#if defined(USE_TIMESTEP_DILATION_FOR_ZOOMS)
    /* freeze this particle's dilation factor for the step we are about to assign. Every later
       conversion of its integer step back to physical time must use this same value, so that the
       physical landing time of the step cannot mutate as the particle moves. Set here at entry so
       it covers every exit from this function. */
    P[p].TimestepDilationFactor = return_timestep_dilation_factor(p, P);
#endif

#ifdef IO_GRADUAL_SNAPSHOT_RESTART // if on the first timestep of a snapshot restart, start at the lowest allowed timestep to minimize any transient effects
    if(RestartFlag == 2 && All.Ti_Current == 0) {return 2;}
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    P[p].SuperTimestepFlag = 0;
    if( (P[p].Type == 5) && P[p].is_in_a_binary ) // candidate: need to decide whether to use super timestepping for binaries
    {
#if (SINGLE_STAR_TIMESTEPPING == 1) // to be conservative, use the semimajor axis, ie. the internal timescale is the orbital period
	    double dt_bin = P[p].Min_Sink_OrbitalTime / (2.*M_PI); // sqrt(a^3/GM) for binary
	    if(0.03*P[p].COM_dt_tidal>dt_bin) {P[p].SuperTimestepFlag=2;
	    } // external timestep is appropriately larger than 'internal' timestep, so use super-timestepping routine
#else // to be more aggressive, use the instantaneous orbital timescale, ie. freefall time from the CURRENT orbital separation. This lets us super step an orbit on the close passages, even when it is affected by tides at apopase
	    double dr = P[p].comp_dx.norm();
	    double dt_bin = sqrt(dr*dr*dr / (All.G * (P[p].Mass + P[p].comp_Mass)));
        if(0.005*P[p].COM_dt_tidal>dt_bin) {P[p].SuperTimestepFlag=2;} // external timestep is appropriately larger than 'internal' timestep, so use super-timestepping routine [constant here stricter for more aggressive routine]
#endif
    }
#endif

    
    
#if defined(SPECIAL_POINT_MOTION)
    {
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
        if(P[p].Type != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
#endif
        {
            Vec3<double> acc = P[p].GravAccel * All.cf_a2inv;
#ifdef PMGRID
            acc += P[p].GravPM * All.cf_a2inv;
#endif
            if(P[p].Type==0) {
                acc += CellP[p].HydroAccel;
#ifdef TURB_DRIVING
                acc += CellP[p].TurbAccel;
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
                acc += CellP[p].Rad_Accel;
#endif
            }
            P[p].Acc_Total_PrevStep = acc;
        }
    }
#endif

    
    if(flag == 0)
    {
        ax = All.cf_a2inv * P[p].GravAccel[0];
        ay = All.cf_a2inv * P[p].GravAccel[1];
        az = All.cf_a2inv * P[p].GravAccel[2];
#ifdef PMGRID
        ax += All.cf_a2inv * P[p].GravPM[0];
        ay += All.cf_a2inv * P[p].GravPM[1];
        az += All.cf_a2inv * P[p].GravPM[2];
#endif
        
#if defined(TIDAL_TIMESTEP_CRITERION)
#if defined(RT_USE_GRAVTREE) && !defined(SINGLE_STAR_FB_RT_HEATING)
        if(P[p].Type>0) // strictly this is better for accuracy, but not necessary
#endif
        ax = ay = az = 0.0; // we're getting our gravitational timestep criterion from the tidal tensor, but still want to do the accel criterion for other forces
#endif

        if(P[p].Type == 0)
        {
            ax += CellP[p].HydroAccel[0];
            ay += CellP[p].HydroAccel[1];
            az += CellP[p].HydroAccel[2];
#ifdef TURB_DRIVING
            ax += CellP[p].TurbAccel[0];
            ay += CellP[p].TurbAccel[1];
            az += CellP[p].TurbAccel[2];
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            ax += CellP[p].Rad_Accel[0];
            ay += CellP[p].Rad_Accel[1];
            az += CellP[p].Rad_Accel[2];
#endif
        }

#if defined(CBE_INTEGRATOR)
        if(CBE_INTEGRATOR_DOES_TYPE(P[p].Type))
        {   /* CBE moment-flux acceleration enters the accel timestep like HydroAccel does for gas;
             * reduces to the hydro acceleration when all bases are identical isotropic Gaussians */
            double a_cbe[3]; cbe_particle_moment_accel(p, a_cbe);
            ax += a_cbe[0]; ay += a_cbe[1]; az += a_cbe[2];
        }
#endif

        ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
        *aphys = ac;
    }
    else
    {ac = *aphys;}

    if(ac == 0) {ac = 1.0e-30;}


    if(flag > 0)
    {
        /* this is the non-standard mode; use timestep to get the maximum acceleration tolerated */
        dt = flag * unit_integertime_in_physical(p, P); /* convert dloga to physical timestep  */
        ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * ForceSoftening_KernelRadius(p) / (dt * dt);
        *aphys = ac;
        return flag;
    }
    {double h_for_accel_dt = KERNEL_CORE_SIZE * ForceSoftening_KernelRadius(p);
#if defined(CBE_INTEGRATOR)
    if(CBE_INTEGRATOR_DOES_TYPE(P[p].Type)) {h_for_accel_dt = KERNEL_CORE_SIZE * Get_Particle_Size_AGS(p);} /* AGS particle size, not force-softening, like gas */
#endif
#ifdef GRAIN_FLUID
    if(((1 << P[p].Type) & (GRAIN_PTYPES)) && (h_for_accel_dt <= 0)) {h_for_accel_dt = P[p].Get_Particle_Size() * All.cf_atime;} /* for grain particles without gravity, use the inter-particle spacing as the characteristic length scale */
#endif
    dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * h_for_accel_dt / ac);}

#if (defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)) && defined(GALSF) && defined(GALSF_FB_MECHANICAL)
    if(is_galsf_stellar_candidate_type(P[p].Type, All.ComovingIntegrationOn) && (P[p].Mass>0))
    {
        if((All.ComovingIntegrationOn)) // sort of a hack here, but acceptable in applications
        {
            double h_min = All.ForceSoftening[P[p].Type], ags_h = DMIN(DMAX(P[p].KernelRadius, h_min), 10.*h_min);
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            ags_h = DMIN(DMAX(P[p].AGS_KernelRadius , DMAX(P[p].KernelRadius, h_min)) , DMAX(100.*h_min, 10.*P[p].AGS_KernelRadius));
#endif
            dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime  * KERNEL_CORE_SIZE * ags_h / ac);
        }
    }
#endif


#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    double tidal_mag = P[p].tidal_tensorps.frobenius_norm(); // can estimate time derivative here, via: dt_ttmag = (tidal_mag-P[p].tidal_tensor_mag_prev) / get_particle_timestep_in_physical(p);
    double dt_tidalsoft = All.CourantFac * NUMDIMS * DMAX(DMAX(get_particle_timestep_in_physical(p, P), dt), All.MinSizeTimestep) * (tidal_mag+P[p].tidal_tensor_mag_prev) / (fabs(tidal_mag-P[p].tidal_tensor_mag_prev) + MIN_REAL_NUMBER);
    if(((1 << P[p].Type) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)) && (P[p].tidal_tensor_mag_prev>0 && All.Time>All.TimeBegin)) {dt = DMIN(dt, dt_tidalsoft);} // use as a timestep criterion for tidal-ags-active particles
    P[p].tidal_tensor_mag_prev = tidal_mag; // save it (overwriting previous value)
    {
        double tt2 = P[p].tidal_tensorps.frobenius_norm_sq();
        double tracett = P[p].tidal_tensorps.trace(); /* compute numbers needed below */
        double H_eff = ForceSoftening_KernelRadius(p); /* get value to calculate H we need to use in the equations below */
        if(tidal_mag > 0) {P[p].tidal_zeta *= -All.G*(H_eff-All.ForceSoftening[P[p].Type])/(2.*NUMDIMS*tt2 - 0 * 8.*M_PI*tracett*(All.G*P[p].Mass/pow(H_eff,NUMDIMS)));} else {P[p].tidal_zeta = 0;}
        P[p].tidal_tensorps_prevstep = P[p].tidal_zeta * P[p].tidal_tensorps; /* save for next iteration in gravtree */
    }
#endif


/* Safety factor on the two-body sink timestep. A Hermite-eligible particle divides it back out: the
   4th-order integrator tolerates the longer step that these 2nd-order-calibrated coefficients shorten. */
#ifndef SINK_TIMESTEP_SAFETY_FACTOR
#define SINK_TIMESTEP_SAFETY_FACTOR (0.3)
#endif

#ifdef TIDAL_TIMESTEP_CRITERION // tidal criterion obtains the same energy error in an optimally-softened Plummer sphere over ~100 crossing times as the Power 2003 criterion
    double tidal_mag_dt = P[p].tidal_tensorps.frobenius_norm_sq();
    double dt_tidal = sqrt(All.ErrTolIntAccuracy / (All.cf_a3inv * sqrt(tidal_mag_dt / 6))); // recovers sqrt(eta) * tdyn for a Keplerian potential
    if(P[p].Type == 0) {dt_tidal = DMIN(sqrt(All.ErrTolIntAccuracy/(All.G*CellP[p].Density*All.cf_a3inv)), dt_tidal);} // gas self-gravity timescale as a bare minimum
    if(All.ComovingIntegrationOn){ // floor to the dynamical time of the universe
        double rho0 = (H0_CGS*H0_CGS*(3./(8.*M_PI*GRAVITY_G_CGS))*All.cf_a3inv / UNIT_DENSITY_IN_CGS);
        dt_tidal = DMIN(dt_tidal, sqrt(All.ErrTolIntAccuracy / (All.G * rho0)));
    } 
#ifdef ADAPTIVE_TREEFORCE_UPDATE
    P[p].tdyn_step_for_treeforce = dt_tidal; // hang onto this to decide how frequently to update the treeforce
#endif
#ifdef HERMITE_INTEGRATION
    /* divide out the second-order margin, as the two-body criterion below does. After the tree-update
       cadence above, which wants the unscaled dynamical estimate rather than the integrator's step. */
    if(eligible_for_hermite(p, P)) {dt_tidal /= SINK_TIMESTEP_SAFETY_FACTOR;}
#endif
    
#if (SINGLE_STAR_TIMESTEPPING > 0)
    if(P[p].SuperTimestepFlag>=2) {dt_tidal = sqrt(2*All.ErrTolIntAccuracy) * P[p].COM_dt_tidal;}
#endif
    dt=DMIN(dt,dt_tidal);
#endif

#ifdef SINGLE_STAR_TIMESTEPPING // this ensures that binaries advance in lock-step, which gives superior conservation
    if(P[p].Type == 5)
    {
        double dt_2body = sqrt(2*All.ErrTolIntAccuracy) * SINK_TIMESTEP_SAFETY_FACTOR / (1./P[p].Min_Sink_Approach_Time + 1./P[p].Min_Sink_Freefall_time); // timestep is harmonic mean of freefall and approach time
#ifdef HERMITE_INTEGRATION
        if(eligible_for_hermite(p, P)) dt_2body /= SINK_TIMESTEP_SAFETY_FACTOR;
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    	if(P[p].is_in_a_binary && (P[p].SuperTimestepFlag >= 2)) //binary candidate or a confirmed binary
	    {    // First we need to construct the same 2-body timescale as above, but from the binary parameters. If this is longer than the above, there is another star that is requiring us to
	         // take a short timestep, so we better not super-timestep otherwise we risk messing up that star's integration. But if it is consistent with the above, then we can safely super-timestep
	        double Mtot=P[p].comp_Mass+P[p].Mass, dr=P[p].comp_dx.norm_sq(), dv=P[p].comp_dv.norm_sq(), dv_dot_dx=dot(P[p].comp_dx,P[p].comp_dv), binary_dt_2body=0;
            double r_effective = KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER * ForceSoftening_KernelRadius(p); // plummer-equivalent softening
	        dr += r_effective*r_effective; // add in quadrature for simple softening estimate
            dr=sqrt(dr); if(dv>0) {dv=sqrt(dv);} else {dv=0;}
            double dt_2body_base = 1/(1./P[p].Min_Sink_Approach_Time + 1./P[p].Min_Sink_Freefall_time); // timestep is harmonic mean of freefall and approach time
	        binary_dt_2body = 1. / (dv / dr + sqrt(All.G * Mtot / (dr*dr*dr)));
	        if(fabs(binary_dt_2body - dt_2body_base)/dt_2body_base < 1e-2)
	        { // If consistent with the binary parameters, we choose a super-timestep that gives ~constant number of timesteps per orbit
                double SUPERTIMESTEPPING_NUM_STEPS_PER_ORBIT = 50;
                dt_2body = 2.*M_PI / SUPERTIMESTEPPING_NUM_STEPS_PER_ORBIT * (binary_dt_2body*2); // orbital frequency is |dr x dv| / r^2, so timestep will be inverse to this
	        } else {P[p].SuperTimestepFlag = 0;}  // we still have to take a proper short N-body integration timestep due to a third body whose approach requires careful integration, so no super timestepping is possible
	    }
#endif
        dt = DMIN(dt, dt_2body);
#ifdef HERMITE_INTEGRATION
        if(eligible_for_hermite(p, P)) dt *= 1.4; // gives 10^-6 energy error per orbit for a 0.9 eccentricity binary
#endif
    }
#if defined(SINGLE_STAR_FB_TIMESTEPLIMIT) && !defined(SELFGRAVITY_OFF)
    if(P[p].Type == 0) {dt = DMIN(dt, 0.5 * All.CourantFac * DMIN(P[p].Min_Sink_FeedbackTime, P[p].Min_Sink_Approach_Time));}
#endif    
#endif // SINGLE_STAR_TIMESTEPPING

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    /* Particles using the adaptive kernel radius need a kernel-size signal-speed
       timestep. Which TYPES require it depends on the active physics, so collect
       the OR of all module triggers, then apply the single AGS Courant criterion.
       Gate on AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE (the SSOT for AGS_KernelRadius/
       AGS_vsig existence), NOT on ADAPTIVE_GRAVSOFT_FORALL alone -- that flag is
       only one of several triggers (SIDM / fuzzy / CBE also need this limiter). */
    int need_agscfl = 0;
#if defined(CBE_INTEGRATOR)
    int need_cbe_agscfl = 0;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(((1 << P[p].Type) & (ADAPTIVE_GRAVSOFT_FORALL)) && (P[p].Type > 0)) {need_agscfl = 1;}
#endif
#ifdef DM_SIDM
    if((1 << P[p].Type) & (DM_SIDM)) {need_agscfl = 1;}
#endif
#if defined(CBE_INTEGRATOR)
    if(CBE_INTEGRATOR_DOES_TYPE(P[p].Type))
    {
        need_agscfl = 1;
        need_cbe_agscfl = 1; /* CBE-moment-flux particle: gets the stricter factor below */
    }
#elif defined(DM_FUZZY)
    if(P[p].Type == 1) {need_agscfl = 1;}
#endif
    if(need_agscfl)
    {
        /* make sure smoothing length of non-gas particles doesn't change too much in one timestep */
        double dt_divv = 0.1 / (MIN_REAL_NUMBER + All.cf_a2inv*fabs(P[p].Particle_DivVel)); // with new integration accuracy in gravtree, we may not need to be super-conservative here. old code used pre-factor 0.25 here, see if we can get away with the larger value which is standard for gas below
        if(dt_divv < dt) {dt = dt_divv;}
        double dt_cour = 2. * All.CourantFac * (Get_Particle_Size_AGS(p)*All.cf_atime) / (MIN_REAL_NUMBER + 0.5*P[p].AGS_vsig); // can be generous here, really the signal velocity isn't that important in the collisionless case, but it is important with some of the physics above //
#if defined(CBE_INTEGRATOR)
#if defined(CBE_INTEGRATOR_RP_GAUSSIAN)
        if(need_cbe_agscfl) {dt_cour *= 0.25;} // stricter criterion for CBE moment fluxes (CBE particles only, not other AGS-CFL types) //
#else
        if(need_cbe_agscfl) {dt_cour *= 0.125;} // as above, tighter still: the top-hat face fluxes also transport second moment, through their gomega and pstress terms, and the compact-support signal speed |u|+c_x is not enough margin for those on its own //
#endif
        if(need_cbe_agscfl)
        {   /* CBE mass-depletion criterion: cap the per-basis fractional mass change per step.
             * m_eff-floor: regularize the softened basis mass as m_eff = max(m_b,
             * CBEMassEffFloor*m_cell), so a near-empty placeholder/free-slot basis
             * (m_b ~ eps*m_cell) receiving a cell-scale flux cannot force an absurdly small step,
             * while an occupied basis (m_b > f_floor*m_cell) still constrains smoothly as it gains
             * mass. Continuous (max, not a hard gate) -> no threshold chatter. Timestep-only (does
             * NOT touch the flux/update). */
            for(int b = 0; b < CBE_INTEGRATOR_NBASIS; b++)
            {
                const double m_b = P[p].CBE_basis_moments[b][0];
                if(!(m_b > 0)) continue;
                const double m_eff = DMAX(m_b, All.CBEMassEffFloor * P[p].Mass);
                const double rate_m = DMAX(fabs(P[p].CBE_basis_out_rate_dt[b][0]), fabs(P[p].CBE_basis_moments_dt[b][0]));
                if(rate_m > MIN_REAL_NUMBER) {double dt_m = 0.1 * m_eff / rate_m; if(dt_m < dt) {dt = dt_m;}}
            }
        }
#endif
        if(dt_cour < dt) {dt = dt_cour;}
    }
#endif


#ifdef DM_FUZZY
    if((P[p].Type > 0) && (P[p].AGS_Density > 0))
    {
        /* fuzzy DM admits longitudinal waves with group velocity =(hbar/m_dm)*k, so need a courant criterion, but because of scaling with k (like diffusion), timestep is quadratic in resolution */
        double L_particle_ags_x = Get_Particle_Size_AGS(p) * All.cf_atime;
        double dt_cour_ags_fuzzy = 0.25 * (L_particle_ags_x*L_particle_ags_x) / All.ScalarField_hbar_over_mass; // wavespeed of resolve-able waves
        if(dt_cour_ags_fuzzy < dt) {dt = dt_cour_ags_fuzzy;}
        dt_cour_ags_fuzzy = 0.25 * L_particle_ags_x / sqrt(MIN_REAL_NUMBER + (10./9.)*P[p].AGS_Numerical_QuantumPotential/P[p].Mass); // wavespeed based on 'stored' sub-grid energy [can get comparable]
        if(dt_cour_ags_fuzzy < dt) {dt = dt_cour_ags_fuzzy;}
    }
#endif


#ifdef DO_FLUID_ALTSPECIES_DRAG_CALCULATION
    if(IS_PARTICLE_DRAGVALID(P[p].Type, P[p].FluidType))
    {
        /* a grain has no cell of its own, so the surrounding gas temperature and ionized fraction are carried onto
           it and the composition estimated from them */
        csnd = GAMMA_DEFAULT * P[p].Gas_Temperature
               / (molecular_weight_estimator_for_gas_around_grain(P[p].Gas_Temperature, P[p].Gas_fion) * U_TO_TEMP_UNITS);
        csnd += (P[p].Gas_Velocity - P[p].Vel).norm_sq();
#if defined(DO_FLUID_DRAG_CALCULATION_WITHBFIELDS)
        csnd += P[p].Gas_B.norm_sq() / (2.0 * P[p].Gas_Density + MIN_REAL_NUMBER);
#endif
        csnd = sqrt(csnd);
        double L_particle = P[p].Get_Particle_Size();
        dt_courant = 0.5 * All.CourantFac * (L_particle*All.cf_atime) / csnd;
#if defined(GRAIN_BACKREACTION)
        if(6.*P[p].Grain_AccelTimeMin < dt_courant) {dt_courant = 6.*P[p].Grain_AccelTimeMin;}
#endif
#if defined(GRAIN_LORENTZFORCE) && defined(GRAIN_RDI_TESTPROBLEM)
        if(All.Grain_Charge_Parameter != 0) {double bmag = P[p].Gas_B.norm_sq();
            if(bmag>0) {double dt_gyro = 1. / ((All.Grain_Charge_Parameter*sqrt(1.)/((All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS)*(All.Grain_Size_Max/UNIT_LENGTH_IN_CGS))) * DMIN(100.,pow(All.Grain_Size_Max/P[p].Grain_Size,2)) * sqrt(bmag)); if(dt_gyro>0 && dt_gyro<dt_courant) {dt_courant=dt_gyro;}}} /* this gives t_Lorentz in code units; sqrt[1] reflects expected unity mean density definition, hard-coded for rdi testproblem options here */
#endif
#ifdef PIC_MHD
        if(P[p].MHD_PIC_SubType>=3)
        {
            double lorentz_units = UNIT_B_IN_GAUSS * UNIT_VEL_IN_CGS * (ELECTRONCHARGE_CGS/(PROTONMASS_CGS*C_LIGHT_CGS)) / (UNIT_VEL_IN_CGS/UNIT_TIME_IN_CGS); // code velocity to CGS and B to Gauss, times base units e/(mp*c), then convert 'back' to code-units acceleration
            double reduced_C = PIC_SPEEDOFLIGHT_REDUCTION * C_LIGHT_CODE, charge_to_mass_ratio_dimensionless = All.PIC_Charge_to_Mass_Ratio;
#ifdef PIC_MHD_NEW_RSOL_METHOD
            lorentz_units *= PIC_SPEEDOFLIGHT_REDUCTION; // the rsol enters by slowing down the forces here, acts as a unit shift for time
#endif
            double B2=(P[p].Gas_B * All.cf_a2inv).norm_sq(), beta2=(P[p].Vel * (1.0/(All.cf_atime*reduced_C))).norm_sq(); /* get magnitude and unit vector for B, and vector beta [-true- beta here] */
            double gamma_lorentz = 1./sqrt(DMAX(1.-DMAX(DMIN(beta2,1.),0.),MIN_REAL_NUMBER)); // calculate lorentz factor (with safety factors included to prevent accidental nan here //
            double dt_courant_pic = 0.5 / ((charge_to_mass_ratio_dimensionless/gamma_lorentz) * sqrt(B2) * lorentz_units); /* dt = 0.5/omega_gyro*/
            if(dt_courant_pic < dt_courant) dt_courant = dt_courant_pic;
        }
#endif
        if(dt_courant < dt) dt = dt_courant;
    }
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
    if(P[p].Type>-1) {double dt_inj = 0.1 * P[p].KernelRadius / C_LIGHT_CODE_REDUCED; if(P[p].Type==4) {dt_inj*=0.25;} if(dt_inj < dt) {dt = dt_inj;}}
#endif
#endif


    if((P[p].Type == 0) && (P[p].Mass > 0))
        {
            csnd = 0.5 * CellP[p].MaxSignalVel ;
            double L_particle = P[p].Get_Particle_Size();
            dt_courant = All.CourantFac * (L_particle*All.cf_atime) / csnd;
#if defined(SINK_WIND_SPAWN) && !defined(SINK_RIAF_SUBEDDINGTON_MODEL)
            if(P[p].ID == All.SpawnedWindCellID) {dt_courant *= 0.5;} // be more careful if this is a spawned-in gas cell
#endif
            if(dt_courant < dt) dt = dt_courant;

#ifdef MHD_BATTERY_MECHANISMS
            /* The battery builds a field out of nothing, at a rate set by the thermodynamic
               gradients rather than by anything the hydro timestep knows about, so the Courant
               condition does not see it. Resolve the time the source takes to change the field
               it is building. The floor on the numerator matters: while the field grows linearly
               from zero, |B|/|dB/dt| IS the elapsed time, so the bare ratio would drive the step
               to zero at the start of a run for a field far too weak to matter. Adding a small
               dimensionless fraction of the thermal pressure stops that, and is negligible once
               the field carries any dynamical weight. */
            if(CellP[p].DtB_battery_magnitude > MIN_REAL_NUMBER)
            {
                const double eta_battery = 0.03;      /* fraction of the field-doubling time per step */
                const double eps_battery = 1.0e-4;    /* dimensionless floor, in units of the thermal pressure */
                const double b_sq = (CellP[p].Bfield() * All.cf_a2inv).norm_sq(); /* physical, matching the stored rate */
                const double p_thermal = CellP[p].Pressure * All.cf_a3inv;
                const double dt_battery = eta_battery * sqrt(b_sq + eps_battery * p_thermal)
                                          / CellP[p].DtB_battery_magnitude;
                if(dt_battery < dt) {dt = dt_battery;}
            }
#endif

            double dt_prefac_diffusion;
            dt_prefac_diffusion = 0.5;
#if (defined(GALSF) || defined(DIFFUSION_OPTIMIZERS)) && !defined(MHD_NON_IDEAL)
            dt_prefac_diffusion = 1.8;
#endif
#ifdef SUPER_TIMESTEP_DIFFUSION
            double dt_superstep_explicit = 1.e10 * dt;
#endif



#ifdef CONDUCTION
            if(CellP[p].Kappa_Conduction > 0) /* no conductivity means no conduction constraint at all; without this test the regularizer below stands in for the zero and returns a finite limit that is not physical. A real conductivity at low density genuinely does demand a short step, so that case is left alone */
            {
                double L_cond_inv = sqrt(CellP[p].Gradients.InternalEnergy[0]*CellP[p].Gradients.InternalEnergy[0] +
                                         CellP[p].Gradients.InternalEnergy[1]*CellP[p].Gradients.InternalEnergy[1] +
                                         CellP[p].Gradients.InternalEnergy[2]*CellP[p].Gradients.InternalEnergy[2]) / CellP[p].InternalEnergy;
                double L_cond = DMAX(L_particle , 1./(L_cond_inv + 1./L_particle)) * All.cf_atime;
                double dt_conduction = dt_prefac_diffusion * L_cond*L_cond / (MIN_REAL_NUMBER + CellP[p].Kappa_Conduction);
                // since we use CONDUCTIVITIES, not DIFFUSIVITIES, we need to add a power of density to get the right units //
                dt_conduction *= CellP[p].Density * All.cf_a3inv;
#ifdef SUPER_TIMESTEP_DIFFUSION
                if(dt_conduction < dt_superstep_explicit) dt_superstep_explicit = dt_conduction; // explicit time-step
                double dt_advective = dt_conduction * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_cond_inv))*All.cf_atime / L_cond);
                if(dt_advective < dt) dt = dt_advective; // 'advective' timestep: needed to limit super-stepping
#else
                if(dt_conduction < dt) dt = dt_conduction; // normal explicit time-step
#endif
            }
#endif


/* No battery-specific CFL: the battery EMF is a SOURCE, not a self-feedback
   growth mode. Pure source dB/dt = const integrates as a linear ramp -- no
   intrinsic stability constraint. Once B grows to dynamically-relevant
   values, the existing MaxSignalVel-driven CFL just below picks up the
   Alfven speed and limits dt naturally. An energy-fraction limiter on the
   battery contribution to dB/dt is applied per-cell in
   hydro_toplevel.cc::out2particle_hydra, which is the right place
   (limits |dE_mag|/dt, not |B|/dt). */

#ifdef MHD_NON_IDEAL
            {
                double b_grad = 0, b_mag = 0;
                for(int k=0;k<3;k++)
                {
                    b_grad += CellP[p].Gradients.B[k].norm_sq();
                    b_mag += CellP[p].Bfield_component(k) * CellP[p].Bfield_component(k);
                }
                double L_cond_inv = MIN_REAL_NUMBER + sqrt(b_grad / (MIN_REAL_NUMBER + b_mag));
                double L_cond = DMAX(0.5*L_particle , DMIN(L_particle , 1./(L_cond_inv + 1./L_particle))) * All.cf_atime;
                L_cond = DMIN( L_particle , DMAX(1./L_cond_inv, 0.5*L_particle) ) * All.cf_atime; // more conservative estimator - may be needed sometimes to deal accurately with steep local gradients //
                double diff_coeff = fabs(CellP[p].Eta_MHD_OhmicResistivity_Coeff) + fabs(CellP[p].Eta_MHD_HallEffect_Coeff) + fabs(CellP[p].Eta_MHD_AmbiPolarDiffusion_Coeff);
                double dt_conduction = dt_prefac_diffusion * L_cond*L_cond / (MIN_REAL_NUMBER + diff_coeff);
#ifdef SUPER_TIMESTEP_DIFFUSION
                if(dt_conduction < dt_superstep_explicit) dt_superstep_explicit = dt_conduction; // explicit time-step
                double dt_advective = dt_conduction * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_cond_inv))*All.cf_atime / L_cond);
                if(dt_advective < dt) dt = dt_advective; // 'advective' timestep: needed to limit super-stepping
#else
                if(dt_conduction < dt) dt = dt_conduction; // normal explicit time-step
#endif
            }
#endif


#ifdef COSMIC_RAY_FLUID
            int k_CRegy;
            for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
            {
                if(Get_Gas_CosmicRayPressure(p, k_CRegy, CellP) > 1.0e-20)
                {
                    int explicit_timestep_on, cr_diffusion_opt = 1;
                    double CRPressureGradScaleLength = Get_CosmicRayGradientLength(p,k_CRegy, P, CellP);
                    double L_cr_weak; L_cr_weak = CRPressureGradScaleLength;
                    double kappa_cr_eff = fabs(CellP[p].CosmicRayDiffusionCoeff[k_CRegy]);
                    kappa_cr_eff *= cosmicrayfluid_rsol_corrfac(k_CRegy); // account for RSOL factor as it actually appears in the flux eqn in code units with this RSOL form
                    double L_cr_strong = DMAX(L_particle*All.cf_atime , 1./(1./CRPressureGradScaleLength + 1./(L_particle*All.cf_atime)));
                    double coeff_inv = 0.67 * L_cr_strong * dt_prefac_diffusion / (1.e-33 + kappa_cr_eff * (GAMMA_COSMICRAY(k_CRegy)-1.));
                    double dt_conduction =  L_cr_strong * coeff_inv; /* true diffusion requires the stronger timestep criterion be applied */
                    explicit_timestep_on = 1;
#if (CRFLUID_DIFFUSION_MODEL < 0)
                    dt_conduction = L_cr_weak * coeff_inv; /* streaming allows weaker timestep criterion because it's really an advection equation */
                    explicit_timestep_on = 0;
#endif
#ifdef GALSF
                    /* for multi-physics problems, we will use a more aggressive timestep criterion
                     based on whether or not the cosmic ray physics are relevant for what we are modeling */
                    if((CellP[p].CosmicRayEnergy[k_CRegy]==0)||(CellP[p].DtCosmicRayEnergy[k_CRegy]==0))
                    {
                        dt_conduction = 10. * dt;
                    } else {
                        double delta_cr = dt_conduction*fabs(CellP[p].DtCosmicRayEnergy[k_CRegy]);
                        double dL_cr = CRPressureGradScaleLength / (L_particle*All.cf_atime);
                        double thres_dL = 2., thres_egy = 1.e-3;
                        if(cr_diffusion_opt==1) {thres_dL = 1.; thres_egy = 1.e-2;}
                        if((dL_cr > thres_dL) || (delta_cr < thres_egy*CellP[p].CosmicRayEnergy[k_CRegy]))
                        {
                            double dt_weak = DMIN(L_cr_weak*coeff_inv , (delta_cr + 1.e-4*CellP[p].CosmicRayEnergy[k_CRegy])/fabs(CellP[p].DtCosmicRayEnergy[k_CRegy]));
                            if((dL_cr > thres_dL+1.) && (delta_cr < 0.1*thres_egy*CellP[p].CosmicRayEnergy[k_CRegy])) {dt_conduction = dt_weak; explicit_timestep_on = 0;}
                        }
                    }
#endif
#ifdef SUPER_TIMESTEP_DIFFUSION
                    if(explicit_timestep_on==1)
                    {
                        if(dt_prefac_diffusion > 1) {dt_conduction *= 0.5;}
                        if(dt_conduction < dt_superstep_explicit) dt_superstep_explicit = dt_conduction; // explicit time-step
                        double dt_advective = dt_conduction * DMAX(1 , DMAX(L_cr_strong,L_cr_weak)/L_cr_strong);
                        if(dt_advective < dt) dt = dt_advective; // 'advective' timestep: needed to limit super-stepping
                    } else {
                        if(dt_conduction < dt) dt = dt_conduction; // this is an advective timestep and super-stepping doesn't apply
                    }
#else
                    double cr_m1_speed = CRFLUID_REDUCED_C_CODE(k_CRegy); // pull for use below
                    if(cr_diffusion_opt==1)
                    {
                        if(CellP[p].CosmicRayEnergy[k_CRegy] > 0)
                        {
                            double cr_speed = cr_m1_speed;
                            //double crv=0; int k; for(k=0;k<3;k++) {crv+=CellP[p].CosmicRayFlux[k_CRegy][k]*CellP[p].CosmicRayFlux[k_CRegy][k];} if(crv > 0) {crv = sqrt(crv) / CellP[p].CosmicRayEnergy[k_CRegy];}
                            cr_speed = DMAX( DMIN(cr_m1_speed , CellP[p].MaxSignalVel) , DMIN(cr_m1_speed , kappa_cr_eff/(P[p].Get_Particle_Size()*All.cf_atime))); // default to min of free-streaming/diffusion speed
                            double dt_courant_CR = 0.4 * (L_particle*All.cf_atime) / cr_speed;
                            dt_conduction = dt_courant_CR; // per TK, strictly enforce this timestep //
                        } else {dt_conduction=10.*dt;}
                    } else {
                        double dt_courant_CR = 0.4 * (L_particle*All.cf_atime) / cr_m1_speed;
                        dt_conduction = dt_courant_CR; // per TK, strictly enforce this timestep //
                    }
                    if(dt_conduction < dt) {
#ifdef TRANSPORT_SUBCYCLE
                        CellP[p].Transport_Dt_Subcycle = DMIN(CellP[p].Transport_Dt_Subcycle, dt_conduction);
                        double dt_max_hydro = TRANSPORT_SUBCYCLE * dt_conduction;
                        if(dt_max_hydro < dt) {dt = dt_max_hydro;}
#else
                        dt = dt_conduction; // normal explicit time-step
#endif
                    }
#endif
                }
            }
#endif


#if defined(RADTRANSFER)
            {
                double dt_rad = 1.e10 * dt; // make some ridiculously large number here
                    
                /* first check if we are using an explicit diffusion-type solver (FLD, OTVET). need to consider the standard diffusive timestep, which we calculate below */
#if (defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION)) && defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_EVOLVE_FLUX) /* for explicit diffusion, we include the usual second-order diffusion timestep */
                int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++)
                {
#if defined(RT_SOLVER_EXPLICIT) // explicit solver -- need diffusion timestep //
                    double gradETmag = CellP[p].Gradients.Rad_E_gamma_ET[kf].norm_sq();
                    double L_ETgrad_inv = sqrt(gradETmag) / (1.e-37 + CellP[p].Rad_E_gamma[kf] * CellP[p].Density/P[p].Mass);
                    double L_RT_diffusion = DMIN(L_particle , 1./(3.*L_ETgrad_inv)) * All.cf_atime;
                    double dt_rt_diffusion = dt_prefac_diffusion * L_RT_diffusion*L_RT_diffusion / (MIN_REAL_NUMBER + rt_diffusion_coefficient(p,kf, CellP));
                    double dt_advective = dt_rt_diffusion * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_ETgrad_inv))*All.cf_atime / L_RT_diffusion);
                    double dt_rt_work = All.CourantFac * DMIN( L_RT_diffusion / csnd , L_particle*All.cf_atime / ((2./3.)*sqrt(CellP[p].Rad_E_gamma[kf]/P[p].Mass)) ); /* time-step related to radiation work, radiation soundspeed, relevant in strongly-coupled limit */
#ifdef RT_FLUXLIMITER /* if we are flux-limited, we can account for the flux limiter making the timestep advective */
                    if(dt_advective > dt_rt_diffusion) {dt_rt_diffusion *= 1. + (1.-CellP[p].Rad_Flux_Limiter[kf]) * DMAX(0,(dt_advective/dt_rt_diffusion-1.));}
                    dt_advective = All.CourantFac * 0.5 * (L_particle*All.cf_atime) / C_LIGHT_CODE_REDUCED;
                    dt_rt_diffusion = DMAX(dt_rt_diffusion, dt_advective);
                    dt_rt_work /= MIN_REAL_NUMBER + CellP[p].Rad_Flux_Limiter[kf];
                    if((CellP[p].Rad_Flux_Limiter[kf] <= 0)||(dt_rt_diffusion<=0)) {dt_rt_diffusion = 1.e9 * dt;}
#endif
                    if((CellP[p].Rad_E_gamma[kf] <= MIN_REAL_NUMBER) || (CellP[p].Rad_E_gamma_Pred[kf] <= MIN_REAL_NUMBER)) {dt_rt_diffusion = dt_advective;} /* if the radiation is totally negligible, just use an advective timestep instead */
#ifdef SUPER_TIMESTEP_DIFFUSION /* if super-timestepping, limit the -super- step with the advective step, since you can still get inaccuracies if this is not respected */
                    if(dt_rt_diffusion < dt_superstep_explicit) dt_superstep_explicit = dt_rt_diffusion; // explicit time-step
                    dt_advective = dt_rt_diffusion * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_ETgrad_inv))*All.cf_atime / L_RT_diffusion);
                    if(dt_advective < dt_rad) dt_rad = dt_advective; // 'advective' timestep: needed to limit super-stepping
#else
                    if(dt_rt_diffusion < dt_rad) dt_rad = dt_rt_diffusion; // normal explicit time-step
                    if(dt_rt_work < dt_rad) {dt_rad = dt_rt_work;} // normal explicit time-step
#endif
#endif // explicit-solver check
#if defined(RT_RAD_PRESSURE_FORCES) // -regardless- of if using an explicit solver, here the acceleration isn't saved to Rad_Accel so we calculate that timestep constraint
                    double gradErad = CellP[p].Gradients.Rad_E_gamma_ET[kf].norm_sq();
                    double radacc = 0; // radiation acceleration for a timestep criterion; needs a density to divide by, and a cell without one imposes no constraint here
                    if(CellP[p].Density > 0) {radacc = CellP[p].flux_limiter(kf) * (sqrt(gradErad) / CellP[p].Density) / All.cf_atime;}
                    if(gradErad > 0 && radacc > 0)
                    {
                        double dt_radacc = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * DMAX(ForceSoftening_KernelRadius(p), P[p].KernelRadius) / radacc);
                        if(dt_radacc < dt_rad) {dt_rad = dt_radacc;}
                    }
#endif
                } // end of loop over frequency bins
#endif // end of conditional to check if we're using FLD or OTVET with an explicit solver

                
                /* now consider the (simpler) CFL-type condition required for advective solvers like M1 or intensity/ray integrators */
#if defined(RT_M1) || defined(RT_LOCALRAYGRID)
                dt_courant = All.CourantFac * (L_particle*All.cf_atime) / C_LIGHT_CODE_REDUCED; /* courant-type criterion, using the reduced speed of light */
#if defined(SINGLE_STAR_STARFORGE_DEFAULTS)
                dt_courant = 0.4 * (L_particle*All.cf_atime) / C_LIGHT_CODE_REDUCED; /* hacked here for starforge, where mike's experimentation suggests we can get away with a slightly larger courant factor. remains experimental. courant-type criterion, using the reduced speed of light - here we hardcode the most aggressive possible Courant factor as an optimization */
#ifdef SINK_WIND_SPAWN
                if((CellP[p].MaxSignalVel > 0.5*C_LIGHT_CODE_REDUCED) || (P[p].ID == All.SpawnedWindCellID && P[p].Type == 0)) {dt_courant *= 0.5;} // be more careful if this is a jet cell or there are transluminal velocities
#endif
#endif                
#if defined(GALSF) && !defined(SINGLE_STAR_SINK_DYNAMICS) && defined(GALSF_FB_FIRE_STELLAREVOLUTION) // custom hacks for FIRE-RT tests; can override CFL condition with diffusion timestep certain limits
                int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++)
                {
                    double dt_rt_diffusion = dt_prefac_diffusion * (L_particle*All.cf_atime)*(L_particle*All.cf_atime) / (MIN_REAL_NUMBER + rt_diffusion_coefficient(p,kf, CellP));
                    if((CellP[p].Rad_E_gamma[kf] <= MIN_REAL_NUMBER) || (CellP[p].Rad_E_gamma_Pred[kf] <= MIN_REAL_NUMBER) || (CellP[p].Rad_E_gamma[kf] < 1.e-5*P[p].Mass*CellP[p].InternalEnergy)) {dt_rt_diffusion = 1.e10 * dt;} /* ignore particles where the radiation energy density is negligible */
                    dt_rad = DMIN(dt_rad, dt_rt_diffusion);
                }
                /* In the optically-thin / near-zero-opacity limit the FLD-style diffusion estimate above
                   diverges (rt_diffusion_coefficient ~ c_reduced/(kappa*rho) grows without bound as
                   kappa->0), which would demand an unphysical sub-free-streaming step. For the M1
                   advective solver the reduced-c Courant limit is the correct floor in that limit, so the
                   diffusion term may only relax the step upward, never push dt_rad below it. */
                dt_rad = DMAX(dt_rad, dt_courant);
                if(All.ComovingIntegrationOn) {dt_courant = DMAX(dt_courant, DMIN(dt_rad, 1.e3*dt_courant));}
#endif
                if(dt_courant < dt_rad) {dt_rad = dt_courant;}
#endif // explicit advective-type solver check

                
                /* one more check - we can optionally limit the timestep for explicit chemical timesteps: implicit solve is fine locally, but gets propagation somewhat wrong if timesteps too large, as that depends on opacity, which depends on ionization step! */
#if defined(RT_CHEM_PHOTOION) && defined(RT_TIMESTEP_LIMIT_RECOMBINATION) /* make sure this doesn't overshoot the recombination time for the opacity to change for ionizing photons */
                double ne_cgs = (CellP[p].Density * All.cf_a3inv * UNIT_DENSITY_IN_NHCGS), dt_recombination = All.CourantFac * (3.3e12/ne_cgs) / UNIT_TIME_IN_CGS;
                double dt_change = 1.e10*dt; if((CellP[p].Rad_E_gamma[RT_FREQ_BIN_H0] > 0)&&(fabs(CellP[p].Dt_Rad_E_gamma[RT_FREQ_BIN_H0])>0)) {dt_change = CellP[p].Rad_E_gamma[RT_FREQ_BIN_H0] / fabs(CellP[p].Dt_Rad_E_gamma[RT_FREQ_BIN_H0]);}
                dt_recombination = DMIN(DMAX(dt_recombination,dt_change), DMAX(dt_courant,dt_rad));
                if(dt_recombination < dt_rad) {dt_rad = dt_recombination;}
#endif

                if(dt_rad < dt) {
#ifdef TRANSPORT_SUBCYCLE
                    CellP[p].Transport_Dt_Subcycle = DMIN(CellP[p].Transport_Dt_Subcycle, dt_rad);
                    /* limit hydro dt so the subcycle cap is never exceeded, while allowing subcycling speedup */
                    double dt_max_hydro = (TRANSPORT_SUBCYCLE - 0.5) * dt_rad; /* -0.5 safety margin for rounding */
                    if(dt_max_hydro < dt) {dt = dt_max_hydro;}
#else
                    dt = dt_rad; // set the actual radiation timestep!
#endif
                }
            }
#endif // RADTRANSFER
            

#ifdef VISCOSITY
            if((CellP[p].Eta_ShearViscosity != 0) || (CellP[p].Zeta_BulkViscosity != 0)) /* no viscosity means no viscous constraint: same stand-in-regularizer trap as the conduction step above */
            {
                int kv1; double dv_mag=0,v_mag=1.0e-33;
                v_mag += P[p].Vel.norm_sq();
                double dv_mag_all = 0.0;
                for(kv1=0;kv1<3;kv1++)
                {
                    double dvmag_tmp = CellP[p].Gradients.Velocity[kv1].norm_sq();
                    dv_mag += dvmag_tmp /DMAX(P[p].Vel[kv1]*P[p].Vel[kv1],0.01*v_mag);
                    dv_mag_all += dvmag_tmp;
                }
                dv_mag = sqrt(DMAX(dv_mag, dv_mag_all/v_mag));
                double L_visc = DMAX(L_particle , 1. / (dv_mag + 1./L_particle)) * All.cf_atime;
                double visc_coeff = sqrt(CellP[p].Eta_ShearViscosity*CellP[p].Eta_ShearViscosity + CellP[p].Zeta_BulkViscosity*CellP[p].Zeta_BulkViscosity);
                double dt_viscosity = 0.25 * L_visc*L_visc / (1.0e-33 + visc_coeff) * CellP[p].Density * All.cf_a3inv;
                // since we use VISCOSITIES, not DIFFUSIVITIES, we need to add a power of density to get the right units //
#ifdef SUPER_TIMESTEP_DIFFUSION
                if(dt_viscosity < dt_superstep_explicit) dt_superstep_explicit = dt_viscosity; // explicit time-step
                double dt_advective = dt_viscosity * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + dv_mag))*All.cf_atime / L_visc);
                if(dt_advective < dt) dt = dt_advective; // 'advective' timestep: needed to limit super-stepping
#else
                if(dt_viscosity < dt) dt = dt_viscosity; // normal explicit time-step
#endif
            }
#endif

#if defined(GRAIN_BACKREACTION)
            if(6.*P[p].Grain_AccelTimeMin < dt) {dt = 6.*P[p].Grain_AccelTimeMin;}
#endif


#ifdef TURB_DIFFUSION
            {
#ifdef TURB_DIFF_METALS
                int k_species; double L_tdiff = L_particle * All.cf_atime; // don't use gradient b/c ill-defined pre-enrichment
                for(k_species=0;k_species<NUM_METAL_SPECIES;k_species++)
                {
                    double dt_tdiff = L_tdiff*L_tdiff / (1.0e-33 + CellP[p].TD_DiffCoeff); // here, we use DIFFUSIVITIES, so there is no extra density power in the equation //
                    if(dt_tdiff < dt) dt = dt_tdiff; // normal explicit time-step
                }
#endif
            }
#endif


#if defined(DIVBCLEANING_DEDNER)
            double fac_magnetic_pressure = 1. / All.cf_atime;
            double phi_b_units = Get_Gas_PhiField(p) / ( All.cf_atime * CellP[p].MaxSignalVel);
            double vA2_for_vsig = 0; // as in the fast-wavespeed scan above: no density to divide by means no magnetic contribution here, not an infinite signal speed
            if(CellP[p].Density > 0) {vA2_for_vsig = fac_magnetic_pressure * (CellP[p].Bfield().norm_sq() + phi_b_units*phi_b_units) / CellP[p].Density;}
            double vsig1 =  sqrt( CellP[p].effective_soundspeed()*CellP[p].effective_soundspeed() + vA2_for_vsig );

            dt_courant = 0.8 * All.CourantFac * (All.cf_atime*L_particle) / vsig1; // 2.0 factor may be added (PFH) //
            if(dt_courant < dt) {dt = dt_courant;}
#endif

            /* make sure that the velocity divergence does not imply a too large change of density or kernel length in the step */
            double divVel = P[p].Particle_DivVel;
            if(divVel != 0)
            {
                dt_divv = 1.5 / fabs(All.cf_a2inv * divVel);
                if(dt_divv < dt) {dt = dt_divv;}
            }


#if defined(TURB_DRIVING) && !defined(TURB_DRIVING_UPDATE_FORCE_ON_TURBUPDATE)
                /* gas cannot step larger than major updates to turbulent driving routine */
                double dt_turb_driving = 1.9 * st_return_dt_between_updates();
                if (dt > dt_turb_driving) {dt = dt_turb_driving;}
#endif
            

#ifdef SUPER_TIMESTEP_DIFFUSION
            /* now use the timestep information above to limit the super-stepping timestep */
            {
                int N_substeps = 5; /*!< number of sub-steps per super-timestep for super-timestepping algorithm */
                double nu_substeps = 0.04; /*!< damping parameter (0<nu<1), optimal behavior around ~1/sqrt[N_substeps] */

                /*!< pre-calculate the multipliers needed for the super-timestep sub-step */
                if(CellP[p].Super_Timestep_j == 0) {CellP[p].Super_Timestep_Dt_Explicit = dt_superstep_explicit;} // reset dt_explicit //
                double j_p_super = (double)(CellP[p].Super_Timestep_j + 1);
                double dt_superstep = CellP[p].Super_Timestep_Dt_Explicit / ((nu_substeps+1) + (nu_substeps-1) * cos(M_PI * (2*j_p_super - 1) / (2*(double)N_substeps)));

                double dt_touse = dt_superstep;
                if((dt <= dt_superstep)||(CellP[p].Super_Timestep_j > 0))
                {
                    /* if(dt <= dt_superstep): other constraints beat our super-step, so it doesn't matter: iterate */
                    /* if(CellP[p].Super_Timestep_j > 0): don't break mid-cycle, so iterate */
                    CellP[p].Super_Timestep_j++; if(CellP[p].Super_Timestep_j>=N_substeps) {CellP[p].Super_Timestep_j=0;} /*!< increment substep 'j' and loop if it cycles fully */
                } else {
                    /* ok, j=0 and dt > dt_superstep [the next super-step matters for starting a new cycle]: think about whether to start */
                    double dt_pred = dt_superstep * All.cf_hubble_a;
                    if(dt_pred > All.MaxSizeTimestep) {dt_pred = All.MaxSizeTimestep;}
                    if(dt_pred < All.MinSizeTimestep) {dt_pred = All.MinSizeTimestep;}
                    /* convert our physical timestep into the dimensionless units of the code */
                    integertime ti_min=TIMEBASE, ti_step = (integertime) (dt_pred / All.Timebase_interval);
                    /* check against valid limits */
                    if(ti_step<=1) {ti_step=2;}
                    if(ti_step>=TIMEBASE) {ti_step=TIMEBASE-1;}
                    while(ti_min > ti_step) {ti_min >>= 1;}  /* make it a power 2 subdivision */
                    ti_step = ti_min;
                    /* now turn it into a timebin */
                    int bin = get_timestep_bin(ti_step);
                    int binold = P[p].TimeBin;
                    if(bin > binold)  /* timestep wants to increase: check whether it wants to move into a valid timebin */
                    {
                        while(TimeBinActive[bin] == 0 && bin > binold) {bin--;} /* make sure the new step is synchronized */
                    }
                    /* now convert this -back- to a physical timestep */
                    double dt_allowed = GET_INTEGERTIME_FROM_TIMEBIN(bin) * unit_integertime_in_physical(-1, P);
                    if(dt_superstep > 1.5*dt_allowed)
                    {
                        /* the next allowed timestep [because of synchronization] is not big enough to fit the 'big step'
                            part of the super-stepping cycle. rather than 'waste' our timestep which will knock us into a
                            lower bin and defeat the super-stepping, we simply take the -safe- explicit timestep and
                            wait until the desired time-bin synchs up, so we can super-step */
                        dt_touse = dt_superstep_explicit; // use the safe [normal explicit] timestep and -do not- cycle j //
                    } else {
                        /* ok, we can jump up in bins to use our super-step; begin the cycle! */
                        CellP[p].Super_Timestep_j++; if(CellP[p].Super_Timestep_j>=N_substeps) {CellP[p].Super_Timestep_j=0;}
                    }
                }
                if(dt < dt_touse) {dt = dt_touse;} // set the actual timestep [now that we've appropriately checked everything above] //
            }
#endif

#ifdef NUCLEAR_NETWORK
            /* nuclear burning timestep limiter: the ODE solver handles stiffness internally
               via subcycling, and the energy cap in nuclear.cc prevents runaway injection.
               only limit the hydro timestep if the actual energy deposited last step was a
               significant fraction of the internal energy (> 10%), to improve operator-split accuracy. */
            if(CellP[p].NuclearEnergyGenerationRate != 0 && CellP[p].InternalEnergy > 0) {
                double de_last = fabs(CellP[p].NuclearEnergyGenerationRate) * dt; /* approx energy change at current dt */
                if(de_last > 0.5 * CellP[p].InternalEnergy) { /* would change u by >50% */
                    double dt_nuclear = 0.5 * CellP[p].InternalEnergy / fabs(CellP[p].NuclearEnergyGenerationRate);
                    if(dt_nuclear > 0 && dt_nuclear < dt) dt = dt_nuclear;
                }
            }
#endif

        } // closes if(P[p].Type == 0) [gas particle check] //


#if defined(DM_SIDM)
    /* Reduce time-step if this particle got interaction probabilities > 0.2 during the last time-step */
    if((1 << P[p].Type) & (DM_SIDM))
    {
        if(P[p].dtime_sidm > 0) {if(P[p].dtime_sidm < dt) {dt = P[p].dtime_sidm;}}
        if(dt > 0)
        {
            double p_target = 0.2; // desired maximum probability per timestep
            double vsig_fac = P[p].AGS_vsig*All.cf_atime/sqrt(3.);
            Vec3<double> dV = {vsig_fac, vsig_fac, vsig_fac}; // convert signal vel to velocity dispersion for estimating rates
#ifdef GRAIN_COLLISIONS
            double p_dt = prob_of_grain_interaction(P[p].Mass, P[p].Grain_Size, 0., P[p].AGS_KernelRadius, P[p].AGS_KernelRadius, dV, dt, p, P); // probability of interacting with another grain super-particle well within kernel, assuming same mass, H, and V~signalvel, for current timestep dt
#else
            double p_dt = prob_of_interaction(P[p].Mass, P[p].Mass, 0., P[p].AGS_KernelRadius, P[p].AGS_KernelRadius, dV, dt); // probability of interacting with another DM particle well within kernel, assuming same mass, H, and V~signalvel, for current timestep dt
#endif
            if(p_dt > p_target) {dt *= p_target / p_dt;}
        }
    }
#endif


    // add a 'stellar evolution timescale' criterion to the timestep, to prevent too-large jumps in feedback //
#if defined(GALSF_FB_FIRE_RT_HIIHEATING) || defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_FIRE_RT_LONGRANGE) || (defined(GALSF) && defined(RADTRANSFER))
    if(is_galsf_stellar_candidate_type(P[p].Type, All.ComovingIntegrationOn) && (P[p].Mass>0))
    {
        double star_age = evaluate_stellar_age_Gyr(p);
        double dt_stellar_evol;
        dt_stellar_evol = DMAX(2.0e-4, star_age/250.); // restrict to small steps for young stars //
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
#if defined(SNE_NONSINK_SPAWN)
        double mcorr = 1.e-4 * (P[p].Mass*UNIT_MASS_IN_SOLAR) / 0.1; // expectation of ()/X SNe per timestep -- here 0.1
        if(star_age > 0.044) {mcorr *= 0.02;} // into Ia regime, lower SNR means we can substantially relax this mass-dependent criterion
        if(mcorr > 1) {dt_stellar_evol /= DMIN(mcorr, 100.);} // don't use - ok to have multiple at low-res, but don't want too-big a jump or miss key stellar evolution
#endif
/* // below not necessary with newer code, can be safely skipped for optimization //
        double mcorr = 1.e-4 * (P[p].Mass*UNIT_MASS_IN_SOLAR) / 0.1; // expectation of ()/X SNe per timestep -- here 0.1
        if(star_age > 0.044) {mcorr *= 0.02;} // into Ia regime, lower SNR means we can substantially relax this mass-dependent criterion
        if(mcorr > 1) {dt_stellar_evol /= DMIN(mcorr, 10.);} // don't use - ok to have multiple at low-res, but don't want too-big a jump or miss key stellar evolution
*/
#else
        double mcorr = 1.e-5 * (P[p].Mass*UNIT_MASS_IN_SOLAR);
        if(mcorr < 1 && mcorr > 0) {dt_stellar_evol /= mcorr;}
#endif
        if(dt_stellar_evol < 1.e-6) {dt_stellar_evol = 1.e-6;}
        dt_stellar_evol /= (UNIT_TIME_IN_GYR); // convert to code units //
        if(dt_stellar_evol>0) {if(dt_stellar_evol<dt) {dt = dt_stellar_evol;}}
    }
#endif


#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    if(is_particle_a_special_zoom_target(p))
    {
        double dt_special_max = 1000./UNIT_TIME_IN_YR; // set a maximum physical timestep to prevent this centering from jumping
        if(dt > dt_special_max) {dt = dt_special_max;}
    }
#endif
    
    
#ifdef SINK_PARTICLES
    if(P[p].Type == 5)
    {
#if !defined(SINGLE_STAR_SINK_DYNAMICS) && defined(GALSF)
      double dt_accr = 4.2e5 / UNIT_TIME_IN_YR; // this is the 1% of Salpeter timescale; not relevant for low radiative efficiency
#else
      double dt_accr = All.MaxSizeTimestep;
#endif
        if(P[p].Sink_Mdot > 0 && P[p].Sink_Mass > 0 && All.Time > All.TimeBegin)
        {
#if (defined(SINK_GRAVCAPTURE_GAS) || defined(SINK_WIND_KICK)) && !defined(SINGLE_STAR_SINK_DYNAMICS)
            /* really want prefactor to be ratio of median gas mass to sink mass */
            dt_accr = 0.001 * DMAX(P[p].Sink_Mass, All.MaxMassForParticleSplit) / P[p].Sink_Mdot;
#if defined(SINK_WIND_KICK)
            dt_accr *= DMAX(0.1, All.Sink_accreted_fraction);
#endif
#else
            dt_accr = 0.1 * DMIN(P[p].Sink_Mass, All.MaxMassForParticleSplit) / P[p].Sink_Mdot;
#endif
#ifdef SINGLE_STAR_FB_JETS	    
            dt_accr = DMIN(dt_accr, target_mass_for_wind_spawning(p) / P[p].Sink_Mdot); 
#endif
        } // if(P[p].Sink_Mdot > 0 && P[p].Sink_Mass > 0)
#if defined(SINK_SEED_GROWTH_TESTS) || defined(FIRE_BHS)
        double dt_evol = 4.2e5 / UNIT_TIME_IN_YR; // totally arbitrary hard-coding here //
#ifdef TURB_DRIVING
        if(dt_evol > 1.e-3*st_return_mode_correlation_time()) {dt_evol=1.e-3*st_return_mode_correlation_time();}
#endif
        if(dt_accr > dt_evol) {dt_accr=dt_evol;}
#endif
        if(dt_accr > 0 && dt_accr < dt) {dt = dt_accr;}

        double dt_ngbs = 4.1 * get_physical_timestep_from_timebin(P[p].Sink_TimeBinGasNeighbor, p, P); /* standard wakeup-type threshold: use this by default here, unless dynamical interaction important (e.g. back-rx term from oscillation of sink c-o-m, which is important for single-sink sims */
        if(dt > dt_ngbs && dt_ngbs > 0) {dt = 1.01 * dt_ngbs; }

#if defined(SINGLE_STAR_TIMESTEPPING)
	    if(P[p].DensityAroundParticle > 0)
	    {
            double eps = DMAX( KERNEL_CORE_SIZE*ForceSoftening_KernelRadius(p), P[p].Sink_dr_to_NearestGasNeighbor);
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
            eps = DMAX(eps, P[p].SinkRadius);
#endif
            if(eps < MAX_REAL_NUMBER) {eps = DMAX(P[p].Get_Particle_Size(), eps);} else {eps = P[p].Get_Particle_Size();}
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
            eps = DMAX(eps, KERNEL_CORE_SIZE*P[p].AGS_KernelRadius);
#endif
            double dt_ff = sqrt(2*All.ErrTolIntAccuracy * pow(eps*All.cf_atime,3) / (All.G * P[p].Mass)); // fraction of the freefall time of the nearest gas particle from rest
            if(dt > dt_ff && dt_ff > 0) {dt = 1.01 * dt_ff;}

            double L_particle = P[p].Get_Particle_Size();
            double vsig = P[p].Sink_SurroundingGasVel;
#if defined(SINGLE_STAR_FB_TIMESTEPLIMIT) && !defined(NOGRAVITY)
            vsig += P[p].MaxFeedbackVel;
#endif                        
            double dt_cour_sink = All.CourantFac * (L_particle*All.cf_atime) / vsig;
            if(dt > dt_cour_sink && dt_cour_sink > 0 && isfinite(dt_cour_sink)) {dt = 1.01 * dt_cour_sink;}
        }
        if(P[p].StellarAge == All.Time)
        {   // want a brand new sink to be on the lowest occupied timebin
            long bin; for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinCount[bin] > 0) break;}
            double dt_min =  get_physical_timestep_from_timebin(bin, p, P);
            if(dt > dt_min && dt_min > 0) dt = 1.01 * dt_min;
        }
#endif // SINGLE_STAR_TIMESTEPPING
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
#ifdef SINGLE_STAR_FB_WINDS
        if(P[p].ProtoStellarStage == 5) {
            double mdot_spawn = single_star_wind_mdot(p,1);
            if(mdot_spawn > 0) {
                double dm_spawn = target_mass_for_wind_spawning(p), dt_spawn = dm_spawn / mdot_spawn;
                if(dt > dt_spawn && dt_spawn > 0) {dt = 1.01 * dt_spawn;}
            }}
#endif
#ifdef SINGLE_STAR_FB_SNE
        if ( (P[p].ProtoStellarStage == 6) && ( (P[p].Sink_Mass > 0) || (P[p].unspawned_wind_mass > 0) ) ) { //Star going supernova, still has mass to eject
            double eps = DMIN(KERNEL_CORE_SIZE*ForceSoftening_KernelRadius(p), P[p].KernelRadius);
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
            eps = DMAX(eps, P[p].SinkRadius);
#endif
            double t_clear=eps/single_star_SN_velocity(p);
            if(t_clear > 0 && dt > 0) {dt=DMIN(dt, DMAX(0.5*t_clear, 1.01*All.MinSizeTimestep));}; // time needed spawned wind particles to clear the sink so that we don't spawn on top of them (leading to progressively smaller timesteps from each spawn until crashing the code)
        }
#endif
#endif
    } // if(P[p].Type == 5)

#if defined(SINK_WIND_SPAWN_SET_BFIELD_POLTOR) /* KYSu: here for de-bugging jet injection model right now */
    if((P[p].Type==5) || (P[p].Type==0 && P[p].ID==All.SpawnedWindCellID && CellP[p].IniDen<0)) {if(dt>All.Sink_spawn_injectionradius/All.Sink_outflow_velocity && All.Sink_spawn_injectionradius>0 && All.Sink_outflow_velocity>0) {dt=All.Sink_spawn_injectionradius/All.Sink_outflow_velocity;}}
#endif
#endif // SINK_PARTICLES
    



    /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected, All.cf_hubble_a=1. */
    dt *= All.cf_hubble_a;

#ifdef ONLY_PM
    dt = All.MaxSizeTimestep;
#endif

    if(dt >= All.MaxSizeTimestep) {dt = All.MaxSizeTimestep;}

    if(dt >= dt_displacement) {dt = dt_displacement;}

    if((dt < All.MinSizeTimestep)||(((integertime) (dt / All.Timebase_interval)) <= 1))
    {
        PRINT_WARNING("Timestep wants to be below the limit `MinSizeTimestep'");
        double agrav_pm=0, agrav = P[p].GravAccel.norm() * All.cf_a2inv;
#ifdef PMGRID
        agrav_pm = P[p].GravPM.norm() * All.cf_a2inv;
#endif
        if(P[p].Type == 0)
        {
            double aturb=0, arad=0, ahydro = CellP[p].HydroAccel.norm();
#ifdef TURB_DRIVING
            aturb = CellP[p].TurbAccel.norm();
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            arad = CellP[p].Rad_Accel.norm();
#endif
            PRINT_WARNING("\n Cell-ID=%llu  dt_desired=%g dt_Courant=%g dt_Accel=%g\n accel_tot=%g accel_gravTree=%g accel_gravPM=%g accel_hydro=%g accel_rad=%g accel_turb=%g Pos_xyz=(%g|%g|%g) Vel_xyz=(%g|%g|%g)\n KernelRadius=%g Density=%g InternalEnergy=%g dtInternalEnergy=%g divV=%g Pressure=%g Cs_Eff=%g vAlfven=%g f_ion=%g\n csnd_for_signalspeed=%g eps_forcesoftening=%g mass=%g type=%d condition_number=%g Nngb=%g\n NVT=%.17g/%.17g/%.17g %.17g/%.17g/%.17g %.17g/%.17g/%.17g\n",
                          (unsigned long long) P[p].ID, dt, dt_courant*All.cf_hubble_a, sqrt(2*All.ErrTolIntAccuracy*All.cf_atime*ForceSoftening_KernelRadius(p) / ac)*All.cf_hubble_a,
                          ac, agrav, agrav_pm, ahydro, arad, aturb, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].Vel[0]/All.cf_atime, P[p].Vel[1]/All.cf_atime, P[p].Vel[2]/All.cf_atime,
                          P[p].KernelRadius*All.cf_atime, CellP[p].Density*All.cf_a3inv, CellP[p].InternalEnergy, CellP[p].DtInternalEnergy, P[p].Particle_DivVel*All.cf_a2inv,
                          CellP[p].Pressure*All.cf_a3inv, CellP[p].effective_soundspeed(), CellP[p].Alfven_speed(), Get_Gas_Ionized_Fraction(p, P, CellP),
                          csnd, ForceSoftening_KernelRadius(p)*All.cf_atime, P[p].Mass, P[p].Type, CellP[p].ConditionNumber, P[p].NumNgb,
                          CellP[p].NV_T[0][0],CellP[p].NV_T[0][1],CellP[p].NV_T[0][2],CellP[p].NV_T[1][0],CellP[p].NV_T[1][1],CellP[p].NV_T[1][2],CellP[p].NV_T[2][0],CellP[p].NV_T[2][1],CellP[p].NV_T[2][2]);
        }
        else // if(P[p].Type == 0)
        {
            PRINT_WARNING("Part-ID=%llu  dt_desired=%g dt_Accel=%g\n accel_tot=%g accel_gravTree=%g accel_gravPM=%g  mass=%g pos_xyz=(%g|%g|%g) vel_xyz=(%g|%g|%g) soft=%g type=%d\n",
                          (unsigned long long) P[p].ID, dt, sqrt(2*All.ErrTolIntAccuracy*All.cf_atime*ForceSoftening_KernelRadius(p) / ac)*All.cf_hubble_a,
                          ac, agrav, agrav_pm, P[p].Mass, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].Vel[0], P[p].Vel[1], P[p].Vel[2], ForceSoftening_KernelRadius(p), P[p].Type);
        }
        fflush(stdout); fprintf(stderr, "\n @ fflush \n");
#ifdef STOP_WHEN_BELOW_MINTIMESTEP
        /* Replaced endrun(888) with a CONTROLLED
         * stop request. Routing dt-floor through MPI_Abort (what
         * endrun(non-zero) does) bypasses Kokkos / CUDA-aware-MPI
         * cleanup and has correlated with Vista jobs stuck in SLURM CG
         * state. find_timesteps() collects the request at its end via
         * MPI_Allreduce; run() handles the global stop after the
         * find_timesteps call. NO MPI here -- this is inside a
         * per-particle loop where rank counts will differ. */
        if(P[p].Mass > 0) {
            gizmo_request_controlled_stop(888,
                "STOP_WHEN_BELOW_MINTIMESTEP: timestep < MinSizeTimestep",
                __FILE__, __LINE__, __FUNCTION__);
        }
#endif
        dt = All.MinSizeTimestep;
    }

    if(dt > 0.5 * TIMEBASE * All.Timebase_interval) {dt = 0.5 * TIMEBASE * All.Timebase_interval;} /* prevent integer timeline overflow */
    ti_step = (integertime) (dt / All.Timebase_interval);
    /* Floor the step on the integer timeline. This is not the MinSizeTimestep policy: a run that
       has asked for a controlled stop still has to reach the next phase boundary to take it, and a
       zero-length step here would instead abort hard from inside a per-particle loop. */
    if(ti_step<=1) ti_step=2;

    if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
        printf("\nError: A timestep of size zero was assigned on the integer timeline. Code must stop.\n"
               "Task=%d Part-ID=%llu dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%lld ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
               ThisTask, (unsigned long long) P[p].ID, dt, dt_courant, dt_divv, dt_displacement,
               All.Timebase_interval, (long long) ti_step, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].GravAccel[0], P[p].GravAccel[1], P[p].GravAccel[2]);
#ifdef PMGRID
        printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
        fflush(stdout); endrun(818);
        /* Soft-stop drains at the find_timesteps poll; until then, hand back a minimal
         * in-bounds integer-timeline step (matches the non-STOP_WHEN_BELOW_MINTIMESTEP
         * floor below) so the caller's get_timestep_bin() stays in range. This is a
         * continuation sentinel for a dying run, NOT a recovered/physical timestep. */
        ti_step = 2;
    }

    return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
    int i, type;
    int count[6];
    long long count_sum[6];
    double v[6], v_sum[6], mim[6], mnm[6], min_mass[6], mean_mass[6];
    double dt, dmean, asmth = 0;

    dt_displacement = All.MaxSizeTimestep;

    if(All.ComovingIntegrationOn)
    {
        for(type = 0; type < 6; type++)
        {
            count[type] = 0;
            v[type] = 0;
            mim[type] = 1.0e30;
            mnm[type] = 0;
        }

        for(i = 0; i < NumPart; i++)
        {
            if(P[i].Mass > 0)
            {
                double v2 = P[i].Vel.norm_sq();
                if(v2 > 0 && isfinite(v2)) {
                    count[P[i].Type]++;
                    if(P[i].Type == 0) {v[P[i].Type] += P[i].Mass * v2;} else {v[P[i].Type] += v2;} /* for gas use a weighted average to deal with extreme cell-mass difference situations */
                    if(mim[P[i].Type] > P[i].Mass) {mim[P[i].Type] = P[i].Mass;}
                    mnm[P[i].Type] += P[i].Mass;
                }
            }
        }

        MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(mnm, mean_mass, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumup_large_ints(6, count, count_sum);
        if(mean_mass[0] > 0) {v_sum[0] /= mean_mass[0];} /* for gas use a weighted average to deal with extreme cell-mass difference situations */

#ifdef GALSF
        /* add star and gas particles together to treat them on equal footing, using the original gas particle spacing. */
        double vsum0_0 = v_sum[0], minmass0_0 = min_mass[0], meanmass0_0 = mean_mass[0]; long long countsum0_0=count_sum[0];
        v_sum[0] += v_sum[4]; count_sum[0] += count_sum[4];
        if(count_sum[0] > 0) {
            if(count_sum[4]<=1 || (v_sum[0]+v_sum[4])*count_sum[4] < v_sum[4]*(count_sum[0]+count_sum[4])) {
                v_sum[4] += vsum0_0; count_sum[4] += countsum0_0; mean_mass[4] += meanmass0_0; min_mass[4] = DMAX(DMAX(min_mass[4],minmass0_0),mean_mass[4]/count_sum[4]);}}
        //v_sum[4] = v_sum[0]; count_sum[4] = count_sum[0]; if(count_sum[0] > 0) {min_mass[0] = min_mass[4] = (mean_mass[0] + mean_mass[4]) / count_sum[0];}
#ifdef SINK_PARTICLES
        vsum0_0 = v_sum[0]; minmass0_0 = min_mass[0]; meanmass0_0 = mean_mass[0]; countsum0_0=count_sum[0];
        v_sum[0] += v_sum[5]; count_sum[0] += count_sum[5];
        if(count_sum[0] > 0) {
            if(count_sum[5]<=1 || (v_sum[0]+v_sum[5])*count_sum[5] < v_sum[5]*(count_sum[0]+count_sum[5])) {
                v_sum[5] += vsum0_0; count_sum[5] += countsum0_0; mean_mass[5] += meanmass0_0; min_mass[5] = DMAX(DMAX(min_mass[5],minmass0_0),mean_mass[5]/count_sum[5]);}}
        //v_sum[5] = v_sum[0]; count_sum[5] = count_sum[0]; min_mass[5] = min_mass[0];
#endif
#ifdef SPECIAL_POINT_MOTION
        v_sum[SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES] = v_sum[0];
        count_sum[SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES] = count_sum[0];
        min_mass[SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES] = min_mass[0];
#endif
#endif

        if(ThisTask == 0) {printf("Global displacement time constraint computation: \n");}
        for(type = 0; type < 6; type++)
        {
            if(count_sum[type] > 0 && v_sum[type] > 0)
            {
#ifdef GALSF
                if(type == 0 || type == 4)
#else
                if(type == 0)
#endif
                    dmean = pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);
                else
                    dmean = pow(min_mass[type] / ((All.OmegaMatter - All.OmegaBaryon) * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);

#ifdef SINK_PARTICLES
                if(type == 5) {dmean = pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);}
#endif
                dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
                asmth = All.Asmth[0];
#ifdef PM_PLACEHIGHRESREGION
                if(((1 << type) & (PM_PLACEHIGHRESREGION))) {asmth = All.Asmth[1];}
#endif
                if(asmth < dmean) {dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);}
#endif
                if(ThisTask == 0) {printf(" ..type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);}
                if(dt < dt_displacement && dt > 0) {dt_displacement = dt;}
            }
        }

        if(ThisTask == 0) {printf(" ..global displacement time constraint: %g  (All.MaxSizeTimestep=%g)\n", dt_displacement, All.MaxSizeTimestep);}
    }
}



int get_timestep_bin(integertime ti_step)
{
    int bin = -1;

    if(ti_step == 0)
        return 0;

    if(ti_step == 1)
    {
        printf("time-step of integer size 1 not allowed (task=%d)\n", ThisTask); fflush(stdout);
        endrun(90001006);
        return 0;   /* graceful: bad-stop set; return bin 0 (the loop below yields 0 for ti_step==1); drains at the find_timesteps poll */
    }

    while(ti_step)
    {
        bin++;
        ti_step >>= 1;
    }

    return bin;
}





/* Authoritative rebuild of the wakeup dirty sidecar from P[]. Runs when the
 * index mapping changed (domain decomp / rearrange / init) or on first use. */
void wakeup_sidecar_rebuild(void)
{
    if(WakeupDirty) {
        for(int i = 0; i < NumPart; i++) { WakeupDirty[i] = (P[i].wakeup != 0) ? 1 : 0; }
    }
    WakeupDirtyValid = 1;
}

/* The particles the last process_wake_ups moved to a shorter step (particles_woken_last). */
static std::vector<int> WokenParticles;

void process_wake_ups(void)
{
    WokenParticles.clear();
#ifdef FORCE_EQUAL_TIMESTEPS
    return; /* no wakeups if all particles are on the same timestep */
#endif
    int i, n, max_time_bin_active, bin, binold, prev, next; long long ntot;
    integertime dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;

    /* find the next kick time */
    for(n = 0, ti_next_kick = TIMEBASE; n < TIMEBINS; n++)
    {
        if(TimeBinCount[n])
        {
            if(n > 0)
            {
                dt_bin = GET_INTEGERTIME_FROM_TIMEBIN(n);
                ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
            }
            else {dt_bin = 0; ti_next_for_bin = All.Ti_Current;}
            if(ti_next_for_bin < ti_next_kick) {ti_next_kick = ti_next_for_bin;}
        }
    }

    MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_TYPE_TIME, MPI_MIN, MPI_COMM_WORLD);

    PRINT_STATUS("Predicting next timestep: %g", (ti_next_kick_global - All.Ti_Current) * All.Timebase_interval);
    max_time_bin_active = 0;
    /* get the highest bin, that is active next time */
    for(n = 0; n < TIMEBINS; n++)
    {
        dt_bin = (((integertime) 1) << n);
        if((ti_next_kick_global % dt_bin) == 0) {max_time_bin_active = n;}
    }

    /* move the particle into the highest bin, that is active in the next timestep and that is lower than its last timebin */
    bin = 0; for(n = 0; n < TIMEBINS; n++) {if(TimeBinCount[n] > 0) {bin = n; break;}}
    n = 0;

    MPI_Allreduce(&NeedToWakeupParticles_local, &NeedToWakeupParticles, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); // if one process processes wakeups then they all should, just in case a woke particle gets swapped to another process before we get here

    int wakeup_bin_offset = 0;
    while(((integertime)1 << wakeup_bin_offset) < (integertime)WAKEUP) wakeup_bin_offset++;

    if(NeedToWakeupParticles){
	/* Floor for positive (relative) wakeup: the lowest currently-OCCUPIED
	 * and ACTIVE bin (global across ranks). Without this floor, repeated
	 * a-wakes-b-wakes-c shells of relative wakeup can drive bins below
	 * anything currently being processed (waker_bin - offset cascades down
	 * each generation). The floor caps relative wakeup at the most-aggressive
	 * bin already in flight, preventing the multiplicative cascade while
	 * preserving legitimate hydro-style subcycle wakeups in [floor, max]. */
	int local_lowest_occupied_active_bin = TIMEBINS;
	for(int nb = 0; nb < TIMEBINS; nb++) {   /* own index: n is the woken-particle counter, zeroed above */
	    if(TimeBinActive[nb] && TimeBinCount[nb] > 0) {
	        local_lowest_occupied_active_bin = nb;
	        break;
	    }
	}
	int lowest_occupied_active_bin = TIMEBINS;
	MPI_Allreduce(&local_lowest_occupied_active_bin, &lowest_occupied_active_bin,
	              1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	/* Dirty-sidecar scan: reject on the contiguous WakeupDirty byte array
	 * instead of striding the fat P[] struct. Rebuild from P[] if the index
	 * mapping changed since the last scan (domain decomp / rearrange / init). */
	if(WakeupDirty && !WakeupDirtyValid) { wakeup_sidecar_rebuild(); }
	for(i = 0; i < NumPart; i++)
	{
	    if(WakeupDirty) {
		if(!WakeupDirty[i]) {continue;}
		if(!P[i].wakeup) {WakeupDirty[i] = 0; continue;}   /* superset false positive — self-clear */
	    } else {
		if(!P[i].wakeup) {continue;}
	    }
#if !defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
	    if(P[i].Type != 0) {continue;} // only gas particles can be awakened
#endif
	    if(P[i].Mass <= 0) {continue;}
	    binold = P[i].TimeBin;
	    if(TimeBinActive[binold]) {continue;}

	    if(P[i].wakeup > 0) {
		/* hydro wakeup: target timestep = dt_waker / WAKEUP */
		int waker_bin = P[i].wakeup - 1;
		bin = IMAX(0, waker_bin - wakeup_bin_offset);
		/* Floor at the lowest-currently-occupied-and-active bin (see comment
		 * above the i-loop). Prevents the multiplicative wakeup cascade. */
		if(bin < lowest_occupied_active_bin) bin = lowest_occupied_active_bin;
#ifdef STOP_WHEN_BELOW_MINTIMESTEP
		/* If wakeup assigns a bin whose physical timestep is below MinSizeTimestep, abort with full
		 * context. Without this, the wakeup-application path silently floors at bin 0 (dt ~ Timebase_interval) bypassing the get_timestep() warning. */
		{
		    double dt_assigned_physical = (double)GET_INTEGERTIME_FROM_TIMEBIN(bin) * All.Timebase_interval;
		    if(dt_assigned_physical < All.MinSizeTimestep) {
			printf("\n[WAKEUP-MINTIMESTEP] ABORT: rank=%d i=%d ID=%llu Type=%d  binold=%d -> bin=%d  waker_bin=%d  wakeup_bin_offset=%d  dt_assigned_physical=%g < MinSizeTimestep=%g  Mass=%g\n",
			    ThisTask, i, (unsigned long long)P[i].ID, P[i].Type, binold, bin,
			    waker_bin, wakeup_bin_offset,
			    dt_assigned_physical, All.MinSizeTimestep, P[i].Mass);
			if(P[i].Type == 0) {
			    printf("[WAKEUP-MINTIMESTEP]   Pressure=%g Density=%g InternalEnergy=%g MaxSignalVel=%g\n",
				CellP[i].Pressure, CellP[i].Density, CellP[i].InternalEnergy, CellP[i].MaxSignalVel);
			}
			fflush(stdout);
			endrun(889);
		    }
		}
#endif
	    } else {
		/* generic wakeup (sinks, merge, etc.): use highest active bin */
		bin = max_time_bin_active;
	    }
	    if(bin > max_time_bin_active) {bin = max_time_bin_active;} /* must be active at next kick */
	    if(bin >= binold) {bin = binold;} /* don't increase timestep */

	    if(bin != binold)
	    {
		integertime tstart = P[i].Ti_begstep + P[i].integertime_step(); /* the step this particle is actually on, which a previous demotion may have truncated below its bin length */
		integertime t_2 = P[i].Ti_current;
		if(t_2 > tstart) {tstart = t_2;}
		integertime tend = All.Ti_Current;

		TimeBinCount[binold]--;
		if(P[i].Type == 0) {TimeBinCountGas[binold]--;}

		prev = PrevInTimeBin[i];
		next = NextInTimeBin[i];

		if(FirstInTimeBin[binold] == i) {FirstInTimeBin[binold] = next;}
		if(LastInTimeBin[binold] == i) {LastInTimeBin[binold] = prev;}
		if(prev >= 0) {NextInTimeBin[prev] = next;}
		if(next >= 0) {PrevInTimeBin[next] = prev;}

		if(TimeBinCount[bin] > 0)
		{
		    PrevInTimeBin[i] = LastInTimeBin[bin];
		    NextInTimeBin[LastInTimeBin[bin]] = i;
		    NextInTimeBin[i] = -1;
		    LastInTimeBin[bin] = i;
		}
		else
		{
		    FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
		    PrevInTimeBin[i] = NextInTimeBin[i] = -1;
		}
		TimeBinCount[bin]++;
		if(P[i].Type == 0) {TimeBinCountGas[bin]++;}
		P[i].TimeBin = bin;
        if(TimeBinActive[bin]) {NumForceUpdate++;}
		n++;
		WokenParticles.push_back(i);

		/* The kick this particle already received covers past the time it is being woken to. Do NOT
		   try to reverse it: re-deriving the increment with reversed bounds would only cancel the
		   original if the acceleration and energy rate were unchanged since, and they are not, so it
		   injects energy instead. Saitoh & Makino (2009) eq (3) avoids integrating the system backwards
		   for exactly this reason and instead sets the new time consistent with the system time, which
		   is what the truncation below does. */
		if(tend < tstart) {set_predicted_quantities_for_extra_physics(i);}
		/* End the step at the current system time, so the interval the applied kick already covered is
		   not integrated twice. dt_step stays authoritative, but is no longer a power-of-two bin length
		   and so deliberately disagrees with TimeBin: a step must be read from integertime_step(), never
		   derived from the bin. */
		{
		    integertime dt_truncated = All.Ti_Current - P[i].Ti_current;
		    if(dt_truncated > 0) {P[i].Ti_begstep = P[i].Ti_current; P[i].dt_step = dt_truncated;}
		    else {P[i].Ti_begstep = All.Ti_Current; P[i].dt_step = GET_INTEGERTIME_FROM_TIMEBIN(bin);}
		}
#if defined(USE_TIMESTEP_DILATION_FOR_ZOOMS)
        /* a wakeup starts a new step for this particle, so freeze its dilation factor at the
           position it now holds, as a normal timestep assignment would. This must follow the
           partial kick reversed just above, which cancels the kick it undoes only while the
           factor still matches the one that kick was applied with. */
        P[i].TimestepDilationFactor = return_timestep_dilation_factor(i, P);
#endif
		if(P[i].Ti_current < All.Ti_Current) {P[i].Ti_current=All.Ti_Current;}
	    }
	}
    }

    sumup_large_ints(1, &n, &ntot);
    if(ThisTask == 0) {if(ntot > 0) {printf("%d%09d particles activated (in wakeup check).\n", (int) (ntot / 1000000000), (int) (ntot % 1000000000));}}
    NeedToWakeupParticles = 0;
    NeedToWakeupParticles_local = 0;
}

/* The particles the last process_wake_ups moved to a shorter step, and how many.  A wake-up changes how
   fast a particle moves without kicking it through the active list: it re-freezes the particle's dilation
   factor, and under the finite-volume kick, which runs by active bin, a particle woken into a bin active
   now is kicked with the active set although ActiveParticleList was built before it woke.  So the motion
   bounds that hold these particles are raised over them once the first half-kick has run. */
int particles_woken_last(const int **idx)
{
    *idx = WokenParticles.data();
    return (int)WokenParticles.size();
}




#ifdef BOX_SHEARING
void calc_shearing_box_pos_offset(void) /* function that calculates the shear-offset between the shear-periodic boundaries in a shearing box */
{
    /* Shearing_Box_*_Offset are macros (allvars.h:118-119) that already
       expand to (All.Shearing_Box_*_Offset). Use the bare macro names — the
       `All.` prefix was double-resolving via the macro and yielding
       "expected a member name" under nvc++. */
    Shearing_Box_Pos_Offset = Shearing_Box_Vel_Offset * All.Time;
    while(Shearing_Box_Pos_Offset > boxSize_Y) {Shearing_Box_Pos_Offset -= boxSize_Y;}
}
#endif


/* timestep_dilation_factor, unit_integertime_in_physical, get_physical_timestep_from_timebin,
   get_particle_timestep_in_physical: definitions now in timestep_functions.h (single source of truth).
   Include with non-inline linkage to provide externally-visible symbols. */
#undef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION
#include "timestep_functions.h"


/* Timestep dilation for zoom-in runs with extreme dynamic range. The dilation factor f = 1/a <= 1
   is the rate at which a particle's clock advances relative to the global integer timeline: the
   assigned integer step is divided by f, and every conversion of an integer step back to physical
   time multiplies by f, so the physical landing time of the step is unchanged while the local
   dynamics are integrated more finely. f is frozen for each particle when its timestep is assigned
   (see get_timestep) and cached in P[].TimestepDilationFactor; read the cached value with
   timestep_dilation_factor(). The two functions below are the live evaluations, needed only where
   no cache exists: at timestep assignment, and for tree nodes. Both are built on the device-callable
   helpers in timestep_functions.h, so they sit below its inclusion. */

/* live dilation factor for particle i. Called at timestep assignment; all other consumers read the
   frozen value via timestep_dilation_factor(). */
double return_timestep_dilation_factor(int i, struct particle_data *pp)
{
#if !defined(USE_TIMESTEP_DILATION_FOR_ZOOMS)
    (void)i; (void)pp; return 1;
#else

    if(All.Time <= All.TimeBegin) {return 1;}
    if(i < 0) {return 1;}
#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    if(pp[i].Type != 4 && pp[i].Type != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) {return 1;} /* only do cosmological 'stars' type -and- the special smoothing-source-types */
#else
    if(pp[i].Type != 4) {return 1;} /* only do cosmological 'stars' type */
#endif
#endif

    /* now specify some dilation factor a(r) or otherwise */
    double a = 1;

#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    double r_to_sink = pp[i].Min_Distance_to_Sink;
    if(pp[i].Type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) {r_to_sink = 0;}
    double wt = weight_function_for_weighted_motion_smoothing(r_to_sink, 0);
    if(wt > 0 && wt < 1) {a = 1. / wt;}
#endif

#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    double r = distance_to_nearest_refinement_center(pp[i].Pos);
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 3)
    r = sqrt(pp[i].Pos[0]*pp[i].Pos[0] + pp[i].Pos[1]*pp[i].Pos[1] + pp[i].Pos[2]*pp[i].Pos[2]);
#endif
    a = nuclear_zoom_dilation_amplitude(r);
#endif

    return 1. / a;
#endif
}


/*! Live dilation factor at the center of mass of tree node 'no'. Host callers pass the
 *  global node array; see return_node_timestep_dilation_factor_P for the rationale. */
double return_node_timestep_dilation_factor(int no)
{
    return return_node_timestep_dilation_factor_P(no, Nodes);
}

double get_particle_feedback_timestep_in_physical(int i, struct particle_data *pp)
{
#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    return pp[i].integertime_step() * (All.Timebase_interval / All.cf_hubble_a); /* no dilation */
#elif defined(GALSF_LIMIT_FBTIMESTEPS_FROM_BELOW)
    return DMAX(get_particle_timestep_in_physical(i, pp), All.Dt_Min_Between_FBCalc_Gyr / UNIT_TIME_IN_GYR);
#else
    return get_particle_timestep_in_physical(i, pp);
#endif
}
