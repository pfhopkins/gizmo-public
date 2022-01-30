#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"


/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part (adding/removing calls, re-ordering some routines, and
 * adding hooks to new elements such as particle splitting, as necessary)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */
void run(void)
{
    CPU_Step[CPU_MISC] += measure_time();

    if(RestartFlag != 1)		/* need to compute forces at initial synchronization time, unless we restarted from restart files */
    {
        output_log_messages();

        domain_Decomposition(0, 0, 0);

        set_non_standard_physics_for_current_time();

        compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */

        compute_hydro_densities_and_forces();	/* densities, gradients, & hydro-accels for synchronous particles */

        calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */
    }

    while(1)			/* main timestep iteration loop */
    {
        compute_statistics();	/* regular statistics outputs (like total energy) */

        write_cpu_log();		/* output some CPU usage log-info (accounts for everything needed up to the current sync-point) */

        if((All.Ti_Current >= TIMEBASE) || (All.Time > All.TimeMax)) /* check whether we reached the final time */
        {
            if(ThisTask == 0) {printf("\nFinal time=%g reached. Simulation ends.\n", All.TimeMax);}
            restart(0); /* write a restart file to allow continuation of the run for a larger value of TimeMax */
            if(All.Ti_lastoutput != All.Ti_Current) {savepositions(All.SnapshotFileCount++);} /* make a snapshot at the final time in case none has produced at this time; this will be overwritten if All.TimeMax is increased and the run is continued */
            break;
        }

        find_timesteps();		/* find-timesteps */
        int TreeReconstructFlag_local = TreeReconstructFlag;
#ifdef HERMITE_INTEGRATION
        HermiteOnlyFlag = 1;
        gravity_tree();	/* re-compute gravitational accelerations for synchronous particles */
        HermiteOnlyFlag = 0;
#endif
        do_first_halfstep_kick();	/* half-step kick at beginning of timestep for synchronous particles */

        find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
                                             * If needed, this function will also write an output file
                                             * at the desired time.
                                             */

        output_log_messages();	/* write some info to log-files */

        set_non_standard_physics_for_current_time();	/* update auxiliary physics for current time */

        int reconstructed_tree = 0;
#if defined(SINGLE_STAR_SINK_DYNAMICS)
        if(All.NumForcesSinceLastDomainDecomp > All.TreeDomainUpdateFrequency * All.TotNumPart) {TreeReconstructFlag_local = 1;}
#endif
        MPI_Allreduce(&TreeReconstructFlag_local, &TreeReconstructFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); // if one process reconstructs the tree then everbody has to
        if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)	/* check whether we have a big step */
        {
            domain_Decomposition(0, 0, 1);      /* do domain decomposition if step is big enough, and set new list of active particles  */
            reconstructed_tree = 1;
        }
        else if(TreeReconstructFlag) {domain_Decomposition(0, 0, 1); reconstructed_tree = 1;}
        else
        {
            force_update_tree();	/* update tree dynamically with kicks of last step so that it can be reused */
            make_list_of_active_particles();	/* now we can set the new chain list of active particles */
        }

        compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)
        // Need to figure out how frequently we calculate this; below is pretty rough //
#ifdef PMGRID
        if(All.Ti_Current == All.PM_Ti_endstep && get_random_number(1+All.Ti_Current) < 0.05)
#else
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
        {
            disp_density(); /* compute the DM velocity dispersion around gas particles every 20 PM steps, should be sufficient */
        }
#endif
#endif

        /* flag particles which will be feedback centers, so kernel lengths can be computed for them */
#ifdef GALSF_FB_MECHANICAL
        determine_where_SNe_occur(); // for mechanical FB models
#endif
#ifdef GALSF_FB_THERMAL
        determine_where_addthermalFB_events_occur(); // (same, but for simple thermal feedback models)
#endif

        compute_hydro_densities_and_forces();	/* densities, gradients, & hydro-accels for synchronous particles */
        
#ifdef PARTICLE_MERGE_SPLIT_EVERY_TIMESTEP // do merge/split routines every single timestep - need to do it here if we didn't do it during domain decomp on a coarse timestep
        if(!reconstructed_tree)
        {
            merge_and_split_particles();
            rearrange_particle_sequence();
        }
#endif
        
        do_second_halfstep_kick();	/* this does the half-step kick at the end of the timestep */

        calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */

#ifdef HERMITE_INTEGRATION // we do a prediction step using the saved "old" pos, accel and jerk from the beginning of the timestep. Then we recompute accel and jerk and do the correction
        do_hermite_prediction();
        HermiteOnlyFlag = 2;
        gravity_tree();	/* re-compute gravitational accelerations for synchronous particles */
        HermiteOnlyFlag = 0;
        do_hermite_correction();
#endif
        /* Check whether we need to interrupt the run */
        int stopflag = 0;
        if(ThisTask == 0)
        {
            FILE *fd;
            char stopfname[1000];
            sprintf(stopfname, "%sstop", All.OutputDir);
            if((fd = fopen(stopfname, "r")))	/* Is the stop-file present? If yes, interrupt the run. */
            {
                fclose(fd);
                stopflag = 1;
                unlink(stopfname);
            }

            if(CPUThisRun > 0.85 * All.TimeLimitCPU)	/* are we running out of CPU-time ? If yes, interrupt run. */
            {
                printf("reaching time-limit. stopping.\n");
                stopflag = 2;
            }
        }

        MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(stopflag)
        {
            restart(0);		/* write restart file */
            MPI_Barrier(MPI_COMM_WORLD);

            if(stopflag == 2 && ThisTask == 0)
            {
                FILE *fd;
                char contfname[1000];
                sprintf(contfname, "%scont", All.OutputDir);
                if((fd = fopen(contfname, "w")))
                    fclose(fd);

                if(All.ResubmitOn)
                    execute_resubmit_command();
            }
            return;
        }

        if(ThisTask == 0)
        {
            /* is it time to write one of the regularly space restart-files? */
            if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
            {
                All.TimeLastRestartFile = CPUThisRun;
                stopflag = 3;
            }
            else
                stopflag = 0;
        }

        MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(stopflag == 3)
        {
            restart(0);		/* write an occasional restart file */
            stopflag = 0;
            All.TimeLastRestartFile += report_time();
        }

        set_random_numbers();	/* draw a new list of random numbers */

        report_memory_usage(&HighMark_run, "RUN");
    }

}



void set_non_standard_physics_for_current_time(void)
{
#if defined(COOLING) && !defined(CHIMES)
    /* set UV background for the current time */
    IonizeParams();
#endif

#if defined(COOL_METAL_LINES_BY_SPECIES) && !defined(CHIMES)
    /* load the metal-line cooling tables appropriate for the UV background */
    if(All.ComovingIntegrationOn) {LoadMultiSpeciesTables();}
#endif

}



void calculate_non_standard_physics(void)
{
#ifdef PARTICLE_EXCISION
    apply_excision();
#endif


#if defined(TURB_DRIVING) && defined(TURB_DRIVING_SPECTRUMGRID)
    if(All.Time >= All.TimeNextTurbSpectrum) {powerspec_turb(All.FileNumberTurbSpectrum++); All.TimeNextTurbSpectrum += All.TimeBetTurbSpectrum;}
#endif


#ifdef GALSF /* PFH set of feedback routines; here because for e.g. strong SNe, obtain better stability if they are coupled discretely just -before- the hydro force is computed */
    //compute_stellar_feedback();
#endif

    
#ifdef BLACK_HOLES /***** black hole accretion and feedback *****/
    CPU_Step[CPU_MISC] += measure_time();
    blackhole_accretion();
#ifdef BH_WIND_SPAWN
    double MaxUnSpanMassBH_global;
    MPI_Allreduce(&MaxUnSpanMassBH, &MaxUnSpanMassBH_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if(MaxUnSpanMassBH_global > (BH_WIND_SPAWN)*All.BAL_wind_particle_mass)
    {
        spawn_bh_wind_feedback();
        rearrange_particle_sequence();
        MaxUnSpanMassBH=MaxUnSpanMassBH_global=0.;
    }
#endif
    MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_BLACKHOLES] += measure_time();
#endif


#if (defined(BLACK_HOLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
    if(All.Time >= All.TimeNextOnTheFlyFoF) {fof_fof(-1); /* this will find new black hole seed halos and/or assign host halo masses for the variable wind model */
        if(All.ComovingIntegrationOn) {All.TimeNextOnTheFlyFoF *= All.TimeBetOnTheFlyFoF;} else {All.TimeNextOnTheFlyFoF += All.TimeBetOnTheFlyFoF;}}
#endif

#ifdef RADTRANSFER
    CPU_Step[CPU_MISC] += measure_time();
#if defined(RT_SOURCE_INJECTION)
    int flag; flag=1;
#if !defined(RT_INJECT_PHOTONS_DISCRETELY)
    flag = Flag_FullStep; /* for continous injection, requires all sources and gas be active synchronously or else 2x-counts */
#endif
#if !defined(GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION)
    if(flag) {rt_source_injection();} /* source injection into neighbor gas particles (only on full timesteps, if using non-discrete scheme) */
#endif
#endif
#if defined(RT_DIFFUSION_CG) /* use the CG method to solve the RT diffusion equation implicitly for all particles; do only on full timesteps, requires synchronous timestepping right now */
    if(Flag_FullStep) {All.Radiation_Ti_endstep = All.Ti_Current; rt_diffusion_cg_solve(); All.Radiation_Ti_begstep = All.Radiation_Ti_endstep;}
#endif
#if defined(RT_CHEM_PHOTOION) && (!defined(COOLING) || defined(RT_COOLING_PHOTOHEATING_OLDFORMAT))
    rt_update_chemistry(); /* chemistry updated at sub-stepping as well */
#ifndef IO_REDUCED_MODE
    if(Flag_FullStep) {rt_write_chemistry_stats();}
#endif
#endif
    MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_RTNONFLUXOPS] += measure_time();
#endif // RADTRANSFER block

#ifdef COOLING	/**** radiative cooling and chemistry  *****/
    cooling_parent_routine(); // top-level cooling and chemistry subroutine //
    MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_COOLINGSFR] += measure_time(); // finish time calc for SFR+cooling
#endif


#ifdef GALSF /**** star/sink particle formation *****/
    star_formation_parent_routine(); // top-level star formation routine (because this involves common particle conversions, want to keep this at end of this subroutine) //
    MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_COOLINGSFR] += measure_time(); // finish time calc for SFR+cooling
#endif

#ifdef BH_INTERACT_ON_GAS_TIMESTEP
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]){if(P[i].Type == 5 && P[i].do_gas_search_this_timestep){P[i].dt_since_last_gas_search = 0;}}
#endif

}



void compute_statistics(void)
{
    if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
    {
#if !defined(EVALPOTENTIAL)          // compute_potential is not defined if EVALPOTENTIAL is on //
#ifdef COMPUTE_POTENTIAL_ENERGY
        compute_potential();
#endif
#endif
        energy_statistics();	/* compute and output energy statistics */
        All.TimeLastStatistics += All.TimeBetStatistics;
    }
}



void execute_resubmit_command(void)
{
    char buf[1000];
    sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
    system(buf);
#endif
}



/*! This function finds the next synchronization point of the system
 * (i.e. the earliest point of time any of the particles needs a force
 * computation), and drifts the system to this point of time.  If the
 * system drifts over the desired time of a snapshot file, the
 * function will drift to this moment, generate an output, and then
 * resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, i, prev;
  integertime dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;
  int highest_active_bin, highest_occupied_bin;
  double timeold;

  timeold = All.Time;

  All.NumCurrentTiStep++;	/* we are now moving to the next sync point */

  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE, highest_occupied_bin = 0; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      highest_occupied_bin = n;
	      dt_bin = GET_INTEGERTIME_FROM_TIMEBIN(n);;
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }
	  else
	    {
	      dt_bin = 0;
	      ti_next_for_bin = All.Ti_Current;
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_TYPE_TIME, MPI_MIN, MPI_COMM_WORLD);

  while(ti_next_kick_global >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
        All.Ti_Current = All.Ti_nextoutput;

        if(All.ComovingIntegrationOn) {All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);}
            else {All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;}

        set_cosmo_factors_for_current_time();

        move_particles(All.Ti_nextoutput);
        MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_DRIFT] += measure_time();

#ifdef OUTPUT_POTENTIAL
#if !defined(EVALPOTENTIAL) || (defined(EVALPOTENTIAL) && defined(OUTPUT_RECOMPUTE_POTENTIAL))
        domain_Decomposition(0, 0, 0);
        compute_potential();
#endif
#endif

        savepositions(All.SnapshotFileCount++);	/* write snapshot file */
        All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }


  All.Previous_Ti_Current = All.Ti_Current;
  All.Ti_Current = ti_next_kick_global;

  if(All.ComovingIntegrationOn) {All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);}
    else {All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;}

  set_cosmo_factors_for_current_time();
#ifdef BOX_SHEARING
    calc_shearing_box_pos_offset();
#endif

#ifdef GR_TABULATED_COSMOLOGY_G
  All.G = All.Gini * dGfak(All.Time);
#endif

  All.TimeStep = All.Time - timeold;

  /* mark the bins that will be active */
  for(n = 1, TimeBinActive[0] = 1, NumForceUpdate = TimeBinCount[0], highest_active_bin = 0; n < TIMEBINS;
      n++)
    {
      dt_bin = GET_INTEGERTIME_FROM_TIMEBIN(n);
      if((ti_next_kick_global % dt_bin) == 0)
	{
	  TimeBinActive[n] = 1;
	  NumForceUpdate += TimeBinCount[n];
	  if(TimeBinCount[n])
	    highest_active_bin = n;
	}
      else
	TimeBinActive[n] = 0;
    }

  sumup_large_ints(1, &NumForceUpdate, &GlobNumForceUpdate);
  All.NumForcesSinceLastDomainDecomp += GlobNumForceUpdate;
  MPI_Allreduce(&highest_active_bin, &All.HighestActiveTimeBin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&highest_occupied_bin, &All.HighestOccupiedTimeBin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(GlobNumForceUpdate == All.TotNumPart)
    {
      Flag_FullStep = 1;
      if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
	terminate("Something is wrong with the time bins.\n");
    }
  else
    Flag_FullStep = 0;




  /* move the new set of active/synchronized particles. Note: We do not yet call make_list_of_active_particles(), since we
   * may still need to old list in the dynamic tree update */
  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {if(TimeBinActive[n]) {for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i]) {drift_particle(i, All.Ti_Current);}}}

}


void make_list_of_active_particles(void)
{
    int i, n, prev;
    /* make a link list with the particles in the active time bins */
    FirstActiveParticle = -1;

    for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
        if(TimeBinActive[n])
        {
            for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
            {
                if(P[i].Mass <= 0) {continue;}
                if(prev == -1) {FirstActiveParticle = i;}
                if(prev >= 0) {NextActiveParticle[prev] = i;}
                prev = i;
            }
        }
    }

    if(prev >= 0) {NextActiveParticle[prev] = -1;}
}





/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
integertime find_next_outputtime(integertime ti_curr)
{
  int i, iter = 0;
  integertime ti, ti_next;
  double next, time;

  DumpFlag = 1;
  ti_next = -1;


  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn) {ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);}
          else {ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);}

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount) {All.SnapshotFileCount = i;}
		    }

		  if(ti_next > ti)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount) {All.SnapshotFileCount = i;}
		    }
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      time = All.TimeOfFirstSnapshot;

      iter = 0;

      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(110);
	    }
	}
      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn) {ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);}
        else {ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);}

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }


  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn) {next = All.TimeBegin * exp(ti_next * All.Timebase_interval);}
      else {next = All.TimeBegin + ti_next * All.Timebase_interval;}

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g  (DumpFlag=%d)\n", next, DumpFlag);

    }

  return ti_next;
}




/*! This routine writes for every synchronisation point in the timeline information to two log-files:
 * In FdInfo, we just list the timesteps that have been done, while in
 * FdTimebins we inform about the distribution of particles over the timebins, and which timebins are active on this step.
 * code is stored.
 */
void output_log_messages(void)
{
  double z;
  int i, j;
  long long tot, tot_sph;
  long long tot_count[TIMEBINS];
  long long tot_count_sph[TIMEBINS];
  long long tot_cumulative[TIMEBINS];
  int weight, corr_weight;
  double sum, avg_CPU_TimeBin[TIMEBINS], frac_CPU_TimeBin[TIMEBINS];

  sumup_large_ints(TIMEBINS, TimeBinCount, tot_count);
  sumup_large_ints(TIMEBINS, TimeBinCountSph, tot_count_sph);

#if defined(IO_SUPPRESS_TIMEBIN_STDOUT)
    if((ThisTask == 0) && (All.HighestActiveTimeBin>=(TIMEBINS-IO_SUPPRESS_TIMEBIN_STDOUT)))
#else
    if(ThisTask == 0)
#endif
    {
        if(All.ComovingIntegrationOn)
        {
            z = 1.0 / (All.Time) - 1;
#ifndef IO_REDUCED_MODE
            fprintf(FdInfo, "Sync-Point %lld, Time: %g, Redshift: %g, Nf = %d%09d, Systemstep: %g, Dloga: %g\n",
                    (long long) All.NumCurrentTiStep, All.Time, z, (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
            fflush(FdInfo);
            fprintf(FdTimebin, "Sync-Point %lld, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", (long long) All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
#endif
            printf("\nSync-Point %lld, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", (long long) All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
        }
        else
        {
#ifndef IO_REDUCED_MODE
            fprintf(FdInfo, "Sync-Point %lld, Time: %g, Nf = %d%09d, Systemstep: %g\n", (long long) All.NumCurrentTiStep,
                    All.Time, (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), All.TimeStep);
            fflush(FdInfo);
            fprintf(FdTimebin, "Sync-Point %lld, Time: %g, Systemstep: %g\n", (long long) All.NumCurrentTiStep, All.Time, All.TimeStep);
#endif
            printf("\nSync-Point %lld, Time: %g, Systemstep: %g\n", (long long) All.NumCurrentTiStep, All.Time, All.TimeStep);
        }

        for(i = 1, tot_cumulative[0] = tot_count[0]; i < TIMEBINS; i++) {tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];}


      for(i = 0; i < TIMEBINS; i++)
	{
	  for(j = 0, sum = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++) {sum += All.CPU_TimeBinMeasurements[i][j];}
	  if(All.CPU_TimeBinCountMeasurements[i]) {avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];} else {avg_CPU_TimeBin[i] = 0;}
	}

      for(i = All.HighestOccupiedTimeBin, weight = 1, sum = 0; i >= 0 && tot_count[i] > 0; i--, weight *= 2)
	{
	  if(weight > 1) {corr_weight = weight / 2;} else {corr_weight = weight;}
	  frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
	  sum += frac_CPU_TimeBin[i];
	}

      for(i = All.HighestOccupiedTimeBin; i >= 0 && tot_count[i] > 0; i--) {if(sum) {frac_CPU_TimeBin[i] /= sum;}}


        printf("Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
#ifndef IO_REDUCED_MODE
        fprintf(FdTimebin,"Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
#endif
        for(i = TIMEBINS - 1, tot = tot_sph = 0; i >= 0; i--)
            if(tot_count_sph[i] > 0 || tot_count[i] > 0)
            {
                printf(" %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n", TimeBinActive[i] ? 'X' : ' ', i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
                       GET_INTEGERTIME_FROM_TIMEBIN(i) * All.Timebase_interval, tot_cumulative[i], (i == All.HighestActiveTimeBin) ? '<' : ' ',
                       (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ', avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);
#ifndef IO_REDUCED_MODE
                fprintf(FdTimebin," %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n", TimeBinActive[i] ? 'X' : ' ', i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
                        GET_INTEGERTIME_FROM_TIMEBIN(i) * All.Timebase_interval, tot_cumulative[i], (i == All.HighestActiveTimeBin) ? '<' : ' ',
                        (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ', avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);
#endif
                if(TimeBinActive[i])
                {
                    tot += tot_count[i];
                    tot_sph += tot_count_sph[i];
                }
            }
        printf("               ------------------------\n");
#ifndef IO_REDUCED_MODE
        fprintf(FdTimebin, "               ------------------------\n");
#endif
#ifdef PMGRID
        if(All.PM_Ti_endstep == All.Ti_Current)
        {
            printf("PM-Step. Total: %10llu  %10llu    Sum: %10llu\n\n", tot - tot_sph, tot_sph, tot);
#ifndef IO_REDUCED_MODE
            fprintf(FdTimebin, "PM-Step. Total: %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
#endif
        }
        else
#endif
        {
            printf("Total active:   %10llu  %10llu    Sum: %10llu\n\n", tot - tot_sph, tot_sph, tot);
#ifndef IO_REDUCED_MODE
            fprintf(FdTimebin, "Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
#endif
        }
#ifndef IO_REDUCED_MODE
        fprintf(FdTimebin, "\n");
        fflush(FdTimebin);
#endif
    }

  output_extra_log_messages();
}




void write_cpu_log(void)
{
  double max_CPU_Step[CPU_PARTS], avg_CPU_Step[CPU_PARTS], t0, t1, tsum; int i; t0=0; t1=0; tsum=0;
  CPU_Step[CPU_MISC] += measure_time();

  for(i = 1, CPU_Step[0] = 0; i < CPU_PARTS; i++) {CPU_Step[0] += CPU_Step[i];}

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++) {avg_CPU_Step[i] /= NTask;}

#ifndef IO_REDUCED_MODE
      put_symbol(0.0, 1.0, '#');
      for(i = 1, tsum = 0.0; i < CPU_PARTS; i++)
      {
            if(max_CPU_Step[i] > 0)
            {
              t0 = tsum; t1 = t0 + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
              put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_Symbol[i]);
              tsum += t1 - t0;

              t0 = tsum; t1 = t0 + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
              put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_SymbolImbalance[i]);
              tsum += t1 - t0;
            }
      }
      put_symbol(tsum / max_CPU_Step[0], 1.0, '-');
      fprintf(FdBalance, "Step=%7lld  sec=%10.3f  Nf=%2d%09d  %s\n", (long long) All.NumCurrentTiStep, max_CPU_Step[0], (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), CPU_String); fflush(FdBalance);
#endif

      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
	{
	  All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
	  memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0], &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1], (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
	}

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]++] = max_CPU_Step[0];
    }

    CPUThisRun += CPU_Step[0];

    for(i = 0; i < CPU_PARTS; i++) {CPU_Step[i] = 0;}
    if(ThisTask == 0)
    {
        for(i = 0; i < CPU_PARTS; i++) {All.CPU_Sum[i] += avg_CPU_Step[i];}
    }

#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) // only do the actual -print- operation on global timesteps
#endif
  if(ThisTask == 0)
    {
      fprintf(FdCPU, "Step %lld, Time: %g, CPUs: %d\n",(long long) All.NumCurrentTiStep, All.Time, NTask);
      fprintf(FdCPU, "Nactive=%lld, Imbal(Max/Mean)=%g \n", (long long) GlobNumForceUpdate, (max_CPU_Step[0]/(MIN_REAL_NUMBER + avg_CPU_Step[0])-1.)*NTask+1.);
      fprintf(FdCPU,
	      "total         %10.2f  %5.1f%%\n"
	      "tree+gravity  %10.2f  %5.1f%%\n"
	      "   treebuild  %10.2f  %5.1f%%\n"
	      "   treewalk   %10.2f  %5.1f%%\n"
	      "   treecomm   %10.2f  %5.1f%%\n"
	      "   treeimbal  %10.2f  %5.1f%%\n"
#ifdef PMGRID
          "pm-gravity    %10.2f  %5.1f%%\n"
#endif
#if !defined(EVALPOTENTIAL) && (defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL))
          "potentialeval %10.2f  %5.1f%%\n"
#endif
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
	      "ags-nongas    %10.2f  %5.1f%%\n"
	      "   agsdensity %10.2f  %5.1f%%\n"
	      "   agscomm    %10.2f  %5.1f%%\n"
	      "   agsimbal   %10.2f  %5.1f%%\n"
          "   agsmisc    %10.2f  %5.1f%%\n"
#endif
#ifdef TURB_DIFF_DYNAMIC
          "dyndiff       %10.2f  %5.1f%%\n"
          "   compute    %10.2f  %5.1f%%\n"
          "   comm       %10.2f  %5.1f%%\n"
          "   wait       %10.2f  %5.1f%%\n"
          "   misc       %10.2f  %5.1f%%\n"
          "velsmooth     %10.2f  %5.1f%%\n"
          "   compute    %10.2f  %5.1f%%\n"
          "   comm       %10.2f  %5.1f%%\n"
          "   wait       %10.2f  %5.1f%%\n"
          "   misc       %10.2f  %5.1f%%\n"
#endif
	      "hydro/fluids  %10.2f  %5.1f%%\n"
	      "   dens+grad  %10.2f  %5.1f%%\n"
	      "   denscomm   %10.2f  %5.1f%%\n"
	      "   densimbal  %10.2f  %5.1f%%\n"
	      "   hydrofrc   %10.2f  %5.1f%%\n"
	      "   hydcomm    %10.2f  %5.1f%%\n"
	      "   hydimbal   %10.2f  %5.1f%%\n"
	      "   hmaxupdate %10.2f  %5.1f%%\n"
          "   hydmisc    %10.2f  %5.1f%%\n"
	      "domain        %10.2f  %5.1f%%\n"
          "peano         %10.2f  %5.1f%%\n"
#ifdef FOF
          "fof/subfind   %10.2f  %5.1f%%\n"
#endif
          "drift/splitmg %10.2f  %5.1f%%\n"
	      "kicks         %10.2f  %5.1f%%\n"
	      "io/snapshots  %10.2f  %5.1f%%\n"
#ifdef COOLING
	      "cooling+chem  %10.2f  %5.1f%%\n"
#endif
#ifdef CHIMES
	      " coolchmimbal %10.2f  %5.1f%%\n"
#endif
#ifdef BLACK_HOLES
	      "blackholes    %10.2f  %5.1f%%\n"
#endif
#ifdef GRAIN_FLUID
          "grains        %10.2f  %5.1f%%\n"
#endif
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
          "mech_fb_loop  %10.2f  %5.1f%%\n"
#endif
#if defined(RADTRANSFER)
          "rt_nonfluxops %10.2f  %5.1f%%\n"
#endif
          "misc          %10.2f  %5.1f%%\n",

    All.CPU_Sum[CPU_ALL], 100.0,
    All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2] + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
              + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2] + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEMISC],
    (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2] + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
              + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2] + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEBUILD], (All.CPU_Sum[CPU_TREEBUILD]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2], (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV], (All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2], (All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef PMGRID
    All.CPU_Sum[CPU_MESH], (All.CPU_Sum[CPU_MESH]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#if !defined(EVALPOTENTIAL) && (defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL))
    All.CPU_Sum[CPU_POTENTIAL], (All.CPU_Sum[CPU_POTENTIAL]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    All.CPU_Sum[CPU_AGSDENSCOMPUTE] + All.CPU_Sum[CPU_AGSDENSWAIT] + All.CPU_Sum[CPU_AGSDENSCOMM] + All.CPU_Sum[CPU_AGSDENSMISC],
              (All.CPU_Sum[CPU_AGSDENSCOMPUTE] + All.CPU_Sum[CPU_AGSDENSWAIT] + All.CPU_Sum[CPU_AGSDENSCOMM] + All.CPU_Sum[CPU_AGSDENSMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_AGSDENSCOMPUTE], (All.CPU_Sum[CPU_AGSDENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_AGSDENSCOMM], (All.CPU_Sum[CPU_AGSDENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_AGSDENSWAIT], (All.CPU_Sum[CPU_AGSDENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_AGSDENSMISC], (All.CPU_Sum[CPU_AGSDENSMISC]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#ifdef TURB_DIFF_DYNAMIC
    (All.CPU_Sum[CPU_DYNDIFFCOMPUTE] + All.CPU_Sum[CPU_DYNDIFFWAIT] + All.CPU_Sum[CPU_DYNDIFFCOMM] + All.CPU_Sum[CPU_DYNDIFFMISC]), (All.CPU_Sum[CPU_DYNDIFFCOMPUTE] + All.CPU_Sum[CPU_DYNDIFFWAIT] + All.CPU_Sum[CPU_DYNDIFFCOMM] + All.CPU_Sum[CPU_DYNDIFFMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DYNDIFFCOMPUTE], (All.CPU_Sum[CPU_DYNDIFFCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DYNDIFFWAIT], (All.CPU_Sum[CPU_DYNDIFFWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DYNDIFFCOMM], (All.CPU_Sum[CPU_DYNDIFFCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DYNDIFFMISC], (All.CPU_Sum[CPU_DYNDIFFMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    (All.CPU_Sum[CPU_IMPROVDIFFCOMPUTE] + All.CPU_Sum[CPU_IMPROVDIFFWAIT] + All.CPU_Sum[CPU_IMPROVDIFFCOMM] + All.CPU_Sum[CPU_IMPROVDIFFMISC]), (All.CPU_Sum[CPU_IMPROVDIFFCOMPUTE] + All.CPU_Sum[CPU_IMPROVDIFFWAIT] + All.CPU_Sum[CPU_IMPROVDIFFCOMM] + All.CPU_Sum[CPU_IMPROVDIFFMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_IMPROVDIFFCOMPUTE], (All.CPU_Sum[CPU_IMPROVDIFFCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_IMPROVDIFFWAIT], (All.CPU_Sum[CPU_IMPROVDIFFWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_IMPROVDIFFCOMM], (All.CPU_Sum[CPU_IMPROVDIFFCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_IMPROVDIFFMISC], (All.CPU_Sum[CPU_IMPROVDIFFMISC]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
    All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSWAIT] + All.CPU_Sum[CPU_DENSMISC]
              + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC]
              + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE],
    (All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSWAIT] + All.CPU_Sum[CPU_DENSMISC]
              + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC]
              + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DENSCOMPUTE], (All.CPU_Sum[CPU_DENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DENSCOMM], (All.CPU_Sum[CPU_DENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DENSWAIT], (All.CPU_Sum[CPU_DENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDCOMPUTE], (All.CPU_Sum[CPU_HYDCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDCOMM], (All.CPU_Sum[CPU_HYDCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDWAIT], (All.CPU_Sum[CPU_HYDWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEHMAXUPDATE], (All.CPU_Sum[CPU_TREEHMAXUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDMISC] + All.CPU_Sum[CPU_DENSMISC], (All.CPU_Sum[CPU_HYDMISC] + All.CPU_Sum[CPU_DENSMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DOMAIN], (All.CPU_Sum[CPU_DOMAIN]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_PEANO], (All.CPU_Sum[CPU_PEANO]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef FOF
    All.CPU_Sum[CPU_FOF], (All.CPU_Sum[CPU_FOF]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
    All.CPU_Sum[CPU_DRIFT], (All.CPU_Sum[CPU_DRIFT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TIMELINE], (All.CPU_Sum[CPU_TIMELINE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_SNAPSHOT], (All.CPU_Sum[CPU_SNAPSHOT]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef COOLING
    All.CPU_Sum[CPU_COOLINGSFR], (All.CPU_Sum[CPU_COOLINGSFR]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#ifdef CHIMES
    All.CPU_Sum[CPU_COOLSFRIMBAL], (All.CPU_Sum[CPU_COOLSFRIMBAL]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#ifdef BLACK_HOLES
    All.CPU_Sum[CPU_BLACKHOLES], (All.CPU_Sum[CPU_BLACKHOLES]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#ifdef GRAIN_FLUID
    All.CPU_Sum[CPU_DRAGFORCE], (All.CPU_Sum[CPU_DRAGFORCE]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
    All.CPU_Sum[CPU_SNIIHEATING], (All.CPU_Sum[CPU_SNIIHEATING]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#if defined(RADTRANSFER)
    All.CPU_Sum[CPU_RTNONFLUXOPS], (All.CPU_Sum[CPU_RTNONFLUXOPS]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
    All.CPU_Sum[CPU_MISC], (All.CPU_Sum[CPU_MISC]) / All.CPU_Sum[CPU_ALL] * 100);

    fprintf(FdCPU, "\n");
    fflush(FdCPU);
    }
}


void put_symbol(double t0, double t1, char c)
{
    int i, j;
    i = (int) (t0 * CPU_STRING_LEN + 0.5);
    j = (int) (t1 * CPU_STRING_LEN);
    if(i < 0) {i = 0;}
    if(j < 0) {j = 0;}
    if(i >= CPU_STRING_LEN) {i = CPU_STRING_LEN;}
    if(j >= CPU_STRING_LEN) {j = CPU_STRING_LEN;}
    while(i <= j) {CPU_String[i++] = c;}
    CPU_String[CPU_STRING_LEN] = 0;
}




/*! This routine first calls a computation of various global
 * quantities of the particle distribution, and then writes some
 * statistics about the energies in the various particle components to
 * the file FdEnergy.
 */
void energy_statistics(void)
{
#ifndef IO_REDUCED_MODE
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

      fprintf(FdEnergy," \n");
      fflush(FdEnergy);
    }
#endif
}



void output_extra_log_messages(void)
{
#if defined(TURB_DRIVING) && !defined(IO_REDUCED_MODE)
    log_turb_temp();
#endif

#if defined(GR_TABULATED_COSMOLOGY) && !defined(IO_REDUCED_MODE)
    if((ThisTask == 0) && (All.ComovingIntegrationOn == 1)
    {
        double hubble_a;

        hubble_a = hubble_function(All.Time);
        fprintf(FdDE, "%lld %g %e ", (long long) All.NumCurrentTiStep, All.Time, hubble_a);
#ifndef GR_TABULATED_COSMOLOGY_W
        fprintf(FdDE, "%e ", All.DarkEnergyConstantW);
#else
        fprintf(FdDE, "%e %e ", get_wa(All.Time), DarkEnergy_a(All.Time));
#endif
#ifdef GR_TABULATED_COSMOLOGY_G
        fprintf(FdDE, "%e %e", dHfak(All.Time), dGfak(All.Time));
#endif
        fprintf(FdDE, "\n");
        fflush(FdDE);
    }
#endif
}
