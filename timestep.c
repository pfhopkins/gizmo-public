#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "kernel.h"

/*! \file timestep.c
 *  routines for assigning new timesteps
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO; these 
 * modifications include the addition of various timestep criteria, the WAKEUP 
 * additions, and various changes of units and variable naming conventions throughout.
 */

static double dt_displacement = 0;

#ifdef BOX_SHEARING
void calc_shearing_box_pos_offset(void)
{
    Shearing_Box_Pos_Offset = Shearing_Box_Vel_Offset * All.Time;
    while(Shearing_Box_Pos_Offset > boxSize_Y) {Shearing_Box_Pos_Offset -= boxSize_Y;}
}
#endif

void set_cosmo_factors_for_current_time(void)
{
    
    /* These are critical factors used throughout for co-moving integrations. Set them here and
       call THESE, rather than trying to come up with the factors throughout, since that makes debugging a nightmare */
    if(All.ComovingIntegrationOn)
    {
        /* All.cf_atime = a = 1/(1+z), the cosmological scale factor */
        All.cf_atime = All.Time;
        /* All.cf_a2inv is just handy */
        All.cf_a2inv = 1 / (All.Time * All.Time);
        /* All.cf_a3inv * Density_code = Density_physical */
        All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
        /* Pressure_code/Density_code = All.cf_afac1 * Pressure_physical/Density_physical */
        All.cf_afac1 = 1;
        /* All.cf_afac2 * Pressure_code/Density_code * 1/r_code = Pressure_physical/Density_physical * 1/r_physical */
        All.cf_afac2 = 1 / (All.Time * All.cf_afac1);
        /* All.cf_afac3 * sqrt(Pressure_code/Density_code) = sqrt(Pressure_phys/Density_phys) = cs_physical */
        All.cf_afac3 = 1 / sqrt(All.cf_afac1);
        /* time units: proper time dt_phys = 1/hubble_function(a) * dz/(1+z) = dlna / hubble_function(a)
            code time unit in comoving is dlna, so dt_phys = dt_code / All.cf_hubble_a   */
        All.cf_hubble_a = hubble_function(All.Time); /* hubble_function(a) = H(a) = H(z) */
        /* dt_code * v_code/r_code = All.cf_hubble_a2 * dt_phys * v_phys/r_phys */
        All.cf_hubble_a2 = All.Time * All.Time * hubble_function(All.Time);
    }
    else
    {
        All.cf_atime = 1;
        All.cf_a2inv = 1;
        All.cf_a3inv = 1;
        All.cf_afac1 = 1;
        All.cf_afac2 = 1;
        All.cf_afac3 = 1;
        All.cf_hubble_a = 1;
        All.cf_hubble_a2 = 1;
    }
}


/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void find_timesteps(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    
    int i, bin, binold, prev, next;
    integertime ti_step, ti_step_old, ti_min;
    double aphys;
    
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin || dt_displacement == 0)
        find_dt_displacement_constraint(All.cf_hubble_a * All.cf_atime * All.cf_atime);
     
#ifdef DIVBCLEANING_DEDNER
    /* need to calculate the global fastest wave speed to manage the damping terms stably */
    if((All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)||(All.FastestWaveSpeed == 0))
    {
        double fastwavespeed = 0.0;
        double fastwavedecay = 0.0;
        double fac_magnetic_pressure = All.cf_afac1 / All.cf_atime;
        for(i=0;i<NumPart;i++)
        {
            if(P[i].Type==0)
            {
                double vsig2 = 0.5 * All.cf_afac3 * fabs(SphP[i].MaxSignalVel); // in v_phys units //
                double vsig1 = All.cf_afac3 * sqrt( Get_Gas_effective_soundspeed_i(i)*Get_Gas_effective_soundspeed_i(i) + fac_magnetic_pressure * (Get_Gas_BField(i,0)*Get_Gas_BField(i,0)+Get_Gas_BField(i,1)*Get_Gas_BField(i,1)+Get_Gas_BField(i,2)*Get_Gas_BField(i,2)) / SphP[i].Density );
                double vsig0 = DMAX(vsig1,vsig2);

                if(vsig0 > fastwavespeed) fastwavespeed = vsig0; // physical unit
                double hsig0 = Get_Particle_Size(i) * All.cf_atime; // physical unit
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

#ifdef FORCE_EQUAL_TIMESTEPS
    for(i = FirstActiveParticle, ti_min = TIMEBASE; i >= 0; i = NextActiveParticle[i])
    {
        ti_step = get_timestep(i, &aphys, 0);
        if(ti_step < ti_min) {ti_min = ti_step;}
    }
    if(ti_min > (dt_displacement / All.Timebase_interval)) {ti_min = (dt_displacement / All.Timebase_interval);}
    
    ti_step = TIMEBASE;
    while(ti_step > ti_min) {ti_step >>= 1;}
    integertime ti_min_glob;
    MPI_Allreduce(&ti_step, &ti_min_glob, 1, MPI_TYPE_TIME, MPI_MIN, MPI_COMM_WORLD);
#endif
    
    
    /* Now assign new timesteps  */
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef FORCE_EQUAL_TIMESTEPS
        ti_step = ti_min_glob;
#else
        ti_step = get_timestep(i, &aphys, 0);
#endif

        /* make it a power 2 subdivision */
        ti_min = TIMEBASE;
        while(ti_min > ti_step)
            ti_min >>= 1;
        ti_step = ti_min;
        
        bin = get_timestep_bin(ti_step);
        binold = P[i].TimeBin;
        
        if(bin > binold)		/* timestep wants to increase */
        {
            while(TimeBinActive[bin] == 0 && bin > binold)	/* make sure the new step is synchronized */
                bin--;
            
            ti_step = bin ? (((integertime) 1) << bin) : 0;
        }
        
        if(All.Ti_Current >= TIMEBASE)	/* we here finish the last timestep. */
        {
            ti_step = 0;
            bin = 0;
        }
        
        if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
        {
            terminate("we are beyond the end of the timeline");	/* should not happen */
            ti_step = TIMEBASE - All.Ti_Current;
            ti_min = TIMEBASE;
            while(ti_min > ti_step)
                ti_min >>= 1;
            ti_step = ti_min;
        }
        
        if(bin != binold)
        {
            TimeBinCount[binold]--;
            if(P[i].Type == 0)
            {
                TimeBinCountSph[binold]--;
#ifdef GALSF
                TimeBinSfr[binold] -= SphP[i].Sfr;
                TimeBinSfr[bin] += SphP[i].Sfr;
#endif
            }
            
#ifdef BLACK_HOLES
            if(P[i].Type == 5)
            {
                TimeBin_BH_mass[binold] -= BPP(i).BH_Mass;
                TimeBin_BH_dynamicalmass[binold] -= P[i].Mass;
                TimeBin_BH_Mdot[binold] -= BPP(i).BH_Mdot;
                if(BPP(i).BH_Mass > 0)
                    TimeBin_BH_Medd[binold] -= BPP(i).BH_Mdot / BPP(i).BH_Mass;
                TimeBin_BH_mass[bin] += BPP(i).BH_Mass;
                TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
                TimeBin_BH_Mdot[bin] += BPP(i).BH_Mdot;
                if(BPP(i).BH_Mass > 0)
                    TimeBin_BH_Medd[bin] += BPP(i).BH_Mdot / BPP(i).BH_Mass;
            }
#endif
            
            prev = PrevInTimeBin[i];
            next = NextInTimeBin[i];
            
            if(FirstInTimeBin[binold] == i)
                FirstInTimeBin[binold] = next;
            if(LastInTimeBin[binold] == i)
                LastInTimeBin[binold] = prev;
            if(prev >= 0)
                NextInTimeBin[prev] = next;
            if(next >= 0)
                PrevInTimeBin[next] = prev;
            
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
            if(P[i].Type == 0)
                TimeBinCountSph[bin]++;
            
            P[i].TimeBin = bin;
        }
        
#ifndef WAKEUP
        ti_step_old = binold ? (((integertime) 1) << binold) : 0;
#else
        ti_step_old = P[i].dt_step;
#endif
        
        P[i].Ti_begstep += ti_step_old;
        
#if defined(WAKEUP)
        P[i].dt_step = ti_step;
#endif
    }
    
    
    
#ifdef PMGRID
    if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = TIMEBASE;
        while(ti_step > (dt_displacement / All.Timebase_interval))
            ti_step >>= 1;
        
        if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
        {
            bin = get_timestep_bin(ti_step);
            binold = get_timestep_bin(All.PM_Ti_endstep - All.PM_Ti_begstep);
            
            while(TimeBinActive[bin] == 0 && bin > binold)	/* make sure the new step is synchronized */
                bin--;
            
            ti_step = bin ? (((integertime) 1) << bin) : 0;
        }
        
        if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
            ti_step = 0;
        
        All.PM_Ti_begstep = All.PM_Ti_endstep;
        All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;
    }
#endif
    
    
#ifdef WAKEUP
    process_wake_ups();
#endif
    
    CPU_Step[CPU_TIMELINE] += measure_time();
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

#if (SINGLE_STAR_TIMESTEPPING > 0)
    P[p].SuperTimestepFlag = 0;
    if( (P[p].Type==5) && P[p].is_in_a_binary ) // candidate: need to decide whether to use super timestepping for binaries
    { 
#if (SINGLE_STAR_TIMESTEPPING == 1) // to be conservative, use the semimajor axis, ie. the internal timescale is the orbital period
	    double dt_bin = P[p].min_bh_t_orbital / (2.*M_PI); // sqrt(a^3/GM) for binary
	    if(0.03*P[p].COM_dt_tidal>dt_bin) {P[p].SuperTimestepFlag=2;
	    } // external timestep is appropriately larger than 'internal' timestep, so use super-timestepping routine
#else // to be more aggressive, use the instantaneous orbital timescale, ie. freefall time from the CURRENT orbital separation. This lets us super step an orbit on the close passages, even when it is affected by tides at apopase
	    double dr = sqrt(P[p].comp_dx[0]*P[p].comp_dx[0] + P[p].comp_dx[1]*P[p].comp_dx[1] + P[p].comp_dx[2]*P[p].comp_dx[2]);
	    double dt_bin = sqrt(dr*dr*dr / (All.G * (P[p].Mass + P[p].comp_Mass)));
        if(0.005*P[p].COM_dt_tidal>dt_bin) {P[p].SuperTimestepFlag=2;} // external timestep is appropriately larger than 'internal' timestep, so use super-timestepping routine [constant here stricter for more aggressive routine]
#endif
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
#if defined(RT_USE_GRAVTREE) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
        if(P[p].Type>0) // strictly this is better for accuracy, but not necessary
#endif
        ax = ay = az = 0.0; // we're getting our gravitational timestep criterion from the tidal tensor, but still want to do the accel criterion for other forces
#endif

        if(P[p].Type == 0)
        {
            ax += SphP[p].HydroAccel[0];
            ay += SphP[p].HydroAccel[1];
            az += SphP[p].HydroAccel[2];
#ifdef TURB_DRIVING
            ax += SphP[p].TurbAccel[0];
            ay += SphP[p].TurbAccel[1];
            az += SphP[p].TurbAccel[2];
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            ax += SphP[p].Rad_Accel[0];
            ay += SphP[p].Rad_Accel[1];
            az += SphP[p].Rad_Accel[2];
#endif
        }

        ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
#ifdef TURB_DRIVING
        if(P[p].Type==0)
        { /* because the turbulent acceleration is a random variable, we dont want it to catch us by surprise if it moves up, so
            we include a safety factor here which (very crudely) approximates the maximum amplitude it could reach */
            double a_max_safety = 1.4 * sqrt(pow(All.StKmax,NUMDIMS) * All.StEnergy / All.StDecay);
            ac = sqrt(ac*ac + a_max_safety*a_max_safety);
        }
#endif
        
        *aphys = ac;
    }
    else
    {ac = *aphys;}
    
    if(ac == 0) ac = 1.0e-30;
    
    
    if(flag > 0)
    {
        /* this is the non-standard mode; use timestep to get the maximum acceleration tolerated */
        dt = flag * All.Timebase_interval;
        dt /= All.cf_hubble_a;	/* convert dloga to physical timestep  */
        
        ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * All.ForceSoftening[P[p].Type] / (dt * dt);
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * PPP[p].AGS_Hsml / (dt * dt);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
        if(P[p].Type==0) {ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * PPP[p].Hsml / (dt * dt);}
#endif
        *aphys = ac;
        return flag;
    }

    dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * All.ForceSoftening[P[p].Type] / ac);

#ifdef ADAPTIVE_GRAVSOFT_FORALL
    dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime  * KERNEL_CORE_SIZE * PPP[p].AGS_Hsml / ac);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
    if(P[p].Type == 0) {dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * PPP[p].Hsml / ac);}
#endif

#if (defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)) && defined(GALSF) && defined(GALSF_FB_MECHANICAL)
    if(((P[p].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[p].Type == 2)||(P[p].Type==3))))&&(P[p].Mass>0))
    {
        if((All.ComovingIntegrationOn))
        {
            double ags_h = DMAX(PPP[p].Hsml, All.ForceSoftening[P[p].Type]);
            ags_h = DMIN(ags_h, 10.*All.ForceSoftening[P[p].Type]);
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            ags_h = DMAX(PPP[p].AGS_Hsml , DMAX(PPP[p].Hsml,All.ForceSoftening[P[p].Type]));
            ags_h = DMIN(ags_h, DMAX(100.*All.ForceSoftening[P[p].Type] , 10.*PPP[p].AGS_Hsml));
#endif
            dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime  * KERNEL_CORE_SIZE * ags_h / ac);
        }
    }
#endif

#ifdef TIDAL_TIMESTEP_CRITERION // tidal criterion obtains the same energy error in an optimally-softened Plummer sphere over ~100 crossing times as the Power 2003 criterion
    double dt_tidal = 0.; {int k; for(k=0; k<3; k++) {dt_tidal += P[p].tidal_tensorps[k][k]*P[p].tidal_tensorps[k][k];}} // this is diagonalized already in the gravity loop
    dt_tidal = sqrt(2. * All.ErrTolIntAccuracy / sqrt(dt_tidal / 6)); // recovers sqrt(eta) * tdyn for a Keplerian potential
#if (SINGLE_STAR_TIMESTEPPING > 0)
    if(P[p].SuperTimestepFlag>=2) {dt_tidal = sqrt(2*All.ErrTolIntAccuracy) * P[p].COM_dt_tidal;}
#endif
    dt=DMIN(dt,dt_tidal);
#endif
    
#ifdef SINGLE_STAR_TIMESTEPPING // this ensures that binaries advance in lock-step, which gives superior conservation
    if(P[p].Type == 5)
    {
        double dt_2body = sqrt(2*All.ErrTolIntAccuracy) * 0.3 / (1./P[p].min_bh_approach_time + 1./P[p].min_bh_freefall_time); // timestep is harmonic mean of freefall and approach time
#ifdef HERMITE_INTEGRATION
        if(eligible_for_hermite(p)) dt_2body /= 0.3;
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    	if(P[p].is_in_a_binary && (P[p].SuperTimestepFlag >= 2)) //binary candidate or a confirmed binary
	    {    // First we need to construct the same 2-body timescale as above, but from the binary parameters. If this is longer than the above, there is another star that is requiring us to
	         // take a short timestep, so we better not super-timestep otherwise we risk messing up that star's integration. But if it is consistent with the above, then we can safely super-timestep
	        double Mtot=P[p].comp_Mass+P[p].Mass, dr=0,dv=0,dv_dot_dx=0, binary_dt_2body=0;
	        for(k=0;k<3;k++) {dr+=P[p].comp_dx[k]*P[p].comp_dx[k]; dv+=P[p].comp_dv[k]*P[p].comp_dv[k]; dv_dot_dx+=P[p].comp_dx[k]*P[p].comp_dv[k];}
	        dr += All.SofteningTable[5]*All.SofteningTable[5]; dr=sqrt(dr); if(dv>0) {dv=sqrt(dv);} else {dv=0;}
            double dt_2body_base = 1/(1./P[p].min_bh_approach_time + 1./P[p].min_bh_freefall_time); // timestep is harmonic mean of freefall and approach time
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
        if(eligible_for_hermite(p)) dt *= 1.4; // gives 10^-6 energy error per orbit for a 0.9 eccentricity binary
#endif
    }
#endif // SINGLE_STAR_TIMESTEPPING
    
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    /* make sure smoothing length of non-gas particles doesn't change too much in one timestep */
    if(((1 << P[p].Type) & (ADAPTIVE_GRAVSOFT_FORALL)) && (P[p].Type > 0))
    {
        double dt_divv = 0.1 / (MIN_REAL_NUMBER + All.cf_a2inv*fabs(P[p].Particle_DivVel)); // with new integration accuracy in gravtree, we may not need to be super-conservative here. old code used pre-factor 0.25 here, see if we can get away with the larger value which is standard for gas below
        if(dt_divv < dt) {dt = dt_divv;}
        double dt_cour = 2. * All.CourantFac * (Get_Particle_Size_AGS(p)*All.cf_atime) / (MIN_REAL_NUMBER + 0.5*P[p].AGS_vsig*All.cf_afac3); // can be generous here, really the signal velocity isn't that important in the collisionless case, but it is important with some of the physics above //
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


#ifdef GRAIN_FLUID
    if((1 << P[p].Type) & (GRAIN_PTYPES))
    {
        csnd = convert_internalenergy_soundspeed2(p, P[p].Gas_InternalEnergy);
        int k; for(k=0;k<3;k++) {csnd += (P[p].Gas_Velocity[k]-P[p].Vel[k])*(P[p].Gas_Velocity[k]-P[p].Vel[k]);}
#if defined(GRAIN_LORENTZFORCE)
        for(k=0;k<3;k++) {csnd += P[p].Gas_B[k]*P[p].Gas_B[k] / (2.0 * P[p].Gas_Density);}
#endif
        csnd = sqrt(csnd);
        double L_particle = Get_Particle_Size(p);
        dt_courant = 0.5 * All.CourantFac * (L_particle*All.cf_atime) / csnd;
#if defined(GRAIN_BACKREACTION)
        if(P[p].Grain_AccelTimeMin < dt_courant) {dt_courant = P[p].Grain_AccelTimeMin;}
#endif
#ifdef PIC_MHD
        if(P[p].Grain_SubType==3)
        {
            double lorentz_units = UNIT_B_IN_GAUSS; // code B to Gauss
            lorentz_units *= UNIT_VEL_IN_CGS * (ELECTRONCHARGE/(PROTONMASS*C_LIGHT)); // code velocity to CGS, times base units e/(mp*c)
            lorentz_units /= UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS; // convert 'back' to code-units acceleration
            double reduced_C = PIC_SPEEDOFLIGHT_REDUCTION * C_LIGHT_CODE;
            double charge_to_mass_ratio_dimensionless = All.PIC_Charge_to_Mass_Ratio;
            double v2=0, B2=0; for(k=0;k<3;k++) {v2=P[p].Vel[k]*P[p].Vel[k]; B2+=P[p].Gas_B[k]*P[p].Gas_B[k];}
            double gamma=1/sqrt(1-v2/(reduced_C*reduced_C));
            double dt_courant_pic = 0.5 * gamma / (charge_to_mass_ratio_dimensionless * sqrt(B2) * lorentz_units); /* dt = 0.5/omega_g */
            if(dt_courant_pic < dt_courant) dt_courant = dt_courant_pic;
        }
#endif
        if(dt_courant < dt) dt = dt_courant;
    }
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
    if(P[p].Type>0) {double dt_inj = 0.1 * PPP[p].Hsml / C_LIGHT_CODE_REDUCED; if(dt_inj < dt) {dt = dt_inj;}}
#endif
#endif
    
    
    if((P[p].Type == 0) && (P[p].Mass > 0))
        {
            csnd = 0.5 * SphP[p].MaxSignalVel * All.cf_afac3;
            double L_particle = Get_Particle_Size(p);
            
            dt_courant = All.CourantFac * (L_particle*All.cf_atime) / csnd;
            if(dt_courant < dt) dt = dt_courant;

            double dt_prefac_diffusion;
            dt_prefac_diffusion = 0.5;
#if defined(GALSF) || defined(DIFFUSION_OPTIMIZERS)
            dt_prefac_diffusion = 1.8;
#endif
#ifdef SUPER_TIMESTEP_DIFFUSION
            double dt_superstep_explicit = 1.e10 * dt;
#endif
            
            
            
#ifdef CONDUCTION
            {
                double L_cond_inv = sqrt(SphP[p].Gradients.InternalEnergy[0]*SphP[p].Gradients.InternalEnergy[0] +
                                         SphP[p].Gradients.InternalEnergy[1]*SphP[p].Gradients.InternalEnergy[1] +
                                         SphP[p].Gradients.InternalEnergy[2]*SphP[p].Gradients.InternalEnergy[2]) / SphP[p].InternalEnergy;
                double L_cond = DMAX(L_particle , 1./(L_cond_inv + 1./L_particle)) * All.cf_atime;
                double dt_conduction = dt_prefac_diffusion * L_cond*L_cond / (MIN_REAL_NUMBER + SphP[p].Kappa_Conduction);
                // since we use CONDUCTIVITIES, not DIFFUSIVITIES, we need to add a power of density to get the right units //
                dt_conduction *= SphP[p].Density * All.cf_a3inv;
#ifdef SUPER_TIMESTEP_DIFFUSION
                if(dt_conduction < dt_superstep_explicit) dt_superstep_explicit = dt_conduction; // explicit time-step
                double dt_advective = dt_conduction * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_cond_inv))*All.cf_atime / L_cond);
                if(dt_advective < dt) dt = dt_advective; // 'advective' timestep: needed to limit super-stepping
#else
                if(dt_conduction < dt) dt = dt_conduction; // normal explicit time-step
#endif
            }
#endif


#ifdef MHD_NON_IDEAL
            {
                double b_grad = 0, b_mag = 0;
                int k; for(k=0;k<3;k++)
                {
                    int k2;
                    for(k2=0;k2<3;k2++)
                    {
                        double tmp_grad = SphP[p].Gradients.B[k][k2];
                        b_grad += tmp_grad * tmp_grad;
                    }
                    double tmp_grad = Get_Gas_BField(p,k);
                    b_mag += tmp_grad * tmp_grad;
                }
                double L_cond_inv = sqrt(b_grad / (1.e-37 + b_mag));
                double L_cond = DMAX(L_particle , 1./(L_cond_inv + 1./L_particle)) * All.cf_atime;
                double diff_coeff = SphP[p].Eta_MHD_OhmicResistivity_Coeff + SphP[p].Eta_MHD_HallEffect_Coeff + SphP[p].Eta_MHD_AmbiPolarDiffusion_Coeff;
                double dt_conduction = dt_prefac_diffusion * L_cond*L_cond / (1.0e-37 + diff_coeff);
#ifdef SUPER_TIMESTEP_DIFFUSION
                if(dt_conduction < dt_superstep_explicit) dt_superstep_explicit = dt_conduction; // explicit time-step
                double dt_advective = dt_conduction * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_cond_inv))*All.cf_atime / L_cond);
                if(dt_advective < dt) dt = dt_advective; // 'advective' timestep: needed to limit super-stepping
#else
                if(dt_conduction < dt) dt = dt_conduction; // normal explicit time-step
#endif
            }
#endif
            
            

            
#if defined(RADTRANSFER)
            {
                int kf;
                double dt_rad = 1.e10 * dt;
                for(kf=0;kf<N_RT_FREQ_BINS;kf++)
                {
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_EVOLVE_FLUX) /* for explicit diffusion, we include the usual second-order diffusion timestep */
                    double gradETmag=0; for(k=0;k<3;k++) {gradETmag += SphP[p].Gradients.Rad_E_gamma_ET[kf][k]*SphP[p].Gradients.Rad_E_gamma_ET[kf][k];}
                    double L_ETgrad_inv = sqrt(gradETmag) / (1.e-37 + SphP[p].Rad_E_gamma[kf] * SphP[p].Density/P[p].Mass);
                    double L_RT_diffusion = DMIN(L_particle , 1./(3.*L_ETgrad_inv)) * All.cf_atime;
                    double dt_rt_diffusion = dt_prefac_diffusion * L_RT_diffusion*L_RT_diffusion / (MIN_REAL_NUMBER + rt_diffusion_coefficient(p,kf));
                    double dt_advective = dt_rt_diffusion * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_ETgrad_inv))*All.cf_atime / L_RT_diffusion);
                    double dt_rt_work = All.CourantFac * DMIN( L_RT_diffusion / csnd , L_particle*All.cf_atime / ((2./3.)*sqrt(SphP[p].Rad_E_gamma[kf]/P[p].Mass)) ); /* time-step related to radiation work, radiation soundspeed, relevant in strongly-coupled limit */
#ifdef RT_FLUXLIMITER
                    if(dt_advective > dt_rt_diffusion) {dt_rt_diffusion *= 1. + (1.-SphP[p].Rad_Flux_Limiter[kf]) * DMAX(0,(dt_advective/dt_rt_diffusion-1.));}
                    dt_advective = All.CourantFac * 0.5 * (L_particle*All.cf_atime) / C_LIGHT_CODE_REDUCED;
                    dt_rt_diffusion = DMAX(dt_rt_diffusion, dt_advective);
                    dt_rt_work /= MIN_REAL_NUMBER + SphP[p].Rad_Flux_Limiter[kf];
                    if((SphP[p].Rad_Flux_Limiter[kf] <= 0)||(dt_rt_diffusion<=0)) {dt_rt_diffusion = 1.e9 * dt;}
#endif
                    if((SphP[p].Rad_E_gamma[kf] <= MIN_REAL_NUMBER) || (SphP[p].Rad_E_gamma_Pred[kf] <= MIN_REAL_NUMBER)) {dt_rt_diffusion = dt_advective;}
#ifdef SUPER_TIMESTEP_DIFFUSION
                    if(dt_rt_diffusion < dt_superstep_explicit) dt_superstep_explicit = dt_rt_diffusion; // explicit time-step
                    dt_advective = dt_rt_diffusion * DMAX(1,DMAX(L_particle , 1/(MIN_REAL_NUMBER + L_ETgrad_inv))*All.cf_atime / L_RT_diffusion);
                    if(dt_advective < dt_rad) dt_rad = dt_advective; // 'advective' timestep: needed to limit super-stepping
#else
                    if(dt_rt_diffusion < dt_rad) dt_rad = dt_rt_diffusion; // normal explicit time-step
                    if(dt_rt_work < dt_rad) {dt_rad = dt_rt_work;} // normal explicit time-step
#endif
#endif
#if defined(RT_RAD_PRESSURE_FORCES) && defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_EVOLVE_FLUX) // here the acceleration isn't saved to Rad_Accel so we calculate that timestep constraint
                    double gradErad=0; for(k=0;k<3;k++) {gradErad+=SphP[p].Gradients.Rad_E_gamma_ET[kf][k]*SphP[p].Gradients.Rad_E_gamma_ET[kf][k];}
                    double radacc = return_flux_limiter(p,kf) * (sqrt(gradErad) / SphP[p].Density) / All.cf_atime; // radiation acceleration for a timestep criterion
                    if(gradErad > 0 && radacc > 0)
                    {
                        double dt_radacc = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * DMIN(All.ForceSoftening[0], PPP[p].Hsml) / radacc);
                        if(dt_radacc < dt_rad) {dt_rad = dt_radacc;}
                    }
#endif
                }
                /* even with a fully-implicit solver, we require a CFL-like criterion on timesteps (much larger steps allowed for stability, but not accuracy) */
                dt_courant = All.CourantFac * (L_particle*All.cf_atime) / C_LIGHT_CODE_REDUCED; /* courant-type criterion, using the reduced speed of light */
#ifdef RT_M1
#ifndef GALSF
                dt_rad = dt_courant;
#endif
                double L_RT_diffusion = L_particle*All.cf_atime;
                for(kf=0;kf<N_RT_FREQ_BINS;kf++)
                {
                    double dt_rt_diffusion = dt_prefac_diffusion * L_RT_diffusion*L_RT_diffusion / (MIN_REAL_NUMBER + rt_diffusion_coefficient(p,kf));
#ifdef GALSF
                    /* ignore particles where the radiation energy density is basically non-existant */
                    if((SphP[p].Rad_E_gamma[kf] <= MIN_REAL_NUMBER) ||
                       (SphP[p].Rad_E_gamma_Pred[kf] <= MIN_REAL_NUMBER) ||
                       (SphP[p].Rad_E_gamma[kf] < 1.e-5*P[p].Mass*SphP[p].InternalEnergy)) {dt_rt_diffusion = 1.e10 * dt;}
#endif
                    if(dt_rt_diffusion < dt_rad) dt_rad = dt_rt_diffusion;
                }
                if(dt_rad > 1.e3*dt_courant) {dt_rad = 1.e3*dt_courant;}
                if(dt_courant > dt_rad) {dt_rad = dt_courant;}
#if defined(RT_CHEM_PHOTOION)
                /* since we're allowing rather large timesteps above in special conditions, make sure this doesn't overshoot the recombination time for the opacity to 
                    change, which can happen particularly for ionizing photons */
                if(kf==RT_FREQ_BIN_H0)
                {
                    double ne_cgs = (SphP[p].Density * All.cf_a3inv * UNIT_DENSITY_IN_NHCGS);
                    double dt_recombination = All.CourantFac * (3.3e12/ne_cgs) / UNIT_TIME_IN_CGS;
                    double dt_change = 1.e10*dt; if((SphP[p].Rad_E_gamma[kf] > 0)&&(fabs(SphP[p].Dt_Rad_E_gamma[kf])>0)) {dt_change = SphP[p].Rad_E_gamma[kf] / fabs(SphP[p].Dt_Rad_E_gamma[kf]);}
                    dt_recombination = DMIN(DMAX(dt_recombination,dt_change), DMAX(dt_courant,dt_rad));
                    if(dt_recombination < dt_rad) {dt_rad = dt_recombination;}
                }
#endif
#else
                if(dt_courant < dt_rad) {dt_rad = dt_courant;}
#endif
                if(dt_rad < dt) dt = dt_rad;
            }
#endif
    
    
#ifdef VISCOSITY
            {
                int kv1,kv2; double dv_mag=0,v_mag=1.0e-33;
                for(kv1=0;kv1<3;kv1++) {v_mag+=P[p].Vel[kv1]*P[p].Vel[kv1];}
                double dv_mag_all = 0.0;
                for(kv1=0;kv1<3;kv1++)
                {
                    double dvmag_tmp = 0;
                    for(kv2=0;kv2<3;kv2++) {dvmag_tmp+=SphP[p].Gradients.Velocity[kv1][kv2]*SphP[p].Gradients.Velocity[kv1][kv2];}
                    dv_mag += dvmag_tmp /DMAX(P[p].Vel[kv1]*P[p].Vel[kv1],0.01*v_mag);
                    dv_mag_all += dvmag_tmp;
                }
                dv_mag = sqrt(DMAX(dv_mag, dv_mag_all/v_mag));
                double L_visc = DMAX(L_particle , 1. / (dv_mag + 1./L_particle)) * All.cf_atime;
                double visc_coeff = sqrt(SphP[p].Eta_ShearViscosity*SphP[p].Eta_ShearViscosity + SphP[p].Zeta_BulkViscosity*SphP[p].Zeta_BulkViscosity);
                double dt_viscosity = 0.25 * L_visc*L_visc / (1.0e-33 + visc_coeff) * SphP[p].Density * All.cf_a3inv;
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
            if(P[p].Grain_AccelTimeMin < dt) {dt = P[p].Grain_AccelTimeMin;}
#endif


#ifdef TURB_DIFFUSION
            {
#ifdef TURB_DIFF_METALS
                int k_species; double L_tdiff = L_particle * All.cf_atime; // don't use gradient b/c ill-defined pre-enrichment
                for(k_species=0;k_species<NUM_METAL_SPECIES;k_species++)
                {
                    double dt_tdiff = L_tdiff*L_tdiff / (1.0e-33 + SphP[p].TD_DiffCoeff); // here, we use DIFFUSIVITIES, so there is no extra density power in the equation //
                    if(dt_tdiff < dt) dt = dt_tdiff; // normal explicit time-step
                }
#endif
            }
#endif
            
            
#if defined(DIVBCLEANING_DEDNER) 
            double fac_magnetic_pressure = All.cf_afac1 / All.cf_atime;
            double phi_b_units = Get_Gas_PhiField(p) / (All.cf_afac3 * All.cf_atime * SphP[p].MaxSignalVel);
            double vsig1 = All.cf_afac3 * sqrt( Get_Gas_effective_soundspeed_i(p)*Get_Gas_effective_soundspeed_i(p) +
                    fac_magnetic_pressure * (Get_Gas_BField(p,0)*Get_Gas_BField(p,0) +
                                             Get_Gas_BField(p,1)*Get_Gas_BField(p,1)+
                                             Get_Gas_BField(p,2)*Get_Gas_BField(p,2) +
                                             phi_b_units*phi_b_units) / SphP[p].Density );

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
	    
	    
            
            
#ifdef SUPER_TIMESTEP_DIFFUSION
            /* now use the timestep information above to limit the super-stepping timestep */
            {
                int N_substeps = 5; /*!< number of sub-steps per super-timestep for super-timestepping algorithm */
                double nu_substeps = 0.04; /*!< damping parameter (0<nu<1), optimal behavior around ~1/sqrt[N_substeps] */
                
                /*!< pre-calculate the multipliers needed for the super-timestep sub-step */
                if(SphP[p].Super_Timestep_j == 0) {SphP[p].Super_Timestep_Dt_Explicit = dt_superstep_explicit;} // reset dt_explicit //
                double j_p_super = (double)(SphP[p].Super_Timestep_j + 1);
                double dt_superstep = SphP[p].Super_Timestep_Dt_Explicit / ((nu_substeps+1) + (nu_substeps-1) * cos(M_PI * (2*j_p_super - 1) / (2*(double)N_substeps)));

                double dt_touse = dt_superstep;
                if((dt <= dt_superstep)||(SphP[p].Super_Timestep_j > 0))
                {
                    /* if(dt <= dt_superstep): other constraints beat our super-step, so it doesn't matter: iterate */
                    /* if(SphP[p].Super_Timestep_j > 0): don't break mid-cycle, so iterate */
                    SphP[p].Super_Timestep_j++; if(SphP[p].Super_Timestep_j>=N_substeps) {SphP[p].Super_Timestep_j=0;} /*!< increment substep 'j' and loop if it cycles fully */
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
                    double dt_allowed = (bin ? (((integertime) 1) << bin) : 0) * All.Timebase_interval / All.cf_hubble_a;
                    if(dt_superstep > 1.5*dt_allowed)
                    {
                        /* the next allowed timestep [because of synchronization] is not big enough to fit the 'big step' 
                            part of the super-stepping cycle. rather than 'waste' our timestep which will knock us into a 
                            lower bin and defeat the super-stepping, we simply take the -safe- explicit timestep and 
                            wait until the desired time-bin synchs up, so we can super-step */
                        dt_touse = dt_superstep_explicit; // use the safe [normal explicit] timestep and -do not- cycle j //
                    } else {
                        /* ok, we can jump up in bins to use our super-step; begin the cycle! */
                        SphP[p].Super_Timestep_j++; if(SphP[p].Super_Timestep_j>=N_substeps) {SphP[p].Super_Timestep_j=0;}
                    }
                }
                if(dt < dt_touse) {dt = dt_touse;} // set the actual timestep [now that we've appropriately checked everything above] //
            }
#endif
            
        } // closes if(P[p].Type == 0) [gas particle check] //
    
    
#if defined(DM_SIDM)
    /* Reduce time-step if this particle got interaction probabilities > 0.2 during the last time-step */
    if((1 << P[p].Type) & (DM_SIDM))
    {
        if(P[p].dt_step_sidm > 0)
        {
            double dt_sidm_physical = P[p].dt_step_sidm * All.Timebase_interval / All.cf_hubble_a;
            if(dt_sidm_physical < dt) {dt = dt_sidm_physical;}
        }
        if(dt > 0)
        {
            double p_target = 0.2; // desired maximum probability per timestep
            double dV[3]; for(k=0;k<3;k++) {dV[k]=P[p].AGS_vsig*All.cf_afac3*All.cf_atime/sqrt(3.);} // convert signal vel to velocity dispersion for estimating rates
#ifdef GRAIN_COLLISIONS
            double p_dt = prob_of_grain_interaction(return_grain_cross_section_per_unit_mass(p),P[p].Mass,0.,PPP[p].AGS_Hsml,dV,dt,p); // probability of interacting with another grain super-particle well within kernel, assuming same mass, H, and V~signalvel, for current timestep dt
#else
            double p_dt = prob_of_interaction(P[p].Mass,0.,PPP[p].AGS_Hsml,dV,dt); // probability of interacting with another DM particle well within kernel, assuming same mass, H, and V~signalvel, for current timestep dt
#endif
            if(p_dt > p_target) {dt = p_target;}
        }
    }
#endif
    
    
    // add a 'stellar evolution timescale' criterion to the timestep, to prevent too-large jumps in feedback //
#if defined(YOUNGSTARWINDDRIVING) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF_FB_MECHANICAL) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    if(((P[p].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[p].Type == 2)||(P[p].Type==3))))&&(P[p].Mass>0))
    {
        double star_age = evaluate_stellar_age_Gyr(P[p].StellarAge);
        double dt_stellar_evol;
        if(star_age<0.1)
        {
            dt_stellar_evol = DMAX(2.0e-4, star_age/250.); // restrict to small steps for young stars //
        } else {
            dt_stellar_evol = star_age/10.;
        }
        // PFH: temporarily modifying the terms above while Marcel studies them: turns out not to be necessary to use as strict a mass-dependent timestep, so faster to comment out //
        double mcorr = 1.e-5 * (P[p].Mass*UNIT_MASS_IN_SOLAR);
        if(mcorr < 1 && mcorr > 0) {dt_stellar_evol /= mcorr;}
        if(dt_stellar_evol < 1.e-6) {dt_stellar_evol = 1.e-6;}
        dt_stellar_evol /= (UNIT_TIME_IN_GYR); // convert to code units //
        if(dt_stellar_evol>0) {if(dt_stellar_evol<dt) {dt = dt_stellar_evol;}}
    }
#endif
    
    
#ifdef BLACK_HOLES

#ifdef BH_WAKEUP_GAS
    if(P[p].Type == 0)
    {
    	double dt_bh = 2.*(P[p].LowestBHTimeBin ? ((integertime)1 <<  P[p].LowestBHTimeBin) : 0) * All.Timebase_interval/All.cf_hubble_a;
    	if(dt>dt_bh) {dt=1.01*dt_bh; P[p].LowestBHTimeBin=TIMEBINS;}
    }
#endif	

    if(P[p].Type == 5)
    {
#if !defined(SINGLE_STAR_SINK_DYNAMICS) && defined(GALSF)
      double dt_accr = 4.2e5 / UNIT_TIME_IN_YR; // this is the 1% of Salpeter timescale; not relevant for low radiative efficiency
#else
      double dt_accr = All.MaxSizeTimestep;
#endif      
        if(BPP(p).BH_Mdot > 0 && BPP(p).BH_Mass > 0)
        {
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
            /* really want prefactor to be ratio of median gas mass to bh mass */
            dt_accr = 0.001 * DMAX(BPP(p).BH_Mass, All.MaxMassForParticleSplit) / BPP(p).BH_Mdot;
#if defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
            dt_accr *= DMAX(0.1, All.BAL_f_accretion);
#endif
#else
            dt_accr = 0.05 * DMAX(BPP(p).BH_Mass , All.MaxMassForParticleSplit) / BPP(p).BH_Mdot;
#endif 	    
#ifdef SINGLE_STAR_SINK_DYNAMICS
            dt_accr = All.MinMassForParticleMerger / BPP(p).BH_Mdot;
#endif
        } // if(BPP(p).BH_Mdot > 0 && BPP(p).BH_Mass > 0)
#ifdef BH_SEED_GROWTH_TESTS
            double dt_evol = 1.e4 / UNIT_TIME_IN_YR;
#ifdef TURB_DRIVING
            if(dt_evol > 1.e-3*All.StDecay) {dt_evol=1.e-3*All.StDecay;}
#endif
            if(dt_accr > dt_evol) {dt_accr=dt_evol;}
#endif
            if(dt_accr > 0 && dt_accr < dt) {dt = dt_accr;}

        double dt_ngbs = (BPP(p).BH_TimeBinGasNeighbor ? (((integertime) 1) << BPP(p).BH_TimeBinGasNeighbor) : 0) * All.Timebase_interval / All.cf_hubble_a;
#ifndef SINGLE_STAR_SINK_DYNAMICS
        dt_ngbs *= 4.; /* standard wakeup-type threshold: use this by default here, unless dynamical interaction important (e.g. back-rx term from oscillation of BH c-o-m, which is important for single-sink sims */
#endif
        if(dt > dt_ngbs && dt_ngbs > 0) {dt = 1.01 * dt_ngbs; }

#if defined(SINGLE_STAR_TIMESTEPPING)
	    if(P[p].DensAroundStar > 0)
	    {
            double eps = DMAX( KERNEL_CORE_SIZE*All.ForceSoftening[5], BPP(p).BH_dr_to_NearestGasNeighbor);
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
            eps = DMAX(eps, BPP(p).SinkRadius);
#endif
            if(eps < MAX_REAL_NUMBER) {eps = DMAX(Get_Particle_Size(p), eps);} else {eps = Get_Particle_Size(p);}
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
            eps = DMAX(eps, KERNEL_CORE_SIZE*P[p].AGS_Hsml);
#endif		
            double dt_ff = sqrt(2*All.ErrTolIntAccuracy * pow(eps*All.cf_atime,3) / (All.G * P[p].Mass)); // fraction of the freefall time of the nearest gas particle from rest
            if(dt > dt_ff && dt_ff > 0) {dt = 1.01 * dt_ff;}
		
            double L_particle = Get_Particle_Size(p);           
            double dt_cour_sink = 0.5 * All.CourantFac * (L_particle*All.cf_atime) / P[p].BH_SurroundingGasVel;
            if(dt > dt_cour_sink && dt_cour_sink > 0) {dt = 1.01 * dt_cour_sink;}
        }
        if(P[p].StellarAge == All.Time)
        {   // want a brand new sink to be on the lowest occupied timebin
            long bin; for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinCount[bin] > 0) break;}
            double dt_min =  ((bin ? (((integertime) 1) << bin) : 0) * All.Timebase_interval / All.cf_hubble_a);
            if(dt > dt_min && dt_min > 0) dt = 1.01 * dt_min;
        }
#endif
    } // if(P[p].Type == 5)

#endif // BLACK_HOLES
    

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
        double agrav = sqrt(P[p].GravAccel[0]*P[p].GravAccel[0] + P[p].GravAccel[1]*P[p].GravAccel[1] + P[p].GravAccel[2]*P[p].GravAccel[2]);
        if(P[p].Type == 0)
        {
            double ahydro = sqrt(SphP[p].HydroAccel[0]*SphP[p].HydroAccel[0] + SphP[p].HydroAccel[1]*SphP[p].HydroAccel[1] + SphP[p].HydroAccel[2]*SphP[p].HydroAccel[2]);
#ifndef LONGIDS
            PRINT_WARNING("Part-ID=%d  dt=%g dtc=%g ac=%g agrav=%g ahydro=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g m=%g type=%d\n",
                          (int) P[p].ID, dt, dt_courant * All.cf_hubble_a, ac, agrav, ahydro, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], PPP[p].Hsml, csnd, sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac) * All.cf_hubble_a, All.SofteningTable[P[p].Type], P[p].Mass,P[p].Type);
#else
            PRINT_WARNING("Part-ID=%llu  dt=%g dtc=%g ac=%g agrav=%g ahydro=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g m=%g type=%d\n",
                          (MyIDType) P[p].ID, dt, dt_courant * All.cf_hubble_a, ac, agrav, ahydro, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], PPP[p].Hsml, csnd, sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac) * All.cf_hubble_a, All.SofteningTable[P[p].Type], P[p].Mass,P[p].Type);
#endif // ndef LONGIDS
        }
        else // if(P[p].Type == 0)
        {
#ifndef LONGIDS
            PRINT_WARNING("Part-ID=%d  dt=%g ac=%g agrav=%g xyz=(%g|%g|%g) type=%d\n", (int) P[p].ID, dt, ac, agrav, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],P[p].Type);
#else
            PRINT_WARNING("Part-ID=%llu  dt=%g ac=%g agrav=%g xyz=(%g|%g|%g) type=%d\n", (MyIDType) P[p].ID, dt, ac, agrav, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],P[p].Type);
#endif // ndef LONGIDS
        }
        fflush(stdout); fprintf(stderr, "\n @ fflush \n");
#ifdef STOP_WHEN_BELOW_MINTIMESTEP
        endrun(888);
#endif
        dt = All.MinSizeTimestep;
    }
    
    ti_step = (integertime) (dt / All.Timebase_interval);
#ifndef STOP_WHEN_BELOW_MINTIMESTEP
    if(ti_step<=1) ti_step=2;
#endif
    
    if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
        printf("\nError: A timestep of size zero was assigned on the integer timeline, no here!!!\n"
               "We better stop.\n"
               "Task=%d Part-ID=%llu dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%lld ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
               ThisTask, (unsigned long long) P[p].ID, dt, dt_courant, dt_divv, dt_displacement,
               All.Timebase_interval, (long long) ti_step, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].GravAccel[0], P[p].GravAccel[1], P[p].GravAccel[2]);
#ifdef PMGRID
        printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
        fflush(stdout); endrun(818);
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
                count[P[i].Type]++;
                v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
                if(mim[P[i].Type] > P[i].Mass) {mim[P[i].Type] = P[i].Mass;}
                mnm[P[i].Type] += P[i].Mass;
            }
        }
        
        MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(mnm, mean_mass, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumup_large_ints(6, count, count_sum);
        
#ifdef GALSF
        /* add star and gas particles together to treat them on equal footing, using the original gas particle spacing. */
        v_sum[0] += v_sum[4];
        count_sum[0] += count_sum[4];
        v_sum[4] = v_sum[0];
        count_sum[4] = count_sum[0];
        min_mass[0] = min_mass[4] = (mean_mass[0] + mean_mass[4]) / count_sum[0];
#ifdef BLACK_HOLES
        v_sum[0] += v_sum[5];
        count_sum[0] += count_sum[5];
        v_sum[5] = v_sum[0];
        count_sum[5] = count_sum[0];
        min_mass[5] = min_mass[0];
#endif
#endif
        
        if(ThisTask == 0)
            printf("Global displacement time constraint computation: \n");
        for(type = 0; type < 6; type++)
        {
            if(count_sum[type] > 0)
            {
#ifdef GALSF
                if(type == 0 || type == 4)
#else
                if(type == 0)
#endif
                    dmean = pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);
                else
                    dmean = pow(min_mass[type] / ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);
                
#ifdef BLACK_HOLES
                if(type == 5) {dmean = pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);}
#endif
                dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);
                
#ifdef PMGRID
                asmth = All.Asmth[0];
#ifdef PM_PLACEHIGHRESREGION
                if(((1 << type) & (PM_PLACEHIGHRESREGION)))
                    asmth = All.Asmth[1];
#endif
                if(asmth < dmean)
                    dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif
                
                if(ThisTask == 0)
                    printf(" ..type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
                           type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);
                
                if(dt < dt_displacement)
                    dt_displacement = dt;
            }
        }
        
        if(ThisTask == 0)
            printf(" ..global displacement time constraint: %g  (All.MaxSizeTimestep=%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}



int get_timestep_bin(integertime ti_step)
{
    int bin = -1;
    
    if(ti_step == 0)
        return 0;
    
    if(ti_step == 1)
        terminate("time-step of integer size 1 not allowed\n");
    
    while(ti_step)
    {
        bin++;
        ti_step >>= 1;
    }
    
    return bin;
}





#ifdef WAKEUP
void process_wake_ups(void)
{   
    int i, n;
    integertime dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;
    int max_time_bin_active;
    int bin, binold, prev, next;
    long long ntot;
    
    /* find the next kick time */
    for(n = 0, ti_next_kick = TIMEBASE; n < TIMEBINS; n++)
    {
        if(TimeBinCount[n])
        {
            if(n > 0)
            {
                dt_bin = (((integertime) 1) << n);
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
    
    PRINT_STATUS("Predicting next timestep: %g", (ti_next_kick_global - All.Ti_Current) * All.Timebase_interval);
    max_time_bin_active = 0;
    /* get the highest bin, that is active next time */
    for(n = 0; n < TIMEBINS; n++)
    {
        dt_bin = (((integertime) 1) << n);
        
        if((ti_next_kick_global % dt_bin) == 0)
            max_time_bin_active = n;
    }
    
    /* move the particle on the highest bin, that is active in the next timestep and that is lower than its last timebin */
    bin = 0;
    for(n = 0; n < TIMEBINS; n++)
    {
        if(TimeBinCount[n] > 0)
        {
            bin = n;
            break;
        }
    }
    n = 0;

    MPI_Allreduce(&NeedToWakeupParticles_local, &NeedToWakeupParticles, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); // if one process processes wakeups then they all should, just in case a woke particle gets swapped to another process before we get here
    
    if(NeedToWakeupParticles){
	for(i = 0; i < NumPart; i++)
	{
	    if(!PPPZ[i].wakeup)
		continue;
	
#if !defined(AGS_HSML_CALCULATION_IS_ACTIVE)
	    if(P[i].Type != 0) {continue;} // only gas particles can be awakened
#endif
        
	    if(P[i].Mass <= 0)
		continue;       
        
	    binold = P[i].TimeBin;
	    if(TimeBinActive[binold])
		continue;
        
	    bin = max_time_bin_active < binold ? max_time_bin_active : binold;
        
	    if(bin != binold)
	    {
		integertime dt_0 = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
		integertime tstart = P[i].Ti_begstep + dt_0;
		integertime t_2 = P[i].Ti_current;
		if(t_2 > tstart) {tstart = t_2;}
		integertime tend = All.Ti_Current;

		TimeBinCount[binold]--;
		if(P[i].Type == 0)
		    TimeBinCountSph[binold]--;
            
		prev = PrevInTimeBin[i];
		next = NextInTimeBin[i];
            
		if(FirstInTimeBin[binold] == i)
		    FirstInTimeBin[binold] = next;
		if(LastInTimeBin[binold] == i)
		    LastInTimeBin[binold] = prev;
		if(prev >= 0)
		    NextInTimeBin[prev] = next;
		if(next >= 0)
		    PrevInTimeBin[next] = prev;
            
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
		if(P[i].Type == 0)
		    TimeBinCountSph[bin]++;
            
		P[i].TimeBin = bin;
            
		if(TimeBinActive[bin])
		    NumForceUpdate++;
                        
		n++;

		/* reverse part of the last second-half kick this particle received 
		   (to correct it back to its new active time) */
		if(tend < tstart)
		{
		    do_the_kick(i, tstart, tend, P[i].Ti_current, 1);
		    set_predicted_sph_quantities_for_extra_physics(i);
		}
		P[i].Ti_begstep = All.Ti_Current;
		P[i].dt_step = bin ? (((integertime) 1) << bin) : 0;
		if(P[i].Ti_current < All.Ti_Current) {P[i].Ti_current=All.Ti_Current;}
            
	    }
	}
    }
    
    sumup_large_ints(1, &n, &ntot);
    if(ThisTask == 0)
    {
        if(ntot > 0) {printf("%d%09d particles woken up.\n", (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
    }
    NeedToWakeupParticles = 0;
    NeedToWakeupParticles_local = 0;
}
#endif
