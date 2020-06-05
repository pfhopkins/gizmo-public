#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../../allvars.h"
#include "../../proto.h"

/*! \file blackhole_util.c
 *  \brief util routines for memory (de)allocation and array setting for black holes
 */
/*
* This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
* see notes in blackhole.c for details on code history.
*/

#ifdef BLACK_HOLES // master flag [needs to be here to prevent compiler breaking when this is not active] //

/* function for allocating temp BH data struc needed for feedback routines*/
void blackhole_start(void)
{
    int i, Nbh;
    
    /* count the num BHs on this task */
    N_active_loc_BHs=0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type ==5)
        {
            P[i].IndexMapToTempStruc = N_active_loc_BHs;         /* allows access via BlackholeTempInfo[P[i].IndexMapToTempStruc] */
            N_active_loc_BHs++;                     /* N_active_loc_BHs now set for BH routines */
        }
    }
    
    /* allocate the blackhole temp struct */
    if(N_active_loc_BHs>0)
    {
        BlackholeTempInfo = (struct blackhole_temp_particle_data *) mymalloc("BlackholeTempInfo", N_active_loc_BHs * sizeof(struct blackhole_temp_particle_data));
    } else {
        BlackholeTempInfo = (struct blackhole_temp_particle_data *) mymalloc("BlackholeTempInfo", 1 * sizeof(struct blackhole_temp_particle_data));
    }
    
    memset( &BlackholeTempInfo[0], 0, N_active_loc_BHs * sizeof(struct blackhole_temp_particle_data) );
    
    Nbh=0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type ==5)
        {
            BlackholeTempInfo[Nbh].index = i;               /* only meaningful field set here */
            Nbh++;
        }
    }
    
    /* all future loops can now take the following form:
     for(i=0; i<N_active_loc_BHs; i++)
     {
     i_old = BlackholeTempInfo[i].index;
     ...
     }
     */
    
}


/* function for freeing temp BH data struc needed for feedback routines*/
void blackhole_end(void)
{
#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    {
        int bin;
        double mdot, mdot_in_msun_per_year;
        double mass_real, total_mass_real, medd, total_mdoteddington;
        double mass_holes, total_mass_holes, total_mdot;
        
        /* sum up numbers to print for summary of the BH step (blackholes.txt) */
        mdot = mass_holes = mass_real = medd = 0;
        for(bin = 0; bin < TIMEBINS; bin++)
        {
            if(TimeBinCount[bin])
            {
                mass_holes += TimeBin_BH_mass[bin];
                mass_real += TimeBin_BH_dynamicalmass[bin];
                mdot += TimeBin_BH_Mdot[bin];
                medd += TimeBin_BH_Medd[bin];
            }
        }
        MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&medd, &total_mdoteddington, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if((ThisTask == 0) && (total_mdot > 0) && (total_mass_real > 0))
        {
            /* convert to solar masses per yr */
            mdot_in_msun_per_year = total_mdot * UNIT_MASS_IN_SOLAR/UNIT_TIME_IN_YR;
            total_mdoteddington *= 1.0 / bh_eddington_mdot(1);
            fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n", All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year, total_mass_real, total_mdoteddington);
        }
        fflush(FdBlackHoles);
#ifdef BH_OUTPUT_GASSWALLOW
        fflush(FdBhSwallowDetails);
#endif
#ifdef BH_OUTPUT_FORMATION_PROPERTIES
        fflush(FdBhFormationDetails);
#endif
#if !defined(IO_REDUCED_MODE) || defined(BH_OUTPUT_MOREINFO)
        fflush(FdBlackHolesDetails);
#ifdef BH_OUTPUT_MOREINFO
        fflush(FdBhMergerDetails);
#ifdef BH_WIND_KICK
        fflush(FdBhWindDetails);
#endif
#endif
#endif
    }
    myfree(BlackholeTempInfo);
}



/* return the eddington accretion-rate = L_edd/(epsilon_r*c*c) */
double bh_eddington_mdot(double bh_mass)
{
    return (4*M_PI * GRAVITY_G * PROTONMASS / (All.BlackHoleRadiativeEfficiency * C_LIGHT * THOMPSON)) * bh_mass * UNIT_TIME_IN_CGS;
}



/* return the bh luminosity given some accretion rate and mass (allows for non-standard models: radiatively inefficient flows, stellar sinks, etc) */
double bh_lum_bol(double mdot, double mass, long id)
{
    double lum = All.BlackHoleRadiativeEfficiency * mdot * C_LIGHT_CODE*C_LIGHT_CODE; // this is automatically in -physical code units-
#ifdef SINGLE_STAR_SINK_DYNAMICS
    lum = calculate_individual_stellar_luminosity(mdot,mass,id);
#endif
    return All.BlackHoleFeedbackFactor * lum;
}



void blackhole_properties_loop(void) /* Note, normalize_temp_info_struct is now done at the end of blackhole_environment_loop(), so that final quantities are available for the second environment loop if needed */
{
    int i, n; double dt;
    for(i=0; i<N_active_loc_BHs; i++)
    {
        n = BlackholeTempInfo[i].index;
#ifndef WAKEUP /* define the timestep */
        dt = (P[n].TimeBin ? (((integertime) 1) << P[n].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[n].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
        BPP(n).BH_Mdot=0;  /* always initialize/default to zero accretion rate */
        set_blackhole_long_range_rp(i, n);
        set_blackhole_mdot(i, n, dt);
#if defined(BH_DRAG) || defined(BH_DYNFRICTION)
        set_blackhole_drag(i, n, dt);
#endif
        set_blackhole_new_mass(i, n, dt);
        /* results dumped to 'blackhole_details' files at the end of blackhole_final_operations so that BH mass is corrected for mass loss to radiation/bal outflows */
    }// for(i=0; i<N_active_loc_BHs; i++)
}



#endif // master flag
