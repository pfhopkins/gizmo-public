#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "allvars.h"
#include "proto.h"
#include "kernel.h"

/*! \file io.c
 *  \brief Output of a snapshot file to disk.
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly to
 * write out new/modified quantities, as needed)
 */

static int n_type[6];
static long long ntot_type_all[6];

static int n_info;

/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */
void savepositions(int num)
{
    size_t bytes;
    char buf[500];
    int n, filenr, gr, ngroups, primaryTask, lastTask;

    CPU_Step[CPU_MISC] += measure_time();

#ifdef CHIMES_REDUCED_OUTPUT
    if (num % N_chimes_full_output_freq == 0) {Chimes_incl_full_output = 1;} else {Chimes_incl_full_output = 0;}
#endif

    rearrange_particle_sequence();
    /* ensures that new tree will be constructed */
    All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TreeDomainUpdateFrequency * All.TotNumPart);


    if(DumpFlag == 1)
    {
        if(ThisTask == 0)
            printf("\nwriting snapshot file #%d... \n", num);

        size_t MyBufferSize = All.BufferSize;
        if(!(CommBuffer = mymalloc("CommBuffer", bytes = MyBufferSize * 1024 * 1024)))
        {
            printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
            endrun(2);
        }


        if(NTask < All.NumFilesPerSnapshot)
        {
            if(ThisTask == 0)
                printf
                ("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
            endrun(0);
        }
        if(All.SnapFormat < 1 || All.SnapFormat > 3)
        {
            if(ThisTask == 0)
                printf("Unsupported File-Format\n");
            endrun(0);
        }
#ifndef  HAVE_HDF5
        if(All.SnapFormat == 3)
        {
            if(ThisTask == 0)
                printf("Code wasn't compiled with HDF5 support enabled!\n");
            endrun(0);
        }
#endif


        /* determine global and local particle numbers */
        for(n = 0; n < 6; n++)
            n_type[n] = 0;

        for(n = 0; n < NumPart; n++)
            n_type[P[n].Type]++;

        sumup_large_ints(6, n_type, ntot_type_all);

        /* assign processors to output files */
        distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &primaryTask, &lastTask);

        if(All.NumFilesPerSnapshot > 1)
        {
            if(ThisTask == 0)
            {
                sprintf(buf, "%s/snapdir_%03d", All.OutputDir, num);
                mkdir(buf, 02755);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        if(All.NumFilesPerSnapshot > 1)
            sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
        else
            sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);


        ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
        if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel)) {ngroups++;}

        for(gr = 0; gr < ngroups; gr++)
        {
            if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
            {
                write_file(buf, primaryTask, lastTask);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        myfree(CommBuffer);

        if(ThisTask == 0)
            printf("done with snapshot.\n");

        All.Ti_lastoutput = All.Ti_Current;

        CPU_Step[CPU_SNAPSHOT] += measure_time();
}

#ifdef FOF
    if(RestartFlag != 4)
    {
        if(ThisTask == 0)
            printf("\ncomputing group catalogue...\n");

        fof_fof(num);

        if(ThisTask == 0)
            printf("done with group catalogue.\n");

        CPU_Step[CPU_FOF] += measure_time();
    }
#endif

#ifdef OUTPUT_POWERSPEC
    if(RestartFlag != 4)
    {
        if(ThisTask == 0)
            printf("\ncomputing power spectra...\n");

        calculate_power_spectra(num, &ntot_type_all[0]);

        if(ThisTask == 0)
            printf("done with power spectra.\n");

        CPU_Step[CPU_MISC] += measure_time();
    }
#endif

}



/*! This function fills the write buffer with particle data. New output blocks can in
 *  principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
    int n, k, pindex;
    MyOutputFloat *fp;
    MyOutputPosFloat *fp_pos;
    MyIDType *ip;
    int *ip_int;
    float *fp_single;
#ifdef OUTPUT_COOLRATE
    double tcool, u;
#endif
#if defined(OUTPUT_GDE_DISTORTIONTENSOR)
    MyBigFloat half_kick_add[6][6];
#endif
#ifdef MAGNETIC /* NOTE: we always work -internally- in code units where MU_0 = 1; hence the 4pi here; [much simpler, but be sure of your conversions!] */
    double gizmo2gauss = UNIT_B_IN_GAUSS / All.UnitMagneticField_in_gauss;
#endif
#ifdef GDE_DISTORTIONTENSOR
    MyBigFloat flde, psde;
#endif
#ifdef PMGRID
    double dt_gravkick_pm = 0;
    if(All.ComovingIntegrationOn)
        {dt_gravkick_pm = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) - get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);}
    else
        {dt_gravkick_pm = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;}
#endif

    fp = (MyOutputFloat *) CommBuffer;
    fp_single = (float *) CommBuffer;
    fp_pos = (MyOutputPosFloat *) CommBuffer;
    ip = (MyIDType *) CommBuffer;
    ip_int = (int *) CommBuffer;
    pindex = *startindex;

    switch (blocknr)
    {
        case IO_POS:		/* positions */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                    {
                        fp_pos[k] = (MyOutputPosFloat) P[pindex].Pos[k];
#ifdef BOX_PERIODIC
                        double box_length_xyz;
                        if(k==0) {box_length_xyz = boxSize_X;}
                        if(k==1) {box_length_xyz = boxSize_Y;}
                        if(k==2) {box_length_xyz = boxSize_Z;}
                        while(fp_pos[k] < 0) {fp_pos[k] += (MyOutputFloat) box_length_xyz;}
                        while(fp_pos[k] >= box_length_xyz) {fp_pos[k] -= (MyOutputFloat) box_length_xyz;}
#endif
                    }
                    fp_pos += 3;
                    n++;
                }
            break;

        case IO_VEL:		/* velocities [we're drifting here to the snapshot, note this is -not- the exact velocity in-code b/c we're alternating drifts and kicks!] */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#if 1
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) (P[pindex].Vel[k] * sqrt(All.cf_a3inv));} // JUST write the conserved velocity here, not the drifted one in this manner //
#else
                    double dt_gravkick, dt_hydrokick;
                    integertime dt_integerstep = GET_PARTICLE_INTEGERTIME(pindex);
                    dt_hydrokick = (All.Ti_Current - (P[pindex].Ti_begstep + dt_integerstep / 2)) * UNIT_INTEGERTIME_IN_PHYSICAL;
                    if(All.ComovingIntegrationOn) {dt_gravkick = get_gravkick_factor(P[pindex].Ti_begstep, All.Ti_Current) - get_gravkick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_integerstep / 2);} else {dt_gravkick = dt_hydrokick;}
                    for(k = 0; k < 3; k++)
                    {
                        fp[k] = (MyOutputFloat) (P[pindex].Vel[k] + P[pindex].GravAccel[k] * dt_gravkick);
#if (SINGLE_STAR_TIMESTEPPING > 0)
			            if((P[pindex].Type == 5) && (P[pindex].SuperTimestepFlag >= 2)) {fp[k] += (MyOutputFloat) ((P[pindex].COM_GravAccel[k]-P[pindex].GravAccel[k]) * dt_gravkick);}
#endif
                        if(P[pindex].Type == 0) {fp[k] += (MyOutputFloat) (SphP[pindex].HydroAccel[k] * dt_hydrokick * All.cf_atime);}
                    }
#ifdef PMGRID
                    for(k = 0; k < 3; k++) {fp[k] += (MyOutputFloat) (P[pindex].GravPM[k] * dt_gravkick_pm);}
#endif
                    for(k = 0; k < 3; k++) {fp[k] *= (MyOutputFloat) sqrt(All.cf_a3inv);}
#endif
                    fp += 3;
                    n++;
                }
            break;

        case IO_ID:		/* particle ID */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip++ = (MyIDType) P[pindex].ID;
                    n++;
                }
            break;

        case IO_CHILD_ID:		/* particle 'child' ID (for splits/mergers) */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip++ = (MyIDType) P[pindex].ID_child_number;
                    n++;
                }
            break;

        case IO_GENERATION_ID:	/* particle ID generation (for splits/mergers) */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip++ = (MyIDType) P[pindex].ID_generation;
                    n++;
                }
            break;

        case IO_MASS:		/* particle mass */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Mass;
                    n++;
                }
            break;

        case IO_U:			/* internal energy */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) DMAX(All.MinEgySpec, SphP[pindex].InternalEnergyPred);
                    n++;
                }
            break;

        case IO_RHO:		/* density */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].Density;
                    n++;
                }
            break;

        case IO_NE:		/* electron abundance */
#if (defined(COOLING) || defined(RT_CHEM_PHOTOION)) && !defined(CHIMES)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].Ne;
                    n++;
                }
#endif
            break;

        case IO_NH:		/* neutral hydrogen fraction */
#if (defined(COOLING) || defined(RT_CHEM_PHOTOION)) && !defined(CHIMES)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#if defined(RT_CHEM_PHOTOION)
                    *fp++ = (MyOutputFloat) SphP[pindex].HI;
#elif (COOL_GRACKLE_CHEMISTRY > 0)
                    *fp++ = (MyOutputFloat) SphP[pindex].grHI;
#else
                    double u, ne, nh0 = 0, mu = 1, temp, nHeII, nhp, nHe0, nHepp; u = DMAX(All.MinEgySpec, SphP[pindex].InternalEnergy); // needs to be in code units
                    temp = ThermalProperties(u, SphP[pindex].Density * All.cf_a3inv, pindex, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);
                    *fp++ = (MyOutputFloat) nh0;
#endif
                    n++;
                }
#endif
            break;

        case IO_HII:		/* ionized hydrogen abundance */
#if defined(RT_CHEM_PHOTOION)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].HII;
                    n++;
                }
#endif
            break;

        case IO_HeI:		/* neutral Helium */
#if defined(RT_CHEM_PHOTOION_HE)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].HeI;
                    n++;
                }
#endif
            break;

        case IO_HeII:		/* ionized Helium */
#if defined(RT_CHEM_PHOTOION_HE)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].HeII;
                    n++;
                }
#endif
            break;
            
        case IO_INIB:
#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) && defined(BH_DEBUG_SPAWN_JET_TEST)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {*fp++ = (MyOutputFloat) SphP[pindex].IniB[k];}
                    n++;
                }
#endif               
            break;
            
        case IO_IDEN:
#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) && defined(BH_DEBUG_SPAWN_JET_TEST)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].IniDen;
                    n++;
                }
#endif                
            break;
            
        case IO_UNSPMASS:
#if defined(BH_WIND_SPAWN) && defined(BH_DEBUG_SPAWN_JET_TEST)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].unspawned_wind_mass;
                    n++;
                }
#endif           
            break;
            
        case IO_CRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].CoolingRate;
                    n++;
                }
#endif
            break;

        case IO_HRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].HeatingRate;
                    n++;
                }
#endif
            break;

        case IO_NHRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].NetHeatingRateQ;
                    n++;
                }
#endif
            break;

        case IO_HHRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].HydroHeatingRate;
                    n++;
                }
#endif
            break;

        case IO_MCRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].MetalCoolingRate;
                    n++;
                }
#endif
            break;

        case IO_HSML:		/* gas kernel length */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) PPP[pindex].Hsml;
                    n++;
                }
            break;

        case IO_SFR:		/* star formation rate */
#ifdef GALSF
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {   /* units convert to solar masses per yr */
                    *fp++ = (MyOutputFloat) (get_starformation_rate(pindex, 1) * UNIT_MASS_IN_SOLAR / UNIT_TIME_IN_YR);
                    n++;
                }
#endif
            break;

        case IO_AGE:		/* stellar formation time */
#ifdef GALSF
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].StellarAge;
                    n++;
                }
#endif
            break;

        case IO_OSTAR:
#ifdef GALSF_SFR_IMF_SAMPLING
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].IMF_NumMassiveStars;
                    n++;
                  }
#endif
            break;

        case IO_GRAINSIZE:		/* grain size */
#ifdef GRAIN_FLUID
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Grain_Size;
                    n++;
                }
#endif
            break;

        case IO_GRAINTYPE:      /* grain type */
#if defined(PIC_MHD)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) P[pindex].MHD_PIC_SubType;
                    n++;
                }
#endif
            break;

        case IO_VSTURB_DISS:
#if defined(TURB_DRIVING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].DuDt_diss;
                    n++;
                }
#endif
            break;

        case IO_VSTURB_DRIVE:
#if defined(TURB_DRIVING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].DuDt_drive;
                    n++;
                }
#endif
            break;

        case IO_HSMS:		/* kernel length for star particles */
#ifdef SUBFIND
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].DM_Hsml;
                    n++;
                }
#endif
            break;

        case IO_Z:			/* gas and star metallicity */
#ifdef METALS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<NUM_METAL_SPECIES;k++) {fp[k] = (MyOutputFloat) P[pindex].Metallicity[k];}
                    fp += NUM_METAL_SPECIES;
                    n++;
                }
#endif
            break;

        case IO_CHIMES_ABUNDANCES:
#ifdef CHIMES
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    for (k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++) {fp[k] = (MyOutputFloat) ChimesGasVars[pindex].abundances[k];}
                    fp += ChimesGlobalVars.totalNumberOfSpecies;
                    n++;
                }
#endif
            break;


        case IO_CHIMES_MU:
#ifdef CHIMES
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) calculate_mean_molecular_weight(&(ChimesGasVars[pindex]), &ChimesGlobalVars);
                    n++;
                }
#endif
            break;

        case IO_CHIMES_REDUCED:
#ifdef CHIMES_REDUCED_OUTPUT
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    fp[0] = (MyOutputFloat) ChimesGasVars[pindex].abundances[ChimesGlobalVars.speciesIndices[sp_elec]];
                    fp[1] = (MyOutputFloat) ChimesGasVars[pindex].abundances[ChimesGlobalVars.speciesIndices[sp_HI]];
                    fp[2] = (MyOutputFloat) ChimesGasVars[pindex].abundances[ChimesGlobalVars.speciesIndices[sp_H2]];
                    fp[3] = (MyOutputFloat) ChimesGasVars[pindex].abundances[ChimesGlobalVars.speciesIndices[sp_CO]];
                    fp += 4;
                    n++;
                }
#endif
            break;

        case IO_CHIMES_NH:
#if defined(CHIMES_NH_OUTPUT)
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (evaluate_NH_from_GradRho(P[pindex].GradRho,PPP[pindex].Hsml,SphP[pindex].Density,PPP[pindex].NumNgb,1,pindex) * UNIT_SURFDEN_IN_CGS * shielding_length_factor * (1.0 - (P[pindex].Metallicity[0] + P[pindex].Metallicity[1])) / PROTONMASS_CGS);
                    n++;
                }
#endif
            break;

        case IO_CHIMES_STAR_SIGMA:
#if defined(CHIMES_NH_OUTPUT) && defined(OUTPUT_DENS_AROUND_STAR)
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (evaluate_NH_from_GradRho(P[pindex].GradRho,PPP[pindex].Hsml,P[pindex].DensAroundStar,PPP[pindex].NumNgb,0,pindex) * UNIT_SURFDEN_IN_CGS);  // g cm^-2
                    n++;
                }
#endif
            break;

        case IO_CHIMES_FLUX_G0:
#ifdef CHIMES_STELLAR_FLUXES
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
#ifdef CHIMES_HII_REGIONS
                    if(SphP[pindex].DelayTimeHII > 0) {for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) (SphP[pindex].Chimes_G0[k] + SphP[pindex].Chimes_G0_HII[k]);}}
                        else {for(k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) SphP[pindex].Chimes_G0[k];}}
#else
                    for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) SphP[pindex].Chimes_G0[k]; }
#endif
                    fp += CHIMES_LOCAL_UV_NBINS;
                    n++;
                }
#endif
            break;

        case IO_CHIMES_FLUX_ION:
#ifdef CHIMES_STELLAR_FLUXES
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
#ifdef CHIMES_HII_REGIONS
                    if(SphP[pindex].DelayTimeHII > 0) {for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) (SphP[pindex].Chimes_fluxPhotIon[k] + SphP[pindex].Chimes_fluxPhotIon_HII[k]);}}
                        else {for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) SphP[pindex].Chimes_fluxPhotIon[k];}}
#else
                    for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) SphP[pindex].Chimes_fluxPhotIon[k];}
#endif
                    fp += CHIMES_LOCAL_UV_NBINS;
                    n++;
                }
#endif
            break;

        case IO_DENS_AROUND_STAR:
#ifdef OUTPUT_DENS_AROUND_STAR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].DensAroundStar;
                    n++;
                }
#endif
            break;

        case IO_DELAY_TIME_HII:
            break;

        case IO_MOLECULARFRACTION:
#if defined(OUTPUT_MOLECULAR_FRACTION)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#if defined(COOL_MOLECFRAC_NONEQM)
                    *fp++ = (MyOutputFloat) SphP[pindex].MolecularMassFraction_perNeutralH; /* more useful to output this particular value, rather than fH2 */
#else
                    double u, ne, nh0 = 0, mu = 1, temp, nHeII, nhp, nHe0, nHepp; u = DMAX(All.MinEgySpec, SphP[pindex].InternalEnergy); // needs to be in code units
                    temp = ThermalProperties(u, SphP[pindex].Density * All.cf_a3inv, pindex, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);
                    *fp++ = (MyOutputFloat) SphP[pindex].MolecularMassFraction; /* we call the subroutine above to make sure this quantity is as up-to-the-moment updated as possible, going into our next routine */
#endif
                    n++;
                }
#endif
            break;

        case IO_POT:		/* gravitational potential */
#if defined(OUTPUT_POTENTIAL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Potential;
                    n++;
                }
#endif
            break;

        case IO_BH_DIST:
#if defined(BH_CALC_DISTANCES) && defined(OUTPUT_BH_DISTANCES)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].min_dist_to_bh;
                    n++;
                }
#endif
            break;

        case IO_ACCEL:		/* acceleration */
#ifdef OUTPUT_ACCELERATION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++) {fp[k] = (MyOutputFloat) (All.cf_a2inv * P[pindex].GravAccel[k]);}
#ifdef PMGRID
                    for(k = 0; k < 3; k++) {fp[k] += (MyOutputFloat) (All.cf_a2inv * P[pindex].GravPM[k]);}
#endif
                    if(P[pindex].Type == 0) {for(k = 0; k < 3; k++) {fp[k] += (MyOutputFloat) SphP[pindex].HydroAccel[k];}}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_DTENTR:		/* rate of change of internal energy */
#ifdef OUTPUT_CHANGEOFENERGY
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].DtInternalEnergy;
                    n++;
                }
#endif
            break;

        case IO_DELAYTIME:
#ifdef GALSF_SUBGRID_WINDS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].DelayTime;
                    n++;
                }
#endif
            break;

        case IO_TSTP:		/* timestep  */
#ifdef OUTPUT_TIMESTEP
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (GET_PARTICLE_TIMESTEP_IN_PHYSICAL(pindex));
                    n++;
                }
#endif
            break;

        case IO_BFLD:		/* magnetic field  */
#ifdef MAGNETIC
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) (Get_Gas_BField(pindex,k) * All.cf_a2inv * gizmo2gauss);}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_VDIV:		/* Divergence of Vel */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Particle_DivVel;
                    n++;
                }
            break;

        case IO_VORT:		/* Vorticity */
#if defined(TURB_DRIVING) || defined(OUTPUT_VORTICITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) SphP[pindex].Vorticity[k];}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_IMF:		/* parameters describing the IMF  */
#ifdef GALSF_SFR_IMF_VARIATION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < N_IMF_FORMPROPS; k++) {fp[k] = (MyOutputFloat) P[pindex].IMF_FormProps[k];}
                    fp += N_IMF_FORMPROPS;
                    n++;
                }
#endif
            break;

        case IO_COSMICRAY_ENERGY:	/* energy in the cosmic ray field  */
#ifdef COSMIC_RAY_FLUID
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<N_CR_PARTICLE_BINS;k++) {fp[k] = (MyOutputFloat) SphP[pindex].CosmicRayEnergyPred[k];}
                    fp += N_CR_PARTICLE_BINS;
                    n++;
                }
#endif
            break;

            
        case IO_COSMICRAY_SLOPES:    /* piecewise spectral slopes in the cosmic ray spectrum  */
            break;

        case IO_COSMICRAY_KAPPA:    /* local CR diffusion constant */
#if defined(COSMIC_RAY_FLUID) && defined(CRFLUID_DIFFUSION_MODEL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<N_CR_PARTICLE_BINS;k++) {fp[k] = (MyOutputFloat) SphP[pindex].CosmicRayDiffusionCoeff[k];}
                    fp += N_CR_PARTICLE_BINS;
                    n++;
                }
#endif
            break;

        case IO_COSMICRAY_ALFVEN:    /* energy in the resonant (~gyro-radii) Alfven modes field, in the +/- (with respect to B) fields  */
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<2;k++) {int k2; for(k2=0;k2<N_CR_PARTICLE_BINS;k2++) {fp[N_CR_PARTICLE_BINS*k + k2] = (MyOutputFloat) SphP[pindex].CosmicRayAlfvenEnergyPred[k2][k];}}
                    fp += 2*N_CR_PARTICLE_BINS;
                    n++;
                }
#endif
            break;


        case IO_DIVB:		/* divergence of magnetic field  */
#if defined(MAGNETIC) && defined(OUTPUT_BFIELD_DIVCLEAN_INFO)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                { /* divB is saved in physical units */
                    *fp++ = (MyOutputFloat) ((SphP[pindex].divB * gizmo2gauss * (SphP[pindex].Density*All.cf_a3inv / P[pindex].Mass)));
                    n++;
                }
#endif
            break;

        case IO_ABVC:		/* artificial viscosity of particle  */
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (SphP[pindex].alpha * SphP[pindex].alpha_limiter);
                    n++;
                }
#endif
            break;


        case IO_AMDC:		/* artificial magnetic dissipation of particle  */
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].Balpha;
                    n++;
                }
#endif
            break;

        case IO_PHI:		/* divBcleaning fuction of particle  */
#if defined(DIVBCLEANING_DEDNER) && defined(OUTPUT_BFIELD_DIVCLEAN_INFO)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (Get_Gas_PhiField(pindex) * All.cf_a3inv * gizmo2gauss);
                    n++;
                }
#endif
            break;

        case IO_GRADPHI:		/* divBcleaning fuction of particle  */
#if defined(DIVBCLEANING_DEDNER) && defined(OUTPUT_BFIELD_DIVCLEAN_INFO)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) (SphP[pindex].Gradients.Phi[k] * All.cf_a2inv*All.cf_a2inv * gizmo2gauss);}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_COOLRATE:		/* current cooling rate of particle  */
#ifdef OUTPUT_COOLRATE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    double ne = SphP[pindex].Ne; u = SphP[pindex].InternalEnergyPred; tcool = GetCoolingTime(u, SphP[pindex].Density * All.cf_a3inv, ne, pindex); /* get cooling time */
                    double coolrate_to_output = 0; if(tcool != 0) {coolrate_to_output = u / tcool;} /* convert cooling time with current thermal energy to du/dt */
                    *fp++ = (MyOutputFloat) coolrate_to_output;
                    n++;
                }
#endif
            break;

        case IO_BHMASS:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) BPP(pindex).BH_Mass;
                    n++;
                }
#endif
            break;

        case IO_BHDUSTMASS:
#if defined(BLACK_HOLES) && defined(GRAIN_FLUID)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) BPP(pindex).BH_Dust_Mass;
                    n++;
                }
#endif
            break;

        case IO_BHMASSALPHA:
#ifdef BH_ALPHADISK_ACCRETION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) BPP(pindex).BH_Mass_AlphaDisk;
                    n++;
                }
#endif
            break;

        case IO_BH_ANGMOM:
#ifdef BH_FOLLOW_ACCRETED_ANGMOM
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++) {fp[k] = (MyOutputFloat) BPP(pindex).BH_Specific_AngMom[k];}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_BHMDOT:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) BPP(pindex).BH_Mdot;
                    n++;
                }
#endif
            break;

        case IO_R_PROTOSTAR:
            break;

        case IO_MASS_D_PROTOSTAR:
            break;

        case IO_ZAMS_MASS:
            break;

        case IO_STAGE_PROTOSTAR:
            break;
            
        case IO_AGE_PROTOSTAR:
            break;

        case IO_LUM_SINGLESTAR:
            break;

        case IO_BHPROGS:
#ifdef BH_COUNTPROGS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) BPP(pindex).BH_CountProgs;
                    n++;
                }
#endif
            break;

        case IO_ACRB:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Hsml;
                    n++;
                }
#endif
            break;

        case IO_SINKRAD:
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].SinkRadius;
                    n++;
                }
#endif
            break;

        case IO_TIDALTENSORPS:   /* 3x3 configuration-space tidal tensor that is driving the GDE */
#ifdef OUTPUT_TIDAL_TENSOR
            for(n = 0; n < pc; pindex++)

                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                    {
                        int l_tt_tmp;
                        for(l_tt_tmp = 0; l_tt_tmp < 3; l_tt_tmp++)
                        {
                            fp[k * 3 + l_tt_tmp] = (MyOutputFloat) P[pindex].tidal_tensorps[k][l_tt_tmp];
#if defined(PMGRID) && !defined(GDE_DISTORTIONTENSOR)
                            fp[k * 3 + l_tt_tmp] += (MyOutputFloat) P[pindex].tidal_tensorpsPM[k][l_tt_tmp]; // in current code (without GDE_DISTORTIONTENSOR) this isn't necessary because of how the tidal tensor terms are added and diagonalized already in the gravtree operations
#endif

                        }
                    }
                    //fflush(stderr);
                    n++;
                    fp += 9;
                }
#endif
            break;

        case IO_GDE_DISTORTIONTENSOR:   /* full 6D phase-space distortion tensor from GDE integration */
#ifdef OUTPUT_GDE_DISTORTIONTENSOR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    get_half_kick_distortion(pindex, half_kick_add);
                    for(k = 0; k < 6; k++)
                    {
                        for(l = 0; l < 6; l++)
                        {
                            fp[k * 6 + l] = (MyOutputFloat) (P[pindex].distortion_tensorps[k][l] + half_kick_add[k][l]);
                        }
                    }
                    n++;
                    fp += 36;

                }
#endif
            break;

        case IO_CAUSTIC_COUNTER:   /* caustic counter */
#ifdef GDE_DISTORTIONTENSOR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].caustic_counter;
                    n++;
                }
#endif
            break;

        case IO_FLOW_DETERMINANT:   /* physical NON-CUTOFF corrected stream determinant = 1.0/normed stream density * 1.0/initial stream density */
#if defined(GDE_DISTORTIONTENSOR) && !defined(GDE_LEAN)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    get_current_ps_info(pindex, &flde, &psde);
                    *fp++ = (MyOutputFloat) flde;
                    n++;
                }
#endif
            break;

        case IO_STREAM_DENSITY:   /* physical CUTOFF corrected stream density = normed stream density * initial stream density */
#ifdef GDE_DISTORTIONTENSOR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (P[pindex].stream_density);
                    n++;
                }
#endif
            break;

        case IO_PHASE_SPACE_DETERMINANT:   /* determinant of phase-space distortion tensor -> should be 1 due to Liouville theorem */
#ifdef GDE_DISTORTIONTENSOR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    get_current_ps_info(pindex, &flde, &psde);
                    *fp++ = (MyOutputFloat) psde;
                    n++;
                }
#endif
            break;

        case IO_ANNIHILATION_RADIATION:   /* time integrated stream density in physical units */
#if defined(GDE_DISTORTIONTENSOR) && !defined(GDE_LEAN)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (P[pindex].annihilation * GDE_INITDENSITY(pindex));
                    *fp++ = (MyOutputFloat) (P[pindex].analytic_caustics);
                    *fp++ = (MyOutputFloat) (P[pindex].analytic_annihilation * GDE_INITDENSITY(pindex));
                    n++;
                }
#endif
            break;

        case IO_LAST_CAUSTIC:   /* extensive information on the last caustic the particle has passed */
#ifdef OUTPUT_GDE_LASTCAUSTIC
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].lc_Time;
                    *fp++ = (MyOutputFloat) P[pindex].lc_Pos[0];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Pos[1];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Pos[2];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Vel[0];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Vel[1];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Vel[2];
                    *fp++ = (MyOutputFloat) P[pindex].lc_rho_normed_cutoff;

                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_x[0];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_x[1];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_x[2];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_y[0];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_y[1];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_y[2];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_z[0];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_z[1];
                    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_z[2];

                    *fp++ = (MyOutputFloat) P[pindex].lc_smear_x;
                    *fp++ = (MyOutputFloat) P[pindex].lc_smear_y;
                    *fp++ = (MyOutputFloat) P[pindex].lc_smear_z;
                    n++;
                }
#endif
            break;

        case IO_SHEET_ORIENTATION:   /* initial orientation of the CDM sheet where the particle started */
#if defined(GDE_DISTORTIONTENSOR) && (!defined(GDE_LEAN) || defined(GDE_READIC))
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,0,0);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,0,1);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,0,2);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,1,0);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,1,1);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,1,2);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,2,0);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,2,1);
                    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex,2,2);
                    n++;
                }
#endif
            break;

        case IO_INIT_DENSITY:   /* initial stream density in physical units  */
#if defined(GDE_DISTORTIONTENSOR) && (!defined(GDE_LEAN) || defined(GDE_READIC))
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    if(All.ComovingIntegrationOn)
                        {*fp++ = (MyOutputFloat) (GDE_INITDENSITY(pindex) / (GDE_TIMEBEGIN(pindex) * GDE_TIMEBEGIN(pindex) * GDE_TIMEBEGIN(pindex)));}
                    else
                        {*fp++ = (MyOutputFloat) GDE_INITDENSITY(pindex);}
                    n++;
                }
#endif
            break;

        case IO_EOSABAR:
#ifdef EOS_CARRIES_ABAR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].Abar;
                    n++;
                }
#endif
            break;

        case IO_TURB_DYNAMIC_COEFF:
#ifdef TURB_DIFF_DYNAMIC
            for (n = 0; n < pc; pindex++) {
                if (P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].TD_DynDiffCoeff;
                    n++;
                }
            }
#endif
            break;

        case IO_EOSYE:
#ifdef EOS_CARRIES_YE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].Ye;
                    n++;
                }
#endif
            break;

        case IO_EOSTEMP:
#ifdef EOS_CARRIES_TEMPERATURE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].Temperature;
                    n++;
                }
#elif defined(OUTPUT_TEMPERATURE)
	    for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
		    double u, ne, nh0 = 0, mu = 1, temp, nHeII, nhp, nHe0, nHepp; u = DMAX(All.MinEgySpec, SphP[pindex].InternalEnergy); // needs to be in code units                                                                                                                                                  
		    temp = ThermalProperties(u, SphP[pindex].Density * All.cf_a3inv, pindex, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);		  
		    *fp++ = (MyOutputFloat) temp;
		    n++;
                }	    
#endif
            break;

        case IO_PRESSURE:
#if defined(EOS_GENERAL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].Pressure;
                    n++;
                }
#endif
            break;

            case IO_EOSCS:
#if defined(EOS_GENERAL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) SphP[pindex].SoundSpeed;
                    n++;
                }
#endif
            break;

        case IO_CBE_MOMENTS:
            break;

        case IO_EOS_STRESS_TENSOR:
#if defined(EOS_ELASTIC)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                    {
                        int kf;
                        for(kf = 0; kf < 3; kf++)
                            fp[3*k + kf] = (MyOutputFloat) SphP[pindex].Elastic_Stress_Tensor[kf][k];
                    }
                    fp += 9;
                    n++;
                }
#endif
            break;

            case IO_EOSCOMP:
#ifdef EOS_TILLOTSON
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) SphP[pindex].CompositionType;
                    n++;
                }
#endif
            break;

        case IO_PARTVEL:
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++) {fp[k] = (MyOutputFloat) SphP[pindex].ParticleVel[k];}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_RADGAMMA:
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<N_RT_FREQ_BINS;k++) {fp[k] = (MyOutputFloat) SphP[pindex].Rad_E_gamma[k];}
                    fp += N_RT_FREQ_BINS;
                    n++;
                }
#endif
            break;

        case IO_RAD_FLUX:
#if defined(OUTPUT_RT_RAD_FLUX) && defined(RT_EVOLVE_FLUX)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {fp[N_RT_FREQ_BINS*k + kf] = (MyOutputFloat) (SphP[pindex].Rad_Flux_Pred[kf][k] * (SphP[pindex].Density*All.cf_a3inv/P[pindex].Mass));}}
                    fp += 3*N_RT_FREQ_BINS;
                    n++;
                }
#endif
            break;

        case IO_RAD_ACCEL:
#ifdef RT_RAD_PRESSURE_OUTPUT
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) SphP[pindex].Rad_Accel[k];}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_EDDINGTON_TENSOR:
#ifdef RADTRANSFER
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<6;k++) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {fp[N_RT_FREQ_BINS*k + kf] = (MyOutputFloat) SphP[pindex].ET[kf][k];}}
                    fp += 6*N_RT_FREQ_BINS;
                    n++;
                }
#endif
            break;

        case IO_AGS_SOFT:		/* Adaptive Gravitational Softening: softening */
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTGRAVSOFT)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) PPP[pindex].AGS_Hsml;
                    n++;
                }
#endif
            break;

        case IO_AGS_RHO:        /* Adaptive Gravitational Softening: density */
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) PPP[pindex].AGS_Density;
                    n++;
                }
#endif
            break;

        case IO_AGS_QPT:        /* quantum potential (Q) */
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    double d2rho = P[pindex].AGS_Gradients2_Density[0][0] + P[pindex].AGS_Gradients2_Density[1][1] + P[pindex].AGS_Gradients2_Density[2][2]; // laplacian
                    double drho2 = P[pindex].AGS_Gradients_Density[0]*P[pindex].AGS_Gradients_Density[0] + P[pindex].AGS_Gradients_Density[1]*P[pindex].AGS_Gradients_Density[1] + P[pindex].AGS_Gradients_Density[2]*P[pindex].AGS_Gradients_Density[2];
                    double AGS_QuantumPotential = (0.25*All.ScalarField_hbar_over_mass*All.ScalarField_hbar_over_mass / P[pindex].AGS_Density) * (d2rho - 0.5*drho2/P[pindex].AGS_Density);
                    *fp++ = (MyOutputFloat) AGS_QuantumPotential;
                    n++;
                }
#endif
            break;

        case IO_AGS_PSI_RE:        /* real part of wavefunction */
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
#if (DM_FUZZY > 0)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (P[pindex].AGS_Psi_Re * P[pindex].AGS_Density / P[pindex].Mass);
                    n++;
                }
#endif
#endif
            break;

        case IO_AGS_PSI_IM:        /* imaginary part of wavefunction */
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
#if (DM_FUZZY > 0)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (P[pindex].AGS_Psi_Im * P[pindex].AGS_Density / P[pindex].Mass);
                    n++;
                }
#endif
#endif
            break;

        case IO_AGS_ZETA:		/* Adaptive Gravitational Softening: zeta */
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTZETA)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) PPPZ[pindex].AGS_zeta;
                    n++;
                }
#endif
            break;

        case IO_grHI:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grHI;
                    n++;
                }
            }
#endif
            break;

        case IO_grHII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grHII;
                    n++;
                }
            }
#endif
            break;

        case IO_grHM:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grHM;
                    n++;
                }
            }
#endif
            break;

        case IO_grHeI:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grHeI;
                    n++;
                }
            }
#endif
            break;

        case IO_grHeII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grHeII;
                    n++;
                }
            }
#endif
            break;

        case IO_grHeIII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grHeIII;
                    n++;
                }
            }
#endif
            break;

        case IO_grH2I:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grH2I;
                    n++;
                }
            }
#endif
            break;

        case IO_grH2II:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grH2II;
                    n++;
                }
            }
#endif
            break;

        case IO_grDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grDI;
                    n++;
                }
            }
#endif
            break;

        case IO_grDII:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grDII;
                    n++;
                }
            }
#endif
            break;

        case IO_grHDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) SphP[pindex].grHDI;
                    n++;
                }
            }
#endif
            break;

    case IO_TURB_DIFF_COEFF:
#ifdef TURB_DIFF_DYNAMIC
        for (n = 0; n < pc; pindex++) {
            if (P[pindex].Type == type) {
                *fp++ = (MyOutputFloat) SphP[pindex].TD_DiffCoeff;
                n++;
            }
        }
#endif

        break;

    case IO_DYNERROR:
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        for (n = 0; n < pc; pindex++) {
            if (P[pindex].Type == type) {
                *fp++ = (MyOutputFloat) SphP[pindex].TD_DynDiffCoeff_error;
                n++;
            }
        }
#endif
        break;

    case IO_DYNERRORDEFAULT:
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        for (n = 0; n < pc; pindex++) {
            if (P[pindex].Type == type) {
                *fp++ = (MyOutputFloat) SphP[pindex].TD_DynDiffCoeff_error_default;
                n++;
            }
        }
#endif
        break;

        case IO_LASTENTRY:
            endrun(213);
            break;
    }

    *startindex = pindex;
}




/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file.
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
    int bytes_per_blockelement = 0;

    switch (blocknr)
    {
        case IO_POS:
            if(mode)
                bytes_per_blockelement = 3 * sizeof(MyInputPosFloat);
            else
                bytes_per_blockelement = 3 * sizeof(MyOutputPosFloat);
            break;

        case IO_VEL:
        case IO_PARTVEL:
        case IO_ACCEL:
        case IO_BFLD:
        case IO_INIB:
        case IO_GRADPHI:
        case IO_RAD_ACCEL:
        case IO_VORT:
        case IO_BH_ANGMOM:
        case IO_ANNIHILATION_RADIATION:
            if(mode)
                bytes_per_blockelement = 3 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
            break;

        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
            bytes_per_blockelement = sizeof(MyIDType);
            break;

        case IO_BHPROGS:
        case IO_GRAINTYPE:
        case IO_EOSCOMP:
        case IO_STAGE_PROTOSTAR:
            bytes_per_blockelement = sizeof(int);
            break;
            
        case IO_AGE_PROTOSTAR:
        case IO_MASS:
        case IO_BH_DIST:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_IDEN:
        case IO_UNSPMASS:
        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_HSML:
        case IO_SFR:
        case IO_AGE:
        case IO_OSTAR:
        case IO_GRAINSIZE:
        case IO_DELAYTIME:
        case IO_HSMS:
        case IO_POT:
        case IO_DTENTR:
        case IO_TSTP:
        case IO_DIVB:
        case IO_VDIV:
        case IO_ABVC:
        case IO_AMDC:
        case IO_PHI:
        case IO_COOLRATE:
        case IO_BHMASS:
        case IO_BHDUSTMASS:
        case IO_BHMASSALPHA:
        case IO_ACRB:
        case IO_SINKRAD:
        case IO_BHMDOT:
        case IO_R_PROTOSTAR:
        case IO_MASS_D_PROTOSTAR:
        case IO_ZAMS_MASS:
        case IO_LUM_SINGLESTAR:
        case IO_CAUSTIC_COUNTER:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_EOSTEMP:
        case IO_EOSABAR:
        case IO_EOSYE:
        case IO_EOSCS:
        case IO_PRESSURE:
        case IO_INIT_DENSITY:
        case IO_AGS_SOFT:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_PSI_RE:
        case IO_AGS_PSI_IM:
        case IO_AGS_ZETA:
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
        case IO_grH2I:
        case IO_grH2II:
        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
        case IO_TURB_DYNAMIC_COEFF:
        case IO_TURB_DIFF_COEFF:
        case IO_DYNERROR:
        case IO_DYNERRORDEFAULT:
        case IO_CHIMES_MU:
        case IO_CHIMES_NH:
        case IO_CHIMES_STAR_SIGMA:
        case IO_DENS_AROUND_STAR:
        case IO_DELAY_TIME_HII:
        case IO_MOLECULARFRACTION:
            if(mode)
                bytes_per_blockelement = sizeof(MyInputFloat);
            else
                bytes_per_blockelement = sizeof(MyOutputFloat);
            break;


        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_SLOPES:
        case IO_COSMICRAY_KAPPA:
#if defined(COSMIC_RAY_FLUID) || defined(FLAG_NOT_IN_PUBLIC_CODE)
            if(mode)
                bytes_per_blockelement = (N_CR_PARTICLE_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (N_CR_PARTICLE_BINS) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_COSMICRAY_ALFVEN:
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            if(mode)
                bytes_per_blockelement = (2 * N_CR_PARTICLE_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (2 * N_CR_PARTICLE_BINS) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_IMF:
#ifdef GALSF_SFR_IMF_VARIATION
            if(mode)
                bytes_per_blockelement = (N_IMF_FORMPROPS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (N_IMF_FORMPROPS) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_RADGAMMA:
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
            if(mode)
                bytes_per_blockelement = (N_RT_FREQ_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (N_RT_FREQ_BINS) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_RAD_FLUX:
#ifdef RADTRANSFER
            if(mode)
                bytes_per_blockelement = (3*N_RT_FREQ_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (3*N_RT_FREQ_BINS) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_EDDINGTON_TENSOR:
#ifdef RADTRANSFER
            if(mode)
                bytes_per_blockelement = (6*N_RT_FREQ_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (6*N_RT_FREQ_BINS) * sizeof(MyOutputFloat);
#endif
            break;


        case IO_Z:
#ifdef METALS
            if(mode)
                bytes_per_blockelement = (NUM_METAL_SPECIES) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (NUM_METAL_SPECIES) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_CHIMES_ABUNDANCES:
#ifdef CHIMES
            if(mode)
                bytes_per_blockelement = ChimesGlobalVars.totalNumberOfSpecies * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = ChimesGlobalVars.totalNumberOfSpecies * sizeof(MyOutputFloat);
#endif
            break;

        case IO_CHIMES_REDUCED:
            if(mode)
                bytes_per_blockelement = 4 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 4 * sizeof(MyOutputFloat);
            break;

        case IO_CHIMES_FLUX_G0:
        case IO_CHIMES_FLUX_ION:
#if defined(CHIMES) && defined(CHIMES_STELLAR_FLUXES)
            if(mode)
                bytes_per_blockelement = CHIMES_LOCAL_UV_NBINS * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = CHIMES_LOCAL_UV_NBINS * sizeof(MyOutputFloat);
#endif
            break;

        case IO_CBE_MOMENTS:

        case IO_EOS_STRESS_TENSOR:
        case IO_TIDALTENSORPS:
        case IO_SHEET_ORIENTATION:
            if(mode)
                bytes_per_blockelement = 9 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
            break;

        case IO_GDE_DISTORTIONTENSOR:
            if(mode)
                bytes_per_blockelement = 36 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 36 * sizeof(MyOutputFloat);
            break;

        case IO_LAST_CAUSTIC:
            if(mode)
                bytes_per_blockelement = 20 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 20 * sizeof(MyOutputFloat);
            break;

        case IO_LASTENTRY:
            endrun(214);
            break;
    }

    return bytes_per_blockelement;
}

int get_datatype_in_block(enum iofields blocknr)
{
    int typekey;
    switch (blocknr)
    {
#if (defined(OUTPUT_POSITIONS_IN_DOUBLE) && !defined(OUTPUT_IN_DOUBLEPRECISION)) || (defined(INPUT_POSITIONS_IN_DOUBLE) && !defined(INPUT_IN_DOUBLEPRECISION))
        case IO_POS:
            typekey = 3; /* pos outputs in HDF5 are double automatically, to prevent overlaps */
            break;
#endif

        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
#ifdef LONGIDS
            typekey = 2;		/* native long long */
#else
            typekey = 0;		/* native int */
#endif
            break;

        case IO_BHPROGS:
        case IO_GRAINTYPE:
        case IO_EOSCOMP:
        case IO_STAGE_PROTOSTAR:
            typekey = 0;		/* native int */
            break;
            
        default:
            typekey = 1;		/* native MyOutputFloat */
            break;
    }

    return typekey;
}



int get_values_per_blockelement(enum iofields blocknr)
{
    int values = 0;
    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_INIB:
        case IO_PARTVEL:
        case IO_ACCEL:
        case IO_BFLD:
        case IO_GRADPHI:
        case IO_RAD_ACCEL:
        case IO_VORT:
        case IO_BH_ANGMOM:
        case IO_ANNIHILATION_RADIATION:
            values = 3;
            break;

        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
        case IO_MASS:
        case IO_BH_DIST:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_IDEN:
        case IO_UNSPMASS:
        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_HSML:
        case IO_SFR:
        case IO_AGE:
        case IO_OSTAR:
        case IO_GRAINSIZE:
        case IO_GRAINTYPE:
        case IO_DELAYTIME:
        case IO_HSMS:
        case IO_POT:
        case IO_DTENTR:
        case IO_TSTP:
        case IO_VDIV:
        case IO_DIVB:
        case IO_ABVC:
        case IO_AMDC:
        case IO_PHI:
        case IO_COOLRATE:
        case IO_BHMASS:
        case IO_BHDUSTMASS:
        case IO_BHMASSALPHA:
        case IO_ACRB:
        case IO_SINKRAD:
        case IO_BHMDOT:
        case IO_R_PROTOSTAR:
        case IO_MASS_D_PROTOSTAR:
        case IO_ZAMS_MASS:
        case IO_STAGE_PROTOSTAR:
        case IO_AGE_PROTOSTAR:
        case IO_LUM_SINGLESTAR:
        case IO_BHPROGS:
        case IO_CAUSTIC_COUNTER:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_EOSTEMP:
        case IO_EOSABAR:
        case IO_EOSCS:
        case IO_EOSCOMP:
        case IO_EOSYE:
        case IO_PRESSURE:
        case IO_INIT_DENSITY:
        case IO_AGS_SOFT:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_PSI_RE:
        case IO_AGS_PSI_IM:
        case IO_AGS_ZETA:
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
        case IO_grH2I:
        case IO_grH2II:
        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
        case IO_TURB_DYNAMIC_COEFF:
        case IO_TURB_DIFF_COEFF:
        case IO_DYNERROR:
        case IO_DYNERRORDEFAULT:
        case IO_CHIMES_MU:
        case IO_CHIMES_NH:
        case IO_CHIMES_STAR_SIGMA:
        case IO_DENS_AROUND_STAR:
        case IO_DELAY_TIME_HII:
        case IO_MOLECULARFRACTION:
            values = 1;
            break;

        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_SLOPES:
        case IO_COSMICRAY_KAPPA:
#if defined(COSMIC_RAY_FLUID) || defined(FLAG_NOT_IN_PUBLIC_CODE)
            values = N_CR_PARTICLE_BINS;
#endif
            break;

        case IO_COSMICRAY_ALFVEN:
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            values = (2*N_CR_PARTICLE_BINS);
#endif
            break;

        case IO_CBE_MOMENTS:
            break;

        case IO_EDDINGTON_TENSOR:
#ifdef RADTRANSFER
            values = (6*N_RT_FREQ_BINS);
#endif
            break;

        case IO_RAD_FLUX:
#ifdef RADTRANSFER
            values = (3*N_RT_FREQ_BINS);
#endif
            break;

        case IO_RADGAMMA:
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
            values = N_RT_FREQ_BINS;
#endif
            break;

        case IO_Z:
#ifdef METALS
            values = NUM_METAL_SPECIES;
#endif
            break;

        case IO_CHIMES_ABUNDANCES:
#ifdef CHIMES
            values = ChimesGlobalVars.totalNumberOfSpecies;
#endif
            break;

        case IO_CHIMES_REDUCED:
            values = 4;
            break;

        case IO_CHIMES_FLUX_G0:
        case IO_CHIMES_FLUX_ION:
#ifdef CHIMES_STELLAR_FLUXES
            values = CHIMES_LOCAL_UV_NBINS;
#endif
            break;

        case IO_IMF:
#ifdef GALSF_SFR_IMF_VARIATION
            values = N_IMF_FORMPROPS;
#endif
            break;

        case IO_TIDALTENSORPS:
        case IO_SHEET_ORIENTATION:
        case IO_EOS_STRESS_TENSOR:
            values = 9;
            break;

        case IO_GDE_DISTORTIONTENSOR:
            values = 36;
            break;

        case IO_LAST_CAUSTIC:
            values = 20;
            break;

        case IO_LASTENTRY:
            endrun(215);
            break;
    }
    return values;
}




/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 */
long get_particles_in_block(enum iofields blocknr, int *typelist)
{
    long i, nall, nsel, ntot_withmasses, ngas, nstars, nngb, nstars_tot;
    nall = 0; nsel = 0; ntot_withmasses = 0; nstars_tot=0;

    int valid_star_types=16; if(All.ComovingIntegrationOn==0) {valid_star_types+=4+8;} /* used in e.g. ages block below, not for all, but easier to define here */
#ifdef BLACK_HOLES
    valid_star_types += 32;
#endif
    for(i = 0; i < 6; i++)
    {
        typelist[i] = 0;
        if(header.npart[i] > 0)
        {
            nall += header.npart[i];
            typelist[i] = 1;
            if((1 << i) & (valid_star_types)) {nstars_tot += header.npart[i];}
        }
        if(All.MassTable[i] == 0) {ntot_withmasses += header.npart[i];}
    }
    ngas = header.npart[0];
    nstars = header.npart[4];


    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_TSTP:
        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
        case IO_POT:
        case IO_AGS_SOFT:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_PSI_RE:
        case IO_AGS_PSI_IM:
        case IO_AGS_ZETA:
        case IO_BH_DIST:
        case IO_CBE_MOMENTS:
            return nall;
            break;

        case IO_MASS:
            for(i = 0; i < 6; i++)
            {
                typelist[i] = 0;
                if(All.MassTable[i] == 0 && header.npart[i] > 0) {typelist[i] = 1;}
            }
            return ntot_withmasses;
            break;

        case IO_PARTVEL:
        case IO_RAD_ACCEL:
        case IO_RADGAMMA:
        case IO_RAD_FLUX:
        case IO_EDDINGTON_TENSOR:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_INIB:
        case IO_IDEN:
        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_HSML:
        case IO_DELAYTIME:
        case IO_SFR:
        case IO_DTENTR:
        case IO_BFLD:
        case IO_VDIV:
        case IO_VORT:
        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_SLOPES:
        case IO_COSMICRAY_KAPPA:
        case IO_COSMICRAY_ALFVEN:
        case IO_DIVB:
        case IO_ABVC:
        case IO_AMDC:
        case IO_PHI:
        case IO_GRADPHI:
        case IO_COOLRATE:
        case IO_EOSTEMP:
        case IO_EOSABAR:
        case IO_EOSYE:
        case IO_EOSCS:
        case IO_EOS_STRESS_TENSOR:
        case IO_EOSCOMP:
        case IO_PRESSURE:
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
        case IO_grH2I:
        case IO_grH2II:
        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
        case IO_TURB_DYNAMIC_COEFF:
        case IO_TURB_DIFF_COEFF:
        case IO_DYNERROR:
        case IO_DYNERRORDEFAULT:
        case IO_DELAY_TIME_HII:
        case IO_MOLECULARFRACTION:
        case IO_CHIMES_ABUNDANCES:
        case IO_CHIMES_MU:
        case IO_CHIMES_REDUCED:
        case IO_CHIMES_NH:
        case IO_CHIMES_FLUX_G0:
        case IO_CHIMES_FLUX_ION:
            for(i = 1; i < 6; i++) {typelist[i] = 0;}
            return ngas;
            break;

        case IO_AGE:
            for(i=0; i<6; i++) {if(!((1 << i) & (valid_star_types))) {typelist[i]=0;}}
            return nstars_tot;
            break;

        case IO_OSTAR:
        case IO_HSMS:
            for(i = 0; i < 6; i++) {if(i != 4) {typelist[i] = 0;}}
            return nstars;
            break;

        case IO_GRAINSIZE:
        case IO_GRAINTYPE:
            nngb=0;
#ifdef GRAIN_FLUID
            for(i=0;i<6;i++) {if((1 << i) & (GRAIN_PTYPES)) {nngb+=header.npart[i];} else {typelist[i]=0;}}
#endif
            return nngb;
            break;

        case IO_IMF:
            for(i = 1; i < 6; i++) {if(i != 4 && i != 5) {typelist[i] = 0;}}
            return nstars + header.npart[5];
            break;

        case IO_Z:
            for(i = 0; i < 6; i++) {if(i != 0 && i != 4 && i != 5) {typelist[i] = 0;}}
            return ngas + nstars + header.npart[5];
            break;

        case IO_CHIMES_STAR_SIGMA:
            nngb = nstars;
            typelist[0]=typelist[1]=typelist[5]=0;
            if(All.ComovingIntegrationOn) {typelist[2]=typelist[3]=0;} else {nngb += header.npart[2] + header.npart[3];}
            return nngb;
            break;

        case IO_DENS_AROUND_STAR:
            typelist[0] = 0;
            return header.npart[1]+header.npart[2]+header.npart[3]+header.npart[4]+header.npart[5];
            break;

        case IO_BHMASS:
        case IO_BHDUSTMASS:
        case IO_BHMASSALPHA:
        case IO_BH_ANGMOM:
        case IO_UNSPMASS:
        case IO_ACRB:
        case IO_SINKRAD:
        case IO_BHMDOT:
        case IO_R_PROTOSTAR:
        case IO_MASS_D_PROTOSTAR:
        case IO_ZAMS_MASS:
        case IO_STAGE_PROTOSTAR:
        case IO_AGE_PROTOSTAR:
        case IO_LUM_SINGLESTAR:
        case IO_BHPROGS:
            for(i = 0; i < 6; i++) {if(i != 5) {typelist[i] = 0;}}
            return header.npart[5];
            break;

        case IO_TIDALTENSORPS:
        case IO_GDE_DISTORTIONTENSOR:
        case IO_CAUSTIC_COUNTER:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_ANNIHILATION_RADIATION:
        case IO_LAST_CAUSTIC:
        case IO_SHEET_ORIENTATION:
        case IO_INIT_DENSITY:
            for(i = 0; i < 6; i++) {if(((1 << i) & (GDE_TYPES))) {nsel += header.npart[i];} else {typelist[i] = 0;}}
            return nsel;
            break;

        case IO_LASTENTRY:
            endrun(216);
            break;
    }

    endrun(212);
    return 0;
}



/*! This function tells whether a block in the output file is present or not.
 */
int blockpresent(enum iofields blocknr)
{
    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
        case IO_MASS:
        case IO_U:
        case IO_RHO:
        case IO_HSML:
            return 1;			/* always present */
            break;

        case IO_NE:
        case IO_NH:
#if (defined(COOLING) || defined(RADTRANSFER)) && !defined(CHIMES)
            return 1;
#endif
            break;

        case IO_RADGAMMA:
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
            return 1;
#endif
            break;

        case IO_RAD_FLUX:
#if defined(OUTPUT_RT_RAD_FLUX) && defined(RT_EVOLVE_FLUX)
            return 1;
#endif
            break;

        case IO_RAD_ACCEL:
#if defined(RT_RAD_PRESSURE_OUTPUT)
            return 1;
#endif
            break;

        case IO_HSMS:
#if defined(SUBFIND)
            return 1;
#endif
            break;

        case IO_SFR:
        case IO_AGE:
#ifdef GALSF
            return 1;
#endif
            break;

        case IO_GRAINSIZE:
#ifdef GRAIN_FLUID
            return 1;
#endif
            break;

        case IO_GRAINTYPE:
#ifdef PIC_MHD
            return 1;
#endif
            break;

        case IO_Z:
#ifdef METALS
            return 1;
#endif
            break;

        case IO_CHIMES_ABUNDANCES:
#if defined(CHIMES_REDUCED_OUTPUT)
            if(Chimes_incl_full_output == 1) {return 1;} else {return 0;}
#elif defined(CHIMES)
            return 1;
#endif
            break;

        case IO_CHIMES_MU:
#ifdef CHIMES
            return 1;
#endif
            break;

        case IO_CHIMES_REDUCED:
#if defined(CHIMES) && defined(CHIMES_REDUCED_OUTPUT)
            if(Chimes_incl_full_output == 0) {return 1;} else {return 0;}
#endif
            break;

        case IO_CHIMES_NH:
#ifdef CHIMES_NH_OUTPUT
            return 1;
#endif
            break;

        case IO_CHIMES_STAR_SIGMA:
#if defined(CHIMES_NH_OUTPUT) && defined(OUTPUT_DENS_AROUND_STAR)
            return 1;
#endif
            break;

        case IO_CHIMES_FLUX_G0:
#ifdef CHIMES_STELLAR_FLUXES
            return 1;
#endif
            break;

        case IO_CHIMES_FLUX_ION:
#ifdef CHIMES_STELLAR_FLUXES
            return 1;
#endif
            break;

        case IO_DENS_AROUND_STAR:
#ifdef OUTPUT_DENS_AROUND_STAR
            return 1;
#endif
            break;

        case IO_DELAY_TIME_HII:
            break;

        case IO_DELAYTIME:
#ifdef GALSF_SUBGRID_WINDS
            return 1;
#endif
            break;

        case IO_MOLECULARFRACTION:
#if defined(OUTPUT_MOLECULAR_FRACTION)
            return 1;
#endif
            break;

        case IO_HII:
#if defined(RT_CHEM_PHOTOION)
            return 1;
#endif
            break;

        case IO_HeI:
        case IO_HeII:
#if defined(RT_CHEM_PHOTOION_HE)
            return 1;
#endif
            break;
            
        case IO_IDEN:
        case IO_INIB:
#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) && defined(BH_DEBUG_SPAWN_JET_TEST)
            return 1;
#endif         
            break;

        case IO_UNSPMASS:
#if defined(BH_WIND_SPAWN) && defined(BH_DEBUG_SPAWN_JET_TEST)
            return 1;
#endif   
            break;

        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            return 1;
#endif
            break;

        case IO_POT:
#if defined(OUTPUT_POTENTIAL)
            return 1;
#endif
            break;

        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
#if defined(TURB_DRIVING)
            return 1;
#endif
            break;

        case IO_ACCEL:
#ifdef OUTPUT_ACCELERATION
            return 1;
#endif
            break;

        case IO_DTENTR:
#ifdef OUTPUT_CHANGEOFENERGY
            return 1;
#endif
            break;

        case IO_TSTP:
#ifdef OUTPUT_TIMESTEP
            return 1;
#endif
            break;

        case IO_BFLD:
#if defined(MAGNETIC)
            return 1;
#endif
            break;
            
        case IO_DIVB:
#if defined(MAGNETIC) && defined(OUTPUT_BFIELD_DIVCLEAN_INFO)
            return 1;
#endif
            break;

        case IO_VDIV:
        case IO_VORT:
#if defined(TURB_DRIVING) || defined(OUTPUT_VORTICITY)
            return 1;
#else
            return 0;
#endif
            break;

        case IO_IMF:
#ifdef GALSF_SFR_IMF_VARIATION
            return 1;
#endif
            break;

        case IO_OSTAR:
#ifdef GALSF_SFR_IMF_SAMPLING
            return 1;
#endif
            break;

        case IO_COSMICRAY_ENERGY:
#if defined(COSMIC_RAY_FLUID) || defined(FLAG_NOT_IN_PUBLIC_CODE)
            return 1;
#endif
            break;

        case IO_COSMICRAY_SLOPES:
            break;

        case IO_COSMICRAY_KAPPA:
#if defined(COSMIC_RAY_FLUID) && defined(CRFLUID_DIFFUSION_MODEL)
            return 1;
#endif
            break;

        case IO_COSMICRAY_ALFVEN:
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            return 1;
#endif
            break;

        case IO_ABVC:
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
            return 1;
#endif
            break;

        case IO_AMDC:
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
            return 1;
#endif
            break;

        case IO_PHI:
        case IO_GRADPHI:
#if defined(DIVBCLEANING_DEDNER) && defined(OUTPUT_BFIELD_DIVCLEAN_INFO)
            return 1;
#endif
            break;

        case IO_COOLRATE:
#ifdef OUTPUT_COOLRATE
            return 1;
#endif
            break;

        case IO_BH_ANGMOM:
#ifdef BH_FOLLOW_ACCRETED_ANGMOM
            return 1;
#endif
            break;

        case IO_ACRB:
        case IO_SINKRAD:
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
            return 1;
#endif
            break;

        case IO_BHMASS:
#ifdef BLACK_HOLES
            return 1;
#endif
            break;

        case IO_BHDUSTMASS:
#if defined(BLACK_HOLES) && defined(GRAIN_FLUID)
            return 1;
#endif
            break;

        case IO_BHMASSALPHA:
#ifdef BH_ALPHADISK_ACCRETION
            return 1;
#endif
            break;

        case IO_BHMDOT:
#ifdef BLACK_HOLES
            return 1;
#endif
            break;

        case IO_R_PROTOSTAR:
            break;

        case IO_MASS_D_PROTOSTAR:
            break;

        case IO_ZAMS_MASS:
            break;

        case IO_STAGE_PROTOSTAR:
            break;
            
        case IO_AGE_PROTOSTAR:
            break;

        case IO_LUM_SINGLESTAR:
            break;

        case IO_BH_DIST:
#if defined(BH_CALC_DISTANCES) && defined(OUTPUT_BH_DISTANCES)
            return 1;
#endif
            break;

        case IO_BHPROGS:
#ifdef BH_COUNTPROGS
            return 1;
#endif
            break;

        case IO_TIDALTENSORPS:
#ifdef OUTPUT_TIDAL_TENSOR
            return 1;
#endif
            break;

        case IO_GDE_DISTORTIONTENSOR:
        case IO_CAUSTIC_COUNTER:
#ifdef GDE_DISTORTIONTENSOR
            return 1;
#endif
            break;

        case IO_FLOW_DETERMINANT:
#if defined(GDE_DISTORTIONTENSOR) && !defined(GDE_LEAN)
            return 1;
#endif
            break;

        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
#ifdef GDE_DISTORTIONTENSOR
            return 1;
#endif
            break;

        case IO_ANNIHILATION_RADIATION:
#if defined(GDE_DISTORTIONTENSOR) && !defined(GDE_LEAN)
            return 1;
#endif
            break;

        case IO_LAST_CAUSTIC:
#ifdef OUTPUT_GDE_LASTCAUSTIC
            return 1;
#endif
            break;

        case IO_SHEET_ORIENTATION:
        case IO_INIT_DENSITY:
#if defined(GDE_DISTORTIONTENSOR) && (!defined(GDE_LEAN) || defined(GDE_READIC))
            return 1;
#endif
            break;

        case IO_EOSTEMP:
#ifdef EOS_CARRIES_TEMPERATURE
            return 1;
#elif defined(OUTPUT_TEMPERATURE)
	    return 1;
#endif
            break;

        case IO_PRESSURE:
        case IO_EOSCS:
#ifdef EOS_GENERAL
            return 1;
#endif
            break;

        case IO_EOS_STRESS_TENSOR:
#ifdef EOS_ELASTIC
            return 1;
#endif
            break;

        case IO_EOSCOMP:
#ifdef EOS_TILLOTSON
            return 1;
#endif
            break;

        case IO_EOSABAR:
#ifdef EOS_CARRIES_ABAR
            return 1;
#endif
            break;

        case IO_EOSYE:
#ifdef EOS_CARRIES_YE
            return 1;
#endif
            break;

        case IO_CBE_MOMENTS:
            break;

        case IO_PARTVEL:
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            return 1;
#endif
            break;

        case IO_EDDINGTON_TENSOR:
#if defined(RADTRANSFER)
            return 1;
#endif
            break;

        case IO_AGS_SOFT:
#if defined (AGS_HSML_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTGRAVSOFT)
            return 1;
#endif
            break;

        case IO_AGS_RHO:
        case IO_AGS_QPT:
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
            return 1;
#endif
            break;

        case IO_AGS_PSI_RE:
        case IO_AGS_PSI_IM:
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY) && (DM_FUZZY > 0)
            return 1;
#endif
            break;

        case IO_AGS_ZETA:
#if defined(AGS_HSML_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTZETA)
            return 1;
#endif
            break;

        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            return 1;
#endif
            break;

        case IO_grH2I:
        case IO_grH2II:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            return 1;
#endif
            break;

        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            return 1;
#endif
            break;

        case IO_TURB_DYNAMIC_COEFF:
        case IO_TURB_DIFF_COEFF:
#ifdef TURB_DIFF_DYNAMIC
            return 1;
#endif
            break;

        case IO_DYNERRORDEFAULT:
        case IO_DYNERROR:
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
            return 1;
#endif
            break;

        case IO_LASTENTRY: /* will not occur */
            break;
    }

    return 0;			/* default: not present */
}




/*! This function associates a short 4-character block name with each block number.
 *  This is stored in front of each block for snapshot FileFormat=2.
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
    switch (blocknr)
    {
        case IO_POS:
            strncpy(label, "POS ", 4);
            break;
        case IO_VEL:
            strncpy(label, "VEL ", 4);
            break;
        case IO_CHILD_ID:
            strncpy(label, "IDch", 4);
            break;
        case IO_GENERATION_ID:
            strncpy(label, "IDgn", 4);
            break;
        case IO_ID:
            strncpy(label, "ID  ", 4);
            break;
        case IO_MASS:
            strncpy(label, "MASS", 4);
            break;
        case IO_U:
            strncpy(label, "U   ", 4);
            break;
        case IO_RHO:
            strncpy(label, "RHO ", 4);
            break;
        case IO_NE:
            strncpy(label, "NE  ", 4);
            break;
        case IO_NH:
            strncpy(label, "NH  ", 4);
            break;
        case IO_HII:
            strncpy(label, "HII ", 4);
            break;
        case IO_HeI:
            strncpy(label, "HeI ", 4);
            break;
        case IO_HeII:
            strncpy(label, "HeII", 4);
            break;
        case IO_INIB:
            strncpy(label, "INIB", 4);
            break;
        case IO_IDEN:
            strncpy(label, "IDEN", 4);
            break; 
        case IO_UNSPMASS:
            strncpy(label, "USPM", 4);
            break;     
        case IO_CRATE:
            strncpy(label, "CRATE", 4);
            break;
        case IO_HRATE:
            strncpy(label, "HRATE", 4);
            break;
        case IO_NHRATE:
            strncpy(label, "NHRATE", 4);
            break;
        case IO_HHRATE:
            strncpy(label, "HHRATE", 4);
            break;
        case IO_MCRATE:
            strncpy(label, "MCRATE", 4);
            break;
        case IO_HSML:
            strncpy(label, "HSML", 4);
            break;
        case IO_SFR:
            strncpy(label, "SFR ", 4);
            break;
        case IO_AGE:
            strncpy(label, "AGE ", 4);
            break;
        case IO_GRAINSIZE:
            strncpy(label, "GRSZ", 4);
            break;
        case IO_GRAINTYPE:
            strncpy(label, "GRTP", 4);
            break;
        case IO_DELAYTIME:
            strncpy(label, "DETI", 4);
            break;
        case IO_HSMS:
            strncpy(label, "HSMS", 4);
            break;
        case IO_Z:
            strncpy(label, "Z   ", 4);
            break;
        case IO_CHIMES_ABUNDANCES:
            strncpy(label, "CHIM", 4);
            break;
        case IO_CHIMES_MU:
            strncpy(label, "CHMU", 4);
            break;
        case IO_CHIMES_REDUCED:
            strncpy(label, "REDU", 4);
            break;
        case IO_CHIMES_NH:
            strncpy(label, "CHNH", 4);
            break;
        case IO_CHIMES_STAR_SIGMA:
            strncpy(label, "CHST", 4);
            break;
        case IO_CHIMES_FLUX_G0:
            strncpy(label, "CHGO", 4);
            break;
        case IO_CHIMES_FLUX_ION:
            strncpy(label, "CHIO", 4);
            break;
        case IO_DENS_AROUND_STAR:
            strncpy(label, "DNST", 4);
            break;
        case IO_DELAY_TIME_HII:
            strncpy(label, "DHII", 4);
            break;
        case IO_MOLECULARFRACTION:
            strncpy(label, "FMOL", 4);
            break;
        case IO_POT:
            strncpy(label, "POT ", 4);
            break;
        case IO_ACCEL:
            strncpy(label, "ACCE", 4);
            break;
        case IO_DTENTR:
            strncpy(label, "ENDT", 4);
            break;
        case IO_TSTP:
            strncpy(label, "TSTP", 4);
            break;
        case IO_BFLD:
            strncpy(label, "BFLD", 4);
            break;
        case IO_VDIV:
            strncpy(label, "VDIV", 4);
            break;
        case IO_VORT:
            strncpy(label, "VORT", 4);
            break;
        case IO_IMF:
            strncpy(label, "IMF ", 4);
            break;
        case IO_OSTAR:
            strncpy(label, "IMF ", 4);
            break;
        case IO_COSMICRAY_ENERGY:
            strncpy(label, "CREG ", 4);
            break;
        case IO_COSMICRAY_SLOPES:
            strncpy(label, "CRSL ", 4);
            break;
        case IO_COSMICRAY_KAPPA:
            strncpy(label, "CRK ", 4);
            break;
        case IO_COSMICRAY_ALFVEN:
            strncpy(label, "CRAV ", 4);
            break;
        case IO_DIVB:
            strncpy(label, "DIVB", 4);
            break;
        case IO_ABVC:
            strncpy(label, "ABVC", 4);
            break;
        case IO_AMDC:
            strncpy(label, "AMDC", 4);
            break;
        case IO_PHI:
            strncpy(label, "PHI ", 4);
            break;
        case IO_GRADPHI:
            strncpy(label, "GPHI", 4);
            break;
        case IO_COOLRATE:
            strncpy(label, "COOR", 4);
            break;
        case IO_BHMASS:
            strncpy(label, "BHMA", 4);
            break;
        case IO_BHDUSTMASS:
            strncpy(label, "BHDM", 4);
            break;
        case IO_BH_DIST:
            strncpy(label, "BHR ", 4);
            break;
        case IO_BHMASSALPHA:
            strncpy(label, "BHMa", 4);
            break;
        case IO_BH_ANGMOM:
            strncpy(label, "BHmJ", 4);
            break;
        case IO_ACRB:
            strncpy(label, "ACRB", 4);
            break;
        case IO_SINKRAD:
            strncpy(label, "SRAD", 4);
            break;
        case IO_BHMDOT:
            strncpy(label, "BHMD", 4);
            break;
        case IO_R_PROTOSTAR:
            strncpy(label, "RPST", 4);
            break;
        case IO_STAGE_PROTOSTAR:
            strncpy(label, "PSST", 4);
            break;
        case IO_AGE_PROTOSTAR:
            strncpy(label, "PSAG", 4);
            break;
        case IO_MASS_D_PROTOSTAR:
            strncpy(label, "PSMD", 4);
            break;
        case IO_ZAMS_MASS:
            strncpy(label, "PSMS", 4);
            break;
        case IO_LUM_SINGLESTAR:
            strncpy(label, "LUMS", 4);
            break;
        case IO_BHPROGS:
            strncpy(label, "BHPC", 4);
            break;
        case IO_TIDALTENSORPS:
            strncpy(label, "TIPS", 4);
            break;
        case IO_GDE_DISTORTIONTENSOR:
            strncpy(label, "DIPS", 4);
            break;
        case IO_CAUSTIC_COUNTER:
            strncpy(label, "CACO", 4);
            break;
        case IO_FLOW_DETERMINANT:
            strncpy(label, "FLDE", 4);
            break;
        case IO_STREAM_DENSITY:
            strncpy(label, "STDE", 4);
            break;
        case IO_PHASE_SPACE_DETERMINANT:
            strncpy(label, "PSDE", 4);
            break;
        case IO_ANNIHILATION_RADIATION:
            strncpy(label, "ANRA", 4);
            break;
        case IO_LAST_CAUSTIC:
            strncpy(label, "LACA", 4);
            break;
        case IO_SHEET_ORIENTATION:
            strncpy(label, "SHOR", 4);
            break;
        case IO_INIT_DENSITY:
            strncpy(label, "INDE", 4);
            break;
        case IO_EOSTEMP:
            strncpy(label, "TEMP", 4);
            break;
        case IO_EOSABAR:
            strncpy(label, "ABAR", 4);
            break;
        case IO_EOSCS:
            strncpy(label, "EQCS", 4);
            break;
        case IO_EOS_STRESS_TENSOR:
            strncpy(label, "ESTT", 4);
            break;
        case IO_CBE_MOMENTS:
            strncpy(label, "VMOM", 4);
            break;
        case IO_EOSCOMP:
            strncpy(label, "COMP", 4);
            break;
        case IO_PARTVEL:
            strncpy(label, "PVEL", 4);
            break;
        case IO_EOSYE:
            strncpy(label, "YE  ", 4);
            break;
        case IO_PRESSURE:
            strncpy(label, "P   ", 4);
            break;
        case IO_RADGAMMA:
            strncpy(label, "RADG", 4);
            break;
        case IO_RAD_FLUX:
            strncpy(label, "RADF", 4);
            break;
        case IO_RAD_ACCEL:
            strncpy(label, "RADA", 4);
            break;
        case IO_EDDINGTON_TENSOR:
            strncpy(label, "ET", 4);
            break;
        case IO_AGS_SOFT:
            strncpy(label, "AGSH", 4);
            break;
        case IO_AGS_RHO:
            strncpy(label, "ARHO", 4);
            break;
        case IO_AGS_QPT:
            strncpy(label, "AQPT", 4);
            break;
        case IO_AGS_PSI_RE:
            strncpy(label, "PSIR", 4);
            break;
        case IO_AGS_PSI_IM:
            strncpy(label, "PSII", 4);
            break;
        case IO_AGS_ZETA:
            strncpy(label, "AGSZ", 4);
            break;
        case IO_VSTURB_DISS:
            strncpy(label, "VSDI", 4);
            break;
        case IO_VSTURB_DRIVE:
            strncpy(label, "VSDR", 4);
            break;
        case IO_grHI:
            strncpy(label, "gHI", 4);
            break;
        case IO_grHII:
            strncpy(label, "gHII", 4);
            break;
        case IO_grHM:
            strncpy(label, "gHM", 4);
            break;
        case IO_grHeI:
            strncpy(label, "gHeI", 4);
            break;
        case IO_grHeII:
            strncpy(label, "gHe2", 4);
            break;
        case IO_grHeIII:
            strncpy(label, "gHe3", 4);
            break;
        case IO_grH2I:
            strncpy(label, "gH2I", 4);
            break;
        case IO_grH2II:
            strncpy(label, "H2II", 4);
            break;
        case IO_grDI:
            strncpy(label, "gDI", 4);
            break;
        case IO_grDII:
            strncpy(label, "gDII", 4);
            break;
        case IO_grHDI:
            strncpy(label, "gHDI", 4);
            break;
        case IO_TURB_DYNAMIC_COEFF:
            strncpy(label, "tdyn", 4);
            break;
        case IO_TURB_DIFF_COEFF:
            strncpy(label, "turb", 4);
            break;
        case IO_DYNERROR:
            strncpy(label, "derr", 4);
            break;
        case IO_DYNERRORDEFAULT:
            strncpy(label, "derd", 4);
            break;

        case IO_LASTENTRY:
            endrun(217);
            break;
    }
}


void get_dataset_name(enum iofields blocknr, char *buf)
{
    strcpy(buf, "default");

    switch (blocknr)
    {
        case IO_POS:
            strcpy(buf, "Coordinates");
            break;
        case IO_VEL:
            strcpy(buf, "Velocities");
            break;
        case IO_ID:
            strcpy(buf, "ParticleIDs");
            break;
        case IO_CHILD_ID:
            strcpy(buf, "ParticleChildIDsNumber");
            break;
        case IO_GENERATION_ID:
            strcpy(buf, "ParticleIDGenerationNumber");
            break;
        case IO_MASS:
            strcpy(buf, "Masses");
            break;
        case IO_U:
            strcpy(buf, "InternalEnergy");
            break;
        case IO_RHO:
            strcpy(buf, "Density");
            break;
        case IO_NE:
            strcpy(buf, "ElectronAbundance");
            break;
        case IO_NH:
            strcpy(buf, "NeutralHydrogenAbundance");
            break;
        case IO_RADGAMMA:
            strcpy(buf, "PhotonEnergy");
            break;
        case IO_RAD_FLUX:
            strcpy(buf, "PhotonFluxDensity");
            break;
        case IO_RAD_ACCEL:
            strcpy(buf, "RadiativeAcceleration");
            break;
        case IO_HII:
            strcpy(buf, "HII");
            break;
        case IO_HeI:
            strcpy(buf, "HeI");
            break;
        case IO_HeII:
            strcpy(buf, "HeII");
            break;
        case IO_IDEN:
            strcpy(buf, "IniDen");
            break;
        case IO_INIB:
            strcpy(buf, "IniB");
            break;    
        case IO_UNSPMASS:
            strcpy(buf, "Unspawned_Wind_Mass");
            break;     
        case IO_CRATE:
            strcpy(buf, "CoolingRate");
            break;
        case IO_HRATE:
            strcpy(buf, "HeatingRate");
            break;
        case IO_NHRATE:
            strcpy(buf, "NetHeatingRateQ");
            break;
        case IO_HHRATE:
            strcpy(buf, "HydroHeatingRate");
            break;
        case IO_MCRATE:
            strcpy(buf, "MetalCoolingRate");
            break;
        case IO_DELAYTIME:
            strcpy(buf, "DelayTime");
            break;
        case IO_HSML:
            strcpy(buf, "SmoothingLength");
            break;
        case IO_SFR:
            strcpy(buf, "StarFormationRate");
            break;
        case IO_AGE:
            strcpy(buf, "StellarFormationTime");
            break;
        case IO_GRAINSIZE:
            strcpy(buf, "GrainSize");
            break;
        case IO_GRAINTYPE:
            strcpy(buf, "PICParticleType");
            break;
        case IO_HSMS:
            strcpy(buf, "StellarSmoothingLength");
            break;
        case IO_Z:
            strcpy(buf, "Metallicity");
            break;
        case IO_CHIMES_ABUNDANCES:
            strcpy(buf, "ChimesAbundances");
            break;
        case IO_CHIMES_MU:
            strcpy(buf, "ChimesMu");
            break;
        case IO_CHIMES_REDUCED:
            strcpy(buf, "ChimesReducedAbundances");
            break;
        case IO_CHIMES_NH:
            strcpy(buf, "ChimesColumnDensity");
            break;
        case IO_CHIMES_STAR_SIGMA:
            strcpy(buf, "SigmaEff");
            break;
        case IO_CHIMES_FLUX_G0:
            strcpy(buf, "ChimesFluxG0");
            break;
        case IO_CHIMES_FLUX_ION:
            strcpy(buf, "ChimesFluxIon");
            break;
        case IO_DENS_AROUND_STAR:
            strcpy(buf, "DensityAtParticleLocation");
            break;
        case IO_DELAY_TIME_HII:
            strcpy(buf, "DelayTime_HIIRegion_Cooling");
            break;
        case IO_MOLECULARFRACTION:
            strcpy(buf, "MolecularMassFraction");
            break;
        case IO_POT:
            strcpy(buf, "Potential");
            break;
        case IO_ACCEL:
            strcpy(buf, "Acceleration");
            break;
        case IO_DTENTR:
            strcpy(buf, "RateOfChangeOfInternalEnergy");
            break;
        case IO_TSTP:
            strcpy(buf, "TimeStep");
            break;
        case IO_BFLD:
            strcpy(buf, "MagneticField");
            break;
        case IO_VDIV:
            strcpy(buf, "VelocityDivergence");
            break;
        case IO_VORT:
            strcpy(buf, "Vorticity");
            break;
        case IO_IMF:
            strcpy(buf, "IMFFormationProperties");
            break;
        case IO_OSTAR:
            strcpy(buf, "OStarNumber");
            break;
        case IO_COSMICRAY_ENERGY:
            strcpy(buf, "CosmicRayEnergy");
            break;
        case IO_COSMICRAY_SLOPES:
            strcpy(buf, "CosmicRayMomentumDistSlope");
            break;
        case IO_COSMICRAY_KAPPA:
            strcpy(buf, "CosmicRayDiffusivity");
            break;
        case IO_COSMICRAY_ALFVEN:
            strcpy(buf, "CosmicRayAlfvenEnergyPM");
            break;
        case IO_DIVB:
            strcpy(buf, "DivergenceOfMagneticField");
            break;
        case IO_ABVC:
            strcpy(buf, "ArtificialViscosity");
            break;
        case IO_AMDC:
            strcpy(buf, "ArtMagneticDissipation");
            break;
        case IO_PHI:
            strcpy(buf, "DivBcleaningFunctionPhi");
            break;
        case IO_GRADPHI:
            strcpy(buf, "DivBcleaningFunctionGradPhi");
            break;
        case IO_COOLRATE:
            strcpy(buf, "CoolingRate");
            break;
        case IO_BHMASS:
            strcpy(buf, "BH_Mass");
            break;
        case IO_BHDUSTMASS:
            strcpy(buf, "BH_Dust_Mass");
            break;
        case IO_BH_DIST:
            strcpy(buf, "BH_Dist");
            break;
        case IO_BHMASSALPHA:
            strcpy(buf, "BH_Mass_AlphaDisk");
            break;
        case IO_BH_ANGMOM:
            strcpy(buf, "BH_Specific_AngMom");
            break;
        case IO_ACRB:
            strcpy(buf, "BH_AccretionLength");
            break;
        case IO_SINKRAD:
            strcpy(buf, "SinkRadius");
            break;
        case IO_BHMDOT:
            strcpy(buf, "BH_Mdot");
            break;
        case IO_R_PROTOSTAR:
            strcpy(buf, "ProtoStellarRadius_inSolar");
            break;
        case IO_MASS_D_PROTOSTAR:
            strcpy(buf, "Mass_D");
            break;
        case IO_ZAMS_MASS:
            strcpy(buf, "ZAMS_Mass");
            break;
        case IO_STAGE_PROTOSTAR:
            strcpy(buf, "ProtoStellarStage");
            break;
        case IO_AGE_PROTOSTAR:
            strcpy(buf, "ProtoStellarAge");
            break;
        case IO_LUM_SINGLESTAR:
            strcpy(buf, "StarLuminosity_Solar");
            break;
        case IO_BHPROGS:
            strcpy(buf, "BH_NProgs");
            break;
        case IO_TIDALTENSORPS:
            strcpy(buf, "TidalTensor");
            break;
        case IO_GDE_DISTORTIONTENSOR:
            strcpy(buf, "DistortionTensorPS");
            break;
        case IO_CAUSTIC_COUNTER:
            strcpy(buf, "CausticCounter");
            break;
        case IO_FLOW_DETERMINANT:
            strcpy(buf, "FlowDeterminant");
            break;
        case IO_STREAM_DENSITY:
            strcpy(buf, "StreamDensity");
            break;
        case IO_PHASE_SPACE_DETERMINANT:
            strcpy(buf, "PhaseSpaceDensity");
            break;
        case IO_ANNIHILATION_RADIATION:
            strcpy(buf, "AnnihilationRadiation");
            break;
        case IO_LAST_CAUSTIC:
            strcpy(buf, "LastCaustic");
            break;
        case IO_SHEET_ORIENTATION:
            strcpy(buf, "SheetOrientation");
            break;
        case IO_INIT_DENSITY:
            strcpy(buf, "InitDensity");
            break;
        case IO_EOSTEMP:
            strcpy(buf, "Temperature");
            break;
        case IO_EOSABAR:
            strcpy(buf, "Abar");
            break;
        case IO_EOSCS:
            strcpy(buf, "SoundSpeed");
            break;
        case IO_EOS_STRESS_TENSOR:
            strcpy(buf, "StressTensor");
            break;
        case IO_CBE_MOMENTS:
            strcpy(buf, "VlasovMoments");
            break;
        case IO_EOSCOMP:
            strcpy(buf, "CompositionType");
            break;
        case IO_PARTVEL:
            strcpy(buf, "ParticleVelocities");
            break;
        case IO_EOSYE:
            strcpy(buf, "Ye");
            break;
        case IO_PRESSURE:
            strcpy(buf, "Pressure");
            break;
        case IO_EDDINGTON_TENSOR:
            strcpy(buf, "EddingtonTensor");
            break;
        case IO_AGS_SOFT:
            strcpy(buf, "AGS-Softening");
            break;
        case IO_AGS_RHO:
            strcpy(buf, "AGS-Density");
            break;
        case IO_AGS_QPT:
            strcpy(buf, "AGS-QuantumPotentialQ");
            break;
        case IO_AGS_PSI_RE:
            strcpy(buf, "WavefunctionPsi-Real");
            break;
        case IO_AGS_PSI_IM:
            strcpy(buf, "WavefunctionPsi-Imag");
            break;
        case IO_AGS_ZETA:
            strcpy(buf, "AGS-Zeta");
            break;
        case IO_VSTURB_DISS:
            strcpy(buf, "TurbulenceDissipation");
            break;
        case IO_VSTURB_DRIVE:
            strcpy(buf, "TurbulenceDriving");
            break;
        case IO_grHI:
            strcpy(buf, "GrackleHI");
            break;
        case IO_grHII:
            strcpy(buf, "GrackleHII");
            break;
        case IO_grHM:
            strcpy(buf, "GrackleHM");
            break;
        case IO_grHeI:
            strcpy(buf, "GrackleHeI");
            break;
        case IO_grHeII:
            strcpy(buf, "GrackleHeII");
            break;
        case IO_grHeIII:
            strcpy(buf, "GrackleHeIII");
            break;
        case IO_grH2I:
            strcpy(buf, "GrackleH2I");
            break;
        case IO_grH2II:
            strcpy(buf, "GrackleH2II");
            break;
        case IO_grDI:
            strcpy(buf, "GrackleDI");
            break;
        case IO_grDII:
            strcpy(buf, "GrackleDII");
            break;
        case IO_grHDI:
            strcpy(buf, "GrackleHDI");
            break;
        case IO_TURB_DYNAMIC_COEFF:
            strcpy(buf, "DynSmagCoeff");
            break;
        case IO_TURB_DIFF_COEFF:
            strcpy(buf, "TurbDiffCoeff");
            break;
        case IO_DYNERROR:
            strcpy(buf, "DynamicError");
            break;
        case IO_DYNERRORDEFAULT:
            strcpy(buf, "DynamicErrorDefault");
            break;
        case IO_LASTENTRY:
            endrun(218);
            break;
    }
}



/*! This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the gas kernel
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 */
void write_file(char *fname, int writeTask, int lastTask)
{
    int type, bytes_per_blockelement, npart, nextblock, typelist[6];
    int n_for_this_task, n, p, pc, offset = 0, task, i;
    size_t blockmaxlen;
    int ntot_type[6], nn[6];
    enum iofields blocknr;
    char label[8];
    int bnr;
    int blksize;
    MPI_Status status;
    FILE *fd = 0;

#ifdef HAVE_HDF5
    hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
    hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
    herr_t hdf5_status;
    hsize_t dims[2], count[2], start[2];
    int rank = 0, pcsum = 0;
    char buf[500];
#endif

#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

    /* determine particle numbers of each type in file */

    if(ThisTask == writeTask)
    {
        for(n = 0; n < 6; n++)
            ntot_type[n] = n_type[n];

        for(task = writeTask + 1; task <= lastTask; task++)
        {
            MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
            for(n = 0; n < 6; n++)
                ntot_type[n] += nn[n];
        }

        for(task = writeTask + 1; task <= lastTask; task++)
            MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
        MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }

    /* fill file header */
    for(n = 0; n < 6; n++)
    {
        header.npart[n] = (int) ntot_type[n];
        header.npartTotal[n] = (unsigned int) ntot_type_all[n];
        header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }

    if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS) {header.flag_ic_info = FLAG_EVOLVED_2LPT;}
    if(header.flag_ic_info == FLAG_ZELDOVICH_ICS) {header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;}
    if(header.flag_ic_info == FLAG_NORMALICS_2LPT) {header.flag_ic_info = FLAG_EVOLVED_2LPT;}
    if(header.flag_ic_info == 0 && All.ComovingIntegrationOn != 0) {header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;}

    for(n = 0; n < 6; n++) {header.mass[n] = All.MassTable[n];}

    header.time = All.Time;
    if(All.ComovingIntegrationOn) {header.redshift = 1.0 / All.Time - 1;} else {header.redshift = 0;}

    header.flag_sfr = 0;
    header.flag_feedback = 0;
    header.flag_cooling = 0;
    header.flag_stellarage = 0;
    header.flag_metals = 0;

#ifdef COOLING
    header.flag_cooling = 1;
#endif

#ifdef GALSF
    header.flag_sfr = 1;
    header.flag_feedback = 1;
    header.flag_stellarage = 1;
#endif

#ifdef METALS
    header.flag_metals = NUM_METAL_SPECIES;
#endif

    header.num_files = All.NumFilesPerSnapshot;
    header.BoxSize = All.BoxSize;
    header.OmegaMatter = All.OmegaMatter;
    header.OmegaLambda = All.OmegaLambda;
    header.HubbleParam = All.HubbleParam;

#ifdef OUTPUT_IN_DOUBLEPRECISION
    header.flag_doubleprecision = 1;
#else
    header.flag_doubleprecision = 0;
#endif

    /* open file and write header */

    if(ThisTask == writeTask)
    {
        if(All.SnapFormat == 3)
        {
#ifdef HAVE_HDF5
            sprintf(buf, "%s.hdf5", fname);
            hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

            hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

            for(type = 0; type < 6; type++)
            {
                if(header.npart[type] > 0)
                {
                    sprintf(buf, "/PartType%d", type);
                    hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
                }
            }

            write_header_attributes_in_hdf5(hdf5_headergrp);
#endif
        }
        else
        {
            if(!(fd = fopen(fname, "w")))
            {
                printf("can't open file `%s' for writing snapshot.\n", fname);
                endrun(123);
            }

            if(All.SnapFormat == 2)
            {
                blksize = sizeof(int) + 4 * sizeof(char);
                SKIP;
                my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
                nextblock = sizeof(header) + 2 * sizeof(int);
                my_fwrite(&nextblock, sizeof(int), 1, fd);
                SKIP;
            }

            blksize = sizeof(header);
            SKIP;
            my_fwrite(&header, sizeof(header), 1, fd);
            SKIP;
        }
    }

    if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
        n_info = 0;
        InfoBlock = (struct info_block *) mymalloc("InfoBlock", sizeof(struct info_block) * 1000);

        for(bnr = 0; bnr < 1000; bnr++)
        {
            blocknr = (enum iofields) bnr;

            if(blocknr == IO_LASTENTRY)
                break;

            if(blockpresent(blocknr))
            {
                bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);
                npart = get_particles_in_block(blocknr, &typelist[0]);

                if(npart > 0)
                {
                    for(type = 0; type < 6; type++) {InfoBlock[n_info].is_present[type] = typelist[type];}
                    InfoBlock[n_info].ndim = get_values_per_blockelement(blocknr);
                    get_Tab_IO_Label(blocknr, label);
                    for(i = 0; i < 4; i++) {InfoBlock[n_info].label[i] = label[i];}
                    switch (get_datatype_in_block(blocknr))
                    {
                        case 0:
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "LONG    ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "LONGN   ", 8);
                            break;
                        case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
#else
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
#endif
                            break;
                        case 2:
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "LLONG   ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "LLONGN  ", 8);
                            break;
                        case 3:
#ifdef OUTPUT_POSITIONS_IN_DOUBLE
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
#else
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
#endif
                            break;
                    }
                    n_info++;
                }
            }
        }
    }

    for(bnr = 0; bnr < 1000; bnr++)
    {
        blocknr = (enum iofields) bnr;
        if(blocknr == IO_LASTENTRY) {break;}

        if(blockpresent(blocknr))
        {
            bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);
            size_t MyBufferSize = All.BufferSize;
            blockmaxlen = (size_t) ((MyBufferSize * 1024 * 1024) / bytes_per_blockelement);
            npart = get_particles_in_block(blocknr, &typelist[0]);

            if(npart > 0)
            {
                if(ThisTask == 0)
                {
                    char buf[1000];
                    get_dataset_name(blocknr, buf);
                    printf("writing block %d (%s)...\n", bnr, buf);
                }

                if(ThisTask == writeTask)
                {

                    if(All.SnapFormat == 1 || All.SnapFormat == 2)
                    {
                        if(All.SnapFormat == 2)
                        {
                            blksize = sizeof(int) + 4 * sizeof(char);
                            SKIP;
                            get_Tab_IO_Label(blocknr, label);
                            my_fwrite(label, sizeof(char), 4, fd);
                            nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
                            my_fwrite(&nextblock, sizeof(int), 1, fd);
                            SKIP;
                        }
                        blksize = npart * bytes_per_blockelement;
                        SKIP;
                    }
                }

                for(type = 0; type < 6; type++)
                {
                    if(typelist[type])
                    {
#ifdef HAVE_HDF5
                        if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
                        {
                            switch(get_datatype_in_block(blocknr))
                            {
                                case 0:
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
                                    break;
                                    
                                case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                    break;
                                    
                                case 2:
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
                                    break;
                                    
                                case 3:
#ifdef OUTPUT_POSITIONS_IN_DOUBLE
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                    break;
                            }

                            dims[0] = header.npart[type];
                            dims[1] = get_values_per_blockelement(blocknr);
                            if(dims[1] == 1) {rank = 1;} else {rank = 2;}

                            get_dataset_name(blocknr, buf);
                            hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
#ifndef IO_COMPRESS_HDF5
                            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
#else
                            if(dims[0] > 10)
                            {
                            	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
                            	hsize_t cdims[2]; cdims[0] = (hsize_t) (dims[0] / 10); cdims[1] = dims[1];
                            	hdf5_status = H5Pset_chunk (plist_id, rank, cdims);
                            	hdf5_status = H5Pset_deflate (plist_id, 4);
                            	hdf5_dataset = H5Dcreate2(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                            } else {
                            	hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
                            }
#endif
                            pcsum = 0;
                        }
#endif

                        for(task = writeTask, offset = 0; task <= lastTask; task++)
                        {
                            if(task == ThisTask)
                            {
                                n_for_this_task = n_type[type];

                                for(p = writeTask; p <= lastTask; p++)
                                    {if(p != ThisTask) {MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);}}
                            }
                            else
                                MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

                            while(n_for_this_task > 0)
                            {
                                pc = n_for_this_task;
                                if(pc > (int)blockmaxlen) {pc = blockmaxlen;}
                                if(ThisTask == task) {fill_write_buffer(blocknr, &offset, pc, type);}

                                if(ThisTask == writeTask && task != writeTask)
                                    {MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD, &status);}

                                if(ThisTask != writeTask && task == ThisTask)
                                    {MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, MPI_COMM_WORLD);}

                                if(ThisTask == writeTask)
                                {
                                    if(All.SnapFormat == 3)
                                    {
#ifdef HAVE_HDF5
                                        start[0] = pcsum;
                                        start[1] = 0;

                                        count[0] = pc;
                                        count[1] = get_values_per_blockelement(blocknr);
                                        pcsum += pc;

                                        H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
                                                            start, NULL, count, NULL);

                                        dims[0] = pc;
                                        dims[1] = get_values_per_blockelement(blocknr);
                                        if(dims[1] == 1) {rank = 1;} else {rank = 2;}
                                        hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

                                        hdf5_status = H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

                                        H5Sclose(hdf5_dataspace_memory);
#endif
                                    }
                                    else
                                    {
                                        my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                                    }
                                }

                                n_for_this_task -= pc;
                            }
                        }

#ifdef HAVE_HDF5
                        if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
                        {
                            if(All.SnapFormat == 3)
                            {
                                H5Dclose(hdf5_dataset);
                                H5Sclose(hdf5_dataspace_in_file);
                                H5Tclose(hdf5_datatype);
                            }
                        }
#endif
                    }
                }

                if(ThisTask == writeTask)
                {
                    if(All.SnapFormat == 1 || All.SnapFormat == 2)
                    {
                        SKIP;
                    }
                }
            }
        }
    }

    if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
        myfree(InfoBlock);
    }

    if(ThisTask == writeTask)
    {
        if(All.SnapFormat == 3)
        {
#ifdef HAVE_HDF5
            for(type = 5; type >= 0; type--) {if(header.npart[type] > 0) {H5Gclose(hdf5_grp[type]);}}
            H5Gclose(hdf5_headergrp);
            H5Fclose(hdf5_file);
#endif
        }
        else
        {
            fclose(fd);
        }
    }
}




#ifdef HAVE_HDF5
void write_header_attributes_in_hdf5(hid_t handle)
{
    hsize_t adim[1] = { 6 }; hid_t hdf5_dataspace, hdf5_attribute;

    {int tmp=GIZMO_VERSION; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "GIZMO_version", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &tmp); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE); H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, header.npart); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SIMPLE); H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SIMPLE); H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SIMPLE); H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "ComovingIntegrationOn", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.ComovingIntegrationOn); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    if(All.ComovingIntegrationOn)
    {
        hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Omega_Matter", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.OmegaMatter); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

        hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Omega_Lambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.OmegaLambda); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

        hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Omega_Baryon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.OmegaBaryon); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

        hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Omega_Radiation", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.OmegaRadiation); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    }

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    {double tmp=UNIT_MASS_IN_CGS; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "UnitMass_In_CGS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &tmp); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    {double tmp=UNIT_VEL_IN_CGS; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "UnitVelocity_In_CGS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &tmp); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    {double tmp=UNIT_LENGTH_IN_CGS; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "UnitLength_In_CGS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &tmp); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#ifdef MAGNETIC
    {double tmp=UNIT_B_IN_GAUSS; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Internal_UnitB_In_Gauss", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &tmp); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#endif
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Gravitational_Constant_In_Code_Inits", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.G); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Minimum_Mass_For_Cell_Merge", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MinMassForParticleMerger); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Maximum_Mass_For_Cell_Split", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MaxMassForParticleSplit); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SIMPLE); H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "Fixed_ForceSoftening_Keplerian_Kernel_Extent", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.ForceSoftening); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_IC_Info", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

    {int ivar=KERNEL_FUNCTION; hdf5_dataspace=H5Screate(H5S_SCALAR); hdf5_attribute=H5Acreate(handle,"Kernel_Function_ID",H5T_NATIVE_INT,hdf5_dataspace,H5P_DEFAULT);
        H5Awrite(hdf5_attribute,H5T_NATIVE_INT,&ivar); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    hdf5_dataspace=H5Screate(H5S_SCALAR); hdf5_attribute=H5Acreate(handle,"Effective_Kernel_NeighborNumber",H5T_NATIVE_DOUBLE,hdf5_dataspace,H5P_DEFAULT);
    H5Awrite(hdf5_attribute,H5T_NATIVE_DOUBLE,&All.DesNumNgb); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);

#ifdef SUBFIND
    hdf5_dataspace=H5Screate(H5S_SCALAR); hdf5_attribute=H5Acreate(handle,"Subfind_FOFLink_NeighborNumber",H5T_NATIVE_INT,hdf5_dataspace,H5P_DEFAULT);
    H5Awrite(hdf5_attribute,H5T_NATIVE_INT,&All.DesLinkNgb); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(COSMIC_RAY_FLUID) && (CRFLUID_DIFFUSION_MODEL == 0)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "CosmicRayDiffusionCoeff_at_GV_CodeUnits", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.CosmicRayDiffusionCoeff); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(DM_SIDM)
#ifdef GRAIN_COLLISIONS
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_InteractionRenormalization", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_InteractionCrossSection); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_DissipationFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_DissipationFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_KickPerCollision", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_KickPerCollision); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_InteractionVelocityScale", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_InteractionVelocityScale); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#else
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "DM_InteractionCrossSection", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_InteractionCrossSection); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "DM_DissipationFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_DissipationFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "DM_KickPerCollision", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_KickPerCollision); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "DM_InteractionVelocityScale", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DM_InteractionVelocityScale); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#endif

#ifdef DM_SCALARFIELD_SCREENING
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "DM_ScalarBeta", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ScalarBeta); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "DM_ScalarScreeningLength", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ScalarScreeningLength); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif


#ifdef RT_EVOLVE_INTENSITIES
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[2]={N_RT_INTENSITY_BINS,3}; H5Sset_extent_simple(hdf5_dataspace, 2, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "Radiation_Intensity_Grid_Direction_Vectors", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.Rad_Intensity_Direction); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#endif

#ifdef CHIMES
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Flag_CHIMES_Thermal_Evolution_On", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.ChimesThermEvolOn); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef RT_LEBRON
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "LEBRON_PhotonMomentum_Coupled_Fraction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.PhotonMomentum_Coupled_Fraction); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef BH_PHOTONMOMENTUM
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "LEBRON_Sink_PhotonMomentum_Coupled_Fraction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BH_Rad_MomentumFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef GRAIN_FLUID
    {int holder=GRAIN_PTYPES; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "GrainCR_Particle_Type_BitFlag", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#if !defined(PIC_MHD) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_Internal_Density", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Grain_Internal_Density); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_Size_Min", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Grain_Size_Min); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_Size_Max", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Grain_Size_Max); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_Size_Spectrum_Powerlaw", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Grain_Size_Spectrum_Powerlaw); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS) && defined(RT_GENERIC_USER_FREQ)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_Absorbed_vs_Total_Extinction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Grain_Absorbed_Fraction_vs_Total_Extinction); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef GRAIN_RDI_TESTPROBLEM
#ifdef RT_OPACITY_FROM_EXPLICIT_GRAINS
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_Q_at_MaxGrainSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Grain_Q_at_MaxGrainSize); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Grain_Charge_Parameter", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Grain_Charge_Parameter); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Dust_to_Gas_Mass_Ratio", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Dust_to_Gas_Mass_Ratio); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Vertical_Gravity_Strength", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Vertical_Gravity_Strength); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Vertical_Grain_Accel", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Vertical_Grain_Accel); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Vertical_Grain_Accel_Angle", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Vertical_Grain_Accel_Angle); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#ifdef BOX_SHEARING
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Pressure_Gradient_Accel", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Pressure_Gradient_Accel); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#endif
#endif

#ifdef PIC_MHD
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "PIC_Charge_to_Mass_Ratio", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.PIC_Charge_to_Mass_Ratio); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef GALSF
    {double holder=DMAX(All.PhysDensThresh,All.OverDensThresh); hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Density_Threshold_For_SF_CodeUnits", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#if !defined(SINGLE_STAR_SINK_DYNAMICS)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SF_Timescale_At_Density_Threshold_CodeUnits", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MaxSfrTimescale); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#endif

#ifdef GALSF_EFFECTIVE_EQS
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_EgySpecSN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.EgySpecSN); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_FactorSN", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.FactorSN); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_EgySpecCold", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.EgySpecCold); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_FactorEVP", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.FactorEVP); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_FeedbackEnergy", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.FeedbackEnergy); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_TempSupernova", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TempSupernova); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_TempClouds", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TempClouds); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_ISM_EOS_Parameter_FactorForSofterEQS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.FactorForSofterEQS); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif


#ifdef GALSF_SUBGRID_WINDS
    {int holder=GALSF_SUBGRID_WIND_SCALING; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SubGrid_Wind_Model_Scaling_Key", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SubGrid_Wind_Model_WindEfficiency", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.WindEfficiency); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SubGrid_Wind_Model_WindEnergyFraction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.WindEnergyFraction); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SubGrid_Wind_Model_WindFreeTravelMaxTimeFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.WindFreeTravelMaxTimeFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SubGrid_Wind_Model_WindFreeTravelDensFac", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.WindFreeTravelDensFac); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#if (GALSF_SUBGRID_WIND_SCALING>0)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SubGrid_Wind_Model_VariableWindVelFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.VariableWindVelFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SubGrid_Wind_Model_VariableWindSpecMomentum", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.VariableWindSpecMomentum); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#endif



#ifdef GALSF_FB_FIRE_AGE_TRACERS
    {int holder=NUM_AGE_TRACERS; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "AgeTracer_NumberOfBins", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "AgeTracerEventsPerTimeBin", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.AgeTracerRateNormalization); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={NUM_AGE_TRACERS+1}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "AgeTracer_CustomTimeBins", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, All.AgeTracerTimeBins); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#else
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "AgeTracerBinStart", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.AgeTracerBinStart); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "AgeTracerBinEnd", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.AgeTracerBinEnd); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#endif

#ifdef METALS
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={NUM_METAL_SPECIES}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "Solar_Abundances_Adopted", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.SolarAbundances); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    /* assign labels for all metal species for reference in outputs */
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={NUM_METAL_SPECIES}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
        hdf5_attribute = H5Acreate(handle, "Metals_Atomic_Number_Or_Key", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
        int zkey[NUM_METAL_SPECIES],k; zkey[0]=0; /* all metals */ for(k=1;k<NUM_METAL_SPECIES;k++) {zkey[k]=-20;}
        if(NUM_LIVE_SPECIES_FOR_COOLTABLES==10) {zkey[1]=2; zkey[2]=6; zkey[3]=7; zkey[4]=8; zkey[5]=10; zkey[6]=12; zkey[7]=14; zkey[8]=16; zkey[9]=20; zkey[10]=26;} /* He,C,N,O,Ne,Mg,Si,S,Ca,Fe */
        for(k=0;k<NUM_RPROCESS_SPECIES;k++) {zkey[1+NUM_LIVE_SPECIES_FOR_COOLTABLES+k]=-1;}
        for(k=0;k<NUM_AGE_TRACERS;k++) {zkey[1+NUM_LIVE_SPECIES_FOR_COOLTABLES+NUM_RPROCESS_SPECIES+k]=-2;}
        for(k=0;k<NUM_STARFORGE_FEEDBACK_TRACERS;k++) {zkey[1+NUM_LIVE_SPECIES_FOR_COOLTABLES+NUM_RPROCESS_SPECIES+NUM_AGE_TRACERS+k]=-3;}
        H5Awrite(hdf5_attribute, H5T_NATIVE_INT, zkey); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#endif

#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
    {
        int k; double numin[N_RT_FREQ_BINS], numax[N_RT_FREQ_BINS]; for(k=0;k<N_RT_FREQ_BINS;k++) {numin[k]=-20; numax[k]=-20;}
#ifdef RT_CHEM_PHOTOION
#if defined(RT_PHOTOION_MULTIFREQUENCY)
        int i_vec[4] = {RT_FREQ_BIN_H0, RT_FREQ_BIN_He0, RT_FREQ_BIN_He1, RT_FREQ_BIN_He2};
        numin[i_vec[3]]=rt_ion_nu_min[i_vec[3]]; numax[i_vec[3]]=500; for(k=0;k<3;k++) {numin[i_vec[k]]=rt_ion_nu_min[i_vec[k]]; numax[i_vec[k]]=rt_ion_nu_min[i_vec[k+1]];}
#else
        k=RT_FREQ_BIN_H0; numin[k]=13.6; numax[k]=500;
#endif
#endif
#ifdef RT_SOFT_XRAY
        k=RT_FREQ_BIN_SOFT_XRAY; numin[k]=500; numax[k]=2000;
#endif
#ifdef RT_HARD_XRAY
        k=RT_FREQ_BIN_HARD_XRAY; numin[k]=2000; numax[k]=10000;
#endif
#ifdef RT_PHOTOELECTRIC
        k=RT_FREQ_BIN_PHOTOELECTRIC; numin[k]=8; numax[k]=13.6;
#endif
#ifdef RT_LYMAN_WERNER
        k=RT_FREQ_BIN_LYMAN_WERNER; numin[k]=11.2; numax[k]=13.6;
#endif
#ifdef RT_NUV
        k=RT_FREQ_BIN_NUV; numin[k]=3.444; numax[k]=8.;
#endif
#ifdef RT_OPTICAL_NIR
        k=RT_FREQ_BIN_OPTICAL_NIR; numin[k]=0.4133; numax[k]=3.444;
#endif
#ifdef RT_GENERIC_USER_FREQ
        k=RT_FREQ_BIN_GENERIC_USER_FREQ; numin[k]=-1; numax[k]=-1;
#endif
#ifdef RT_INFRARED
        k=RT_FREQ_BIN_INFRARED; numin[k]=0.001; numax[k]=0.4133;
#endif
#ifdef RT_FREEFREE
        k=RT_FREQ_BIN_FREEFREE; numin[k]=-2; numax[k]=-2;
#endif
        
        {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_RT_FREQ_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
            hdf5_attribute = H5Acreate(handle, "Radiation_RHD_Min_Bin_Freq_in_eV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
            H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, numin); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
        {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_RT_FREQ_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
            hdf5_attribute = H5Acreate(handle, "Radiation_RHD_Max_Bin_Freq_in_eV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
            H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, numax); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    }
#endif

    
#if defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK) || defined(BH_WIND_SPAWN)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BAL_f_accretion", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BAL_f_accretion); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BAL_v_outflow", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BAL_v_outflow); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(SINGLE_STAR_FB_JETS)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BAL_f_launch_v", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BAL_f_launch_v); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(BH_COSMIC_RAYS)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BH_CosmicRay_Injection_Efficiency", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BH_CosmicRay_Injection_Efficiency); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef GR_TABULATED_COSMOLOGY
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "DarkEnergyConstantW", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.DarkEnergyConstantW); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef TURB_DIFFUSION
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDiffusion_Coefficient", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDiffusion_Coefficient); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#ifdef TURB_DIFF_DYNAMIC
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDynamicDiffFac", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDynamicDiffFac); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDynamicDiffSmoothing", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDynamicDiffSmoothing); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#endif

#if defined(CONDUCTION) && !defined(CONDUCTION_SPITZER)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Thermal_ConductionCoeff", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ConductionCoeff); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(VISCOSITY) && !defined(VISCOSITY_BRAGINSKII)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "ShearViscosityCoeff", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ShearViscosityCoeff); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BulkViscosityCoeff", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BulkViscosityCoeff); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Effective_Kernel_NeighborNumber_CollisionlessParticles", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.AGS_DesNumNgb); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef DM_FUZZY
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "ScalarField_hbar_over_mass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ScalarField_hbar_over_mass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef TURB_DRIVING
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_DecayTime", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDriving_Global_DecayTime); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_AccelerationPowerVariable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDriving_Global_AccelerationPowerVariable); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_DtTurbUpdates", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDriving_Global_DtTurbUpdates); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_DrivingScaleKMinVar", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDriving_Global_DrivingScaleKMinVar); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_DrivingScaleKMaxVar", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDriving_Global_DrivingScaleKMaxVar); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_SolenoidalFraction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.TurbDriving_Global_SolenoidalFraction); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_DrivingSpectrumKey", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.TurbDriving_Global_DrivingSpectrumKey); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "TurbDriving_Global_DrivingRandomNumberKey", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &All.TurbDriving_Global_DrivingRandomNumberKey); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[2]={7,12}; H5Sset_extent_simple(hdf5_dataspace, 2, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "Tillotson_EOS_parameters_for_all_species", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.Tillotson_EOS_params); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#endif

#ifdef BLACK_HOLES
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BlackHoleAccretionFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BlackHoleAccretionFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BlackHoleFeedbackFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BlackHoleFeedbackFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedBlackHoleMass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedBlackHoleMass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BlackHoleNgbFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BlackHoleNgbFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BlackHoleMaxAccretionRadius", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BlackHoleMaxAccretionRadius); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BlackHoleEddingtonFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BlackHoleEddingtonFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BlackHoleRadiativeEfficiency", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BlackHoleRadiativeEfficiency); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#if defined(BH_SEED_FROM_FOF) || defined(BH_SEED_FROM_LOCALGAS)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedBlackHoleMassSigma", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedBlackHoleMassSigma); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedBlackHoleMinRedshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedBlackHoleMinRedshift); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef BH_SEED_FROM_LOCALGAS
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedBlackHolePerUnitMass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedBlackHolePerUnitMass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef BH_ALPHADISK_ACCRETION
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedAlphaDiskMass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedAlphaDiskMass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef BH_SEED_FROM_FOF
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "MinFoFMassForNewSeed", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MinFoFMassForNewSeed); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef BH_WIND_SPAWN
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BAL_wind_particle_mass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BAL_wind_particle_mass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "BAL_internal_temperature", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.BAL_internal_temperature); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    {unsigned long long holder = (unsigned long long) All.AGNWindID; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Spawned_Cell_ID", H5T_NATIVE_ULLONG, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_ULLONG, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#endif
#endif

}
#endif





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
    size_t nwritten;

    if(size * nmemb > 0)
    {
        if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
        {
            printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
            fflush(stdout);
            endrun(777);
        }
    }
    else {nwritten = 0;}

    return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
    size_t nread;

    if(size * nmemb == 0) {return 0;}

    if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
        if(feof(stream)) {printf("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);}
            else {printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));}
        fflush(stdout);
        endrun(778);
    }
    return nread;
}

/** This is so we don't have to preface every printf call with the if statement to only make it print once. */
void mpi_printf(const char *fmt, ...)
{
    if(ThisTask == 0)
    {
        va_list l;
        va_start(l, fmt);
        vprintf(fmt, l);
        va_end(l);
    }
}

#if defined(IO_SUBFIND_READFOF_FROMIC)
int io_compare_P_ID(const void *a, const void *b)
{
    if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID)) {return -1;}
    if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID)) {return +1;}
    return 0;
}

int io_compare_P_GrNr_SubNr(const void *a, const void *b)
{
    if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr)) {return -1;}
    if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr)) {return +1;}
    if(((struct particle_data *) a)->SubNr < (((struct particle_data *) b)->SubNr)) {return -1;}
    if(((struct particle_data *) a)->SubNr > (((struct particle_data *) b)->SubNr)) {return +1;}
    return 0;
}

int io_compare_P_GrNr_ID(const void *a, const void *b)
{
    if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr)) {return -1;}
    if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr)) {return +1;}
    if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID)) {return -1;}
    if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID)) {return +1;}
    return 0;
}

#endif
