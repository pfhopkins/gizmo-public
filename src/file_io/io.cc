#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "../core/timestep_functions.h"   /* the Hermite predicted state written by IO_HERMITE_SYNC */
#if defined(GALSF_ISMDUSTCHEM_MODEL)
/* get_ISMDustChemEvo_bin_mass is KOKKOS_INLINE_FUNCTION (Chunk 1, 2026-05-26). */
#include "../solids/ism_dust_chemistry_functions.h"
#endif

/*! \file io.c
 *  \brief Output of a snapshot file to disk.
 */
/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified heavily
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO. Mostly to
 * write out new/modified quantities, but also changing the units of
 * snapshots, adding additional basic functionality for types of input/output,
 * changing read conditions and dynamic dependence on run-time and compile-time
 * variables, changing flexibility for certain dataset layouts and architectures,
 * adding compression options, adding extensive header options and
 * output of various parameterfile and run-time options for later user reference,
 * and rewriting some subroutines to be compatible with modern libraries.
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
    char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE], outputdir[DEFAULT_PATH_BUFFERSIZE_TOUSE];
    int n, filenr, gr, ngroups, primaryTask, lastTask;

    CPU_Step[CPU_MISC] += measure_time();

#ifdef CHIMES_REDUCED_OUTPUT
    if (num % N_chimes_full_output_freq == 0) {Chimes_incl_full_output = 1;} else {Chimes_incl_full_output = 0;}
#endif

    rearrange_particle_sequence();
    gizmo_exit_bad_stop_if_requested("io:after_rearrange"); /* a refused rearrange must not be followed by a write */
    All.NumForcesSinceLastTreeBuild = (long long) (1 + All.TreeRebuild_ActiveFraction * All.TotNumPart);    /* ensures that new tree will be constructed */



    if(DumpFlag == 1)
    {
        if(ThisTask == 0)
            printf("\nwriting snapshot file #%d... \n", num);

        size_t comm_chunk_mb = All.CommChunkSize;
        bytes = comm_chunk_mb * 1024 * 1024;
        /* Collective OOM preflight (mirrors allocate_memory): mymalloc() does
         * NOT return NULL on capacity failure -- it routes to the TEMP-hard
         * path -- so check fit BEFORE the call and bad-stop gracefully. All-rank
         * here (inside the collective DumpFlag==1 snapshot-write block). */
        if(!gizmo_alloc_fits_all_ranks(gizmo_mymalloc_rounded_size(bytes), 1))
        {
            gizmo_request_controlled_stop(2, "savepositions: CommBuffer arena allocation won't fit on >=1 rank", __FILE__, __LINE__, __FUNCTION__);
            return;
        }
        if(!(CommBuffer = mymalloc("CommBuffer", bytes)))
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

        snprintf(outputdir, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s", All.OutputDir);
#if IO_SINKS_ONLY_SNAPSHOT_FREQUENCY > 0
        if( (All.SnapshotFileCount-1) % (IO_SINKS_ONLY_SNAPSHOT_FREQUENCY+1) )
        {
            snprintf(outputdir, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/stars_only", All.OutputDir);
            mkdir(outputdir, 02755);
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
                snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/snapdir_%03d", outputdir, num);
                mkdir(buf, 02755);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        if(All.NumFilesPerSnapshot > 1)
            snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/snapdir_%03d/%s_%03d.%d", outputdir, num, All.SnapshotFileBase, num, filenr);
        else
            snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/%s_%03d", outputdir, All.SnapshotFileBase, num);


        ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
        if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel)) {ngroups++;}

        for(gr = 0; gr < ngroups; gr++)
        {
            if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
            {
                write_file(buf, primaryTask, lastTask);
            }
            /* All-rank per-turn boundary: a write_file fault (failed open / I/O error) set a soft bad-stop on the
             * writeTask; this poll (Allreduce = barrier-equivalent sync) collects it and finalizes cleanly before
             * the next turn, in place of the closing MPI_Barrier. */
            gizmo_exit_bad_stop_if_requested("io:write_file_turn");
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

#if defined(OUTPUT_POWERSPEC) && defined(PMGRID) && defined(BOX_PERIODIC)
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



/*! This function fills the write buffer with particle data. New output blocks can in principle be added here.
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
#ifdef MAGNETIC /* NOTE: we always work -internally- in code units where MU_0 = 1; hence the 4pi here; [much simpler, but be sure of your conversions!] */
    double gizmo2gauss = UNIT_B_IN_GAUSS / All.UnitMagneticField_in_gauss;
#endif
#ifdef PMGRID
    double dt_gravkick_pm = get_gravkick_factor((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2, All.Ti_Current, -1, 0);
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
                    Vec3<double> pos_reported = gizmo_reported_position(P[pindex].Pos);
                    for(k = 0; k < 3; k++)
                    {
                        fp_pos[k] = (MyOutputPosFloat) pos_reported[k];
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
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) (P[pindex].Vel[k] * sqrt(All.cf_a3inv));} // JUST write the conserved velocity here, not the drifted one in this manner //
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

        case IO_REFINE_FLAG:	/* IC tag defining the nuclear-zoom refinement anchor */
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_TAG_ANCHOR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip++ = (MyIDType) P[pindex].Refinement_Flag;
                    n++;
                }
#endif
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
                    *fp++ = (MyOutputFloat) DMAX(All.MinEgySpec, CellP[pindex].InternalEnergyPred);
                    n++;
                }
            break;

        case IO_RHO:		/* density */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Density;
                    n++;
                }
            break;

        case IO_NE:		/* electron abundance */
#if (defined(COOLING) || defined(RT_CHEM_PHOTOION)) && !defined(CHIMES)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Ne;
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
                    *fp++ = (MyOutputFloat) CellP[pindex].HI;
#elif (COOL_GRACKLE_CHEMISTRY > 0)
                    *fp++ = (MyOutputFloat) CellP[pindex].grHI;
#else
                    *fp++ = (MyOutputFloat) CellP[pindex].HI; /* saved by the cooling solver, and read back by read_ic */
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
                    *fp++ = (MyOutputFloat) CellP[pindex].HII;
                    n++;
                }
#endif
            break;

        case IO_HeI:		/* neutral Helium */
#if defined(RT_CHEM_PHOTOION_HE)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].HeI;
                    n++;
                }
#endif
            break;

        case IO_HeII:		/* ionized Helium */
#if defined(RT_CHEM_PHOTOION_HE)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].HeII;
                    n++;
                }
#endif
            break;
                        
        case IO_UNSPMASS:
#if defined(SINK_WIND_SPAWN)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].unspawned_wind_mass;
                    n++;
                }
#endif           
            break;

        case IO_UNSPJETMASS:
#if defined(SINK_WIND_SPAWN) && defined(SINGLE_STAR_FB_JETS)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].unspawned_jet_mass;
                    n++;
                }
#endif
            break;
            
        case IO_CRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].CoolingRate;
                    n++;
                }
#endif
            break;

        case IO_HRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].HeatingRate;
                    n++;
                }
#endif
            break;

        case IO_NHRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].NetHeatingRateQ;
                    n++;
                }
#endif
            break;

        case IO_HHRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].HydroHeatingRate;
                    n++;
                }
#endif
            break;

        case IO_MCRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].MetalCoolingRate;
                    n++;
                }
#endif
            break;

        case IO_DCRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING) && defined(GALSF_ISMDUSTCHEM_MODEL) && defined(GALSF_ISMDUSTCHEM_HIGHTEMPDUSTCOOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].DustCoolingRate;
                    n++;
                }
#endif
            break;

        case IO_PHRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].PElecHeatingRate;
                    n++;
                }
#endif
            break;

        case IO_KERNELRADIUS:		/* gas kernel length */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].KernelRadius;
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
            
        case IO_DTOSTAR:
#ifdef GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].TimeDistribOfStarFormation;
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

        case IO_DUST_TO_GAS:        /* dust to gas mass ratio */
#ifdef OUTPUT_DUST_TO_GAS_RATIO
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Metallicity[0]*return_dust_to_metals_ratio_vs_solar(pindex,0, P, CellP);
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

        case IO_FLUIDTYPE:      /* Lagrangian partition ID */
#ifdef HYDRO_MULTIFLUID
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) P[pindex].FluidType;
                    n++;
                }
#endif
            break;

        case IO_VSTURB_DISS:
#if defined(TURB_DRIVING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].DuDt_diss;
                    n++;
                }
#endif
            break;

        case IO_VSTURB_DRIVE:
#if defined(TURB_DRIVING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].DuDt_drive;
                    n++;
                }
#endif
            break;

        case IO_HSMS:		/* kernel length for star particles */
#ifdef SUBFIND
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].DM_KernelRadius;
                    n++;
                }
#endif
            break;

        case IO_Z:			/* gas and star metallicity */
#ifdef METALS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<NUM_METAL_SPECIES;k++)
                    {
                        MyOutputFloat val = (MyOutputFloat) P[pindex].Metallicity[k];
#ifndef OUTPUT_IN_DOUBLEPRECISION
                        /* Quantize to ~3.3 significant decimal digits (10-bit mantissa) by zeroing
                           the low 13 mantissa bits of the float32. This is lossless at the precision
                           level meaningful for metal abundances, and makes the data ~2x more
                           compressible with shuffle+deflate (IO_COMPRESS_HDF5). */
                        {unsigned int *ibits = (unsigned int *)&val; *ibits &= 0xFFFFE000u;}
#endif
                        fp[k] = val;
                    }
                    fp += NUM_METAL_SPECIES;
                    n++;
                }
#endif
            break;

        case IO_DUSTCHEMZMET:         /* gas dust metallicity and sources of dust formation */
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < NUM_ISMDUSTCHEM_ELEMENTS; k++) {fp[k] = (MyOutputFloat) CellP[pindex].ISMDustChem_Dust_Metal[k];}
                    double source_total = 0.; // Convert dust source mass fractions to fraction of dust
                    for(k = 0; k < NUM_ISMDUSTCHEM_SOURCES; k++) source_total += CellP[pindex].ISMDustChem_Dust_Source[k];
                    for(k = 0; k < NUM_ISMDUSTCHEM_SOURCES; k++) 
                    {
                        if(source_total<=0 || CellP[pindex].ISMDustChem_Dust_Metal[0]<=0) {fp[k+NUM_ISMDUSTCHEM_ELEMENTS] = (MyOutputFloat) 0.;}
                        else {fp[k+NUM_ISMDUSTCHEM_ELEMENTS] = (MyOutputFloat) CellP[pindex].ISMDustChem_Dust_Source[k]/source_total;}
                    }
                    fp += NUM_ISMDUSTCHEM_ELEMENTS + NUM_ISMDUSTCHEM_SOURCES;
                    n++;
                }
#endif
        break;
            
        case IO_DUSTCHEMSPECIESMET:    /* gas dust species following Species routines */
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < NUM_ISMDUSTCHEM_SPECIES; k++) {fp[k] = (MyOutputFloat) CellP[pindex].ISMDustChem_Dust_Species[k];}
                    fp += NUM_ISMDUSTCHEM_SPECIES;
                    n++;
                }
#endif
        break;

        case IO_ISMDUSTCHEMMOL:    /* sub-resolved molecular gas properties (fraction of C in CO and fraction of gas that is in dense molecular phase) used in dust routines wiuthout grain size evolution */
#if defined(GALSF_ISMDUSTCHEM_MODEL) && !defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    fp[0] = (MyOutputFloat) CellP[pindex].ISMDustChem_MassFractionInDenseMolecular;
                    fp[1] = (MyOutputFloat) CellP[pindex].ISMDustChem_C_in_CO;
                    fp += 2;
                    n++;
                }
#endif
        break;

        case IO_MACHNUM:    /* sub-resolved molecular gas properties (fraction of C in CO and fraction of gas that is in dense molecular phase) used in dust routines */
#if defined(OUTPUT_MACH_NUMBER)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].ISMDustChem_MachNumber;
                    n++;
                }
#endif
        break;

        case IO_SHOCKMACHNUM:    /* local shock Mach number from pairwise Riemann problem */
#if defined(OUTPUT_SHOCK_MACH_NUMBER)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].ShockMachNumber;
                    n++;
                }
#endif
        break;

        case IO_DUSTCHEMGRAINBINNUMBERS:    /* number of grains for each grain size bin for each dust species */
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    int k1, k2;
                    for(k1=0;k1<NUM_ISMDUSTCHEM_SPECIES;k1++) {
                        for(k2=0;k2<NUM_ISMDUSTCHEM_SIZE_BINS;k2++) {
                            fp[NUM_ISMDUSTCHEM_SIZE_BINS*k1 + k2] = (MyOutputFloat) (CellP[pindex].ISMDustChem_Dust_NumberInBin[k1][k2]/UNIT_GRAIN_NUMBER);
                        }
                    }
                    fp += (NUM_ISMDUSTCHEM_SPECIES*NUM_ISMDUSTCHEM_SIZE_BINS);
                    n++;
                }
#endif
        break;
        
        case IO_DUSTCHEMGRAINBINMASS:    /* mass for each grain size bin for each dust species */
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    int k1, k2;
                    for(k1=0;k1<NUM_ISMDUSTCHEM_SPECIES;k1++) {
                        for(k2=0;k2<NUM_ISMDUSTCHEM_SIZE_BINS;k2++) {     
                            fp[NUM_ISMDUSTCHEM_SIZE_BINS*k1 + k2] = (MyOutputFloat) (get_ISMDustChemEvo_bin_mass(pindex, k1, k2, CellP) / UNIT_MASS_IN_CGS);
                        }
                    }
                    fp += (NUM_ISMDUSTCHEM_SPECIES*NUM_ISMDUSTCHEM_SIZE_BINS);
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
                    *fp++ = (MyOutputFloat) (evaluate_NH_from_GradRho(P[pindex].GradRho,P[pindex].KernelRadius,CellP[pindex].Density,P[pindex].NumNgb,1,pindex) * UNIT_SURFDEN_IN_CGS * shielding_length_factor * (1.0 - (P[pindex].Metallicity[0] + P[pindex].Metallicity[1])) / PROTONMASS_CGS);
                    n++;
                }
#endif
            break;

        case IO_CHIMES_STAR_SIGMA:
#if defined(CHIMES_NH_OUTPUT) && defined(OUTPUT_DENS_AROUND_NONGAS)
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (evaluate_NH_from_GradRho(P[pindex].GradRho,P[pindex].KernelRadius,P[pindex].DensityAroundParticle,P[pindex].NumNgb,0,pindex) * UNIT_SURFDEN_IN_CGS);  // g cm^-2
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
                    if(CellP[pindex].DelayTimeHII > 0) {for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) (CellP[pindex].Chimes_G0[k] + CellP[pindex].Chimes_G0_HII[k]);}}
                        else {for(k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) CellP[pindex].Chimes_G0[k];}}
#else
                    for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) CellP[pindex].Chimes_G0[k]; }
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
                    if(CellP[pindex].DelayTimeHII > 0) {for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) (CellP[pindex].Chimes_fluxPhotIon[k] + CellP[pindex].Chimes_fluxPhotIon_HII[k]);}}
                        else {for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) CellP[pindex].Chimes_fluxPhotIon[k];}}
#else
                    for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {fp[k] = (MyOutputFloat) CellP[pindex].Chimes_fluxPhotIon[k];}
#endif
                    fp += CHIMES_LOCAL_UV_NBINS;
                    n++;
                }
#endif
            break;

        case IO_DENS_AROUND_STAR:
#ifdef OUTPUT_DENS_AROUND_NONGAS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].DensityAroundParticle;
                    n++;
                }
#endif
            break;

        case IO_DELAY_TIME_HII:
#if defined(GALSF_FB_FIRE_RT_HIIHEATING) && defined(OUTPUT_DELAY_TIME_HII)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].DelayTimeHII;
                    n++;
                }
#endif
            break;

        case IO_MOLECULARFRACTION:
#if defined(OUTPUT_MOLECULAR_FRACTION)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#if defined(COOL_MOLECFRAC_NONEQM)
                    *fp++ = (MyOutputFloat) CellP[pindex].MolecularMassFraction_perNeutralH; /* more useful to output this particular value, rather than fH2 */
#else
                    *fp++ = (MyOutputFloat) CellP[pindex].MolecularMassFraction; /* as the cooling solver last left it */
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

        case IO_HERMITE_POS:
        case IO_HERMITE_VEL:
#if defined(IO_HERMITE_SYNC) && defined(HERMITE_INTEGRATION)
            /* A mutually consistent position-velocity pair at the output time.
             *
             * The ordinary Coordinates/Velocities blocks are not such a pair: positions are
             * drifted to the output time while velocities are whatever the last kick left
             * (see the IO_VEL comment). Vis-viva, orbital elements and kinetic energy are all
             * meaningless on that mixture -- in the e=0.9 binary test it makes |a/a0 - 1| swing
             * by three orders of magnitude within a single orbit.
             *
             * For a particle the Hermite integrator owns, the start-of-step state is retained,
             * so the consistent state at any time within the step is the same series the
             * integrator predicts with. It is the predicted state, not the corrected one -- the
             * corrector needs an end-of-step acceleration that does not exist at an arbitrary
             * output time -- making it 4th order in position and 3rd in velocity, far above the
             * error being diagnosed, and unlike the mixed state it is a state the system
             * actually passes through. Anything else falls back to the ordinary values, so the
             * datasets are always well defined. */
            {
                struct HermiteWalkState hermite_out = hermite_walk_state_snapshot();
                struct DriftKickTableView hermite_out_tables = drift_kick_table_view_host();
                for(n = 0; n < pc; pindex++)
                    if(P[pindex].Type == type)
                    {
                        Vec3<double> r_out = P[pindex].Pos, v_out = P[pindex].Vel;
                        /* the eligibility predicate, not the type bitmask alone: a star being
                           super-timestepped, or one that just accreted, is advanced by something
                           other than this series, so its retained base does not describe where it
                           goes and the plain values are the honest answer for it */
                        if((P[pindex].Mass > 0) && eligible_for_hermite(pindex, P)) {
                            hermite_predict_source_state(pindex, P, hermite_out, &hermite_out_tables, r_out, v_out);
                        }
                        /* written exactly as IO_POS and IO_VEL write theirs: positions through the
                           reported-frame un-shift, wrapped into the box, at position precision;
                           velocities scaled to the conserved variable. A coordinate output that
                           skipped the un-shift would sit in a different frame from Coordinates
                           whenever the gravity solver has randomized the box. */
                        if(blocknr == IO_HERMITE_POS)
                        {
                            Vec3<double> r_reported = gizmo_reported_position(r_out);
                            for(k = 0; k < 3; k++)
                            {
                                fp_pos[k] = (MyOutputPosFloat) r_reported[k];
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
                        }
                        else
                        {
                            for(k = 0; k < 3; k++) {fp[k] = (MyOutputFloat) (v_out[k] * sqrt(All.cf_a3inv));}
                            fp += 3;
                        }
                        n++;
                    }
            }
#endif
            break;

        case IO_SINK_DIST:
#if defined(SINK_CALC_DISTANCES) && defined(OUTPUT_SINK_DISTANCES)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Min_Distance_to_Sink;
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
#ifndef OUTPUT_HYDROACCELERATION
                    if(P[pindex].Type == 0) {for(k = 0; k < 3; k++) {fp[k] += (MyOutputFloat) CellP[pindex].HydroAccel[k];}}
#endif
                    fp += 3;
                    n++;
                }
#endif
            break;

            
        case IO_HYDROACCEL:        /* 'hydro' acceleration */
#ifdef OUTPUT_HYDROACCELERATION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    if(P[pindex].Type == 0) {for(k = 0; k < 3; k++) {fp[k] = (MyOutputFloat) CellP[pindex].HydroAccel[k];}}
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
                    *fp++ = (MyOutputFloat) CellP[pindex].DtInternalEnergy;
                    n++;
                }
#endif
            break;

        case IO_DELAYTIME:
#ifdef GALSF_SUBGRID_WINDS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].DelayTime;
                    n++;
                }
#endif
            break;

        case IO_TSTP:		/* timestep  */
#ifdef OUTPUT_TIMESTEP
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (get_particle_timestep_in_physical(pindex, P));
                    n++;
                }
#endif
            break;

        case IO_BFLD:		/* magnetic field  */
#ifdef MAGNETIC
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    {Vec3<double> Btmp = CellP[pindex].Bfield() * All.cf_a2inv * gizmo2gauss; for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) Btmp[k];}}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_AMBIPOLAR:   /* Ambipolar resistivity */
#if defined(MHD_NON_IDEAL) && defined(OUTPUT_MHD_RESISTIVITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Eta_MHD_AmbiPolarDiffusion_Coeff;
                    n++;
                }
#endif            	    
            break;

        case IO_OHMIC:   /* Ohmic resistivity */
#if defined(MHD_NON_IDEAL) && defined(OUTPUT_MHD_RESISTIVITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Eta_MHD_OhmicResistivity_Coeff;
                    n++;
                }
#endif	    
            break;

        case IO_HALL:   /* Hall resistivity */
#if defined(MHD_NON_IDEAL) && defined(OUTPUT_MHD_RESISTIVITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Eta_MHD_HallEffect_Coeff;
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
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) CellP[pindex].Vorticity[k];}
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
                    for(k=0;k<N_CR_PARTICLE_BINS;k++) {fp[k] = (MyOutputFloat) CellP[pindex].CosmicRayEnergyPred[k];}
                    fp += N_CR_PARTICLE_BINS;
                    n++;
                }
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (CellP[pindex].SubGrid_CosmicRayEnergyDensity * P[pindex].Mass / CellP[pindex].Density);
                    n++;
                }
#endif
            break;

            
        case IO_COSMICRAY_SLOPES:    /* piecewise spectral slopes in the cosmic ray spectrum  */
#if defined(COSMIC_RAY_FLUID) && defined(CRFLUID_EVOLVE_SPECTRUM)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    //for(k=0;k<N_CR_PARTICLE_BINS;k++) {fp[k] = CellP[pindex].CosmicRay_PwrLaw_Slopes_in_Bin[k];} // write out the saved slope
                    for(k=0;k<N_CR_PARTICLE_BINS;k++) {fp[k] = (MyOutputFloat) CR_return_spectral_slope_target(pindex,k, CellP);} // calculate the slope to write out
                    fp += N_CR_PARTICLE_BINS;
                    n++;
                }
#endif
            break;

        case IO_COSMICRAY_KAPPA:    /* local CR diffusion constant */
#if defined(COSMIC_RAY_FLUID) && defined(CRFLUID_DIFFUSION_MODEL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<N_CR_PARTICLE_BINS;k++) {fp[k] = (MyOutputFloat) CellP[pindex].CosmicRayDiffusionCoeff[k];}
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
                    for(k=0;k<2;k++) {int k2; for(k2=0;k2<N_CR_PARTICLE_BINS;k2++) {fp[N_CR_PARTICLE_BINS*k + k2] = (MyOutputFloat) CellP[pindex].CosmicRayAlfvenEnergyPred[k2][k];}}
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
                    *fp++ = (MyOutputFloat) ((CellP[pindex].divB * gizmo2gauss * (CellP[pindex].Density*All.cf_a3inv / P[pindex].Mass)));
                    n++;
                }
#endif
            break;

        case IO_ABVC:		/* artificial viscosity of particle  */
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (CellP[pindex].alpha * CellP[pindex].alpha_limiter);
                    n++;
                }
#endif
            break;


        case IO_AMDC:		/* artificial magnetic dissipation of particle  */
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Balpha;
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
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) (CellP[pindex].Gradients.Phi[k] * All.cf_a2inv*All.cf_a2inv * gizmo2gauss);}
                    fp += 3;
                    n++;
                }
#endif
            break;

            
        case IO_GRADRHO:        /* density gradient */
#if defined(OUTPUT_GRADIENT_RHO)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) (CellP[pindex].Gradients.Density[k] * All.cf_a2inv*All.cf_a2inv);} // output in physical units
                    fp += 3;
                    n++;
                }
#endif
            break;

                                                               
        case IO_GRADVEL:        /* velocity gradients */
#if defined(OUTPUT_GRADIENT_VEL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    int k1,k2;
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {fp[k1*3 + k2] = (MyOutputFloat) (CellP[pindex].Gradients.Velocity[k1][k2] * All.cf_a2inv);}} // output in physical units
                    fp += 9;
                    n++;
                }
#endif
            break;

            
        case IO_GRADMAG:        /* magnetic field gradients */
#if defined(OUTPUT_GRADIENT_MAG)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    int k1,k2;
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {fp[k1*3 + k2] = (MyOutputFloat) (CellP[pindex].Gradients.B[k1][k2] * UNIT_B_IN_GAUSS * All.cf_a2inv/All.cf_atime);}} // output in physical units
                    fp += 9;
                    n++;
                }
#endif
            break;

            
        case IO_COOLRATE:		/* current cooling rate of particle  */
#ifdef OUTPUT_COOLRATE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    double ne = CellP[pindex].Ne, ne_out=ne; u = CellP[pindex].InternalEnergyPred; tcool = GetCoolingTime(u, CellP[pindex].Density * All.cf_a3inv, ne, &ne_out, pindex, P, CellP); /* get cooling time */
                    double coolrate_to_output = 0; if(tcool != 0) {coolrate_to_output = u / tcool;} /* convert cooling time with current thermal energy to du/dt */
                    *fp++ = (MyOutputFloat) coolrate_to_output;
                    n++;
                }
#endif
            break;

        case IO_SINKMASS:
#ifdef SINK_PARTICLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Sink_Mass;
                    n++;
                }
#endif
            break;

        case IO_SINKDUSTMASSACC:
#if defined(SINK_PARTICLES) && defined(GRAIN_FLUID)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Sink_Dust_Mass;
                    n++;
                }
#endif
            break;

        case IO_SINKMASSALPHA:
#ifdef SINK_ALPHADISK_ACCRETION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Sink_Mass_Reservoir;
                    n++;
                }
#endif
            break;

        case IO_SINK_ANGMOM:
#ifdef SINK_FOLLOW_ACCRETED_ANGMOM
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++) {fp[k] = (MyOutputFloat) P[pindex].Sink_Specific_AngMom[k];}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_SINKMDOT:
#ifdef SINK_PARTICLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Sink_Mdot;
                    n++;
                }
#endif
            break;

        case IO_R_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].ProtoStellarRadius_inSolar;
                    n++;
                }
#endif
            break;

        case IO_MASS_D_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Mass_D;
                    n++;
                }
#endif
            break;

        case IO_ZAMS_MASS:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].ZAMS_Mass;
                    n++;
                }
#endif
            break;

        case IO_STAGE_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) P[pindex].ProtoStellarStage;
                    n++;
                }
#endif
            break;
            
        case IO_AGE_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].ProtoStellarAge;
                    n++;
                }
#endif
            break;

        case IO_LUM_SINGLESTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].StarLuminosity_Solar;
                    n++;
                }
#endif
            break;

        case IO_SINKPROGS:
#ifdef SINK_COUNTPROGS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) P[pindex].Sink_CountProgs;
                    n++;
                }
#endif
            break;

        case IO_ACRB:
#ifdef SINK_PARTICLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].KernelRadius;
                    n++;
                }
#endif
            break;

        case IO_SINKRAD:
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].SinkRadius;
                    n++;
                }
#endif
            break;

        case IO_SINK_FORM_MASS:
#ifdef SINK_PARTICLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].Sink_Formation_Mass;
                    n++;
                }
#endif
            break;	    

        case IO_TIDALTENSORPS:   /* 3x3 configuration-space tidal tensor */
#ifdef OUTPUT_TIDAL_TENSOR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++) {int l_tt_tmp; for(l_tt_tmp = 0; l_tt_tmp < 3; l_tt_tmp++) {fp[k * 3 + l_tt_tmp] = (MyOutputFloat) P[pindex].tidal_tensorps[k][l_tt_tmp];}}
                    n++;
                    fp += 9;
                }
#endif
            break;


        case IO_EOSABAR:
#ifdef EOS_CARRIES_ABAR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Abar;
                    n++;
                }
#endif
            break;

        case IO_TURB_DYNAMIC_COEFF:
#ifdef TURB_DIFF_DYNAMIC
            for (n = 0; n < pc; pindex++) {
                if (P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].TD_DynDiffCoeff;
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
                    *fp++ = (MyOutputFloat) CellP[pindex].Ye;
                    n++;
                }
#endif
            break;

#ifdef NUCLEAR_NETWORK
        case IO_NUCLEAR_COMPOSITION:
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < NUM_NUCLEAR_SPECIES; k++) {*fp++ = (MyOutputFloat) P[pindex].Metallicity[NUCLEAR_SPECIES_OFFSET_IN_METALLICITY + k];}
                    n++;
                }
            break;
#endif

        case IO_EOSTEMP:
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Temperature;
                    n++;
                }
            break;

#ifdef TWO_TEMPERATURE_PLASMA
        case IO_TWOTEMP_TE:
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].T_e_cell;
                    n++;
                }
            break;
#endif

        case IO_PRESSURE:
#if defined(EOS_GENERAL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Pressure;
                    n++;
                }
#endif
            break;

            case IO_EOSCS:
#if defined(EOS_GENERAL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].SoundSpeed;
                    n++;
                }
#endif
            break;

        case IO_CBE_MOMENTS:
#if defined(CBE_INTEGRATOR)
            for(n = 0; n < pc; pindex++)
            {
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < CBE_INTEGRATOR_NBASIS; k++)
                    {
                        int kf;
                        for(kf = 0; kf < CBE_INTEGRATOR_NMOMENTS; kf++)
                        {
                            fp[CBE_INTEGRATOR_NMOMENTS*k + kf] = (MyOutputFloat) P[pindex].CBE_basis_moments[k][kf];
                        }
                    }
                    fp += (CBE_INTEGRATOR_NMOMENTS*CBE_INTEGRATOR_NBASIS);
                    n++;
                }
            }
#endif
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
                            fp[3*k + kf] = (MyOutputFloat) CellP[pindex].Elastic_Stress_Tensor[kf][k];
                    }
                    fp += 9;
                    n++;
                }
#endif
            break;

            case IO_EOSCOMP:
#if defined(EOS_TILLOTSON) || defined(EOS_ANEOS)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) CellP[pindex].CompositionType;
                    n++;
                }
#endif
            break;

        case IO_DAMAGE_POROSITY_DAMAGE:
#ifdef EOS_DAMAGE_POROSITY
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Damage;
                    n++;
                }
#endif
            break;

        case IO_DAMAGE_POROSITY_DISTENTION:
#ifdef EOS_DAMAGE_POROSITY
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Distention;
                    n++;
                }
#endif
            break;

        case IO_DAMAGE_POROSITY_ACTVCRACKS:
#ifdef EOS_DAMAGE_POROSITY
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].ActiveCracks;
                    n++;
                }
#endif
            break;

            case IO_EOSPHASE:
#ifdef EOS_ANEOS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = (int) CellP[pindex].PhaseID;
                    n++;
                }
#endif
            break;

        case IO_PARTVEL:
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++) {fp[k] = (MyOutputFloat) CellP[pindex].ParticleVel[k];}
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
                    for(k=0;k<N_RT_FREQ_BINS;k++) {fp[k] = (MyOutputFloat) CellP[pindex].Rad_E_gamma[k];}
                    fp += N_RT_FREQ_BINS;
                    n++;
                }
#endif
            break;

        case IO_RAD_OPACITY:
#if defined(OUTPUT_RT_RAD_OPACITY) && defined(RADTRANSFER)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<N_RT_FREQ_BINS;k++) {fp[k] = (MyOutputFloat) CellP[pindex].Rad_Kappa[k];}
                    fp += N_RT_FREQ_BINS;
                    n++;
                }
#endif
            break;
            
        case IO_RAD_TEMP:
#if defined(RADTRANSFER) && defined(RT_INFRARED)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Radiation_Temperature;
                    n++;
                }
#endif
            break;

        case IO_DUST_TEMP:
#if (defined(RADTRANSFER) && defined(RT_INFRARED)) || (defined(OUTPUT_DUST_TEMPERATURE) && GALSF_FB_FIRE_STELLAREVOLUTION > 2)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) CellP[pindex].Dust_Temperature;
                    n++;
                }
#endif
            break;
            
        case IO_RAD_FLUX:
#if defined(OUTPUT_RT_RAD_FLUX) && defined(RT_EVOLVE_FLUX)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<3;k++) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {fp[N_RT_FREQ_BINS*k + kf] = (MyOutputFloat) (CellP[pindex].Rad_Flux_Pred[kf][k] * (CellP[pindex].Density*All.cf_a3inv/P[pindex].Mass));}}
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
                    for(k=0;k<3;k++) {fp[k] = (MyOutputFloat) CellP[pindex].Rad_Accel[k];}
                    fp += 3;
                    n++;
                }
#endif
            break;

        case IO_EDDINGTON_TENSOR:
#if defined(OUTPUT_EDDINGTON_TENSOR)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k=0;k<6;k++) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {fp[N_RT_FREQ_BINS*k + kf] = (MyOutputFloat) CellP[pindex].ET[kf].data[k];}}
                    fp += 6*N_RT_FREQ_BINS;
                    n++;
                }
#endif
            break;

        case IO_SOFT:		/* Adaptive Gravitational Softening: softening */
#if defined(OUTPUT_SOFTENING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) (ForceSoftening_KernelRadius(pindex));
                    n++;
                }
#endif
            break;

        case IO_AGS_HKERN:        /* neighbor search radius for various particle types */
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].AGS_KernelRadius;
                    n++;
                }
#endif
            break;

        case IO_AGS_RHO:        /* Adaptive Gravitational Softening: density */
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#if defined(DM_FUZZY)
                    *fp++ = (MyOutputFloat) P[pindex].AGS_Density;
#else
                    /* code density = Mass / AGS effective volume (same estimator the
                     * CBE collision operator uses; get_particle_volume_ags) */
                    double V_ags_out = get_particle_volume_ags(pindex);
                    *fp++ = (MyOutputFloat) ((V_ags_out > 0) ? (P[pindex].Mass / V_ags_out) : 0);
#endif
                    n++;
                }
#endif
            break;

        case IO_AGS_QPT:        /* quantum potential (Q) */
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
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
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
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
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
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
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTZETA)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (MyOutputFloat) P[pindex].AGS_zeta;
                    n++;
                }
#endif
            break;

        case IO_grHI:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grHI;
                    n++;
                }
            }
#endif
            break;

        case IO_grHII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grHII;
                    n++;
                }
            }
#endif
            break;

        case IO_grHM:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grHM;
                    n++;
                }
            }
#endif
            break;

        case IO_grHeI:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grHeI;
                    n++;
                }
            }
#endif
            break;

        case IO_grHeII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grHeII;
                    n++;
                }
            }
#endif
            break;

        case IO_grHeIII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grHeIII;
                    n++;
                }
            }
#endif
            break;

        case IO_grH2I:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grH2I;
                    n++;
                }
            }
#endif
            break;

        case IO_grH2II:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grH2II;
                    n++;
                }
            }
#endif
            break;

        case IO_grDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grDI;
                    n++;
                }
            }
#endif
            break;

        case IO_grDII:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grDII;
                    n++;
                }
            }
#endif
            break;

        case IO_grHDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = (MyOutputFloat) CellP[pindex].grHDI;
                    n++;
                }
            }
#endif
            break;

    case IO_TURB_DIFF_COEFF:
#ifdef TURB_DIFF_DYNAMIC
        for (n = 0; n < pc; pindex++) {
            if (P[pindex].Type == type) {
                *fp++ = (MyOutputFloat) CellP[pindex].TD_DiffCoeff;
                n++;
            }
        }
#endif

        break;

    case IO_DYNERROR:
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        for (n = 0; n < pc; pindex++) {
            if (P[pindex].Type == type) {
                *fp++ = (MyOutputFloat) CellP[pindex].TD_DynDiffCoeff_error;
                n++;
            }
        }
#endif
        break;

    case IO_DYNERRORDEFAULT:
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        for (n = 0; n < pc; pindex++) {
            if (P[pindex].Type == type) {
                *fp++ = (MyOutputFloat) CellP[pindex].TD_DynDiffCoeff_error_default;
                n++;
            }
        }
#endif
        break;

        case IO_LASTENTRY:
            endrun(213);   /* soft bad-stop: unreachable IO_LASTENTRY (internal: unknown iofield); buffer left unfilled, drains at write-turn poll */
            break;
    }

    *startindex = pindex;
}




/*! Width of one float-valued element of an input block, as staged in CommBuffer. HDF5 inputs are
 *  read at the file's own precision -- H5Dread converts the stored dtype into the requested memory
 *  type, so these blocks stage as double whatever the file holds. The MyInputFloat/MyInputPosFloat
 *  typedefs describe the unformatted-binary layout, which must be known to the bit, and apply only
 *  to that path. Single source of truth for the read-side element width.
 */
size_t get_input_float_bytes(enum iofields blocknr)
{
    if(All.ICFormat == 3) {return sizeof(double);}
    return (blocknr == IO_POS) ? sizeof(MyInputPosFloat) : sizeof(MyInputFloat);
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
        case IO_HERMITE_POS:
            if(mode)
                bytes_per_blockelement = 3 * sizeof(MyInputPosFloat);
            else
                bytes_per_blockelement = 3 * sizeof(MyOutputPosFloat);
            break;

        case IO_VEL:
        case IO_PARTVEL:
        case IO_HERMITE_VEL:
        case IO_ACCEL:
        case IO_HYDROACCEL:
        case IO_BFLD:
        case IO_GRADPHI:
        case IO_GRADRHO:
        case IO_RAD_ACCEL:
        case IO_VORT:
        case IO_SINK_ANGMOM:
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

        case IO_REFINE_FLAG:
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_TAG_ANCHOR
            bytes_per_blockelement = sizeof(MyIDType);
#endif
            break;

        case IO_SINKPROGS:
        case IO_GRAINTYPE:
        case IO_EOSCOMP:
        case IO_EOSPHASE:
        case IO_STAGE_PROTOSTAR:
        case IO_FLUIDTYPE:
            bytes_per_blockelement = sizeof(int);
            break;

        case IO_DAMAGE_POROSITY_DAMAGE:
        case IO_DAMAGE_POROSITY_DISTENTION:
        case IO_DAMAGE_POROSITY_ACTVCRACKS:
            if(mode)
                bytes_per_blockelement = sizeof(MyInputFloat);
            else
                bytes_per_blockelement = sizeof(MyOutputFloat);
            break;
            
        case IO_AGE_PROTOSTAR:
        case IO_MASS:
        case IO_SINK_DIST:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_RAD_TEMP:
        case IO_DUST_TEMP:
        case IO_UNSPMASS:
        case IO_UNSPJETMASS:
        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_DCRATE:
        case IO_PHRATE:
        case IO_KERNELRADIUS:
        case IO_SFR:
        case IO_AGE:
        case IO_OSTAR:
        case IO_DTOSTAR:
        case IO_GRAINSIZE:
        case IO_DUST_TO_GAS:
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
        case IO_SINKMASS:
        case IO_SINKDUSTMASSACC:
        case IO_SINKMASSALPHA:
        case IO_ACRB:
        case IO_SINKRAD:
        case IO_SINK_FORM_MASS:
        case IO_SINKMDOT:
        case IO_R_PROTOSTAR:
        case IO_MASS_D_PROTOSTAR:
        case IO_ZAMS_MASS:
        case IO_LUM_SINGLESTAR:
        case IO_EOSTEMP:
#ifdef TWO_TEMPERATURE_PLASMA
        case IO_TWOTEMP_TE:
#endif
        case IO_EOSABAR:
        case IO_EOSYE:
        case IO_EOSCS:
        case IO_PRESSURE:
        case IO_SOFT:
        case IO_AGS_HKERN:
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
        case IO_AMBIPOLAR:   /* the resistivities are scalar coefficients, one per particle */
        case IO_HALL:
        case IO_OHMIC:
            if(mode)
                bytes_per_blockelement = sizeof(MyInputFloat);
            else
                bytes_per_blockelement = sizeof(MyOutputFloat);
            break;


        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_SLOPES:
        case IO_COSMICRAY_KAPPA:
#if defined(COSMIC_RAY_FLUID) || defined(COSMIC_RAY_SUBGRID_LEBRON)
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

        case IO_RAD_OPACITY:
#if defined(OUTPUT_RT_RAD_OPACITY) && defined(RADTRANSFER)
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
#if defined(OUTPUT_EDDINGTON_TENSOR)
            if(mode)
                bytes_per_blockelement = (6*N_RT_FREQ_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (6*N_RT_FREQ_BINS) * sizeof(MyOutputFloat);
#endif
            break;


#ifdef NUCLEAR_NETWORK
        case IO_NUCLEAR_COMPOSITION:
            if(mode)
                bytes_per_blockelement = (NUM_NUCLEAR_SPECIES) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (NUM_NUCLEAR_SPECIES) * sizeof(MyOutputFloat);
            break;
#endif

        case IO_Z:
#ifdef METALS
            if(mode)
                bytes_per_blockelement = (NUM_METAL_SPECIES) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (NUM_METAL_SPECIES) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_DUSTCHEMZMET:
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            if(mode)
                bytes_per_blockelement = (NUM_ISMDUSTCHEM_ELEMENTS + NUM_ISMDUSTCHEM_SOURCES) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (NUM_ISMDUSTCHEM_ELEMENTS + NUM_ISMDUSTCHEM_SOURCES) * sizeof(MyOutputFloat);
#endif
            break;

        case IO_DUSTCHEMSPECIESMET:
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            if(mode)
                bytes_per_blockelement = (NUM_ISMDUSTCHEM_SPECIES) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (NUM_ISMDUSTCHEM_SPECIES) * sizeof(MyOutputFloat);
#endif
            break;    

        case IO_ISMDUSTCHEMMOL:
#if defined(GALSF_ISMDUSTCHEM_MODEL) && !defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            if(mode)
                bytes_per_blockelement = 2 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 2 * sizeof(MyOutputFloat);
#endif
            break;

            case IO_MACHNUM:
#if defined(OUTPUT_MACH_NUMBER)
            if(mode)
                bytes_per_blockelement = sizeof(MyInputFloat);
            else
                bytes_per_blockelement = sizeof(MyOutputFloat);
#endif
            break;

            case IO_SHOCKMACHNUM:
#if defined(OUTPUT_SHOCK_MACH_NUMBER)
            if(mode)
                bytes_per_blockelement = sizeof(MyInputFloat);
            else
                bytes_per_blockelement = sizeof(MyOutputFloat);
#endif
            break;

        case IO_DUSTCHEMGRAINBINNUMBERS:
        case IO_DUSTCHEMGRAINBINMASS:
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            if(mode)
                bytes_per_blockelement = (NUM_ISMDUSTCHEM_SPECIES * NUM_ISMDUSTCHEM_SIZE_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (NUM_ISMDUSTCHEM_SPECIES * NUM_ISMDUSTCHEM_SIZE_BINS) * sizeof(MyOutputFloat);
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
#ifdef CBE_INTEGRATOR
            if(mode)
                bytes_per_blockelement = (CBE_INTEGRATOR_NBASIS*CBE_INTEGRATOR_NMOMENTS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (CBE_INTEGRATOR_NBASIS*CBE_INTEGRATOR_NMOMENTS) * sizeof(MyOutputFloat);
            break;
#endif

        case IO_GRADVEL:
        case IO_GRADMAG:
        case IO_EOS_STRESS_TENSOR:
        case IO_TIDALTENSORPS:
            if(mode)
                bytes_per_blockelement = 9 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
            break;

        case IO_LASTENTRY:
            endrun(214);   /* soft bad-stop: unreachable IO_LASTENTRY; nonzero sentinel so caller's bufsize/bytes does not divide by zero */
            bytes_per_blockelement = 1;
            break;
    }

    /* rescale float-valued input blocks from the compile-time binary layout to the width they are
     * actually staged at (see get_input_float_bytes) -- for HDF5 that is double, always. Keyed on
     * All.ICFormat, the same variable that selects the H5Dread branch in read_ic(), so the sizing
     * here and the reader cannot disagree about which format is being read. */
    if(mode && All.ICFormat == 3)
    {
        int typekey_read = get_datatype_in_block(blocknr);
        if(typekey_read == 1 || typekey_read == 3)
        {
            size_t width_binary = (blocknr == IO_POS) ? sizeof(MyInputPosFloat) : sizeof(MyInputFloat);
            /* every float-valued case above is a whole number of input elements wide; if some future
             * block is not, the rescale would silently mis-size it against empty_read_buffer's unpack */
            if(bytes_per_blockelement % ((int) width_binary)) {endrun(219);}
            bytes_per_blockelement = (int) ((bytes_per_blockelement / width_binary) * get_input_float_bytes(blocknr));
        }
    }

    return bytes_per_blockelement;
}

int get_datatype_in_block(enum iofields blocknr)
{
    int typekey;
    switch (blocknr)
    {
#if !defined(OUTPUT_IN_DOUBLEPRECISION) || (defined(INPUT_POSITIONS_IN_DOUBLE) && !defined(INPUT_IN_DOUBLEPRECISION))
        case IO_POS:
        case IO_HERMITE_POS:
            typekey = 3; /* pos outputs in HDF5 are double automatically, to prevent overlaps */
            break;
#endif

        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
            typekey = 2;		/* native long long */
            break;

#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_TAG_ANCHOR
        case IO_REFINE_FLAG:
            typekey = 2;		/* native long long */
            break;
#endif

        case IO_SINKPROGS:
        case IO_GRAINTYPE:
        case IO_EOSCOMP:
        case IO_EOSPHASE:
        case IO_STAGE_PROTOSTAR:
        case IO_FLUIDTYPE:
            typekey = 0;		/* native int */
            break;

        case IO_DAMAGE_POROSITY_DAMAGE:
        case IO_DAMAGE_POROSITY_DISTENTION:
        case IO_DAMAGE_POROSITY_ACTVCRACKS:
            typekey = 1;		/* native MyOutputFloat */
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
        case IO_PARTVEL:
        case IO_HERMITE_POS:
        case IO_HERMITE_VEL:
        case IO_ACCEL:
        case IO_HYDROACCEL:
        case IO_BFLD:
        case IO_GRADPHI:
        case IO_GRADRHO:
        case IO_RAD_ACCEL:
        case IO_VORT:
        case IO_SINK_ANGMOM:
            values = 3;
            break;

        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
        case IO_REFINE_FLAG:
        case IO_MASS:
        case IO_SINK_DIST:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_RAD_TEMP:
        case IO_DUST_TEMP:
        case IO_UNSPMASS:
        case IO_UNSPJETMASS:
        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_DCRATE:
        case IO_PHRATE:
        case IO_KERNELRADIUS:
        case IO_SFR:
        case IO_AGE:
        case IO_OSTAR:
        case IO_DTOSTAR:
        case IO_GRAINSIZE:
        case IO_DUST_TO_GAS:
        case IO_GRAINTYPE:
        case IO_FLUIDTYPE:
        case IO_DELAYTIME:
        case IO_HSMS:
        case IO_POT:
        case IO_DTENTR:
        case IO_TSTP:
        case IO_VDIV:
        case IO_DIVB:
	    case IO_AMBIPOLAR:
        case IO_HALL:
        case IO_OHMIC:
        case IO_ABVC:
        case IO_AMDC:
        case IO_PHI:
        case IO_COOLRATE:
        case IO_SINKMASS:
        case IO_SINKDUSTMASSACC:
        case IO_SINKMASSALPHA:
        case IO_ACRB:
        case IO_SINKRAD:
        case IO_SINK_FORM_MASS:
        case IO_SINKMDOT:
        case IO_R_PROTOSTAR:
        case IO_MASS_D_PROTOSTAR:
        case IO_ZAMS_MASS:
        case IO_STAGE_PROTOSTAR:
        case IO_AGE_PROTOSTAR:
        case IO_LUM_SINGLESTAR:
        case IO_SINKPROGS:
        case IO_EOSTEMP:
#ifdef TWO_TEMPERATURE_PLASMA
        case IO_TWOTEMP_TE:
#endif
        case IO_EOSABAR:
        case IO_EOSCS:
        case IO_EOSCOMP:
        case IO_EOSPHASE:
        case IO_EOSYE:
        case IO_PRESSURE:
        case IO_SOFT:
        case IO_AGS_HKERN:
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
        case IO_DAMAGE_POROSITY_DAMAGE:
        case IO_DAMAGE_POROSITY_DISTENTION:
        case IO_DAMAGE_POROSITY_ACTVCRACKS:
            values = 1;
            break;

#ifdef NUCLEAR_NETWORK
        case IO_NUCLEAR_COMPOSITION:
            values = NUM_NUCLEAR_SPECIES;
            break;
#endif

        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_SLOPES:
        case IO_COSMICRAY_KAPPA:
#if defined(COSMIC_RAY_FLUID) || defined(COSMIC_RAY_SUBGRID_LEBRON)
            values = N_CR_PARTICLE_BINS;
#endif
            break;

        case IO_COSMICRAY_ALFVEN:
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            values = (2*N_CR_PARTICLE_BINS);
#endif
            break;

        case IO_CBE_MOMENTS:
#ifdef CBE_INTEGRATOR
            values = (CBE_INTEGRATOR_NBASIS*CBE_INTEGRATOR_NMOMENTS);
#endif
            break;

        case IO_EDDINGTON_TENSOR:
#if defined(OUTPUT_EDDINGTON_TENSOR)
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
            
        case IO_RAD_OPACITY:
#if defined(OUTPUT_RT_RAD_OPACITY) && defined(RADTRANSFER)
            values = N_RT_FREQ_BINS;
#endif
            break;

        case IO_Z:
#ifdef METALS
            values = NUM_METAL_SPECIES;
#endif
            break;

        case IO_DUSTCHEMZMET:
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            values = (NUM_ISMDUSTCHEM_ELEMENTS + NUM_ISMDUSTCHEM_SOURCES);
#endif
            break;

        case IO_DUSTCHEMSPECIESMET:
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            values = NUM_ISMDUSTCHEM_SPECIES;
#endif
            break;

        case IO_ISMDUSTCHEMMOL:
#if defined(GALSF_ISMDUSTCHEM_MODEL) && !defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            values = 2;
#endif
            break;

        case IO_MACHNUM:
#if defined(OUTPUT_MACH_NUMBER)
            values = 1;
#endif
            break;

        case IO_SHOCKMACHNUM:
#if defined(OUTPUT_SHOCK_MACH_NUMBER)
            values = 1;
#endif
            break;

        case IO_DUSTCHEMGRAINBINNUMBERS:
        case IO_DUSTCHEMGRAINBINMASS:
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            values = (NUM_ISMDUSTCHEM_SPECIES*NUM_ISMDUSTCHEM_SIZE_BINS);
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

        case IO_GRADVEL:
        case IO_GRADMAG:
        case IO_TIDALTENSORPS:
        case IO_EOS_STRESS_TENSOR:
            values = 9;
            break;

        case IO_LASTENTRY:
            endrun(215);   /* soft bad-stop: unreachable IO_LASTENTRY; nonzero sentinel so downstream sizing does not divide by zero */
            values = 1;
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

    /* NOTE: this is an I/O-block membership bitmask (which types write e.g. the ages block),
       a DIFFERENT concept from is_galsf_stellar_candidate_type() (feedback/source activation):
       it is a type bitmask, includes sinks (type 5), and needs no zoom/grain exclusion since
       writing an extra harmless output field for a special particle is not a physics change. */
    int valid_star_types=16; if(All.ComovingIntegrationOn==0) {valid_star_types+=4+8;} /* used in e.g. ages block below, not for all, but easier to define here */
#ifdef SINK_PARTICLES
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
        case IO_REFINE_FLAG:
        case IO_POT:
        case IO_SOFT:
        case IO_AGS_HKERN:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_PSI_RE:
        case IO_AGS_PSI_IM:
        case IO_AGS_ZETA:
        case IO_SINK_DIST:
        case IO_CBE_MOMENTS:
        case IO_TIDALTENSORPS:
        case IO_DUST_TO_GAS:
        case IO_FLUIDTYPE:
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

        case IO_HERMITE_POS:
        case IO_HERMITE_VEL:
#if defined(IO_HERMITE_SYNC) && defined(HERMITE_INTEGRATION)
            /* only the types the integrator owns, and only those: for anything else the values
               would be a byte-for-byte duplicate of Coordinates and Velocities.  Both flags are
               tested because they are independent -- the output flag can be set on a build whose
               integrator is disabled, where the blocks are absent altogether. */
            {
                int n_hermite = 0;
                for(i = 0; i < 6; i++)
                {
                    typelist[i] = ((1 << i) & (HERMITE_INTEGRATION)) ? 1 : 0;
                    if(typelist[i]) {n_hermite += header.npart[i];}
                }
                return n_hermite;
            }
#endif
            break;

        case IO_PARTVEL:
        case IO_RAD_ACCEL:
        case IO_HYDROACCEL:
        case IO_RADGAMMA:
        case IO_RAD_OPACITY:
        case IO_RAD_TEMP:
        case IO_DUST_TEMP:
        case IO_RAD_FLUX:
        case IO_EDDINGTON_TENSOR:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_DCRATE:
        case IO_PHRATE:
        case IO_DELAYTIME:
        case IO_SFR:
        case IO_DTENTR:
        case IO_BFLD:
        case IO_AMBIPOLAR:
        case IO_HALL:
        case IO_OHMIC:        
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
        case IO_GRADRHO:
        case IO_GRADVEL:
        case IO_GRADMAG:
        case IO_COOLRATE:
        case IO_EOSTEMP:
#ifdef TWO_TEMPERATURE_PLASMA
        case IO_TWOTEMP_TE:
#endif
        case IO_EOSABAR:
        case IO_EOSYE:
#ifdef NUCLEAR_NETWORK
        case IO_NUCLEAR_COMPOSITION:
#endif
        case IO_EOSCS:
        case IO_EOS_STRESS_TENSOR:
        case IO_EOSCOMP:
        case IO_EOSPHASE:
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
        case IO_DUSTCHEMZMET:
        case IO_DUSTCHEMSPECIESMET:
        case IO_ISMDUSTCHEMMOL:
        case IO_MACHNUM:
        case IO_SHOCKMACHNUM:
        case IO_DUSTCHEMGRAINBINNUMBERS:
        case IO_DUSTCHEMGRAINBINMASS:
        case IO_DAMAGE_POROSITY_DAMAGE:
        case IO_DAMAGE_POROSITY_DISTENTION:
        case IO_DAMAGE_POROSITY_ACTVCRACKS:
            for(i = 1; i < 6; i++) {typelist[i] = 0;}
            return ngas;
            break;

        case IO_KERNELRADIUS:
#if defined(GRAIN_FLUID)
            return nall;
#endif
            for(i = 1; i < 6; i++) {typelist[i] = 0;}
            return ngas;
            
        case IO_AGE:
            for(i=0; i<6; i++) {if(!((1 << i) & (valid_star_types))) {typelist[i]=0;}}
            return nstars_tot;
            break;

        case IO_OSTAR:
        case IO_DTOSTAR:
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

        case IO_SINKMASS:
        case IO_SINKDUSTMASSACC:
        case IO_SINKMASSALPHA:
        case IO_SINK_ANGMOM:
        case IO_UNSPMASS:
        case IO_UNSPJETMASS:
        case IO_ACRB:
        case IO_SINKRAD:
        case IO_SINK_FORM_MASS:
        case IO_SINKMDOT:
        case IO_R_PROTOSTAR:
        case IO_MASS_D_PROTOSTAR:
        case IO_ZAMS_MASS:
        case IO_STAGE_PROTOSTAR:
        case IO_AGE_PROTOSTAR:
        case IO_LUM_SINGLESTAR:
        case IO_SINKPROGS:
            for(i = 0; i < 6; i++) {if(i != 5) {typelist[i] = 0;}}
            return header.npart[5];
            break;

        case IO_LASTENTRY:
            endrun(216);   /* soft bad-stop: unreachable IO_LASTENTRY; return 0 (no particles in unknown block), drains at next poll */
            return 0;
    }

    endrun(212);   /* soft bad-stop: unreachable end of get_particles_in_block (internal: unknown iofield) */
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
        case IO_KERNELRADIUS:
            return 1;			/* always present */
            break;

        case IO_REFINE_FLAG:
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_TAG_ANCHOR
            return 1;			/* present iff the tag-anchor mode is compiled in */
#endif
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

        case IO_RAD_OPACITY:
#if defined(OUTPUT_RT_RAD_OPACITY) && defined(RADTRANSFER)
            return 1;
#endif
            break;

        case IO_RAD_TEMP:
#if defined(RADTRANSFER) && defined(RT_INFRARED)
            return 1;
#endif
            break;

        case IO_DUST_TEMP:
#if (defined(RADTRANSFER) && defined(RT_INFRARED)) || (defined(OUTPUT_DUST_TEMPERATURE) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2))
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

        case IO_DUST_TO_GAS:
#ifdef OUTPUT_DUST_TO_GAS_RATIO
            return 1;
#endif
            break;

        case IO_GRAINTYPE:
#ifdef PIC_MHD
            return 1;
#endif
            break;

        case IO_FLUIDTYPE:
#ifdef HYDRO_MULTIFLUID
            return 1;
#endif
            break;

        case IO_Z:
#ifdef METALS
            return 1;
#endif
            break;

        case IO_DUSTCHEMZMET:
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            return 1;
#endif
            break;               
          
        case IO_DUSTCHEMSPECIESMET:
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            return 1;
#endif
            break;

        case IO_ISMDUSTCHEMMOL:
#if defined(GALSF_ISMDUSTCHEM_MODEL) && !defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            return 1;
#endif
            break; 

        case IO_DUSTCHEMGRAINBINNUMBERS:
        case IO_DUSTCHEMGRAINBINMASS:
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            return 1;
#endif
            break;

        case IO_MACHNUM:
#if defined(OUTPUT_MACH_NUMBER)
            return 1;
#endif
            break;

        case IO_SHOCKMACHNUM:
#if defined(OUTPUT_SHOCK_MACH_NUMBER)
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
#if defined(CHIMES_NH_OUTPUT) && defined(OUTPUT_DENS_AROUND_NONGAS)
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
#ifdef OUTPUT_DENS_AROUND_NONGAS
            return 1;
#endif
            break;

        case IO_DELAY_TIME_HII:
#if defined(GALSF_FB_FIRE_RT_HIIHEATING) && defined(OUTPUT_DELAY_TIME_HII)
            return 1;
#endif
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

        case IO_DAMAGE_POROSITY_DAMAGE:
        case IO_DAMAGE_POROSITY_DISTENTION:
        case IO_DAMAGE_POROSITY_ACTVCRACKS:
#ifdef EOS_DAMAGE_POROSITY
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
            
        case IO_UNSPMASS:
#if defined(SINK_WIND_SPAWN)
            return 1; /* mass already removed from the sink and awaiting spawning: evolving state, so written unconditionally rather than as an optional diagnostic -- otherwise a snapshot restart silently loses it */
#endif   
            break;

        case IO_UNSPJETMASS:
#if defined(SINK_WIND_SPAWN) && defined(SINGLE_STAR_FB_JETS)
            return 1;
#endif
            break;

        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_PHRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
            return 1;
#endif
            break;

        case IO_DCRATE:
#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING) && defined(GALSF_ISMDUSTCHEM_MODEL) && defined(GALSF_ISMDUSTCHEM_HIGHTEMPDUSTCOOLING)
            return 1;
#endif
            break;

        case IO_POT:
#if defined(OUTPUT_POTENTIAL)
            return 1;
#endif
            break;

        case IO_HERMITE_POS:
        case IO_HERMITE_VEL:
#if defined(IO_HERMITE_SYNC) && defined(HERMITE_INTEGRATION)
            return 1;   /* only where a synchronized pair can actually be built; with the integrator
                           disabled the blocks are absent rather than carrying the ordinary mixed
                           state under a name that promises otherwise */
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

        case IO_HYDROACCEL:
#ifdef OUTPUT_HYDROACCELERATION
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

        case IO_AMBIPOLAR:
#if defined(MHD_NON_IDEAL) && defined(OUTPUT_MHD_RESISTIVITY)
	    return 1;
#endif
            break;
        case IO_HALL:
#if defined(MHD_NON_IDEAL) && defined(OUTPUT_MHD_RESISTIVITY)
            return 1;
#endif
            break;        
        case IO_OHMIC:            
#if defined(MHD_NON_IDEAL) && defined(OUTPUT_MHD_RESISTIVITY)
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

        case IO_DTOSTAR:
#ifdef GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF
            return 1;
#endif
            break;

        case IO_COSMICRAY_ENERGY:
#if defined(COSMIC_RAY_FLUID) || defined(COSMIC_RAY_SUBGRID_LEBRON)
            return 1;
#endif
            break;

        case IO_COSMICRAY_SLOPES:
#if defined(COSMIC_RAY_FLUID) && defined(CRFLUID_EVOLVE_SPECTRUM)
            return 1;
#endif
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

        case IO_GRADRHO:
#if defined(OUTPUT_GRADIENT_RHO)
            return 1;
#endif
            break;

        case IO_GRADVEL:
#if defined(OUTPUT_GRADIENT_VEL)
            return 1;
#endif
            break;

        case IO_GRADMAG:
#if defined(OUTPUT_GRADIENT_MAG)
            return 1;
#endif
            break;

        case IO_COOLRATE:
#ifdef OUTPUT_COOLRATE
            return 1;
#endif
            break;

        case IO_SINK_ANGMOM:
#ifdef SINK_FOLLOW_ACCRETED_ANGMOM
            return 1;
#endif
            break;

        case IO_ACRB:
        case IO_SINKRAD:
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
            return 1;
#endif
            break;

        case IO_SINK_FORM_MASS:
#ifdef SINK_PARTICLES
            return 1;
#endif
            break;
	    

        case IO_SINKMASS:
#ifdef SINK_PARTICLES
            return 1;
#endif
            break;

        case IO_SINKDUSTMASSACC:
#if defined(SINK_PARTICLES) && defined(GRAIN_FLUID)
            return 1;
#endif
            break;

        case IO_SINKMASSALPHA:
#ifdef SINK_ALPHADISK_ACCRETION
            return 1;
#endif
            break;

        case IO_SINKMDOT:
#ifdef SINK_PARTICLES
            return 1;
#endif
            break;

        case IO_R_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            return 1;
#endif
            break;

        case IO_MASS_D_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            return 1;
#endif
            break;

        case IO_ZAMS_MASS:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            return 1;
#endif
            break;

        case IO_STAGE_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            return 1;
#endif
            break;
            
        case IO_AGE_PROTOSTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            return 1;
#endif
            break;

        case IO_LUM_SINGLESTAR:
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
            return 1;
#endif
            break;

        case IO_SINK_DIST:
#if defined(SINK_CALC_DISTANCES) && defined(OUTPUT_SINK_DISTANCES)
            return 1;
#endif
            break;

        case IO_SINKPROGS:
#ifdef SINK_COUNTPROGS
            return 1;
#endif
            break;

        case IO_TIDALTENSORPS:
#ifdef OUTPUT_TIDAL_TENSOR
            return 1;
#endif
            break;

        case IO_EOSTEMP:
            return 1;
            break;

#ifdef TWO_TEMPERATURE_PLASMA
        case IO_TWOTEMP_TE:
            return 1;
            break;
#endif

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
#if defined(EOS_TILLOTSON) || defined(EOS_ANEOS)
            return 1;
#endif
            break;

        case IO_EOSPHASE:
#ifdef EOS_ANEOS
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

#ifdef NUCLEAR_NETWORK
        case IO_NUCLEAR_COMPOSITION:
            return 1;
#endif

        case IO_CBE_MOMENTS:
#ifdef CBE_INTEGRATOR
            return 1;
#endif
            break;

        case IO_PARTVEL:
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            return 1;
#endif
            break;

        case IO_EDDINGTON_TENSOR:
#if defined(OUTPUT_EDDINGTON_TENSOR)
            return 1;
#endif
            break;

        case IO_SOFT:
#if defined(OUTPUT_SOFTENING)
            return 1;
#endif
            break;

        case IO_AGS_HKERN:
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
            return 1;
#endif
            break;

        case IO_AGS_RHO:
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && (defined(DM_FUZZY) || defined(CBE_INTEGRATOR))
            return 1;
#endif
            break;

        case IO_AGS_QPT:
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY)
            return 1;
#endif
            break;

        case IO_AGS_PSI_RE:
        case IO_AGS_PSI_IM:
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && defined(DM_FUZZY) && (DM_FUZZY > 0)
            return 1;
#endif
            break;

        case IO_AGS_ZETA:
#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTZETA)
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
        case IO_REFINE_FLAG:
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_TAG_ANCHOR
            strncpy(label, "RefF", 4);
#endif
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
        case IO_UNSPMASS:
            strncpy(label, "USPM", 4);
            break;
        case IO_UNSPJETMASS:
            strncpy(label, "USPJ", 4);
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
        case IO_DCRATE:
            strncpy(label, "DCRATE", 4);
            break;        
        case IO_PHRATE:
            strncpy(label, "PHRATE", 4);
            break;
        case IO_KERNELRADIUS:
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
        case IO_DUST_TO_GAS:
            strncpy(label, "GRDG", 4);
            break;
        case IO_GRAINTYPE:
            strncpy(label, "GRTP", 4);
            break;
        case IO_FLUIDTYPE:
            strncpy(label, "FLTP", 4);
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
        case IO_DUSTCHEMZMET:
            strncpy(label, "DZ  ", 4);
            break;
        case IO_DUSTCHEMSPECIESMET:
            strncpy(label, "SPEZ", 4);
            break;
        case IO_ISMDUSTCHEMMOL:
            strncpy(label, "DMOL", 4);
            break;
        case IO_MACHNUM:
            strncpy(label, "MACH", 4);
            break;
        case IO_SHOCKMACHNUM:
            strncpy(label, "SMAC", 4);
            break;
        case IO_DAMAGE_POROSITY_DAMAGE:
            strncpy(label, "DAMG", 4);
            break;
        case IO_DAMAGE_POROSITY_DISTENTION:
            strncpy(label, "DIST", 4);
            break;
        case IO_DAMAGE_POROSITY_ACTVCRACKS:
            strncpy(label, "ACTV", 4);
            break;
        case IO_DUSTCHEMGRAINBINNUMBERS:
            strncpy(label, "DBNU", 4);
            break;
        case IO_DUSTCHEMGRAINBINMASS:
            strncpy(label, "DBMA", 4);
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
        case IO_HERMITE_POS:
            strncpy(label, "HPOS", 4);
            break;
        case IO_HERMITE_VEL:
            strncpy(label, "HVEL", 4);
            break;
        case IO_ACCEL:
            strncpy(label, "ACCE", 4);
            break;
        case IO_HYDROACCEL:
            strncpy(label, "ACCH", 4);
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
        case IO_AMBIPOLAR:
            strncpy(label, "AMBI", 4);
            break;
        case IO_OHMIC:
            strncpy(label, "OHM ", 4);
            break;
        case IO_HALL:
            strncpy(label, "HALL", 4);
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
            strncpy(label, "NOST", 4);
            break;
        case IO_DTOSTAR:
            strncpy(label, "DTST", 4);
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
        case IO_GRADRHO:
            strncpy(label, "GRHO", 4);
            break;
        case IO_GRADVEL:
            strncpy(label, "GVEL", 4);
            break;
        case IO_GRADMAG:
            strncpy(label, "GRDB", 4);
            break;
        case IO_COOLRATE:
            strncpy(label, "COOR", 4);
            break;
        case IO_SINKMASS:
            strncpy(label, "BHMA", 4);
            break;
        case IO_SINKDUSTMASSACC:
            strncpy(label, "BHDM", 4);
            break;
        case IO_SINK_DIST:
            strncpy(label, "BHR ", 4);
            break;
        case IO_SINKMASSALPHA:
            strncpy(label, "BHMa", 4);
            break;
        case IO_SINK_ANGMOM:
            strncpy(label, "BHmJ", 4);
            break;
        case IO_ACRB:
            strncpy(label, "ACRB", 4);
            break;
        case IO_SINKRAD:
            strncpy(label, "SRAD", 4);
            break;
        case IO_SINK_FORM_MASS:
            strncpy(label, "SMAS", 4);
            break;	    
        case IO_SINKMDOT:
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
        case IO_SINKPROGS:
            strncpy(label, "BHPC", 4);
            break;
        case IO_TIDALTENSORPS:
            strncpy(label, "TIPS", 4);
            break;
        case IO_EOSTEMP:
            strncpy(label, "TEMP", 4);
            break;
#ifdef TWO_TEMPERATURE_PLASMA
        case IO_TWOTEMP_TE:
            strncpy(label, "TELE", 4);
            break;
#endif
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
        case IO_EOSPHASE:
            strncpy(label, "PHID", 4);
            break;
        case IO_PARTVEL:
            strncpy(label, "PVEL", 4);
            break;
        case IO_EOSYE:
            strncpy(label, "YE  ", 4);
            break;
#ifdef NUCLEAR_NETWORK
        case IO_NUCLEAR_COMPOSITION:
            strncpy(label, "NUCX", 4);
            break;
#endif
        case IO_PRESSURE:
            strncpy(label, "P   ", 4);
            break;
        case IO_RADGAMMA:
            strncpy(label, "RADG", 4);
            break;
        case IO_RAD_OPACITY:
            strncpy(label, "RADO", 4);
            break;
        case IO_RAD_TEMP:
            strncpy(label, "RADT", 4);
            break;
        case IO_DUST_TEMP:
            strncpy(label, "DTMP", 4);
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
        case IO_SOFT:
            strncpy(label, "SOFT", 4);
            break;
        case IO_AGS_HKERN:
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
            endrun(217);   /* soft bad-stop: unreachable IO_LASTENTRY (internal: unknown iofield) */
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
        case IO_REFINE_FLAG:
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_TAG_ANCHOR
            strcpy(buf, "RefinementFlag");
#endif
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
        case IO_RAD_OPACITY:
            strcpy(buf, "PhotonOpacity");
            break;
        case IO_RAD_TEMP:
            strcpy(buf, "IRBand_Radiation_Temperature");
            break;
        case IO_DUST_TEMP:
            strcpy(buf, "Dust_Temperature");
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
        case IO_UNSPMASS:
            strcpy(buf, "Unspawned_Wind_Mass");
            break;
        case IO_UNSPJETMASS:
            strcpy(buf, "Unspawned_Jet_Mass");
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
        case IO_DCRATE:
            strcpy(buf, "DustCoolingRate");
            break;        
        case IO_PHRATE:
            strcpy(buf, "PElecHeatingRate");
            break;
        case IO_DELAYTIME:
            strcpy(buf, "DelayTime");
            break;
        case IO_KERNELRADIUS:
            strcpy(buf, "KernelMaxRadius");
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
        case IO_DUST_TO_GAS:
            strcpy(buf, "DustToGasRatio_Local");
            break;
        case IO_GRAINTYPE:
            strcpy(buf, "PICParticleType");
            break;
        case IO_FLUIDTYPE:
            strcpy(buf, "FluidType");
            break;
        case IO_HSMS:
            strcpy(buf, "StellarKernelMaxRadius");
            break;
        case IO_Z:
            strcpy(buf, "Metallicity");
            break;
        case IO_DUSTCHEMZMET:
            strcpy(buf, "DustMetallicity");
            break;
        case IO_DUSTCHEMSPECIESMET:
            strcpy(buf, "DustSpeciesAbundance");
            break;
        case IO_ISMDUSTCHEMMOL:
            strcpy(buf, "DustMolecularSpeciesFractions");
            break;
        case IO_MACHNUM:
            strcpy(buf, "MachNumber");
            break;
        case IO_SHOCKMACHNUM:
            strcpy(buf, "ShockMachNumber");
            break;
        case IO_DAMAGE_POROSITY_DAMAGE:
            strcpy(buf, "Damage");
            break;
        case IO_DAMAGE_POROSITY_DISTENTION:
            strcpy(buf, "Distention");
            break;
        case IO_DAMAGE_POROSITY_ACTVCRACKS:
            strcpy(buf, "ActiveCracks");
            break;
        case IO_DUSTCHEMGRAINBINNUMBERS:
            strcpy(buf, "DustBinNumbers");
            break;
        case IO_DUSTCHEMGRAINBINMASS:
            strcpy(buf, "DustBinMasses");
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
        case IO_HERMITE_POS:
            strcpy(buf, "HermiteSyncCoordinates");
            break;
        case IO_HERMITE_VEL:
            strcpy(buf, "HermiteSyncVelocities");
            break;
        case IO_ACCEL:
            strcpy(buf, "Acceleration");
            break;
        case IO_HYDROACCEL:
            strcpy(buf, "HydroAcceleration");
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
        case IO_AMBIPOLAR:
            strcpy(buf, "AmbipolarResistivity");
            break;            
        case IO_OHMIC:
            strcpy(buf, "OhmicResistivity");
            break;            
        case IO_HALL:
            strcpy(buf, "HallResistivity");
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
        case IO_DTOSTAR:
            strcpy(buf, "TimeDistribOfStarFormation");
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
        case IO_GRADRHO:
            strcpy(buf, "DensityGradient");
            break;
        case IO_GRADVEL:
            strcpy(buf, "VelocityGradient");
            break;
        case IO_GRADMAG:
            strcpy(buf, "MagneticFieldGradient");
            break;
        case IO_COOLRATE:
            strcpy(buf, "CoolingRate");
            break;
        case IO_SINKMASS:
            strcpy(buf, "Sink_Mass");
            break;
        case IO_SINKDUSTMASSACC:
            strcpy(buf, "Sink_Dust_Mass");
            break;
        case IO_SINK_DIST:
            strcpy(buf, "Sink_Distance");
            break;
        case IO_SINKMASSALPHA:
            strcpy(buf, "Sink_Mass_Reservoir");
            break;
        case IO_SINK_ANGMOM:
            strcpy(buf, "Sink_Specific_AngMom");
            break;
        case IO_ACRB:
            strcpy(buf, "Sink_AccretionLength");
            break;
        case IO_SINKRAD:
            strcpy(buf, "Sink_Radius");
            break;
        case IO_SINK_FORM_MASS:
            strcpy(buf, "Sink_InitialMass");
            break;	    
        case IO_SINKMDOT:
            strcpy(buf, "Sink_Mdot");
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
        case IO_SINKPROGS:
            strcpy(buf, "Sink_NProgenitors");
            break;
        case IO_TIDALTENSORPS:
            strcpy(buf, "TidalTensor");
            break;
        case IO_EOSTEMP:
            strcpy(buf, "Temperature");
            break;
#ifdef TWO_TEMPERATURE_PLASMA
        case IO_TWOTEMP_TE:
            strcpy(buf, "ElectronTemperature");
            break;
#endif
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
        case IO_EOSPHASE:
            strcpy(buf, "PhaseID");
            break;
        case IO_PARTVEL:
            strcpy(buf, "ParticleVelocities");
            break;
        case IO_EOSYE:
            strcpy(buf, "Ye");
            break;
#ifdef NUCLEAR_NETWORK
        case IO_NUCLEAR_COMPOSITION:
            strcpy(buf, "NuclearComposition");
            break;
#endif
        case IO_PRESSURE:
            strcpy(buf, "Pressure");
            break;
        case IO_EDDINGTON_TENSOR:
            strcpy(buf, "EddingtonTensor");
            break;
        case IO_SOFT:
            strcpy(buf, "Softening_KernelRadius");
            break;
        case IO_AGS_HKERN:
            strcpy(buf, "AGS-KernelRadius");
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
            endrun(218);   /* soft bad-stop: unreachable IO_LASTENTRY (internal: unknown iofield); buf already "default" */
            break;
    }
}


/*! Legacy-compatibility dataset name for HDF5 snapshot restarts. Several IO fields were
 *  renamed in the modern code (gas kernel length KernelMaxRadius <- SmoothingLength; the
 *  sink fields Sink_... <- legacy BH_... or Sink... names). When restarting from an older snapshot that still
 *  carries the pre-rename name, the primary H5Dopen on the modern name fails; read_ic then
 *  retries with the alternate name returned here. Fills buf with an empty string for any
 *  field whose name is unchanged (no fallback attempted). HDF5-only. */
void get_dataset_name_legacy_alias(enum iofields blocknr, char *buf)
{
    buf[0] = '\0';
    switch (blocknr)
    {
        case IO_KERNELRADIUS:   strcpy(buf, "SmoothingLength");        break;   /* gas kernel length */
        case IO_HSMS:           strcpy(buf, "StellarSmoothingLength"); break;   /* non-gas kernel length */
        case IO_SINKMASS:       strcpy(buf, "BH_Mass");               break;
        case IO_SINKMASSALPHA:  strcpy(buf, "BH_Mass_AlphaDisk");     break;
        case IO_SINK_ANGMOM:    strcpy(buf, "BH_Specific_AngMom");    break;
        case IO_ACRB:           strcpy(buf, "BH_AccretionLength");    break;
        case IO_SINKRAD:        strcpy(buf, "SinkRadius");            break;
        case IO_SINK_FORM_MASS: strcpy(buf, "SinkInitialMass");       break;
        case IO_SINKMDOT:       strcpy(buf, "BH_Mdot");               break;
        case IO_SINKPROGS:      strcpy(buf, "BH_NProgs");             break;
        case IO_SINK_DIST:      strcpy(buf, "BH_Dist");               break;
        case IO_SINKDUSTMASSACC:strcpy(buf, "BH_Dust_Mass");          break;
        default:                                                      break;
    }
}



/*! This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the gas cells.  If
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
    int write_status = 0;   /* soft bad-stop flag for this (writeTask-local) file write; on a failed open we skip the
                             * fd/my_fwrite calls but KEEP participating in every MPI_Recv from peers (the send/recv
                             * choreography is particle-count-driven, not open-status-driven, so peers never deadlock).
                             * The bad-stop drains at the savepositions per-turn poll. */

    hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
    hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
    herr_t hdf5_status;
    hsize_t dims[2], count[2], start[2];
    int rank = 0, pcsum = 0;
    char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];

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

            snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s.hdf5", fname);
            hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

            hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

            for(type = 0; type < 6; type++)
            {
                if(header.npart[type] > 0)
                {
                    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "/PartType%d", type);
                    hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
                }
            }

            write_header_attributes_in_hdf5(hdf5_headergrp);

        }
        else
        {
            if(!(fd = fopen(fname, "w")))
            {
                printf("can't open file `%s' for writing snapshot.\n", fname);
                endrun(123); write_status = 123;   /* soft: cannot open snapshot file; skip local fd writes below, drain at the savepositions per-turn poll */
            }

            if(!write_status)
            {
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
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
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
            size_t comm_chunk_mb = All.CommChunkSize;
            blockmaxlen = (size_t) ((comm_chunk_mb * 1024 * 1024) / bytes_per_blockelement);
            npart = get_particles_in_block(blocknr, &typelist[0]);

            if(npart > 0)
            {
                if(ThisTask == 0)
                {
                    char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
                    get_dataset_name(blocknr, buf);
                    printf("writing block %d (%s)...\n", bnr, buf);
                }

                if(ThisTask == writeTask && !write_status)
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
#if IO_SINKS_ONLY_SNAPSHOT_FREQUENCY > 0
                    if ( typelist[type] && ( (type!=0) || !( (All.SnapshotFileCount-1) % (IO_SINKS_ONLY_SNAPSHOT_FREQUENCY+1) ) ) ) //we skip type 0 (gas) data for the reduced snapshots
#else
                    if(typelist[type])
#endif
                    {

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
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
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
                            	/* Use many small chunks (~4096 rows) so incremental hyperslab writes
                            	   don't trigger repeated decompress-recompress of huge chunks. The old
                            	   code used dims[0]/10 which could be millions of rows, causing massive
                            	   slowdowns with MPI parallel writes. */
                            	hsize_t target_chunk_rows = 4096;
                            	if(target_chunk_rows > dims[0]) {target_chunk_rows = dims[0];}
                            	hsize_t cdims[2]; cdims[0] = target_chunk_rows; cdims[1] = dims[1];
                            	hdf5_status = H5Pset_chunk(plist_id, rank, cdims);
                            	/* Shuffle reorders bytes for better compression (nearly free CPU cost).
                            	   Deflate/gzip level 1 is fast and gives ~1.5-2x compression on float data.
                            	   The combination shuffle+deflate(1) is much faster than deflate(4) alone
                            	   and achieves similar or better compression ratios. */
                            	hdf5_status = H5Pset_shuffle(plist_id);
                            	hdf5_status = H5Pset_deflate(plist_id, 1);
                            	hdf5_dataset = H5Dcreate2(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                            	H5Pclose(plist_id);
                            } else {
                            	hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
                            }
#endif
                            pcsum = 0;
                        }

                        
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

                                    }
                                    else if(!write_status)   /* binary write skipped on failed open; the MPI_Recv above still ran so peers do not deadlock */
                                    {
                                        my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                                    }
                                }

                                n_for_this_task -= pc;
                            }
                        }


                        if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
                        {
                            if(All.SnapFormat == 3)
                            {
                                H5Dclose(hdf5_dataset);
                                H5Sclose(hdf5_dataspace_in_file);
                                H5Tclose(hdf5_datatype);
                            }
                        }

                    }
                }

                if(ThisTask == writeTask && !write_status)
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

            for(type = 5; type >= 0; type--) {if(header.npart[type] > 0) {H5Gclose(hdf5_grp[type]);}}
            H5Gclose(hdf5_headergrp);
            H5Fclose(hdf5_file);

        }
        else
        {
            if(!write_status) {fclose(fd);}   /* fd is NULL on a failed open; do not fclose(NULL) */
        }
    }
}





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

#if defined(CRFLUID_EVOLVE_SPECTRUM)
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_SPECIES}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Active_Species_ID_List", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, All.CR_species_ID_active_list); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Species_ID", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, All.CR_species_ID_in_bin); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Charge", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_global_charge_in_bin); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Min_Rigidity_GV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_global_min_rigidity_in_bin); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Med_Rigidity_GV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_global_rigidity_at_bin_center); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Max_Rigidity_GV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_global_max_rigidity_in_bin); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Max_Rigidity_GV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_global_max_rigidity_in_bin); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[2]={N_CR_PARTICLE_BINS,N_CR_SPECTRUM_LUT}; H5Sset_extent_simple(hdf5_dataspace, 2, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Energy_vs_Number_GlobalLookupTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_global_slope_lut); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[2]={N_CR_PARTICLE_SPECIES,N_CR_PARTICLE_SPECIES}; H5Sset_extent_simple(hdf5_dataspace, 2, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Primary_to_Secondary_Products_List", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, All.CR_secondary_species_listref); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Radioactive_Decay_Rate_persecond", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_rad_decay_coeff); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_CR_PARTICLE_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Total_Fragmentation_Rate_cmcubedpersecond", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_frag_coeff); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[2]={N_CR_PARTICLE_BINS,N_CR_PARTICLE_SPECIES}; H5Sset_extent_simple(hdf5_dataspace, 2, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Secondary_Product_TargetBins", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, All.CR_secondary_target_bin); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}

    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[2]={N_CR_PARTICLE_BINS,N_CR_PARTICLE_SPECIES}; H5Sset_extent_simple(hdf5_dataspace, 2, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "CR_Bin_Fragmentation_Coefficients_To_Secondaries", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.CR_frag_secondary_coeff); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
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

#ifdef SINK_PHOTONMOMENTUM
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "LEBRON_Sink_PhotonMomentum_Coupled_Fraction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Sink_Rad_MomentumFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef GRAIN_FLUID
    {int holder=GRAIN_PTYPES; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "GrainCR_Particle_Type_BitFlag", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#if !defined(PIC_MHD) || defined(GRAIN_FLUID_AND_PIC_BOTH_DEFINED)
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

#ifdef GALSF_FB_FIRE_RT_LOCALRP
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "LEBRON_Radiation_Pressure_Local_Momentum_Normalization_FromSSP", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.RP_Local_Momentum_Renormalization); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
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

#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SNe_Energy_Normalization_FromSSP", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SNe_Energy_Renormalization); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "StellarMassLoss_Rate_Normalization_FromSSP", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.StellarMassLoss_Rate_Renormalization); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "StellarMassLoss_Energy_Normalization_FromSSP", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.StellarMassLoss_Energy_Renormalization); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(CR_DYNAMICAL_INJECTION_IN_SNE)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "CosmicRay_SNeFraction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.CosmicRay_SNeFraction); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#ifdef GALSF_FB_FIRE_RT_HIIHEATING
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "HIIRegion_fLum_Coupled", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.HIIRegion_fLum_Coupled); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
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

#ifdef GALSF_ISMDUSTCHEM_MODEL
    {int holder=NUM_ISMDUSTCHEM_SPECIES; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "ISMDustChem_NumberOfSpecies", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    /* list out composition (element key and number of atoms per element) of silicates since this changes between versions */
    // key for each element
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES}; 
    H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "Silicates_Element_Key", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, All.ISMDustChem_SilicateMetallicityFieldIndexTable); 
    H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    // number of atoms per element
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES}; 
    H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
    hdf5_attribute = H5Acreate(handle, "Silicates_Element_Number", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.ISMDustChem_SilicateNumberOfAtomsTable);
    H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#ifdef GALSF_ISMDUSTCHEM_GRAINSIZEEVO
    {int holder=NUM_ISMDUSTCHEM_SIZE_BINS; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "ISMDustChem_Num_Grain_Size_Bins", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    {double tmp=UNIT_GRAIN_NUMBER; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "UnitGrainNumber", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &tmp); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    {double tmp=UNIT_GRAIN_LENGTH; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "UnitGrainLength_in_cm", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &tmp); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "ISMDustChem_Grain_Size_Min", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ISMDustChem_Grain_Size_Min); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "ISMDustChem_Grain_Size_Max", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.ISMDustChem_Grain_Size_Max); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#endif
#endif // METALS

#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_RT_FREQ_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
        hdf5_attribute = H5Acreate(handle, "Radiation_RHD_Min_Bin_Freq_in_eV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.RHD_bins_nu_min_ev); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
    {hdf5_dataspace = H5Screate(H5S_SIMPLE); hsize_t tmp_dim[1]={N_RT_FREQ_BINS}; H5Sset_extent_simple(hdf5_dataspace, 1, tmp_dim, NULL);
        hdf5_attribute = H5Acreate(handle, "Radiation_RHD_Max_Bin_Freq_in_eV", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, All.RHD_bins_nu_max_ev); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#endif

#if defined(RT_ISRF_BACKGROUND)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "InterstellarRadiationFieldStrength", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.InterstellarRadiationFieldStrength); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
    
#if defined(SINK_WIND_KICK) || defined(SINK_WIND_SPAWN)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Sink_accreted_fraction", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Sink_accreted_fraction); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Sink_outflow_velocity", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Sink_outflow_velocity); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(SINGLE_STAR_FB_JETS)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Sink_outflow_jetlaunchvelscaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Sink_outflow_jetlaunchvelscaling); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif

#if defined(SINK_COSMIC_RAYS)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Sink_CosmicRay_Injection_Efficiency", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Sink_CosmicRay_Injection_Efficiency); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
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

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
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

#ifdef SINK_PARTICLES
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SinkAccretionFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SinkAccretionFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SinkFeedbackFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SinkFeedbackFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedSinkMass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedSinkMass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SinkNgbFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SinkNgbFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SinkMaxAccretionRadius", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SinkMaxAccretionRadius); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SinkEddingtonFactor", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SinkEddingtonFactor); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SinkRadiativeEfficiency", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SinkRadiativeEfficiency); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#if defined(SINK_SEED_FROM_FOF) || defined(SINK_SEED_FROM_LOCALGAS)
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedSinkMassSigma", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedSinkMassSigma); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedSinkMinRedshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedSinkMinRedshift); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef SINK_SEED_FROM_LOCALGAS
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedSinkPerUnitMass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedSinkPerUnitMass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef SINK_ALPHADISK_ACCRETION
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "SeedReservoirMass", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.SeedReservoirMass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef SINK_SEED_FROM_FOF
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "MinFoFMassForNewSeed", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.MinFoFMassForNewSeed); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
#ifdef SINK_WIND_SPAWN
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Cell_Spawn_Mass_ratio", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Sink_outflow_particlemass); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#ifdef SINGLE_STAR_FB_WINDS
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Cell_Spawn_Mass_ratio_MS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Cell_Spawn_Mass_ratio_MS); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
#endif
    hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Sink_outflow_temperature", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Sink_outflow_temperature); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);
    {unsigned long long holder = (unsigned long long) All.SpawnedWindCellID; hdf5_dataspace = H5Screate(H5S_SCALAR); hdf5_attribute = H5Acreate(handle, "Spawned_Cell_ID", H5T_NATIVE_ULLONG, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_ULLONG, &holder); H5Aclose(hdf5_attribute); H5Sclose(hdf5_dataspace);}
#endif
#endif

}






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
            endrun(777);   /* soft bad-stop: my_fwrite I/O error; return the actual short count and drain at the next caller poll. */
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
        endrun(778);   /* soft bad-stop: my_fread I/O error; return the actual short count (reads <= requested into a caller-sized buffer, no OOB), drain at the read-turn poll. Caller-side partial-read guard is a tracked follow-up. */
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
