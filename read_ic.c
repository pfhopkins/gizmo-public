#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>


#include "allvars.h"
#include "proto.h"

/* This function reads initial conditions that are in the default file format
 * of Gadget, i.e. snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with restartflag==2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameterfile (except for changing the name of the IC file to the snapshot,
 * and ensuring the format tag matches it).  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasTemp>0 is given, the gas temperature will be initialzed to this
 * value assuming a mean colecular weight either corresponding to complete
 * neutrality, or full ionization.
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part (adding/removing read items and changing variable units as necessary)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void read_ic(char *fname)
{
    long i; int num_files, rest_files, ngroups, gr, filenr, primaryTask, lastTask, groupTaskIterator;
    double u_init, molecular_weight; char buf[500];

    CPU_Step[CPU_MISC] += measure_time();

#ifdef RESCALEVINI
    if(ThisTask == 0 && RestartFlag == 0) {fprintf(stdout, "Rescaling v_ini !\n"); fflush(stdout);}
#endif

    NumPart = 0;
    N_gas = 0;
    All.TotNumPart = 0;

    num_files = find_files(fname);

    rest_files = num_files;

    while(rest_files > NTask)
    {
        sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
        if(All.ICFormat == 3) {sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));}

        ngroups = NTask / All.NumFilesWrittenInParallel;
        if((NTask % All.NumFilesWrittenInParallel)) {ngroups++;}
        groupTaskIterator = (ThisTask / ngroups) * ngroups;

        for(gr = 0; gr < ngroups; gr++)
        {
            if(ThisTask == (groupTaskIterator + gr)) {read_file(buf, ThisTask, ThisTask);}	/* ok, it's this processor's turn */
            MPI_Barrier(MPI_COMM_WORLD);
        }
        rest_files -= NTask;
    }


    if(rest_files > 0)
    {
        distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &primaryTask, &lastTask);

        if(num_files > 1)
        {
            sprintf(buf, "%s.%d", fname, filenr);
            if(All.ICFormat == 3)
                sprintf(buf, "%s.%d.hdf5", fname, filenr);
        }
        else
        {
            sprintf(buf, "%s", fname);
            if(All.ICFormat == 3)
                sprintf(buf, "%s.hdf5", fname);
        }

        ngroups = rest_files / All.NumFilesWrittenInParallel;
        if((rest_files % All.NumFilesWrittenInParallel))
            ngroups++;

        for(gr = 0; gr < ngroups; gr++)
        {
            if((filenr / All.NumFilesWrittenInParallel) == gr) {read_file(buf, primaryTask, lastTask);}	/* ok, it's this processor's turn */
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }


    myfree(CommBuffer);


    if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
    {
        /* this makes sure that masses are initialized in the case that the mass-block is empty for this particle type */
        for(i = 0; i < NumPart; i++) {if(All.MassTable[P[i].Type] != 0) {P[i].Mass = All.MassTable[P[i].Type];}}
    }

    /* zero this out, since various operations in the code will want to change particle
     masses and keeping MassTable fixed won't allow that to happen */
    for(i=0;i<6;i++) All.MassTable[i]=0;


#if defined(BLACK_HOLES)
#if defined(BH_SWALLOWGAS) || defined(BH_BONDI) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK) || defined(BH_GRAVACCRETION) || defined(BH_GRAVCAPTURE_GAS) || defined(BH_SEED_FROM_LOCALGAS)
    if(RestartFlag == 0) {All.MassTable[5] = 0;}
#endif
#endif

#ifdef GALSF
    if(RestartFlag == 0)
    {
        if(All.MassTable[4] == 0 && All.MassTable[0] > 0)
        {
            All.MassTable[0] = 0;
            All.MassTable[4] = 0;
        }
    }
#endif

#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
    if(RestartFlag == 0)
    {
        All.MassTable[2] = 0;
        All.MassTable[3] = 0;
        All.MassTable[4] = 0;
    }
#endif

    u_init = All.InitGasTemp / ((GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS);

    molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* assume full ionization */
    if(All.InitGasTemp < 1.0e4) {molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);} /* assume neutral gas */
#if defined(COOL_LOW_TEMPERATURES) || defined(COOL_GRACKLE)
    if(All.InitGasTemp < 1.0e3 && All.ComovingIntegrationOn==0) {molecular_weight =  1. / ( HYDROGEN_MASSFRAC*0.5 + (1-HYDROGEN_MASSFRAC)/4. + 1./(16.+12.));} /* assume fully molecular [self-consistency requires cooling can handle this, and that this is intended to represent dense gas, not e.g. neutral early-universe gas] */
#endif

    u_init /= molecular_weight;

    All.InitGasU = u_init;

    if(RestartFlag == 0)
    {
        if(All.InitGasTemp > 0)
        {
            for(i = 0; i < N_gas; i++)
            {
                if(ThisTask == 0 && i == 0) // && SphP[i].InternalEnergy == 0)
                    {printf("Initializing u from InitGasTemp : InitGasTemp=%g InitGasU=%g MinEgySpec=%g SphP[0].InternalEnergy=%g\n",
                           All.InitGasTemp,All.InitGasU,All.MinEgySpec,SphP[i].InternalEnergy);}

                SphP[i].InternalEnergy = All.InitGasU;
            }
        }
    }

    for(i = 0; i < N_gas; i++) {SphP[i].InternalEnergyPred = SphP[i].InternalEnergy = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);}
    MPI_Barrier(MPI_COMM_WORLD);
    if(ThisTask == 0) {printf("Reading done. Total number of particles :  %d%09d\n\n", (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000)); fflush(stdout);}

    CPU_Step[CPU_SNAPSHOT] += measure_time();
}


/*! This function reads out the buffer that was filled with particle data.
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
    long n, k; MyInputFloat *fp; MyInputPosFloat *fp_pos; MyIDType *ip; int *ip_int; float *fp_single;
    fp = (MyInputFloat *) CommBuffer;
    fp_pos = (MyInputPosFloat *) CommBuffer;
    fp_single = (float *) CommBuffer;
    ip = (MyIDType *) CommBuffer;
    ip_int = (int *) CommBuffer;

    switch(blocknr)
    {
        case IO_POS:		/* positions */
            for(n = 0; n < pc; n++)
                for(k = 0; k < 3; k++)
                {
                    P[offset + n].Pos[k] = *fp_pos++;
                    // P[offset + n].Pos[k] += 0.5*All.BoxSize; /* manually turn on for some ICs */
                }

            for(n = 0; n < pc; n++) {P[offset + n].Type = type;}	/* initialize type here as well */
            break;

        case IO_VEL:		/* velocities */
            for(n = 0; n < pc; n++)
            {
                for(k = 0; k < 3; k++)
                {
#ifdef RESCALEVINI
                    /* scaling v to use same IC's for different cosmologies */
                    if(RestartFlag == 0) {P[offset + n].Vel[k] = (*fp++) * All.VelIniScale;} else {P[offset + n].Vel[k] = *fp++;}
#else
                    P[offset + n].Vel[k] = *fp++;
#endif
                }
            }
            break;

        case IO_ID:		/* particle ID */
            for(n = 0; n < pc; n++) {P[offset + n].ID = *ip++;}
            break;


        case IO_CHILD_ID:		// particle child ID //
            if(RestartFlag == 2)
            {
                for(n = 0; n < pc; n++) {P[offset + n].ID_child_number = *ip++;}
            }
            break;

        case IO_GENERATION_ID:		// particle generation ID //
            if(RestartFlag == 2) {for(n = 0; n < pc; n++) {P[offset + n].ID_generation = *ip++;}}
            break;

        case IO_MASS:		/* particle mass */
            for(n = 0; n < pc; n++) {P[offset + n].Mass = *fp++;}
            break;


        case IO_SHEET_ORIENTATION:	/* initial particle sheet orientation */
#if defined(GDE_DISTORTIONTENSOR) && defined(GDE_READIC)
            for(n = 0; n < pc; n++)
            {
#ifndef GDE_LEAN
                P[offset + n].V_matrix[0][0] = *fp++;
                P[offset + n].V_matrix[0][1] = *fp++;
                P[offset + n].V_matrix[0][2] = *fp++;
                P[offset + n].V_matrix[1][0] = *fp++;
                P[offset + n].V_matrix[1][1] = *fp++;
                P[offset + n].V_matrix[1][2] = *fp++;
                P[offset + n].V_matrix[2][0] = *fp++;
                P[offset + n].V_matrix[2][1] = *fp++;
                P[offset + n].V_matrix[2][2] = *fp++;
#else
                *fp+=8;
#endif
            }
#endif
            break;

        case IO_INIT_DENSITY:	/* initial stream density */
#if defined(GDE_DISTORTIONTENSOR) && defined(GDE_READIC)
            for(n = 0; n < pc; n++) {GDE_INITDENSITY(offset + n) = *fp++;}
            break;
#endif

        case IO_CAUSTIC_COUNTER:	/* initial caustic counter */
#if defined(GDE_DISTORTIONTENSOR) && defined(GDE_READIC)
            for(n = 0; n < pc; n++) {P[offset + n].caustic_counter = *fp++;}
            break;
#endif

        case IO_U:			/* temperature */
            for(n = 0; n < pc; n++) {SphP[offset + n].InternalEnergy = *fp++;}
            break;

        case IO_RHO:		/* density */
            for(n = 0; n < pc; n++) {SphP[offset + n].Density = *fp++;}
            break;

        case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(RT_CHEM_PHOTOION)
#ifndef CHIMES
            for(n = 0; n < pc; n++) {SphP[offset + n].Ne = *fp++;}
#endif
#endif
            break;


        case IO_HSML:		/* gas kernel length */
            for(n = 0; n < pc; n++) {PPP[offset + n].Hsml = *fp++;}
            break;

        case IO_DELAYTIME:
#ifdef GALSF_SUBGRID_WINDS
            for(n = 0; n < pc; n++) {SphP[offset + n].DelayTime = *fp++;}
#endif
            break;

        case IO_AGE:		/* Age of stars */
#ifdef GALSF
            for(n = 0; n < pc; n++) {P[offset + n].StellarAge = *fp++;}
#endif
            break;

        case IO_GRAINSIZE:
#ifdef GRAIN_FLUID
            for(n = 0; n < pc; n++) {P[offset + n].Grain_Size = *fp++;}
#endif
            break;

        case IO_GRAINTYPE:
#if defined(PIC_MHD)
            for(n = 0; n < pc; n++) {P[offset + n].MHD_PIC_SubType = *ip_int++;}
#endif
            break;

        case IO_Z:			/* Gas and star metallicity */
#ifdef METALS
            for(n = 0; n < pc; n++) {
                int nmax=NUM_METAL_SPECIES;
                if(RestartFlag==2 && All.ICFormat==3 && header.flag_metals<NUM_METAL_SPECIES && header.flag_metals>0) {nmax=header.flag_metals;} // special clause to catch cases where read-in snapshot did not use all the metals fields we want to read now
                for(k=0;k<nmax;k++) {P[offset + n].Metallicity[k] = *fp++;} // normal read-in
                if(nmax<NUM_METAL_SPECIES) {for(k=nmax;k<NUM_METAL_SPECIES;k++) {P[offset + n].Metallicity[k]=0;}} // any extra fields zero'd
            }
#endif
            break;

        case IO_BFLD:		/* Magnetic field */
#ifdef MAGNETIC
            for(n = 0; n < pc; n++)
            {
                for(k = 0; k < 3; k++) {SphP[offset + n].BPred[k] = *fp++;}
                SphP[offset + n].divB = 0;
#ifdef DIVBCLEANING_DEDNER
                SphP[offset + n].Phi = 0;
                SphP[offset + n].PhiPred = 0;
#endif
            }
#endif
            break;

        case IO_BHMASS:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; n++) {P[offset + n].BH_Mass = *fp++;}
#endif
            break;

        case IO_BHDUSTMASS:
#if defined(BLACK_HOLES) && defined(GRAIN_FLUID)
            for(n = 0; n < pc; n++) {P[offset + n].BH_Dust_Mass = *fp++;}
#endif
            break;

        case IO_BH_DIST:
            break;

        case IO_BH_ANGMOM:
#ifdef BH_FOLLOW_ACCRETED_ANGMOM
            for(n = 0; n < pc; n++) {for(k = 0; k < 3; k++) {P[offset + n].BH_Specific_AngMom[k] = *fp++;}}
#endif
            break;

        case IO_BHMASSALPHA:
#ifdef BH_ALPHADISK_ACCRETION
            for(n = 0; n < pc; n++) {P[offset + n].BH_Mass_AlphaDisk = *fp++;}
#endif
            break;

        case IO_BHMDOT:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; n++) {P[offset + n].BH_Mdot = *fp++;}
#endif
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
            for(n = 0; n < pc; n++) {P[offset + n].BH_CountProgs = *ip_int++;}
#endif
            break;

        case IO_EOSTEMP:
#ifdef EOS_CARRIES_TEMPERATURE
            for(n = 0; n < pc; n++) {SphP[offset + n].Temperature = *fp++;}
#endif
            break;

        case IO_EOSABAR:
#ifdef EOS_CARRIES_ABAR
            for(n = 0; n < pc; n++) {SphP[offset + n].Abar = *fp++;}
#endif
            break;

        case IO_EOSCOMP:
#ifdef EOS_TILLOTSON
            for(n = 0; n < pc; n++) {SphP[offset + n].CompositionType = *ip_int++;}
#endif
            break;

        case IO_EOSYE:
#ifdef EOS_CARRIES_YE
            for(n = 0; n < pc; n++) {SphP[offset + n].Ye = *fp++;}
#endif
            break;

        case IO_PARTVEL:
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==1)||(HYDRO_FIX_MESH_MOTION==2)||(HYDRO_FIX_MESH_MOTION==3))
            for(n = 0; n < pc; n++) {for(k = 0; k < 3; k++) {SphP[offset + n].ParticleVel[k] = *fp++;}}
#endif
            break;


        case IO_RADGAMMA:
#ifdef RADTRANSFER
            for(n = 0; n < pc; n++) {for(k = 0; k < N_RT_FREQ_BINS; k++) {SphP[offset + n].Rad_E_gamma[k] = *fp++;}}
#endif
            break;

            /* adaptive softening parameters */
        case IO_AGS_SOFT:
#if defined (AGS_HSML_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTGRAVSOFT)
            for(n = 0; n < pc; n++) {PPP[offset + n].AGS_Hsml = *fp++;}
#endif
            break;

        case IO_AGS_ZETA:
#if defined (AGS_HSML_CALCULATION_IS_ACTIVE) && defined(AGS_OUTPUTZETA)
            for(n = 0; n < pc; n++) {PPPZ[offset + n].AGS_Zeta = *fp++;}
#endif
            break;

        case IO_CHIMES_ABUNDANCES:
#if defined(CHIMES) && !defined(CHIMES_INITIALISE_IN_EQM)
            for (n = 0; n < pc; n++)
            {
    	        allocate_gas_abundances_memory(&(ChimesGasVars[offset + n]), &ChimesGlobalVars);
	            for (k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++) {ChimesGasVars[offset + n].abundances[k] = (ChimesFloat) (*fp++);}
#ifdef CHIMES_TURB_DIFF_IONS
                chimes_update_turbulent_abundances(n, 1);
#endif
            }
#endif
            break;

        case IO_COSMICRAY_ENERGY:
#ifdef COSMIC_RAY_FLUID
#ifdef CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART
            for(n = 0; n < pc; n++) {SphP[offset + n].CosmicRayEnergy[0] = *fp++;}
#else
            for(n = 0; n < pc; n++) {for(k=0; k<N_CR_PARTICLE_BINS; k++) {SphP[offset + n].CosmicRayEnergy[k] = *fp++;}}
#endif
#endif
            break;

        case IO_COSMICRAY_SLOPES:
            break;

        case IO_COSMICRAY_ALFVEN:
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(n = 0; n < pc; n++) {
                int k2; for(k=0;k<2;k++) {for(k2=0;k2<N_CR_PARTICLE_BINS;k2++) {
                        SphP[offset + n].CosmicRayAlfvenEnergy[k2][k] = fp[N_CR_PARTICLE_BINS*k + k2];}}
                fp += 2*N_CR_PARTICLE_BINS;
            }
#endif
            break;

        case IO_OSTAR:
#ifdef GALSF_SFR_IMF_SAMPLING
             for(n = 0; n < pc; n++) {P[offset + n].IMF_NumMassiveStars = *fp++;}
#endif
            break;
            
        case IO_UNSPMASS:
#if defined(BH_WIND_SPAWN) && defined(BH_DEBUG_SPAWN_JET_TEST)
             for(n = 0; n < pc; n++) {P[offset + n].unspawned_wind_mass = *fp++;}
#endif
            break; 
            
        case IO_IDEN:
#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) && defined(BH_DEBUG_SPAWN_JET_TEST)
             for(n = 0; n < pc; n++) {SphP[offset + n].IniDen = *fp++;}
#endif
            break;
            
        case IO_INIB:        
#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) && defined(BH_DEBUG_SPAWN_JET_TEST)
             for(n = 0; n < pc; n++) {for(k=0;k<3;k++) {SphP[offset + n].IniB[k] = *fp++;}}
#endif      
            break;  
            
        case IO_TURB_DYNAMIC_COEFF:
#ifdef TURB_DIFF_DYNAMIC
            for (n = 0; n < pc; n++) {SphP[offset + n].TD_DynDiffCoeff = *fp++;}
#endif
            break;

        case IO_SINKRAD:
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
            for(n = 0; n < pc; n++) {P[offset + n].SinkRadius = *fp++;}
#endif
            break;
            
        case IO_MOLECULARFRACTION:
#if defined(COOL_MOLECFRAC_NONEQM) & !defined(IO_MOLECFRAC_NOT_IN_ICFILE)
            for (n = 0; n < pc; n++) {SphP[offset + n].MolecularMassFraction_perNeutralH = *fp++;}
#endif
            break;


        /* the other input fields (if present) are not needed to define the
             initial conditions of the code */

        case IO_COSMICRAY_KAPPA:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_PSI_RE:
        case IO_AGS_PSI_IM:
        case IO_EOSCS:
        case IO_EOS_STRESS_TENSOR:
        case IO_CBE_MOMENTS:
        case IO_SFR:
        case IO_POT:
        case IO_ACCEL:
        case IO_DTENTR:
        case IO_RAD_ACCEL:
        case IO_GDE_DISTORTIONTENSOR:
        case IO_CRATE:
        case IO_HRATE:
        case IO_NHRATE:
        case IO_HHRATE:
        case IO_MCRATE:
        case IO_HeII:
        case IO_HeI:
        case IO_HII:
        case IO_NH:
        case IO_TSTP:
        case IO_IMF:
        case IO_DIVB:
        case IO_ABVC:
        case IO_COOLRATE:
        case IO_AMDC:
        case IO_PHI:
        case IO_GRADPHI:
        case IO_TIDALTENSORPS:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_ANNIHILATION_RADIATION:
        case IO_PRESSURE:
        case IO_EDDINGTON_TENSOR:
        case IO_LAST_CAUSTIC:
        case IO_HSMS:
        case IO_ACRB:
        case IO_RAD_FLUX:
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
        case IO_TURB_DIFF_COEFF:
        case IO_DYNERROR:
        case IO_DYNERRORDEFAULT:
        case IO_VDIV:
        case IO_VORT:
        case IO_CHIMES_MU:
        case IO_CHIMES_REDUCED:
        case IO_CHIMES_NH:
        case IO_CHIMES_STAR_SIGMA:
        case IO_DENS_AROUND_STAR:
        case IO_DELAY_TIME_HII:
        case IO_CHIMES_FLUX_G0:
        case IO_CHIMES_FLUX_ION:
            break;

        case IO_LASTENTRY:
            endrun(220);
            break;
    }
}



/*! This function reads a snapshot file and distributes the data it contains
 *  to tasks 'readTask' to 'lastTask'.
 */
void read_file(char *fname, int readTask, int lastTask)
{
    size_t blockmaxlen;
    long long i, n_in_file, n_for_this_task, ntask, pc, offset = 0, task, nall, nread, nstart, npart;
    int blksize1, blksize2, type, bnr, bytes_per_blockelement, nextblock, typelist[6];
    MPI_Status status;
    FILE *fd = 0;
    char label[4], buf[500];
    enum iofields blocknr;
    size_t bytes;

#ifdef HAVE_HDF5
    int rank, pcsum;
    hid_t hdf5_file = 0, hdf5_grp[6], hdf5_dataspace_in_file;
    hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory, hdf5_dataset;
    hsize_t dims[2], count[2], start[2];
#endif

#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}

    if(ThisTask == readTask)
    {
        if(All.ICFormat == 1 || All.ICFormat == 2)
        {
            if(!(fd = fopen(fname, "r")))
            {
                printf("can't open file `%s' for reading initial conditions.\n", fname);
                endrun(123);
            }


            if(All.ICFormat == 2)
            {
                SKIP;
                my_fread(&label, sizeof(char), 4, fd);
                my_fread(&nextblock, sizeof(int), 1, fd);
                printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3], nextblock);
                SKIP2;
            }

            SKIP;
            my_fread(&header, sizeof(header), 1, fd);
            SKIP2;

            if(blksize1 != 256 || blksize2 != 256)
            {
                printf("incorrect header format\n");
                fflush(stdout);
                endrun(890);
                /* Probable error is wrong size of fill[] in header file. Needs to be 256 bytes in total. */
            }
        }


#ifdef HAVE_HDF5
        if(All.ICFormat == 3)
        {
            read_header_attributes_in_hdf5(fname);
            hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
            for(type = 0; type < 6; type++)
            {
                if(header.npart[type] > 0)
                {
                    sprintf(buf, "/PartType%d", type);
                    hdf5_grp[type] = H5Gopen(hdf5_file, buf);
                }
            }
        }
#endif

        for(task = readTask + 1; task <= lastTask; task++)
        {
            MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
        }

    }
    else
    {
        MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);
    }

#ifdef INPUT_IN_DOUBLEPRECISION
    if(header.flag_doubleprecision == 0)
    {
        if(ThisTask == 0) {printf("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");}
        endrun(11);
    }
#else
    if(header.flag_doubleprecision)
    {
        if(ThisTask == 0) {printf("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");}
        endrun(10);
    }
#endif


    if(All.TotNumPart == 0)
    {
        if(header.num_files <= 1)
            for(i = 0; i < 6; i++)
            {
                header.npartTotal[i] = header.npart[i];
                header.npartTotalHighWord[i] = 0;
            }

        All.TotN_gas = header.npartTotal[0] + (((long long) header.npartTotalHighWord[0]) << 32);

        for(i = 0, All.TotNumPart = 0; i < 6; i++)
        {
            All.TotNumPart += header.npartTotal[i];
            All.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
        }


        for(i = 0; i < 6; i++) {All.MassTable[i] = header.mass[i];}

        All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));
        All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));	/* sets the maximum number of particles that may reside on a processor */
#ifdef ALLOW_IMBALANCED_GASPARTICLELOAD
        All.MaxPartSph = All.MaxPart; // PFH: increasing All.MaxPartSph according to this line can allow better load-balancing in some cases. however it leads to more memory problems
        // (PFH: needed to revert the change -- i.e. INCLUDE the line above: commenting it out, while it improved memory useage, causes some instability in the domain decomposition for
        //   sufficiently irregular trees. overall more stable behavior with the 'buffer', albeit at the expense of memory )
#endif


        allocate_memory();

        size_t MyBufferSize = All.BufferSize;
        if(!(CommBuffer = mymalloc("CommBuffer", bytes = MyBufferSize * 1024 * 1024)))
        {
            printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
            endrun(2);
        }

        if(RestartFlag >= 2)
        {
            All.Time = All.TimeBegin = header.time;
            set_cosmo_factors_for_current_time();
        }

    }

    if(ThisTask == readTask)
    {
        for(i = 0, n_in_file = 0; i < 6; i++) {n_in_file += header.npart[i];}

        printf("\nReading file `%s' on task=%d (contains %lld particles.)\n"
               " ..distributing this file to tasks %d-%d\n"
               "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 1 (halo):  %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 2 (disk):  %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 3 (bulge): %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 4 (stars): %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 5 (bndry): %8d  (tot=%6d%09d) masstab=%g\n\n", fname, ThisTask, n_in_file, readTask,
               lastTask, header.npart[0], (int) (header.npartTotal[0] / 1000000000),
               (int) (header.npartTotal[0] % 1000000000), All.MassTable[0], header.npart[1],
               (int) (header.npartTotal[1] / 1000000000), (int) (header.npartTotal[1] % 1000000000),
               All.MassTable[1], header.npart[2], (int) (header.npartTotal[2] / 1000000000),
               (int) (header.npartTotal[2] % 1000000000), All.MassTable[2], header.npart[3],
               (int) (header.npartTotal[3] / 1000000000), (int) (header.npartTotal[3] % 1000000000),
               All.MassTable[3], header.npart[4], (int) (header.npartTotal[4] / 1000000000),
               (int) (header.npartTotal[4] % 1000000000), All.MassTable[4], header.npart[5],
               (int) (header.npartTotal[5] / 1000000000), (int) (header.npartTotal[5] % 1000000000),
               All.MassTable[5]);
        fflush(stdout);
    }


    ntask = lastTask - readTask + 1;


    /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

    for(type = 0, nall = 0; type < 6; type++)
    {
        n_in_file = header.npart[type];
        n_for_this_task = n_in_file / ntask;
        if((ThisTask - readTask) < (n_in_file % ntask)) {n_for_this_task++;}


        if(type == 0)
        {
            if(N_gas + n_for_this_task > All.MaxPartSph)
            {
                printf("Not enough space on task=%d for SPH particles (space for %d, need at least %lld)\n", ThisTask, All.MaxPartSph, N_gas + n_for_this_task);
                fflush(stdout);
                endrun(172);
            }
        }

        nall += n_for_this_task;
    }

    if(NumPart + nall > All.MaxPart)
    {
        printf("Not enough space on task=%d (space for %d, need at least %lld)\n", ThisTask, All.MaxPart, NumPart + nall);
        fflush(stdout);
        endrun(173);
    }

    memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(struct particle_data));
    nstart = N_gas;


    for(bnr = 0; bnr < 1000; bnr++)
    {
        blocknr = (enum iofields) bnr;
        if(blocknr == IO_LASTENTRY) {break;}
        if(RestartFlag == 5 && blocknr > IO_MASS) {continue;}	/* if we only do power spectra, we don't need to read other blocks beyond the mass */

        if(blockpresent(blocknr))
        {
                /* blocks only for restartflag == 0 */
                if(RestartFlag == 0 && blocknr > IO_U
                   && blocknr != IO_BFLD
#ifdef INPUT_READ_HSML
                   && blocknr != IO_HSML
#endif
#ifdef EOS_CARRIES_TEMPERATURE
                   && blocknr != IO_EOSTEMP
#endif
#ifdef EOS_CARRIES_ABAR
                   && blocknr != IO_EOSABAR
#endif
#ifdef EOS_CARRIES_YE
                   && blocknr != IO_EOSYE
#endif
#ifdef EOS_TILLOTSON
                   && blocknr != IO_EOSCOMP
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==1)||(HYDRO_FIX_MESH_MOTION==2)||(HYDRO_FIX_MESH_MOTION==3))
                   && blocknr != IO_PARTVEL
#endif
#if defined(BH_GRAVCAPTURE_FIXEDSINKRADIUS)
                   && blocknr != IO_SINKRAD
#endif
#if defined(CHIMES) && !defined(CHIMES_INITIALISE_IN_EQM)
                   && blocknr != IO_CHIMES_ABUNDANCES
#endif
#ifdef PIC_MHD
                   && blocknr != IO_GRAINTYPE
#endif
                   )
#if defined(GDE_DISTORTIONTENSOR) && defined(GDE_READIC)
                    if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_SHEET_ORIENTATION))
                        if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_INIT_DENSITY))
                            if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_CAUSTIC_COUNTER))
#endif
                                continue;	/* ignore all other blocks in initial conditions */


            if(RestartFlag == 0 && (blocknr == IO_GENERATION_ID || blocknr == IO_CHILD_ID)) {continue;}
#if defined(NO_CHILD_IDS_IN_ICS) || defined(ASSIGN_NEW_IDS)
            if(blocknr == IO_GENERATION_ID || blocknr == IO_CHILD_ID) {continue;}
#endif
            if((RestartFlag == 0) && (All.InitGasTemp > 0) && (blocknr == IO_U)) {continue;}


#ifdef MHD_B_SET_IN_PARAMS
            if(RestartFlag == 0 && blocknr == IO_BFLD) {continue;}
#endif

#ifdef SUBFIND
            if(RestartFlag == 2 && blocknr == IO_HSMS) {continue;}
#endif

#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
#ifndef AGS_OUTPUTGRAVSOFT
            if(blocknr == IO_AGS_SOFT) {continue;}
#endif
#ifndef AGS_OUTPUTZETA
            if(blocknr == IO_AGS_ZETA) {continue;}
#endif
#endif
            
#ifdef CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART
            if(RestartFlag == 2 && blocknr == IO_COSMICRAY_SLOPES) {continue;}
#if (CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART==2)
            if(RestartFlag == 2 && blocknr == IO_COSMICRAY_ENERGY) {continue;}
#endif
#endif

#if !defined(RADTRANSFER)
            if(RestartFlag == 2 && blocknr == IO_RADGAMMA) {continue;}
#endif

#if defined(IO_MOLECFRAC_NOT_IN_ICFILE)
            if(RestartFlag == 2 && blocknr == IO_MOLECULARFRACTION) {continue;}
#endif

            
            if(blocknr == IO_HSMS) {continue;}

#ifdef TURB_DIFF_DYNAMIC
            if (RestartFlag == 0 && blocknr == IO_TURB_DYNAMIC_COEFF) {continue;}
#endif

            if(ThisTask == readTask)
            {
                get_dataset_name(blocknr, buf);
                printf("reading block %d (%s)...\n", bnr, buf);
            }

            bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 1);
            
#if (CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART==1)
            if(RestartFlag == 2 && blocknr == IO_COSMICRAY_ENERGY) {bytes_per_blockelement = (1) * sizeof(MyInputFloat);}
#endif
#ifdef METALS /* some trickery here to enable snapshot-restarts from runs with different numbers of metal species */
            if(blocknr==IO_Z && RestartFlag==2 && All.ICFormat==3 && header.flag_metals<NUM_METAL_SPECIES && header.flag_metals>0) {bytes_per_blockelement = (header.flag_metals) * sizeof(MyInputFloat);}
#endif

            size_t MyBufferSize = All.BufferSize;
            blockmaxlen = (size_t) ((MyBufferSize * 1024 * 1024) / bytes_per_blockelement);
            npart = get_particles_in_block(blocknr, &typelist[0]);

            if(npart > 0)
            {
                    if(ThisTask == readTask)
                    {
                        if(All.ICFormat == 2)
                        {
                            get_Tab_IO_Label(blocknr, label);
                            find_block(label, fd);
                        }

                        if(All.ICFormat == 1 || All.ICFormat == 2) {
                            SKIP;
                            if (blksize1 == 0) { /* workaround for MUSIC ICs */
                              SKIP2;
                              SKIP;
                            }
                        }
                    }

                for(type = 0, offset = 0, nread = 0; type < 6; type++)
                {
                    n_in_file = header.npart[type];
#ifdef HAVE_HDF5
                    pcsum = 0;
#endif
                    if(typelist[type] == 0)
                    {
                        n_for_this_task = n_in_file / ntask;
                        if((ThisTask - readTask) < (n_in_file % ntask)) {n_for_this_task++;}

                        offset += n_for_this_task;
                    }
                    else
                    {
                        for(task = readTask; task <= lastTask; task++)
                        {
                            n_for_this_task = n_in_file / ntask;
                            if((task - readTask) < (n_in_file % ntask)) {n_for_this_task++;}

                            if(task == ThisTask)
                                if(NumPart + n_for_this_task > All.MaxPart)
                                {
                                    printf("too many particles. %d %lld %d\n", NumPart, n_for_this_task, All.MaxPart);
                                    endrun(1313);
                                }


                            do
                            {
                                pc = n_for_this_task;
                                if(pc > (int)blockmaxlen) {pc = blockmaxlen;}

                                if(ThisTask == readTask)
                                {
                                    if(All.ICFormat == 1 || All.ICFormat == 2)
                                    {
                                            my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
                                            nread += pc;
                                    }

#ifdef HAVE_HDF5
                                    if(All.ICFormat == 3 && pc > 0)
                                    {
                                        get_dataset_name(blocknr, buf);
                                        hdf5_dataset = H5Dopen(hdf5_grp[type], buf);

                                        dims[0] = header.npart[type];
                                        dims[1] = get_values_per_blockelement(blocknr);
#if (CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART==1)
                                        if(RestartFlag == 2 && blocknr == IO_COSMICRAY_ENERGY) {dims[1] = 1;}
#endif
#ifdef METALS /* some trickery here to enable snapshot-restarts from runs with different numbers of metal species */
                                        if(blocknr==IO_Z && RestartFlag==2 && All.ICFormat==3 && header.flag_metals<NUM_METAL_SPECIES && header.flag_metals>0) {dims[1] = header.flag_metals;}
#endif
                                        if(dims[1] == 1) {rank = 1;} else {rank = 2;}
                                        hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

                                        dims[0] = pc;
                                        hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

                                        start[0] = pcsum;
                                        start[1] = 0;

                                        count[0] = pc;
                                        count[1] = get_values_per_blockelement(blocknr);
#if (CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART==1)
                                        if(RestartFlag == 2 && blocknr == IO_COSMICRAY_ENERGY) {count[1] = 1;}
#endif
#ifdef METALS /* some trickery here to enable snapshot-restarts from runs with different numbers of metal species */
                                        if(blocknr==IO_Z && RestartFlag==2 && All.ICFormat==3 && header.flag_metals<NUM_METAL_SPECIES && header.flag_metals>0) {count[1] = header.flag_metals;}
#endif
                                        pcsum += pc;

                                        H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                        switch(get_datatype_in_block(blocknr))
                                        {
                                            case 0:
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
                                                break;
                                                
                                            case 1:
#ifdef INPUT_IN_DOUBLEPRECISION
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                                break;
                                                
                                            case 2:
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
                                                break;

                                            case 3:
#if defined(INPUT_POSITIONS_IN_DOUBLE)
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                                break;
                                        }

                                        H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);
                                        H5Tclose(hdf5_datatype);
                                        H5Sclose(hdf5_dataspace_in_memory);
                                        H5Sclose(hdf5_dataspace_in_file);
                                        H5Dclose(hdf5_dataset);
                                    }
#endif
                                }

                                if(ThisTask == readTask && task != readTask && pc > 0) {MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD);}

                                if(ThisTask != readTask && task == ThisTask && pc > 0) {MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask, TAG_PDATA, MPI_COMM_WORLD, &status);}

                                if(ThisTask == task)
                                {
                                    empty_read_buffer(blocknr, nstart + offset, pc, type);
                                    offset += pc;
                                }

                                n_for_this_task -= pc;
                            }
                            while(n_for_this_task > 0);
                        }
                    }
                }

                if(ThisTask == readTask)
                {
                        if(All.ICFormat == 1 || All.ICFormat == 2)
                        {
                            SKIP2;

                            if(blksize1 != blksize2)
                            {
                                printf("incorrect block-sizes detected!\n");
                                printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, bnr, blksize1, blksize2);
                                if(blocknr == IO_ID)
                                {
                                    printf
                                    ("Possible mismatch of 32bit and 64bit ID's in IC file and GIZMO compilation !\n");
                                }
                                fflush(stdout);
                                endrun(1889);
                            }
                        }
                }
            }
        }
    }

    for(type = 0; type < 6; type++)
    {
        n_in_file = header.npart[type];
        n_for_this_task = n_in_file / ntask;
        if((ThisTask - readTask) < (n_in_file % ntask)) {n_for_this_task++;}
        NumPart += n_for_this_task;
        if(type == 0) {N_gas += n_for_this_task;}
    }

    if(ThisTask == readTask)
    {
        if(All.ICFormat == 1 || All.ICFormat == 2) {fclose(fd);}
#ifdef HAVE_HDF5
        if(All.ICFormat == 3)
        {
            for(type = 5; type >= 0; type--) {if(header.npart[type] > 0) {H5Gclose(hdf5_grp[type]);}}
            H5Fclose(hdf5_file);
        }
#endif
    }

}



/*! This function determines on how many files a given snapshot is distributed.
 */
int find_files(char *fname)
{
    FILE *fd; char buf[200], buf1[200]; int dummy;

    sprintf(buf, "%s.%d", fname, 0);
    sprintf(buf1, "%s", fname);

    if(All.ICFormat == 3)
    {
        sprintf(buf, "%s.%d.hdf5", fname, 0);
        sprintf(buf1, "%s.hdf5", fname);
    }

#ifndef  HAVE_HDF5
    if(All.ICFormat == 3)
    {
        if(ThisTask == 0) {printf("Code wasn't compiled with HDF5 support enabled!\n");}
        endrun(0);
    }
#endif

    header.num_files = 0;

    if(ThisTask == 0)
    {
        if((fd = fopen(buf, "r")))
        {
            if(All.ICFormat == 1 || All.ICFormat == 2)
            {
                if(All.ICFormat == 2)
                {
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                }

                my_fread(&dummy, sizeof(dummy), 1, fd);
                my_fread(&header, sizeof(header), 1, fd);
                my_fread(&dummy, sizeof(dummy), 1, fd);
            }
            fclose(fd);
#ifdef HAVE_HDF5
            if(All.ICFormat == 3) {read_header_attributes_in_hdf5(buf);}
#endif
        }
    }

    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

    if(header.num_files > 0) {return header.num_files;}

    if(ThisTask == 0)
    {
        if((fd = fopen(buf1, "r")))
        {
            if(All.ICFormat == 1 || All.ICFormat == 2)
            {
                if(All.ICFormat == 2)
                {
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                }

                my_fread(&dummy, sizeof(dummy), 1, fd);
                my_fread(&header, sizeof(header), 1, fd);
                my_fread(&dummy, sizeof(dummy), 1, fd);
            }
            fclose(fd);
#ifdef HAVE_HDF5
            if(All.ICFormat == 3) {read_header_attributes_in_hdf5(buf1);}
#endif
            header.num_files = 1;
        }
    }

    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

    if(header.num_files > 0) {return header.num_files;}

    if(ThisTask == 0)
    {
        printf("\nCan't find initial conditions file.");
        printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
        fflush(stdout);
    }

    endrun(0);
    return 0;
}



/*! This function assigns a certain number of files to processors, such that
 *  each processor is exactly assigned to one file, and the number of cpus per
 *  file is as homogenous as possible. The number of files may at most be
 *  equal to the number of processors.
 */
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *primary_taskID, int *last)
{
    int ntask, filesleft, filesright, tasksleft;

    if(nfiles > 1)
    {
        ntask = lasttask - firsttask + 1;

        filesleft = (int) ((((double) (ntask / 2)) / ntask) * nfiles);
        if(filesleft <= 0)
            filesleft = 1;
        if(filesleft >= nfiles)
            filesleft = nfiles - 1;

        filesright = nfiles - filesleft;

        tasksleft = ntask / 2;

        distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, primary_taskID, last);
        distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr, primary_taskID, last);
    }
    else
    {
        if(ThisTask >= firsttask && ThisTask <= lasttask)
        {
            *filenr = firstfile;
            *primary_taskID = firsttask;
            *last = lasttask;
        }
    }
}



#ifdef HAVE_HDF5
void read_header_attributes_in_hdf5(char *fname)
{
    hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;

    hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
    H5Aclose(hdf5_attribute);

#ifdef METALS /* some trickery here to enable snapshot-restarts from runs with different numbers of metal species */
    if(RestartFlag==2)
    {
        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Metals");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
        H5Aclose(hdf5_attribute);
    }
#endif
    
    H5Gclose(hdf5_headergrp);
    H5Fclose(hdf5_file);
}
#endif






/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
void find_block(char *label, FILE * fd)
{
    unsigned int blocksize = 0, blksize;
    char blocklabel[5] = { "    " };

#define FBSKIP  {my_fread(&blksize,sizeof(int),1,fd);}

    rewind(fd);

    while(!feof(fd) && blocksize == 0)
    {
        FBSKIP;
        if(blksize != 8)
        {
            printf("Incorrect Format (blksize=%u)!\n", blksize);
            endrun(1891);
        }
        else
        {
            my_fread(blocklabel, 4 * sizeof(char), 1, fd);
            my_fread(&blocksize, sizeof(int), 1, fd);
            /*
             printf("Searching <%c%c%c%c>, found Block <%s> with %d bytes\n",
             label[0],label[1],label[2],label[3],blocklabel,blocksize);
             */
            FBSKIP;
            if(strncmp(label, blocklabel, 4) != 0)
            {
                fseek(fd, blocksize, 1);
                blocksize = 0;
            }
        }
    }
    if(feof(fd))
    {
        printf("Block '%c%c%c%c' not found !\n", label[0], label[1], label[2], label[3]);
        fflush(stdout);
        endrun(1890);
    }
}
