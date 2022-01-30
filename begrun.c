#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"
#include "kernel.h"


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */

/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO. The modifications
 * mostly center on added functionality for new modules, elimination of unnecessary
 * variables, implementing the DEVELOPER_MODE options, and re-organizing the read order
 * to allow easier manipulation on restarts.
 */



void begrun(void)
{
  struct global_data_all_processes all;
#ifdef _OPENMP
  int tid;
#endif
  if(ThisTask == 0)
    {
     printf("Running on %d MPI tasks.\n", NTask);
#ifdef _OPENMP
#pragma omp parallel private(tid)
      {
#pragma omp master
          printf("Using %d OpenMP threads\n", omp_get_num_threads());
          tid = omp_get_thread_num();
      }
#endif

      printf("\nSize of particle structure       %d  [bytes]\n", (int) sizeof(struct particle_data));
      printf("Size of hydro-cell structure   %d  [bytes]\n\n", (int) sizeof(struct sph_particle_data));

    }

#ifdef CHIMES_TURB_DIFF_IONS
  // Check that TURB_DIFF_METALS and TURB_DIFF_METALS_LOWORDER
  // have also been switched on.
#ifndef TURB_DIFF_METALS
  if (ThisTask == 0)
    {
      printf("ERROR: CHIMES_TURB_DIFF_IONS requires TURB_DIFF_METALS, but this is missing. Aborting.\n");
      endrun(6572);
    }
#endif // !(TURB_DIFF_METALS)
#ifndef TURB_DIFF_METALS_LOWORDER
  if (ThisTask == 0)
    {
      printf("ERROR: CHIMES_TURB_DIFF_IONS requires TURB_DIFF_METALS_LOWORDER, but this is missing. Aborting.\n");
      endrun(6573);
    }
#endif // !(TURB_DIFF_METALS_LOWORDER)
#endif // CHIMES_TURB_DIFF_IONS

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

  mymalloc_init();

#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

#ifdef GR_TABULATED_COSMOLOGY
#ifdef GR_TABULATED_COSMOLOGY_W
  fwa_init();
#endif
#endif

  set_units();
  set_cosmo_factors_for_current_time();
  All.Time = All.TimeBegin;

#ifdef COOLING
  InitCool();
#endif

#ifdef BOX_PERIODIC
  ewald_init();
#endif

#ifdef BOX_PERIODIC
    boxSize = All.BoxSize;
    boxHalf = 0.5 * All.BoxSize;
#endif
#ifdef BOX_LONG_X
    boxSize_X = All.BoxSize * BOX_LONG_X;
    boxHalf_X = 0.5 * boxSize_X;
#endif
#ifdef BOX_LONG_Y
    boxSize_Y = All.BoxSize * BOX_LONG_Y;
    boxHalf_Y = 0.5 * boxSize_Y;
#endif
#ifdef BOX_LONG_Z
    boxSize_Z = All.BoxSize * BOX_LONG_Z;
    boxHalf_Z = 0.5 * boxSize_Z;
#endif

#ifdef BOX_SHEARING
#ifdef BOX_LONG_X
    Shearing_Box_Vel_Offset = BOX_SHEARING_Q * BOX_SHEARING_OMEGA_BOX_CENTER * All.BoxSize * BOX_LONG_X;
#else
    Shearing_Box_Vel_Offset = BOX_SHEARING_Q * BOX_SHEARING_OMEGA_BOX_CENTER * All.BoxSize;
#endif
    calc_shearing_box_pos_offset();
#endif

    /* begin pre-definitions for special boundaries */
#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z)
    special_boundary_condition_xyz_def_reflect[0]=special_boundary_condition_xyz_def_reflect[1]=special_boundary_condition_xyz_def_reflect[2]=BOX_VALUE_FOR_NOTHING_SPECIAL_BOUNDARY_; /* sets arbitrary value code for 'nothing special' */
    special_boundary_condition_xyz_def_outflow[0]=special_boundary_condition_xyz_def_outflow[1]=special_boundary_condition_xyz_def_outflow[2]=BOX_VALUE_FOR_NOTHING_SPECIAL_BOUNDARY_; /* sets arbitrary value code for 'nothing special' */

#if defined(BOX_REFLECT_X)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(BOX_REFLECT_X)
    special_boundary_condition_xyz_def_reflect[0] = BOX_REFLECT_X; /* set to user definition */
#else
    special_boundary_condition_xyz_def_reflect[0] = 0; /* assume special boundary applies to both 'ends' of box, if not specified by user */
#endif
#endif
#if defined(BOX_REFLECT_Y)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(BOX_REFLECT_Y)
    special_boundary_condition_xyz_def_reflect[1] = BOX_REFLECT_Y; /* set to user definition */
#else
    special_boundary_condition_xyz_def_reflect[1] = 0; /* assume special boundary applies to both 'ends' of box, if not specified by user */
#endif
#endif
#if defined(BOX_REFLECT_Z)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(BOX_REFLECT_Z)
    special_boundary_condition_xyz_def_reflect[2] = BOX_REFLECT_Z; /* set to user definition */
#else
    special_boundary_condition_xyz_def_reflect[2] = 0; /* assume special boundary applies to both 'ends' of box, if not specified by user */
#endif
#endif

#if defined(BOX_OUTFLOW_X)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(BOX_OUTFLOW_X)
    special_boundary_condition_xyz_def_outflow[0] = BOX_OUTFLOW_X; /* set to user definition */
#else
    special_boundary_condition_xyz_def_outflow[0] = 0; /* assume special boundary applies to both 'ends' of box, if not specified by user */
#endif
#endif
#if defined(BOX_OUTFLOW_Y)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(BOX_OUTFLOW_Y)
    special_boundary_condition_xyz_def_outflow[1] = BOX_OUTFLOW_Y; /* set to user definition */
#else
    special_boundary_condition_xyz_def_outflow[1] = 0; /* assume special boundary applies to both 'ends' of box, if not specified by user */
#endif
#endif
#if defined(BOX_OUTFLOW_Z)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(BOX_OUTFLOW_Z)
    special_boundary_condition_xyz_def_outflow[2] = BOX_OUTFLOW_Z; /* set to user definition */
#else
    special_boundary_condition_xyz_def_outflow[2] = 0; /* assume special boundary applies to both 'ends' of box, if not specified by user */
#endif
#endif

#endif /* end set of clauses to deal with causal flags for special boundary conditions */




  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42 + ThisTask);	/* start-up seed */

  set_random_numbers();

#ifdef PMGRID
#ifndef ADAPTIVE_GRAVSOFT_FORALL
  if(RestartFlag != 3 && RestartFlag != 4)
#endif
    long_range_init();
#endif

#ifdef SUBFIND
  GrNr = -1;
#endif

#ifdef EOS_TABULATED
    int ierr = eos_init(All.EosTable);
    if(ierr) {
        printf("error initializing the eos");
        endrun(1);
    }
#endif

#ifdef EOS_TILLOTSON
    tillotson_eos_init();
#endif

#ifdef NUCLEAR_NETWORK
    network_init(All.EosSpecies, All.NetworkRates, All.NetworkPartFunc, All.NetworkMasses, All.NetworkWeakrates, &All.nd);
    network_workspace_init(&All.nd, &All.nw);
#endif

#ifdef TURB_DRIVING
    init_turb();
#endif

#if defined(DM_SIDM)
    init_geofactor_table();
#endif

    
  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2 || RestartFlag == 3 || RestartFlag == 4 || RestartFlag == 5 || RestartFlag == 6)
    {
      init();			/* ... read in initial model */
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets
				   all variables in the struct `All'.
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.
				 */

      set_random_numbers();

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.SnapFormat = all.SnapFormat;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MinGasHsmlFractional = all.MinGasHsmlFractional;
      All.MinGasTemp = all.MinGasTemp;
#ifdef CHIMES
      All.ChimesThermEvolOn = all.ChimesThermEvolOn;
#endif

        /* allow softenings to be modified during the run */
        if(All.ComovingIntegrationOn)
        {
            All.SofteningGasMaxPhys = all.SofteningGasMaxPhys;
            All.SofteningHaloMaxPhys = all.SofteningHaloMaxPhys;
            All.SofteningDiskMaxPhys = all.SofteningDiskMaxPhys;
            All.SofteningBulgeMaxPhys = all.SofteningBulgeMaxPhys;
            All.SofteningStarsMaxPhys = all.SofteningStarsMaxPhys;
            All.SofteningBndryMaxPhys = all.SofteningBndryMaxPhys;
        }
        All.SofteningGas = all.SofteningGas;
        All.SofteningHalo = all.SofteningHalo;
        All.SofteningDisk = all.SofteningDisk;
        All.SofteningBulge = all.SofteningBulge;
        All.SofteningStars = all.SofteningStars;
        All.SofteningBndry = all.SofteningBndry;

        All.MaxHsml = all.MaxHsml;
        All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

        All.ErrTolForceAcc = all.ErrTolForceAcc;
        All.NumFilesPerSnapshot = all.NumFilesPerSnapshot;
        All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
        All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

        All.OutputListOn = all.OutputListOn;
        All.CourantFac = all.CourantFac;
        
        All.OutputListLength = all.OutputListLength;
        memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);
        memcpy(All.OutputListFlag, all.OutputListFlag, sizeof(char) * All.OutputListLength);

#ifdef GALSF
        All.CritPhysDensity = all.CritPhysDensity;
        All.MaxSfrTimescale = all.MaxSfrTimescale;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        All.ArtBulkViscConst = all.ArtBulkViscConst;
        All.ViscosityAMin = all.ViscosityAMin;
        All.ViscosityAMax = all.ViscosityAMax;
#endif
#ifdef TURB_DIFFUSION
        All.TurbDiffusion_Coefficient = all.TurbDiffusion_Coefficient;
#endif
#ifdef SPHAV_ARTIFICIAL_CONDUCTIVITY
        All.ArtCondConstant = all.ArtCondConstant;
#endif
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
        All.ArtMagDispConst = all.ArtMagDispConst;
#endif
#ifdef DIVBCLEANING_DEDNER
        All.DivBcleanParabolicSigma = all.DivBcleanParabolicSigma;
        All.DivBcleanHyperbolicSigma = all.DivBcleanHyperbolicSigma;
        All.FastestWaveSpeed = 0.0;
        All.FastestWaveDecay = 0.0;
#endif
#ifdef BLACK_HOLES
        All.BlackHoleEddingtonFactor = all.BlackHoleEddingtonFactor;
        All.SeedBlackHoleMass = all.SeedBlackHoleMass;
        All.BlackHoleNgbFactor = all.BlackHoleNgbFactor;
        All.BlackHoleMaxAccretionRadius = all.BlackHoleMaxAccretionRadius;
        All.BlackHoleRadiativeEfficiency = all.BlackHoleRadiativeEfficiency;
        All.BlackHoleFeedbackFactor = all.BlackHoleFeedbackFactor;
#if defined(BH_SEED_FROM_FOF) || defined(BH_SEED_FROM_LOCALGAS)
        All.SeedBlackHoleMassSigma = all.SeedBlackHoleMassSigma;
        All.SeedBlackHoleMinRedshift = all.SeedBlackHoleMinRedshift;
#ifdef BH_SEED_FROM_LOCALGAS
        All.SeedBlackHolePerUnitMass = all.SeedBlackHolePerUnitMass;
#endif
#endif
#ifdef BH_ALPHADISK_ACCRETION
        All.SeedAlphaDiskMass = all.SeedAlphaDiskMass;
#endif
#ifdef BH_SEED_FROM_FOF
        All.MinFoFMassForNewSeed = all.MinFoFMassForNewSeed;
#endif
#if defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK) || defined(BH_WIND_SPAWN)
        All.BAL_f_accretion = all.BAL_f_accretion;
        All.BAL_v_outflow = all.BAL_v_outflow;
#endif
#if defined(SINGLE_STAR_FB_JETS)
        All.BAL_f_launch_v = all.BAL_f_launch_v;
#endif
#if defined(BH_COSMIC_RAYS)
        All.BH_CosmicRay_Injection_Efficiency = all.BH_CosmicRay_Injection_Efficiency;
#endif
#ifdef BH_WIND_SPAWN
        All.BAL_internal_temperature = all.BAL_internal_temperature;
        All.BAL_wind_particle_mass = all.BAL_wind_particle_mass; // dangeous to change this, as it is also part of the merger criterion!
#endif
#ifdef BH_PHOTONMOMENTUM
        All.BH_Rad_MomentumFactor = all.BH_Rad_MomentumFactor;
#endif
#endif // blackholes
#ifdef RT_LEBRON
        All.PhotonMomentum_Coupled_Fraction = all.PhotonMomentum_Coupled_Fraction;
#endif
#ifdef COSMIC_RAY_FLUID
#if (CRFLUID_DIFFUSION_MODEL == 0)
        All.CosmicRayDiffusionCoeff = all.CosmicRayDiffusionCoeff;
#endif
#endif

#ifdef GALSF_FB_FIRE_AGE_TRACERS
      All.AgeTracerRateNormalization = all.AgeTracerRateNormalization;
#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
      strcpy(All.AgeTracerListFilename, all.AgeTracerListFilename);
#else
      All.AgeTracerBinStart = all.AgeTracerBinStart;
      All.AgeTracerBinEnd = all.AgeTracerBinEnd;
#endif
#endif

#ifdef GR_TABULATED_COSMOLOGY
      All.DarkEnergyConstantW = all.DarkEnergyConstantW;
#endif

      All.MaxNumNgbDeviation = all.MaxNumNgbDeviation;
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
      /* Allow the tolerance over the number of neighbours to vary during the run:
        If it was initially set to a very strict value, convergence in ngb-iteration may at some point fail */
      All.AGS_MaxNumNgbDeviation = all.AGS_MaxNumNgbDeviation;
#endif

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      /*
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.TimebinFile, all.TimebinFile);
      */
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

#ifdef COOL_GRACKLE
      strcpy(All.GrackleDataFile, all.GrackleDataFile);
#endif

#ifdef EOS_TABULATED
        strcpy(All.EosTable, all.EosTable);
#endif

#ifdef NUCLEAR_NETWORK
      strcpy(All.EosSpecies, all.EosSpecies);
      strcpy(All.NetworkRates, all.NetworkRates);
      strcpy(All.NetworkPartFunc, all.NetworkPartFunc);
      strcpy(All.NetworkMasses, all.NetworkMasses);
      strcpy(All.NetworkWeakrates, all.NetworkWeakrates);
      All.nd = all.nd;
      All.nw = all.nw;
      All.NetworkTempThreshold = all.NetworkTempThreshold;
#endif

      if(All.TimeMax != all.TimeMax) {readjust_timebase(All.TimeMax, all.TimeMax);}
    }

#ifdef GALSF_EFFECTIVE_EQS
  init_clouds();
#endif

  char contfname[1000];
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);
  open_outputfiles();

#ifdef PMGRID
  long_range_init_regionsize();
#endif
  reconstruct_timebins();


#ifndef BOX_SHEARING
#if (NUMDIMS==2)
    int i;
    for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[2] = P[i].Vel[2] = P[i].GravAccel[2] = 0;
      if(P[i].Type == 0) {SphP[i].VelPred[2] = SphP[i].HydroAccel[2] = 0;}
    }
#endif
#endif

#ifdef RADTRANSFER
#if defined(RT_EVOLVE_INTENSITIES)
    rt_init_intensity_directions();
#endif
#if defined(RT_DIFFUSION_CG)
    All.Radiation_Ti_begstep = 0;
#endif
#ifdef RT_CHEM_PHOTOION
    rt_get_sigma();
#endif
    rt_set_simple_inits(RestartFlag);
#endif


  if(All.ComovingIntegrationOn) {init_drift_table();}

  if(RestartFlag == 2)
    {All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);}
  else if(RestartFlag == 1)
    {All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);}
  else
    {All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);}

  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the cgs-system
 */
void set_units(void)
{
  /* convert some physical input parameters to internal units */
  if(All.G <= 0) {All.G = GRAVITY_G_CGS * UNIT_MASS_IN_CGS / (UNIT_LENGTH_IN_CGS * UNIT_VEL_IN_CGS*UNIT_VEL_IN_CGS);}
#ifdef GR_TABULATED_COSMOLOGY_G
  All.Gini = All.G;
  All.G = All.Gini * dGfak(All.TimeBegin);
#endif
  All.Hubble_H0_CodeUnits = H0_CGS * UNIT_TIME_IN_CGS;
  if(ThisTask == 0)
    {
      printf("\nCode units to be used: make sure you check these are correct! \n");
      printf("  Hubble H0 (internal units) = %g \n", All.Hubble_H0_CodeUnits);
      printf("  Gravity G (internal units) = %g \n", All.G);
      printf("  unit Mass in g             = %g \n", UNIT_MASS_IN_CGS);
      printf("  unit Length in cm          = %g \n", UNIT_LENGTH_IN_CGS);
      printf("  unit Time in s             = %g \n", UNIT_TIME_IN_CGS);
      printf("  unit Velocity in cm/s      = %g \n", UNIT_VEL_IN_CGS);
      printf("  unit Energy in erg         = %g \n", UNIT_ENERGY_IN_CGS);
      printf("  unit Density in g/cm^3     = %g \n", UNIT_DENSITY_IN_CGS);
      printf("  unit Pressure in erg/cm^3  = %g \n", UNIT_PRESSURE_IN_CGS);
      printf("  unit Luminosity in erg/s   = %g \n", UNIT_LUM_IN_CGS);
      printf("  unit Flux in erg/s/cm^2    = %g \n", UNIT_FLUX_IN_CGS);
      printf("  unit B[internal] in gauss  = %g \n", UNIT_B_IN_GAUSS);
      printf("\n");
    }

    double meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC); /* assumes fully-atomic otherwise */
#ifdef COOL_LOW_TEMPERATURES
    meanweight = 1. / ( HYDROGEN_MASSFRAC*0.5 + (1-HYDROGEN_MASSFRAC)/4. + 1./(16.+12.)); /* assumes fully-molecular if low-temp cooling enabled */
#endif
    All.MinEgySpec = All.MinGasTemp / (meanweight * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS);


#if defined(GALSF)
  /* for historical reasons, we need to convert to "All.MaxSfrTimescale", defined as the SF timescale in code units at the critical physical
     density given above. use the dimensionless SfEffPerFreeFall (which has been read in) to calculate this. This must be done -BEFORE- calling set_units_sfr) */
#ifndef GALSF_EFFECTIVE_EQS
    All.MaxSfrTimescale = (1/All.MaxSfrTimescale) * sqrt(3.*M_PI / (32. * All.G * (All.CritPhysDensity / UNIT_DENSITY_IN_NHCGS)));
#ifdef PROTECT_FROZEN_FIRE
    All.MaxSfrTimescale /= sqrt(meanweight);
#endif
#endif
    set_units_sfr();
#endif


#ifdef DM_FUZZY
    /* For Schroedinger equation: this encodes the coefficient with the mass of the particle: units vel*L = hbar / particle_mass. This is the key variable used throughout */
    All.ScalarField_hbar_over_mass = 591569.0 / ((double)All.ScalarField_hbar_over_mass * UNIT_VEL_IN_CGS * UNIT_LENGTH_IN_CGS);
#endif


#if defined(CONDUCTION_SPITZER) || defined(VISCOSITY_BRAGINSKII)
    /* Note: Because we replace \nabla(T) in the conduction equation with \nabla(u), our conduction coefficient is not the usual kappa, but
     * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with another factor of (meanweight_ion / k_B * (gamma-1)) */
    double meanweight_ion =  4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* mean weight in code units, assuming full ionization */
    double u_to_temp = meanweight_ion * (GAMMA_DEFAULT-1.) * U_TO_TEMP_UNITS; /* for full ionization, assume gas has a monatomic ideal eos gamma=5/3 */
    /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003 ( ApJ 582:162-169, Eq. (5) ) */
    double coulomb_log = 37.8; // Sarazin value (recommendation from PIC calculations) //
    double coefficient = (1.84e-5/coulomb_log) * pow(u_to_temp,3.5) * ((UNIT_TIME_IN_CGS*UNIT_TIME_IN_CGS*UNIT_TIME_IN_CGS) / (UNIT_LENGTH_IN_CGS*UNIT_MASS_IN_CGS)); // ok, this multiplied by the specific energy (u_code)^(3/2) gives the diffusity of u_code, as needed (density term is included in said diffusivity)
#ifdef CONDUCTION_SPITZER
    All.ConductionCoeff *= coefficient;
#endif
#ifdef VISCOSITY_BRAGINSKII
    All.ShearViscosityCoeff *= coefficient * 0.636396*sqrt(ELECTRONMASS_CGS/(PROTONMASS_CGS*meanweight_ion)); // the viscosity coefficient eta is identical in these units up to the order-unity constant, and multiplied by sqrt[m_electron/m_ion] //
    All.BulkViscosityCoeff = 0; // no bulk viscosity in the Braginskii-Spitzer formulation //
#endif
    /* factor used for determining saturation */
    All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow((GAMMA_DEFAULT-1), 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
        / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS_CGS, 3) / pow(ELECTRONCHARGE_CGS, 4) / (UNIT_DENSITY_IN_CGS) * pow(UNIT_SPECEGY_IN_CGS, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electron mean free path in centimeters. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units. */
  All.ElectronFreePathFactor /= UNIT_LENGTH_IN_CGS;
#endif


}




/*!  This function opens various log-files that report on the status and performance of the simulation.
        On restart from restart-files, (start-option 1), the code will append to these files. */
void open_outputfiles(void)
{
  char mode[2], buf[200];
  if(RestartFlag == 0) {strcpy(mode, "w");} else {strcpy(mode, "a");}
  if(ThisTask == 0) {mkdir(All.OutputDir, 02755);}
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef BLACK_HOLES /* Note: This is done by everyone [all tasks can write to these log-files], even if it might be empty */
  if(ThisTask == 0) {sprintf(buf, "%sblackhole_details", All.OutputDir); mkdir(buf, 02755);}
  MPI_Barrier(MPI_COMM_WORLD);
#if !defined(IO_REDUCED_MODE) || defined(BH_OUTPUT_MOREINFO)
  sprintf(buf, "%sblackhole_details/blackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif
#ifdef OUTPUT_SINK_ACCRETION_HIST
  sprintf(buf, "%sblackhole_details/bhswallow_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBhSwallowDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif
#ifdef OUTPUT_SINK_FORMATION_PROPS
  sprintf(buf, "%sblackhole_details/bhformation_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBhFormationDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif
#ifdef BH_OUTPUT_MOREINFO
  sprintf(buf, "%sblackhole_details/bhmergers_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBhMergerDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#ifdef BH_WIND_KICK
  sprintf(buf, "%sblackhole_details/bhwinds_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBhWindDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif
#endif // bh-output-more-info if
#endif // black-holes if

    if(ThisTask != 0) {return;}	/* only the root processors writes to the log files listed below */

    sprintf(buf, "%s%s", All.OutputDir, "cpu.txt");
    if(!(FdCPU = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}

#ifndef IO_REDUCED_MODE
    sprintf(buf, "%s%s", All.OutputDir, "timebin.txt");
    if(!(FdTimebin = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    sprintf(buf, "%s%s", All.OutputDir, "info.txt");
    if(!(FdInfo = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    sprintf(buf, "%s%s", All.OutputDir, "energy.txt");
    if(!(FdEnergy = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    sprintf(buf, "%s%s", All.OutputDir, "timings.txt");
    if(!(FdTimings = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}

    sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
    if(!(FdBalance = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    fprintf(FdBalance, "\n");
    fprintf(FdBalance, "Treewalk1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK1], CPU_SymbolImbalance[CPU_TREEWALK1]);
    fprintf(FdBalance, "Treewalk2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK2], CPU_SymbolImbalance[CPU_TREEWALK2]);
    fprintf(FdBalance, "Treewait1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT1], CPU_SymbolImbalance[CPU_TREEWAIT1]);
    fprintf(FdBalance, "Treewait2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT2], CPU_SymbolImbalance[CPU_TREEWAIT2]);
    fprintf(FdBalance, "Treesend       = '%c' / '%c'\n", CPU_Symbol[CPU_TREESEND], CPU_SymbolImbalance[CPU_TREESEND]);
    fprintf(FdBalance, "Treerecv       = '%c' / '%c'\n", CPU_Symbol[CPU_TREERECV], CPU_SymbolImbalance[CPU_TREERECV]);
    fprintf(FdBalance, "Treebuild      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEBUILD], CPU_SymbolImbalance[CPU_TREEBUILD]);
    fprintf(FdBalance, "Treehmaxupdate = '%c' / '%c'\n", CPU_Symbol[CPU_TREEHMAXUPDATE], CPU_SymbolImbalance[CPU_TREEHMAXUPDATE]);
    fprintf(FdBalance, "Treemisc =       '%c' / '%c'\n", CPU_Symbol[CPU_TREEMISC], CPU_SymbolImbalance[CPU_TREEMISC]);
    fprintf(FdBalance, "Domain decomp  = '%c' / '%c'\n", CPU_Symbol[CPU_DOMAIN], CPU_SymbolImbalance[CPU_DOMAIN]);
    fprintf(FdBalance, "Peano-Hilbert  = '%c' / '%c'\n", CPU_Symbol[CPU_PEANO], CPU_SymbolImbalance[CPU_PEANO]);
    fprintf(FdBalance, "Density compute= '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMPUTE], CPU_SymbolImbalance[CPU_DENSCOMPUTE]);
    fprintf(FdBalance, "Density imbal  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSWAIT], CPU_SymbolImbalance[CPU_DENSWAIT]);
    fprintf(FdBalance, "Density commu  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMM], CPU_SymbolImbalance[CPU_DENSCOMM]);
    fprintf(FdBalance, "Density misc   = '%c' / '%c'\n", CPU_Symbol[CPU_DENSMISC], CPU_SymbolImbalance[CPU_DENSMISC]);
    fprintf(FdBalance, "Hydro compute  = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMPUTE], CPU_SymbolImbalance[CPU_HYDCOMPUTE]);
    fprintf(FdBalance, "Hydro imbalance= '%c' / '%c'\n", CPU_Symbol[CPU_HYDWAIT], CPU_SymbolImbalance[CPU_HYDWAIT]);
    fprintf(FdBalance, "Hydro comm     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMM], CPU_SymbolImbalance[CPU_HYDCOMM]);
    fprintf(FdBalance, "Hydro misc     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDMISC], CPU_SymbolImbalance[CPU_HYDMISC]);
    fprintf(FdBalance, "Drifts         = '%c' / '%c'\n", CPU_Symbol[CPU_DRIFT], CPU_SymbolImbalance[CPU_DRIFT]);
    fprintf(FdBalance, "Kicks          = '%c' / '%c'\n", CPU_Symbol[CPU_TIMELINE], CPU_SymbolImbalance[CPU_TIMELINE]);
    fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL], CPU_SymbolImbalance[CPU_POTENTIAL]);
    fprintf(FdBalance, "PM-gravity     = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
    fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT], CPU_SymbolImbalance[CPU_SNAPSHOT]);
    fprintf(FdBalance, "Blackhole      = '%c' / '%c'\n", CPU_Symbol[CPU_BLACKHOLES], CPU_SymbolImbalance[CPU_BLACKHOLES]);
    fprintf(FdBalance, "Cooling & SFR  = '%c' / '%c'\n", CPU_Symbol[CPU_COOLINGSFR], CPU_SymbolImbalance[CPU_COOLINGSFR]);
    fprintf(FdBalance, "Coolimbal check= '%c' / '%c'\n", CPU_Symbol[CPU_COOLSFRIMBAL], CPU_SymbolImbalance[CPU_COOLSFRIMBAL]);
    fprintf(FdBalance, "FoF & subfind  = '%c' / '%c'\n", CPU_Symbol[CPU_FOF], CPU_SymbolImbalance[CPU_FOF]);
    fprintf(FdBalance, "Grain/PIC part = '%c' / '%c'\n", CPU_Symbol[CPU_DRAGFORCE], CPU_SymbolImbalance[CPU_DRAGFORCE]);
    fprintf(FdBalance, "Mech/Thermal FB= '%c' / '%c'\n", CPU_Symbol[CPU_SNIIHEATING], CPU_SymbolImbalance[CPU_SNIIHEATING]);
    fprintf(FdBalance, "HII-module     = '%c' / '%c'\n", CPU_Symbol[CPU_HIIHEATING], CPU_SymbolImbalance[CPU_HIIHEATING]);
    fprintf(FdBalance, "Local wind kick= '%c' / '%c'\n", CPU_Symbol[CPU_LOCALWIND], CPU_SymbolImbalance[CPU_LOCALWIND]);
    fprintf(FdBalance, "RHD-nonfluxops = '%c' / '%c'\n", CPU_Symbol[CPU_RTNONFLUXOPS], CPU_SymbolImbalance[CPU_RTNONFLUXOPS]);
    fprintf(FdBalance, "AGS-nongas-comp= '%c' / '%c'\n", CPU_Symbol[CPU_AGSDENSCOMPUTE], CPU_SymbolImbalance[CPU_AGSDENSCOMPUTE]);
    fprintf(FdBalance, "AGS-imbal      = '%c' / '%c'\n", CPU_Symbol[CPU_AGSDENSWAIT], CPU_SymbolImbalance[CPU_AGSDENSWAIT]);
    fprintf(FdBalance, "AGS-comm       = '%c' / '%c'\n", CPU_Symbol[CPU_AGSDENSCOMM], CPU_SymbolImbalance[CPU_AGSDENSCOMM]);
    fprintf(FdBalance, "AGS-misc       = '%c' / '%c'\n", CPU_Symbol[CPU_AGSDENSMISC], CPU_SymbolImbalance[CPU_AGSDENSMISC]);
    fprintf(FdBalance, "DynDiffusn-comp= '%c' / '%c'\n", CPU_Symbol[CPU_DYNDIFFCOMPUTE], CPU_SymbolImbalance[CPU_DYNDIFFCOMPUTE]);
    fprintf(FdBalance, "DynDiffusn-imbl= '%c' / '%c'\n", CPU_Symbol[CPU_DYNDIFFWAIT], CPU_SymbolImbalance[CPU_DYNDIFFWAIT]);
    fprintf(FdBalance, "DynDiffusn-comm= '%c' / '%c'\n", CPU_Symbol[CPU_DYNDIFFCOMM], CPU_SymbolImbalance[CPU_DYNDIFFCOMM]);
    fprintf(FdBalance, "DynDiffusn-misc= '%c' / '%c'\n", CPU_Symbol[CPU_DYNDIFFMISC], CPU_SymbolImbalance[CPU_DYNDIFFMISC]);
    fprintf(FdBalance, "MultiDiff-comp = '%c' / '%c'\n", CPU_Symbol[CPU_IMPROVDIFFCOMPUTE], CPU_SymbolImbalance[CPU_IMPROVDIFFCOMPUTE]);
    fprintf(FdBalance, "MultiDiff-imbl = '%c' / '%c'\n", CPU_Symbol[CPU_IMPROVDIFFWAIT], CPU_SymbolImbalance[CPU_IMPROVDIFFWAIT]);
    fprintf(FdBalance, "MultiDiff-comm = '%c' / '%c'\n", CPU_Symbol[CPU_IMPROVDIFFCOMM], CPU_SymbolImbalance[CPU_IMPROVDIFFCOMM]);
    fprintf(FdBalance, "MultiDiff-misc = '%c' / '%c'\n", CPU_Symbol[CPU_IMPROVDIFFMISC], CPU_SymbolImbalance[CPU_IMPROVDIFFMISC]);
    fprintf(FdBalance, "Miscellaneous  = '%c' / '%c'\n", CPU_Symbol[CPU_MISC], CPU_SymbolImbalance[CPU_MISC]);
    fprintf(FdBalance, "\n");
#endif

#ifdef GALSF
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif



#ifdef GALSF_FB_MECHANICAL
    sprintf(buf, "%s%s", All.OutputDir, "SNeIIheating.txt");
    if(!(FdSneIIHeating = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif

#if defined(RT_CHEM_PHOTOION) && !defined(IO_REDUCED_MODE)
  sprintf(buf, "%s%s", All.OutputDir, "rt_photoion_chem.txt");
  if(!(FdRad = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif


#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif

#if defined(TURB_DRIVING) && !defined(IO_REDUCED_MODE)
  sprintf(buf, "%s%s", All.OutputDir, "turb.txt");
  if(!(FdTurb = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif

#if defined(GR_TABULATED_COSMOLOGY) && !defined(IO_REDUCED_MODE)
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0)
  {
	  fprintf(FdDE, "nstep time H(a) ");
#ifndef GR_TABULATED_COSMOLOGY_W
	  fprintf(FdDE, "w0 Omega_L ");
#else
	  fprintf(FdDE, "w(a) Omega_L ");
#endif
#ifdef GR_TABULATED_COSMOLOGY_G
	  fprintf(FdDE, "dH dG ");
#endif
      fprintf(FdDE, "\n"); fflush(FdDE);
  }
#endif

}






/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  char alternate_tag[MAXTAGS][50];
  int pnum, errorFlag = 0;

#ifdef CHIMES
  double Tdust_buf, Tmol_buf, relTol_buf, absTol_buf, expTol_buf, z_reion_buf;
#endif

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;
      for(j=0;j<MAXTAGS;j++) {strcpy(alternate_tag[nt], "-null[invalid_tag_name]-");}

      strcpy(tag[nt], "InitCondFile");
      strcpy(alternate_tag[nt], "Initial_Conditions_File");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      strcpy(alternate_tag[nt], "Output_Directory");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      strcpy(alternate_tag[nt], "Snapshot_Filename_Base");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      strcpy(alternate_tag[nt], "Restart_Filename_Base");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

#ifdef DEVELOPER_MODE
      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitCommand");
      strcpy(alternate_tag[nt], "Shell_Resubmission_Command");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;
#endif

      strcpy(tag[nt], "OutputListFilename");
      strcpy(alternate_tag[nt], "Snapshot_Times_Table_Filename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      strcpy(alternate_tag[nt], "Use_Tabulated_Snapshot_Times");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      strcpy(alternate_tag[nt], "Omega_Matter");
      addr[nt] = &All.OmegaMatter;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      strcpy(alternate_tag[nt], "Omega_Baryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      strcpy(alternate_tag[nt], "Omega_Lambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaRadiation");
      strcpy(alternate_tag[nt], "Omega_Radiation");
      addr[nt] = &All.OmegaRadiation;
      id[nt++] = REAL;
        
      strcpy(tag[nt], "HubbleParam");
      strcpy(alternate_tag[nt], "Hubble_Param_Little_h");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      strcpy(alternate_tag[nt], "Box_Size_In_Code_Units");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxMemSize");
      strcpy(alternate_tag[nt], "Max_Memory_Per_MPI_Task_in_MB");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      strcpy(alternate_tag[nt], "Simulation_Time_of_First_Snapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      strcpy(alternate_tag[nt], "Walltime_in_Seconds_Between_Restartfiles");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

#ifdef DEVELOPER_MODE
      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "TimeBegin");
      strcpy(alternate_tag[nt], "Time_at_ICs_Begin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      strcpy(alternate_tag[nt], "Time_at_End_of_Simulation");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      strcpy(alternate_tag[nt], "Time_Between_Snapshots");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      strcpy(alternate_tag[nt], "UnitVelocity_in_cm_per_seconds");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      strcpy(alternate_tag[nt], "UnitLength_in_cms");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      strcpy(alternate_tag[nt], "UnitMass_in_grams");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

#ifdef MAGNETIC
      strcpy(tag[nt], "UnitMagneticField_in_gauss");
      strcpy(alternate_tag[nt], "UnitMagneticField_in_Gauss_for_ICSnapshotIO");
      addr[nt] = &All.UnitMagneticField_in_gauss;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      strcpy(alternate_tag[nt], "DomainTreeRebuild_ActiveFractionThreshold");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = REAL;

#ifdef DEVELOPER_MODE
        strcpy(tag[nt], "ErrTolIntAccuracy");
        addr[nt] = &All.ErrTolIntAccuracy;
        id[nt++] = REAL;

        strcpy(tag[nt], "ErrTolTheta");
        addr[nt] = &All.ErrTolTheta;
        id[nt++] = REAL;

        strcpy(tag[nt], "CourantFac");
        addr[nt] = &All.CourantFac;
        id[nt++] = REAL;

        strcpy(tag[nt], "ErrTolForceAcc");
        addr[nt] = &All.ErrTolForceAcc;
        id[nt++] = REAL;

        strcpy(tag[nt], "MaxRMSDisplacementFac");
        addr[nt] = &All.MaxRMSDisplacementFac;
        id[nt++] = REAL;

#ifdef HYDRO_SPH
        strcpy(tag[nt], "ArtBulkViscConst");
        addr[nt] = &All.ArtBulkViscConst;
        id[nt++] = REAL;
#ifdef SPHAV_ARTIFICIAL_CONDUCTIVITY
        strcpy(tag[nt], "ArtCondConstant");
        addr[nt] = &All.ArtCondConstant;
        id[nt++] = REAL;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        strcpy(tag[nt], "ViscosityAMin");
        addr[nt] = &All.ViscosityAMin;
        id[nt++] = REAL;

        strcpy(tag[nt], "ViscosityAMax");
        addr[nt] = &All.ViscosityAMax;
        id[nt++] = REAL;
#endif
#ifdef SPH_TP12_ARTIFICIAL_RESISTIVITY
        strcpy(tag[nt], "ArtificialResistivityMax");
        addr[nt] = &All.ArtMagDispConst;
        id[nt++] = REAL;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
        strcpy(tag[nt], "DivBcleaningParabolicSigma");
        addr[nt] = &All.DivBcleanParabolicSigma;
        id[nt++] = REAL;

        strcpy(tag[nt], "DivBcleaningHyperbolicSigma");
        addr[nt] = &All.DivBcleanHyperbolicSigma;
        id[nt++] = REAL;
#endif
#endif // closes DEVELOPER_MODE check


#ifdef GRAIN_FLUID
#ifdef GRAIN_RDI_TESTPROBLEM
        strcpy(tag[nt],"Grain_Charge_Parameter");
        addr[nt] = &All.Grain_Charge_Parameter;
        id[nt++] = REAL;

        strcpy(tag[nt],"Dust_to_Gas_Mass_Ratio");
        addr[nt] = &All.Dust_to_Gas_Mass_Ratio;
        id[nt++] = REAL;

        strcpy(tag[nt],"Vertical_Gravity_Strength");
        addr[nt] = &All.Vertical_Gravity_Strength;
        id[nt++] = REAL;

        strcpy(tag[nt],"Vertical_Grain_Accel");
        addr[nt] = &All.Vertical_Grain_Accel;
        id[nt++] = REAL;

        strcpy(tag[nt],"Vertical_Grain_Accel_Angle");
        addr[nt] = &All.Vertical_Grain_Accel_Angle;
        id[nt++] = REAL;
        
#ifdef BOX_SHEARING
        strcpy(tag[nt],"Pressure_Gradient_Accel");
        addr[nt] = &All.Pressure_Gradient_Accel;
        id[nt++] = REAL;
#endif
        
#ifdef RT_OPACITY_FROM_EXPLICIT_GRAINS
        strcpy(tag[nt],"Grain_Q_at_MaxGrainSize");
        addr[nt] = &All.Grain_Q_at_MaxGrainSize;
        id[nt++] = REAL;
#endif

#endif
#if !defined(PIC_MHD) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        strcpy(tag[nt],"Grain_Internal_Density");
        addr[nt] = &All.Grain_Internal_Density;
        id[nt++] = REAL;

        strcpy(tag[nt],"Grain_Size_Min");
        addr[nt] = &All.Grain_Size_Min;
        id[nt++] = REAL;

        strcpy(tag[nt],"Grain_Size_Max");
        addr[nt] = &All.Grain_Size_Max;
        id[nt++] = REAL;

        strcpy(tag[nt],"Grain_Size_Spectrum_Powerlaw");
        addr[nt] = &All.Grain_Size_Spectrum_Powerlaw;
        id[nt++] = REAL;
#endif
#endif

#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS) && defined(RT_GENERIC_USER_FREQ)
        strcpy(tag[nt],"Grain_Absorbed_vs_Total_Extinction");
        addr[nt] = &All.Grain_Absorbed_Fraction_vs_Total_Extinction;
        id[nt++] = REAL;
#endif

#ifdef PIC_MHD
        strcpy(tag[nt],"PIC_Charge_to_Mass_Ratio");
        addr[nt] = &All.PIC_Charge_to_Mass_Ratio;
        id[nt++] = REAL;
#endif

#if defined(INIT_STELLAR_METALS_AGES_DEFINED)
        strcpy(tag[nt],"InitMetallicity");
        strcpy(alternate_tag[nt],"Initial_Metallicity");
        addr[nt] = &All.InitMetallicityinSolar;
        id[nt++] = REAL;

        strcpy(tag[nt],"InitStellarAge");
        strcpy(alternate_tag[nt],"Initial_StellarAge_NonTypeFourStars");
        addr[nt] = &All.InitStellarAgeinGyr;
        id[nt++] = REAL;
#endif


        
        
        



#ifdef GALSF_FB_FIRE_AGE_TRACERS
        strcpy(tag[nt], "AgeTracerEventsPerTimeBin");
        strcpy(alternate_tag[nt], "AgeTracerActiveTimestepFraction");
        addr[nt] = &All.AgeTracerRateNormalization;
        id[nt++] = REAL;

#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
        strcpy(tag[nt], "AgeTracerListFilename");
        addr[nt] = &All.AgeTracerListFilename;
        id[nt++] = STRING;
#else
        strcpy(tag[nt], "AgeTracerBinStart");
        addr[nt] = &All.AgeTracerBinStart;
        id[nt++] = REAL;

        strcpy(tag[nt], "AgeTracerBinEnd");
        addr[nt] = &All.AgeTracerBinEnd;
        id[nt++] = REAL;
#endif
#endif

#ifdef RT_LEBRON
        strcpy(tag[nt], "PhotonMomentum_Coupled_Fraction");
        addr[nt] = &All.PhotonMomentum_Coupled_Fraction;
        id[nt++] = REAL;
#endif




#ifdef DM_SIDM
#ifdef GRAIN_COLLISIONS
        strcpy(tag[nt], "Grain_InteractionRenormalization");
        addr[nt] = &All.DM_InteractionCrossSection;
        id[nt++] = REAL;

        strcpy(tag[nt], "Grain_DissipationFactor");
        addr[nt] = &All.DM_DissipationFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "Grain_KickPerCollision");
        addr[nt] = &All.DM_KickPerCollision;
        id[nt++] = REAL;

        strcpy(tag[nt], "Grain_InteractionVelocityScale");
        addr[nt] = &All.DM_InteractionVelocityScale;
        id[nt++] = REAL;
#else
        strcpy(tag[nt], "DM_InteractionCrossSection");
        addr[nt] = &All.DM_InteractionCrossSection;
        id[nt++] = REAL;

        strcpy(tag[nt], "DM_DissipationFactor");
        addr[nt] = &All.DM_DissipationFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "DM_KickPerCollision");
        addr[nt] = &All.DM_KickPerCollision;
        id[nt++] = REAL;

        strcpy(tag[nt], "DM_InteractionVelocityScale");
        addr[nt] = &All.DM_InteractionVelocityScale;
        id[nt++] = REAL;
#endif
#endif


        strcpy(tag[nt], "MinGasHsmlFractional");
        strcpy(alternate_tag[nt], "Minimum_Gas_KernelSize_RelativetoSoftening");
        addr[nt] = &All.MinGasHsmlFractional;
        id[nt++] = REAL;

        strcpy(tag[nt], "MaxHsml");
        strcpy(alternate_tag[nt], "Maximum_KernelSize_CodeUnits");
        addr[nt] = &All.MaxHsml;
        id[nt++] = REAL;

        strcpy(tag[nt], "MaxSizeTimestep");
        strcpy(alternate_tag[nt], "Maximum_Timestep_Allowed");
        addr[nt] = &All.MaxSizeTimestep;
        id[nt++] = REAL;

        strcpy(tag[nt], "MinSizeTimestep");
        strcpy(alternate_tag[nt], "Minimum_Timestep_Allowed");
        addr[nt] = &All.MinSizeTimestep;
        id[nt++] = REAL;


        strcpy(tag[nt], "DesNumNgb");
        strcpy(alternate_tag[nt], "Effective_Kernel_NeighborNumber");
        addr[nt] = &All.DesNumNgb;
        id[nt++] = REAL;


#ifdef SUBFIND
      strcpy(tag[nt], "DesLinkNgb");
      strcpy(alternate_tag[nt], "Subfind_FOFLink_NeighborNumber");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = INT;
#endif

#ifdef DEVELOPER_MODE
        strcpy(tag[nt], "MaxNumNgbDeviation");
        addr[nt] = &All.MaxNumNgbDeviation;
        id[nt++] = REAL;
#endif

      strcpy(tag[nt], "ComovingIntegrationOn");
      strcpy(alternate_tag[nt], "Cosmological_Simulation_On");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      strcpy(alternate_tag[nt], "Initial_Conditions_Format");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      strcpy(alternate_tag[nt], "Snapshot_Format");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      strcpy(alternate_tag[nt], "Number_of_Files_per_Snapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      strcpy(alternate_tag[nt], "Number_of_Files_Written_in_Parallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

#ifdef COOL_GRACKLE
        strcpy(tag[nt], "GrackleDataFile");
        addr[nt] = All.GrackleDataFile;
        id[nt++] = STRING;
#endif

      strcpy(tag[nt], "TimeLimitCPU");
      strcpy(alternate_tag[nt], "MaxSimulationWallTime_in_Seconds");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHalo");
      strcpy(alternate_tag[nt], "Softening_Type1");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDisk");
      strcpy(alternate_tag[nt], "Softening_Type2");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulge");
      strcpy(alternate_tag[nt], "Softening_Type3");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGas");
      strcpy(alternate_tag[nt], "Softening_Type0");
      addr[nt] = &All.SofteningGas;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStars");
      strcpy(alternate_tag[nt], "Softening_Type4");
      addr[nt] = &All.SofteningStars;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndry");
      strcpy(alternate_tag[nt], "Softening_Type5");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      strcpy(alternate_tag[nt], "Softening_Type1_MaxPhysLimit");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      strcpy(alternate_tag[nt], "Softening_Type2_MaxPhysLimit");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      strcpy(alternate_tag[nt], "Softening_Type3_MaxPhysLimit");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      strcpy(alternate_tag[nt], "Softening_Type0_MaxPhysLimit");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      strcpy(alternate_tag[nt], "Softening_Type4_MaxPhysLimit");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      strcpy(alternate_tag[nt], "Softening_Type5_MaxPhysLimit");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "BufferSize");
      strcpy(alternate_tag[nt], "MPI_Buffersize_in_MB");
      addr[nt] = &All.BufferSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PartAllocFactor");
      strcpy(alternate_tag[nt], "ParticleNumberMemoryImbalance_Limit");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "GravityConstantInternal");
      strcpy(alternate_tag[nt], "GravityConstant_SetByHand_inCodeUnits");
      addr[nt] = &All.G;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitGasTemp");
      strcpy(alternate_tag[nt], "Initial_Gas_Temperature");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasTemp");
      strcpy(alternate_tag[nt], "Minimum_Gas_Temperature");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = REAL;

#ifdef GDE_DISTORTIONTENSOR
      strcpy(tag[nt], "TidalCorrection");
      addr[nt] = &All.TidalCorrection;
      id[nt++] = REAL;

      strcpy(tag[nt], "DM_velocity_dispersion");
      addr[nt] = &All.DM_velocity_dispersion;
      id[nt++] = REAL;
#endif
#ifdef DM_SCALARFIELD_SCREENING
      strcpy(tag[nt], "ScalarBeta");
      addr[nt] = &All.ScalarBeta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ScalarScreeningLength");
      addr[nt] = &All.ScalarScreeningLength;
      id[nt++] = REAL;
#endif

#ifdef OUTPUT_LINEOFSIGHT
      strcpy(tag[nt], "TimeFirstLineOfSight");
      addr[nt] = &All.TimeFirstLineOfSight;
      id[nt++] = REAL;
#endif


#ifdef COSMIC_RAY_FLUID
#if (CRFLUID_DIFFUSION_MODEL == 0)
        strcpy(tag[nt], "CosmicRayDiffusionCoeff");
        addr[nt] = &All.CosmicRayDiffusionCoeff;
        id[nt++] = REAL;
#endif
#endif


#if (defined(BLACK_HOLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
      strcpy(tag[nt], "TimeBetOnTheFlyFoF");
      addr[nt] = &All.TimeBetOnTheFlyFoF;
      id[nt++] = REAL;
#endif

#ifdef BLACK_HOLES
        strcpy(tag[nt], "BlackHoleAccretionFactor");
        addr[nt] = &All.BlackHoleAccretionFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "BlackHoleEddingtonFactor");
        addr[nt] = &All.BlackHoleEddingtonFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "SeedBlackHoleMass");
        addr[nt] = &All.SeedBlackHoleMass;
        id[nt++] = REAL;

        strcpy(tag[nt], "BlackHoleNgbFactor");
        addr[nt] = &All.BlackHoleNgbFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "BlackHoleMaxAccretionRadius");
        addr[nt] = &All.BlackHoleMaxAccretionRadius;
        id[nt++] = REAL;

        strcpy(tag[nt], "BlackHoleRadiativeEfficiency");
        addr[nt] = &All.BlackHoleRadiativeEfficiency;
        id[nt++] = REAL;

        strcpy(tag[nt], "BlackHoleFeedbackFactor");
        addr[nt] = &All.BlackHoleFeedbackFactor;
        id[nt++] = REAL;

#if defined(BH_SEED_FROM_FOF) || defined(BH_SEED_FROM_LOCALGAS)
        strcpy(tag[nt], "SeedBlackHoleMassSigma");
        addr[nt] = &All.SeedBlackHoleMassSigma;
        id[nt++] = REAL;

        strcpy(tag[nt], "SeedBlackHoleMinRedshift");
        addr[nt] = &All.SeedBlackHoleMinRedshift;
        id[nt++] = REAL;

#ifdef BH_SEED_FROM_LOCALGAS
        strcpy(tag[nt], "SeedBlackHolePerUnitMass");
        addr[nt] = &All.SeedBlackHolePerUnitMass;
        id[nt++] = REAL;
#endif
#endif

#ifdef BH_ALPHADISK_ACCRETION
        strcpy(tag[nt], "SeedAlphaDiskMass");
        addr[nt] = &All.SeedAlphaDiskMass;
        id[nt++] = REAL;
#endif

#ifdef BH_SEED_FROM_FOF
        strcpy(tag[nt], "MinFoFMassForNewSeed");
        addr[nt] = &All.MinFoFMassForNewSeed;
        id[nt++] = REAL;
#endif

#if defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK) || defined(BH_WIND_SPAWN)
        strcpy(tag[nt],"BAL_f_accretion");
        addr[nt] = &All.BAL_f_accretion;
        id[nt++] = REAL;

        strcpy(tag[nt],"BAL_v_outflow");
        addr[nt] = &All.BAL_v_outflow;
        id[nt++] = REAL;
#endif

#if defined(SINGLE_STAR_FB_JETS)
        strcpy(tag[nt],"BAL_f_launch_v");
        addr[nt] = &All.BAL_f_launch_v;
        id[nt++] = REAL;
#endif

#if defined(BH_COSMIC_RAYS)
        strcpy(tag[nt],"BH_CosmicRay_Injection_Efficiency");
        addr[nt] = &All.BH_CosmicRay_Injection_Efficiency;
        id[nt++] = REAL;
#endif


#ifdef BH_WIND_SPAWN
        strcpy(tag[nt], "BAL_internal_temperature");
        addr[nt] = &All.BAL_internal_temperature;
        id[nt++] = REAL;

        strcpy(tag[nt], "BAL_wind_particle_mass");
        addr[nt] = &All.BAL_wind_particle_mass;
        id[nt++] = REAL;
#endif

#ifdef BH_PHOTONMOMENTUM
        strcpy(tag[nt],"BH_FluxMomentumFactor");
        strcpy(alternate_tag[nt], "BH_Rad_MomentumFactor");
        addr[nt] = &All.BH_Rad_MomentumFactor;
        id[nt++] = REAL;
#endif

#endif /* BLACK_HOLES */


#ifdef GALSF
#ifndef GALSF_EFFECTIVE_EQS
      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "SfEffPerFreeFall");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;
      /* for historical reasons, we need to convert to "MaxSfrTimescale",
            defined as the SF timescale in code units at the critical physical
            density given above. use the dimensionless SfEffPerFreeFall
            to calculate this */
#endif

#ifdef GALSF_EFFECTIVE_EQS
      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = REAL;
#endif


#ifdef GALSF_SUBGRID_WINDS
      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelMaxTime");
      addr[nt] = &All.WindFreeTravelMaxTimeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = REAL;

#if (GALSF_SUBGRID_WIND_SCALING>0)
      strcpy(tag[nt], "VariableWindVelFactor");
      addr[nt] = &All.VariableWindVelFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "VariableWindSpecMomentum");
      addr[nt] = &All.VariableWindSpecMomentum;
      id[nt++] = REAL;
#endif
#endif // GALSF_SUBGRID_WINDS

#endif

#ifdef GALSF_EFFECTIVE_EQS
      strcpy(tag[nt], "FactorForSofterEQS");
      addr[nt] = &All.FactorForSofterEQS;
      id[nt++] = REAL;
#endif
#ifdef GR_TABULATED_COSMOLOGY
#ifndef GR_TABULATED_COSMOLOGY_W
      strcpy(tag[nt], "DarkEnergyConstantW");
      addr[nt] = &All.DarkEnergyConstantW;
      id[nt++] = REAL;
#endif
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = REAL;
#endif

#ifdef GR_TABULATED_COSMOLOGY
#if defined(GR_TABULATED_COSMOLOGY_W) || defined(GR_TABULATED_COSMOLOGY_G) || defined(GR_TABULATED_COSMOLOGY_H)
      strcpy(tag[nt], "TabulatedCosmologyFile");
      addr[nt] = All.TabulatedCosmologyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef TURB_DIFFUSION
      strcpy(tag[nt], "TurbDiffusionCoefficient");
      addr[nt] = &All.TurbDiffusion_Coefficient;
      id[nt++] = REAL;

#ifdef TURB_DIFF_DYNAMIC
      strcpy(tag[nt], "TurbDynamicDiffFac");
      addr[nt] = &All.TurbDynamicDiffFac;
      id[nt++] = REAL;
        /*
      strcpy(tag[nt], "TurbDynamicDiffIterations");
      addr[nt] = &All.TurbDynamicDiffIterations;
      id[nt++] = INT;
         */
      strcpy(tag[nt], "TurbDynamicDiffSmoothing");
      addr[nt] = &All.TurbDynamicDiffSmoothing;
      id[nt++] = REAL;

      strcpy(tag[nt], "TurbDynamicDiffMax");
      addr[nt] = &All.TurbDynamicDiffMax;
      id[nt++] = REAL;
#endif
#endif


#if defined(CONDUCTION)
        strcpy(tag[nt], "ConductionCoeff");
        addr[nt] = &All.ConductionCoeff;
        id[nt++] = REAL;
#endif

#if defined(VISCOSITY)
        strcpy(tag[nt], "ShearViscosityCoeff");
        addr[nt] = &All.ShearViscosityCoeff;
        id[nt++] = REAL;

        strcpy(tag[nt], "BulkViscosityCoeff");
        addr[nt] = &All.BulkViscosityCoeff;
        id[nt++] = REAL;
#endif


#ifdef MAGNETIC
#ifdef MHD_B_SET_IN_PARAMS
      strcpy(tag[nt], "BiniX");
      strcpy(alternate_tag[nt], "B_initialvalue_x");
      addr[nt] = &All.BiniX;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniY");
      strcpy(alternate_tag[nt], "B_initialvalue_y");
      addr[nt] = &All.BiniY;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniZ");
      strcpy(alternate_tag[nt], "B_initialvalue_z");
      addr[nt] = &All.BiniZ;
      id[nt++] = REAL;
#endif
#endif /* MAGNETIC */

#ifdef BH_WIND_SPAWN_SET_BFIELD_POLTOR
      strcpy(tag[nt], "BH_spawn_injection_radius");
      addr[nt] = &All.BH_spawn_rinj;
      id[nt++] = REAL;

      strcpy(tag[nt], "BH_spawn_poloidal_B");
      addr[nt] = &All.B_spawn_pol;
      id[nt++] = REAL;

      strcpy(tag[nt], "BH_spawn_toroidal_B");
      addr[nt] = &All.B_spawn_tor;
      id[nt++] = REAL;
#endif
#ifdef BH_WIND_SPAWN_SET_JET_PRECESSION
      strcpy(tag[nt], "BH_jet_precession_degree");
      addr[nt] = &All.BH_jet_precess_degree;
      id[nt++] = REAL;

      strcpy(tag[nt], "BH_jet_precession_period");
      addr[nt] = &All.BH_jet_precess_period;
      id[nt++] = REAL;
#endif
#ifdef BH_DEBUG_FIX_MDOT_MBH
      strcpy(tag[nt], "BH_fb_duty_cycle");
      addr[nt] = &All.BH_fb_duty_cycle;
      id[nt++] = REAL;

      strcpy(tag[nt], "BH_fb_period");
      addr[nt] = &All.BH_fb_period;
      id[nt++] = REAL;
#endif

#ifdef EOS_TABULATED
        strcpy(tag[nt], "EosTable");
        addr[nt] = All.EosTable;
        id[nt++] = STRING;
#endif

#ifdef EOS_TILLOTSON
        strcpy(tag[nt], "Tillotson_EOS_params_a");
        addr[nt] = &All.Tillotson_EOS_params[0][0];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_b");
        addr[nt] = &All.Tillotson_EOS_params[0][1];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_u_0");
        addr[nt] = &All.Tillotson_EOS_params[0][2];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_rho_0");
        addr[nt] = &All.Tillotson_EOS_params[0][3];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_A");
        addr[nt] = &All.Tillotson_EOS_params[0][4];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_B");
        addr[nt] = &All.Tillotson_EOS_params[0][5];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_u_s");
        addr[nt] = &All.Tillotson_EOS_params[0][6];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_u_s_prime");
        addr[nt] = &All.Tillotson_EOS_params[0][7];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_alpha");
        addr[nt] = &All.Tillotson_EOS_params[0][8];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_beta");
        addr[nt] = &All.Tillotson_EOS_params[0][9];
        id[nt++] = REAL;
#endif

#ifdef EOS_ELASTIC
        strcpy(tag[nt], "Tillotson_EOS_params_mu");
        addr[nt] = &All.Tillotson_EOS_params[0][10];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_Y0");
        addr[nt] = &All.Tillotson_EOS_params[0][11];
        id[nt++] = REAL;
#endif


#ifdef NUCLEAR_NETWORK
      strcpy(tag[nt], "EosSpecies");
      addr[nt] = All.EosSpecies;
      id[nt++] = STRING;

        strcpy(tag[nt], "NetworkRates");
      addr[nt] = All.NetworkRates;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkPartFunc");
      addr[nt] = All.NetworkPartFunc;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkMasses");
      addr[nt] = All.NetworkMasses;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkWeakrates");
      addr[nt] = All.NetworkWeakrates;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkTempThreshold");
      addr[nt] = &All.NetworkTempThreshold;
      id[nt++] = REAL;
#endif

#if defined(RT_CHEM_PHOTOION) && !(defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF))
        strcpy(tag[nt], "IonizingLuminosityPerSolarMass_cgs");
        addr[nt] = &All.IonizingLuminosityPerSolarMass_cgs;
        id[nt++] = REAL;

        strcpy(tag[nt], "star_Teff");
        addr[nt] = &All.star_Teff;
        id[nt++] = REAL;
#endif

#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
        strcpy(tag[nt], "AGS_DesNumNgb");
        strcpy(alternate_tag[nt], "AdaptGravSoft_Effective_NeighborNumber");
        addr[nt] = &All.AGS_DesNumNgb;
        id[nt++] = REAL;

#ifdef DEVELOPER_MODE
        strcpy(tag[nt], "AGS_MaxNumNgbDeviation");
        addr[nt] = &All.AGS_MaxNumNgbDeviation;
        id[nt++] = REAL;
#endif
#endif

#ifdef DM_FUZZY
        strcpy(tag[nt], "FuzzyDM_Mass_in_eV");
        addr[nt] = &All.ScalarField_hbar_over_mass;
        id[nt++] = REAL;
#endif

#ifdef TURB_DRIVING

#if defined(TURB_DRIVING_SPECTRUMGRID)
        strcpy(tag[nt], "TimeBetTurbSpectrum"); // time (code) between evaluations of turb pwrspec
        strcpy(alternate_tag[nt], "TurbDrive_TimeBetTurbSpectrum");
        addr[nt] = &All.TimeBetTurbSpectrum;
        id[nt++] = REAL;
#endif

        strcpy(tag[nt], "ST_decay"); // decay time for driving-mode phase correlations
        strcpy(alternate_tag[nt], "TurbDrive_CoherenceTime");
        addr[nt] = &All.TurbDriving_Global_DecayTime;
        id[nt++] = REAL;

        strcpy(tag[nt], "ST_energy"); // energy of driving-scale modes: sets norm of turb (?)
        strcpy(alternate_tag[nt], "TurbDrive_ApproxRMSVturb");
        addr[nt] = &All.TurbDriving_Global_AccelerationPowerVariable;
        id[nt++] = REAL;

        strcpy(tag[nt], "ST_DtFreq"); // time interval for driving updates (set by hand)
        strcpy(alternate_tag[nt], "TurbDrive_TimeBetweenTurbUpdates");
        addr[nt] = &All.TurbDriving_Global_DtTurbUpdates;
        id[nt++] = REAL;

        strcpy(tag[nt], "ST_Kmin"); // minimum driving-k: should be ~2.*M_PI/All.BoxSize
        strcpy(alternate_tag[nt], "TurbDrive_MaxWavelength"); // should be <= BoxSize
        addr[nt] = &All.TurbDriving_Global_DrivingScaleKMinVar;
        id[nt++] = REAL;

        strcpy(tag[nt], "ST_Kmax"); // maximum driving-k: set to couple times Kmin or more if more cascade desired
        strcpy(alternate_tag[nt], "TurbDrive_MinWavelength"); // should be < MaxWavelength
        addr[nt] = &All.TurbDriving_Global_DrivingScaleKMaxVar;
        id[nt++] = REAL;

        strcpy(tag[nt], "ST_SolWeight"); // fractional wt of solenoidal modes (wt*curl + (1-wt)*div)
        strcpy(alternate_tag[nt], "TurbDrive_SolenoidalFraction");
        addr[nt] = &All.TurbDriving_Global_SolenoidalFraction;
        id[nt++] = REAL;

        strcpy(tag[nt], "ST_SpectForm"); // driving pwr-spec: 0=Ek~const; 1=sharp-peak at kc; 2=Ek~k^(-5/3); 3=Ek~k^-2
        strcpy(alternate_tag[nt], "TurbDrive_DrivingSpectrum");
        addr[nt] = &All.TurbDriving_Global_DrivingSpectrumKey;
        id[nt++] = INT;

        strcpy(tag[nt], "ST_Seed"); // random number seed for modes
        strcpy(alternate_tag[nt], "TurbDrive_RandomNumberSeed");
        addr[nt] = &All.TurbDriving_Global_DrivingRandomNumberKey;
        id[nt++] = INT;

        /* Andreas Bauer's paper on turbulence:
         // sub-sonic (Mach~0.3) test: //
         ST_decay        1.
         ST_energy       0.0002 (sigma=0.014)
         ST_DtFreq       0.005
         ST_Kmin         6.27
         ST_Kmax         12.57
         ST_SolWeight    1.
         ST_AmplFac      1.
         ST_Seed         42
         ST_SpectForm    2

         // trans-sonic (Mach~1.2/3.5) test: //
         ST_decay        0.5
         ST_energy       0.21 (sigma=0.21-3.0)
         ST_DtFreq       0.005
         ST_Kmin         6.27
         ST_Kmax         12.57
         ST_SolWeight    1.
         ST_AmplFac      1.
         ST_Seed         42
         ST_SpectForm    2

         // super-sonic (Mach~8.4) test: //
         ST_decay        0.05
         ST_energy       25.0 (sigma=12.247)
         ST_DtFreq       0.005
         ST_Kmin         6.27
         ST_Kmax         18.85
         ST_SolWeight    1.
         ST_AmplFac      1.
         ST_Seed         42
         ST_SpectForm    1
         */
#endif

#ifdef CHIMES
      strcpy(tag[nt], "Chimes_data_path");
      addr[nt] = ChimesDataPath;
      id[nt++] = STRING;

      strcpy(tag[nt], "PhotoIonTable");
      addr[nt] = ChimesPhotoIonTable;
      id[nt++] = STRING;

      strcpy(tag[nt], "EqAbundanceTable");
      addr[nt] = ChimesEqAbundanceTable;
      id[nt++] = STRING;

      strcpy(tag[nt], "Thermal_Evolution_On");
      addr[nt] = &All.ChimesThermEvolOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Chemistry_eqm");
      addr[nt] = &ChimesEqmMode;
      id[nt++] = INT;

      strcpy(tag[nt], "redshift_dependent_UVB_mode");
      addr[nt] = &ChimesUVBMode;
      id[nt++] = INT;

      strcpy(tag[nt], "InitIonState");
      addr[nt] = &ChimesInitIonState;
      id[nt++] = INT;

      strcpy(tag[nt], "StaticMolCooling");
      addr[nt] = &ChimesGlobalVars.StaticMolCooling;
      id[nt++] = INT;

      strcpy(tag[nt], "CellSelfShielding_On");
      addr[nt] = &ChimesGlobalVars.cellSelfShieldingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Shielding_length_factor");
      addr[nt] = &shielding_length_factor;
      id[nt++] = REAL;

      strcpy(tag[nt], "Grain_Temperature");
      addr[nt] = &Tdust_buf;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_rate");
      addr[nt] = &cr_rate;
      id[nt++] = REAL;

      strcpy(tag[nt], "max_mol_temperature");
      addr[nt] = &Tmol_buf;
      id[nt++] = REAL;

      strcpy(tag[nt], "rad_field_norm_factor");
      addr[nt] = &chimes_rad_field_norm_factor;
      id[nt++] = REAL;

      strcpy(tag[nt], "relativeTolerance");
      addr[nt] = &relTol_buf;
      id[nt++] = REAL;

      strcpy(tag[nt], "absoluteTolerance");
      addr[nt] = &absTol_buf;
      id[nt++] = REAL;

      strcpy(tag[nt], "explicitTolerance");
      addr[nt] = &expTol_buf;
      id[nt++] = REAL;

      strcpy(tag[nt], "reionisation_redshift");
      addr[nt] = &z_reion_buf;
      id[nt++] = REAL;

      strcpy(tag[nt], "scale_metal_tolerances");
      addr[nt] = &ChimesGlobalVars.scale_metal_tolerances;
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeCarbon");
      addr[nt] = &ChimesGlobalVars.element_included[0];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeNitrogen");
      addr[nt] = &ChimesGlobalVars.element_included[1];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeOxygen");
      addr[nt] = &ChimesGlobalVars.element_included[2];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeNeon");
      addr[nt] = &ChimesGlobalVars.element_included[3];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeMagnesium");
      addr[nt] = &ChimesGlobalVars.element_included[4];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeSilicon");
      addr[nt] = &ChimesGlobalVars.element_included[5];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeSulphur");
      addr[nt] = &ChimesGlobalVars.element_included[6];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeCalcium");
      addr[nt] = &ChimesGlobalVars.element_included[7];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeIron");
      addr[nt] = &ChimesGlobalVars.element_included[8];
      id[nt++] = INT;

      strcpy(tag[nt], "N_chimes_full_output_freq");
      addr[nt] = &N_chimes_full_output_freq;
      id[nt++] = INT;

      strcpy(tag[nt], "chimes_debug");
      addr[nt] = &ChimesGlobalVars.chimes_debug;
      id[nt++] = INT;

#ifdef CHIMES_STELLAR_FLUXES
      strcpy(tag[nt], "Chimes_f_esc_ion");
      addr[nt] = &All.Chimes_f_esc_ion;
      id[nt++] = REAL;

      strcpy(tag[nt], "Chimes_f_esc_G0");
      addr[nt] = &All.Chimes_f_esc_G0;
      id[nt++] = REAL;
#endif
#endif  // CHIMES

        if((fd = fopen(fname, "r")))
        {
            sprintf(buf, "%s%s", fname, "-usedvalues");
            if(!(fdout = fopen(buf, "w")))
            {
                printf("error opening file '%s' \n", buf);
                errorFlag = 1;
            }
            else
            {
                printf("Obtaining parameters from file '%s':\n", fname);
                while(!feof(fd))
                {

                    *buf = 0;
                    fgets(buf, 200, fd);
                    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                        continue;

                    if(buf1[0] == '%')
                        continue;

                    for(i = 0, j = -1; i < nt; i++)
                        if((strcmp(buf1, tag[i]) == 0) || (strcmp(buf1, alternate_tag[i]) == 0))
                        {
                            j = i;
                            tag[i][0] = 0;
                            break;
                        }

                    if(j >= 0)
                    {
                        switch (id[j])
                        {
                            case REAL:
                                *((double *) addr[j]) = atof(buf2);
                                fprintf(fdout, "%-50s%g\n", buf1, *((double *) addr[j]));
                                fprintf(stdout, "%-50s%g\n", buf1, *((double *) addr[j]));
                                break;
                            case STRING:
                                strcpy((char *) addr[j], buf2);
                                fprintf(fdout, "%-50s%s\n", buf1, buf2);
                                fprintf(stdout, "%-50s%s\n", buf1, buf2);
                                break;
                            case INT:
                                *((int *) addr[j]) = atoi(buf2);
                                fprintf(fdout, "%-50s%d\n", buf1, *((int *) addr[j]));
                                fprintf(stdout, "%-50s%d\n", buf1, *((int *) addr[j]));
                                break;
                        }
                    }
                    else
                    {
#ifdef ALLOWEXTRAPARAMS
                        fprintf(stdout, "Possible warning to be aware of from file %s:   Tag '%s' was specified, but it is being ignored -- make sure this is intended!\n", fname, buf1);
#else
                        fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
                        errorFlag = 1;
#endif
                    }
                }
                fclose(fd);
                fclose(fdout);
                printf("\n");

                i = strlen(All.OutputDir);
                if(i > 0) {if(All.OutputDir[i - 1] != '/') {strcat(All.OutputDir, "/");}}

                sprintf(buf1, "%s%s", fname, "-usedvalues");
                sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
                sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
                int ret; ret = system(buf3);
#endif
            }
        }
        else
        {
            printf("Parameter file %s not found.\n", fname);
            errorFlag = 1;
        }


        for(i = 0; i < nt; i++)
        {
            if(*tag[i])
            {
                if(strcmp("RestartFile",tag[i])==0) {strcpy((char *)addr[i],"restart"); printf("Tag %s (%s) not set in parameter file: defaulting to value = 'restart' \n",tag[i],alternate_tag[i]); continue;}
                if(strcmp("SnapshotFileBase",tag[i])==0) {strcpy((char *)addr[i],"snapshot"); printf("Tag %s (%s) not set in parameter file: defaulting to value = 'snapshot' \n",tag[i],alternate_tag[i]); continue;}
                if(All.ComovingIntegrationOn==0)
                {
                    if(strcmp("Omega0",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to [non-cosmological] value (assuming integration in flat non-expanding space with physical units) = %g \n",tag[i],alternate_tag[i],All.OmegaMatter); continue;}
                    if(strcmp("OmegaLambda",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to [non-cosmological] value (assuming integration in flat non-expanding space with physical units) = %g \n",tag[i],alternate_tag[i],All.OmegaLambda); continue;}
                    if(strcmp("HubbleParam",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to [non-cosmological] value (assuming integration in flat non-expanding space with physical units) = %g \n",tag[i],alternate_tag[i],All.HubbleParam); continue;}
                    if(strcmp("OmegaRadiation",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("OmegaBaryon",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("SofteningGasMaxPhys",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("SofteningHaloMaxPhys",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("SofteningDiskMaxPhys",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("SofteningBulgeMaxPhys",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("SofteningStarsMaxPhys",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("SofteningBndryMaxPhys",tag[i])==0) {*((double *)addr[i])=0; continue;}
                } else {
                    if(strcmp("OmegaRadiation",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to ignore radiation in expansion (Omega_r = %g) \n",tag[i],alternate_tag[i],All.OmegaRadiation); continue;}
                }
                if(All.OutputListOn==0) {
                    if(strcmp("OutputListFilename",tag[i])==0) {strcpy((char *)addr[i],"output_times_dummy.txt"); continue;}
                } else {
                    if(strcmp("TimeOfFirstSnapshot",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("TimeBetSnapshot",tag[i])==0) {*((double *)addr[i])=1.1; continue;}
                }
                if(strcmp("InitGasTemp",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume temperatures defined in ICs (=%g) \n",tag[i],alternate_tag[i],All.InitGasTemp); continue;}
                if(strcmp("MinGasTemp",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume no mininum (=%g) \n",tag[i],alternate_tag[i],All.MinGasTemp); continue;}
                if(strcmp("MinGasHsmlFractional",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume no mininum (=%g) \n",tag[i],alternate_tag[i],All.MinGasHsmlFractional); continue;}
                if(strcmp("MaxHsml",tag[i])==0) {*((double *)addr[i])=MAX_REAL_NUMBER; printf("Tag %s (%s) not set in parameter file: defaulting to assume no maximum (=%g) \n",tag[i],alternate_tag[i],All.MaxHsml); continue;}
                if(strcmp("GravityConstantInternal",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to calculating in terms of other specified units if needed (=%g) \n",tag[i],alternate_tag[i],All.G); continue;}
                if(strcmp("MinSizeTimestep",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to minimum allowed by memory table-size (=%g) \n",tag[i],alternate_tag[i],All.MinSizeTimestep); continue;}
                if(strcmp("NumFilesWrittenInParallel",tag[i])==0) {*((int *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to only main-task writes (=%d) \n",tag[i],alternate_tag[i],All.NumFilesWrittenInParallel); continue;}
                if(strcmp("NumFilesPerSnapshot",tag[i])==0) {*((int *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to single-file snapshots (=%d) \n",tag[i],alternate_tag[i],All.NumFilesPerSnapshot); continue;}
                if(strcmp("SnapFormat",tag[i])==0) {*((int *)addr[i])=3; printf("Tag %s (%s) not set in parameter file: defaulting to standard hdf5 snapshot format (=%d) \n",tag[i],alternate_tag[i],All.SnapFormat); continue;}
                if(strcmp("TimeLimitCPU",tag[i])==0) {*((double *)addr[i])=8.6e4; printf("Tag %s (%s) not set in parameter file: defaulting to 24-hours before auto-shutdown (=%g) \n",tag[i],alternate_tag[i],All.TimeLimitCPU); continue;}
                if(strcmp("CpuTimeBetRestartFile",tag[i])==0) {*((double *)addr[i])=3450.; printf("Tag %s (%s) not set in parameter file: defaulting to write restart checkpoints just under every hour (=%g) \n",tag[i],alternate_tag[i],All.CpuTimeBetRestartFile); continue;}
#if !defined(COOLING) && !defined(GALSF) && !defined(EOS_HELMHOLTZ) && !defined(EOS_ELASTIC) && !defined(EOS_TILLOTSON)
                if(strcmp("UnitLength_in_cm",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitLength_in_cm); continue;}
                if(strcmp("UnitMass_in_g",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitMass_in_g); continue;}
                if(strcmp("UnitVelocity_in_cm_per_s",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitVelocity_in_cm_per_s); continue;}
#ifdef MAGNETIC
                if(strcmp("UnitMagneticField_in_gauss",tag[i])==0) {*((double *)addr[i])=3.5449077018110318; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitMagneticField_in_gauss); continue;}
#endif
#endif
#ifdef CONDUCTION_SPITZER
                if(strcmp("ConductionCoeff",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: code was compiled with Spitzer-Braginski conductivity, so will default to calculating the physical coefficient without arbitrary re-normalization (i.e. user-specified additional coefficient/multipler=%g) \n",tag[i],alternate_tag[i],All.ConductionCoeff); continue;}
#endif
#ifdef VISCOSITY_BRAGINSKII
                if(strcmp("ShearViscosityCoeff",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: code was compiled with Spitzer-Braginski viscosity, so will default to calculating the physical coefficient without arbitrary re-normalization (i.e. user-specified additional coefficient/multipler=%g) \n",tag[i],alternate_tag[i],All.ShearViscosityCoeff); continue;}
                if(strcmp("BulkViscosityCoeff",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: code was compiled with Spitzer-Braginski viscosity, so will default to 0 bulk viscosity as defined by those physics (=%g) \n",tag[i],alternate_tag[i],All.BulkViscosityCoeff); continue;}
#endif
#ifdef TURB_DIFFUSION
                if(strcmp("TurbDiffusionCoefficient",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: code was compiled with turbulent diffusion, so will default to calculating the coefficients without arbitrary re-normalization (i.e. user-specified additional coefficient/multipler=%g) \n",tag[i],alternate_tag[i],All.TurbDiffusion_Coefficient); continue;}
#endif
#if defined(INIT_STELLAR_METALS_AGES_DEFINED)
                if(strcmp("InitMetallicity",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to zero (Z=%g) \n",tag[i],alternate_tag[i],All.InitMetallicityinSolar); continue;}
                if(strcmp("InitStellarAge",tag[i])==0) {*((double *)addr[i])=10.; printf("Tag %s (%s) not set in parameter file: defaulting to very old pre-existing stars [if any exist, otherwise this is irrelevant] (=%g Gyr) \n",tag[i],alternate_tag[i],All.InitStellarAgeinGyr); continue;}
#endif
#ifdef RT_LEBRON
                if(strcmp("PhotonMomentum_Coupled_Fraction",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to use the explicitly-resolved absorption (=%g) \n",tag[i],alternate_tag[i],All.PhotonMomentum_Coupled_Fraction); continue;}
#endif
#if defined(TURB_DRIVING)
                if(strcmp("ST_DtFreq",tag[i])==0) {*((double *)addr[i])=-1; printf("Tag %s (%s) not set in parameter file: defaulting to update turbulent driving fields every 0.01 coherence times (=%g) \n",tag[i],alternate_tag[i],All.TurbDriving_Global_DtTurbUpdates); continue;}
                if(strcmp("ST_decay",tag[i])==0) {*((double *)addr[i])=-1; printf("Tag %s (%s) not set in parameter file: defaulting to assume driving-scale mode coherence time is given by expected rms eddy turnover time ~ L_drive / rms v_turb (=%g) \n",tag[i],alternate_tag[i],All.TurbDriving_Global_DecayTime); continue;}
                if(strcmp("ST_SpectForm",tag[i])==0) {*((int *)addr[i])=2; printf("Tag %s (%s) not set in parameter file: defaulting to assume driving follows a Kolmogorov spectrum (=%d) \n",tag[i],alternate_tag[i],All.TurbDriving_Global_DrivingSpectrumKey); continue;}
                if(strcmp("ST_Seed",tag[i])==0) {*((int *)addr[i])=42; printf("Tag %s (%s) not set in parameter file: defaulting to the answer to everything (=%d) \n",tag[i],alternate_tag[i],All.TurbDriving_Global_DrivingRandomNumberKey); continue;}
                if(strcmp("ST_SolWeight",tag[i])==0) {*((double *)addr[i])=0.5; printf("Tag %s (%s) not set in parameter file: defaulting to assume the so-called natural mix of modes for pressure-free turbulence (=%g) \n",tag[i],alternate_tag[i],All.TurbDriving_Global_SolenoidalFraction); continue;}
#endif
#ifdef GALSF_FB_FIRE_AGE_TRACERS
                if(strcmp("AgeTracerEventsPerTimeBin",tag[i])==0) {*((double *)addr[i])=10; printf("Tag %s (%s) not set in parameter file: defaulting to aim for ~10 age-tracer deposition events per timebin (=%g) \n",tag[i],alternate_tag[i],All.AgeTracerRateNormalization); continue;}
#if !defined(GALSF_FB_FIRE_AGE_TRACERS_CUSTOM)
                if(strcmp("AgeTracerBinStart",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: left-edge of first age-tracer bin is early in stellar evolution (=%g Myr) \n",tag[i],alternate_tag[i],All.AgeTracerBinStart); continue;}
                if(strcmp("AgeTracerBinEnd",tag[i])==0) {*((double *)addr[i])=14000.; printf("Tag %s (%s) not set in parameter file: right-edge of last age-tracer bin is at ~t_Hubble (=%g Myr) \n",tag[i],alternate_tag[i],All.AgeTracerBinEnd); continue;}
#endif
#endif
                printf("ERROR. I miss a required value for tag '%s' (or alternate name '%s') in parameter file '%s'.\n", tag[i], alternate_tag[i], fname);
                errorFlag = 1;
            }
        }

#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
        errorFlag += read_agetracerlist(All.AgeTracerListFilename);
#endif

        if(All.OutputListOn && errorFlag == 0) {errorFlag += read_outputlist(All.OutputListFilename);} else {All.OutputListLength = 0;}
    }

    MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(errorFlag)
    {
        MPI_Finalize();
        exit(0);
    }


    /* now communicate the relevant parameters to the other processes */
    MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
#ifdef CHIMES
    if(ThisTask == 0)
    {
        ChimesGlobalVars.grain_temperature = (ChimesFloat) Tdust_buf;
        ChimesGlobalVars.T_mol = (ChimesFloat) Tmol_buf;
        ChimesGlobalVars.relativeTolerance = (ChimesFloat) relTol_buf;
        ChimesGlobalVars.absoluteTolerance = (ChimesFloat) absTol_buf;
        ChimesGlobalVars.explicitTolerance = (ChimesFloat) expTol_buf;
        ChimesGlobalVars.reionisation_redshift = (ChimesFloat) z_reion_buf;
    }
    MPI_Bcast(&ChimesGlobalVars, sizeof(struct globalVariables), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ChimesDataPath, 256 * sizeof(char), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ChimesEqAbundanceTable, 196 * sizeof(char), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ChimesPhotoIonTable, 196 * sizeof(char), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&chimes_rad_field_norm_factor, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&shielding_length_factor, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cr_rate, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N_chimes_full_output_freq, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ChimesEqmMode, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ChimesUVBMode, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ChimesInitIonState, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif


    /* ok, -NOW- we can properly read the "All" variables; we should do any if/then checks on
     them at this point. if any all variable depends on another, it must be set AFTER this point! */

#ifndef DEVELOPER_MODE
    /*
     %- PFH: these are generally not parameters that should be freely-varied. we're
     %- going to default to hard-coding them, instead, so that only development-level
     %- users are modifying them. However, if you want to set them, here are some
     %- reasonable values that you will need to insert into your parameterfile

     %---- Accuracy of time integration
     ErrTolIntAccuracy       0.010   % <0.02
     CourantFac              0.2 	% <0.40
     MaxRMSDisplacementFac   0.125	% <0.25

     %---- Tree algorithm, force accuracy, domain update frequency
     ErrTolTheta                 0.7	    % 0.7=standard
     ErrTolForceAcc              0.0025	% 0.0025=standard
     %---- Convergence error for evaluating particle volumes
     MaxNumNgbDeviation      0.05    % <<DesNumNgb (values<1 are fine)
     AGS_MaxNumNgbDeviation  2   % same, for adaptive gravsoft: can be much larger

     %--- Dedner Divergence-cleaning Parameters (for MHD)
     DivBcleaningParabolicSigma      0.2  % <1, ~0.2-0.5 needed for stability
     DivBcleaningHyperbolicSigma     1.0  % ~1

     %---------- SPH-Specific Parameters ---------------------------------
     %---- Artificial viscosity
     ArtBulkViscConst    1.0     % multiplies 'standard' AV (use 1.0)
     %---- P&M artificial conductivity (if present); normalized to Alpha_Visc:
     ArtCondConstant     0.25    % multiplies 'standard' (use 0.25-0.5)
     %---- Cullen & Dehnen viscosity suppression
     ViscosityAMin       0.05    % minimum viscosity away from shocks (>0.025)
     ViscosityAMax       2.00    % maximum viscosity in shocks (>1)
     %---- Artificial resistivity (for MHD runs)
     ArtificialResistivityMax    1.  % maximum alpha_B (~1-2) for art. res. (like art. visc)
     */

    All.CourantFac = 0.4;
    All.ErrTolIntAccuracy = 0.02;
    All.ErrTolTheta = 0.7;
    All.ErrTolForceAcc = 0.0025;
    All.MaxRMSDisplacementFac = 0.25;
    All.TimeBetStatistics = 1.0e10;
    strcpy(All.ResubmitCommand,"none");
    All.ResubmitOn = 0;
#ifdef HYDRO_SPH
    All.ArtBulkViscConst = 1.0;
#ifdef SPHAV_ARTIFICIAL_CONDUCTIVITY
    All.ArtCondConstant = 0.25;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
    All.ViscosityAMin = 0.05;
    All.ViscosityAMax = 2.00;
#endif
#ifdef SPH_TP12_ARTIFICIAL_RESISTIVITY
    All.ArtMagDispConst = 1.0;
#endif
#endif // sph
#ifdef DIVBCLEANING_DEDNER
    All.DivBcleanParabolicSigma = 0.2;
    All.DivBcleanHyperbolicSigma = 1.0;
#endif

#ifdef TURB_DIFF_DYNAMIC
    All.TurbDynamicDiffIterations = 0; /* D. Rennehan: This has NOT been tested above 0 */
#endif
    if(All.ComovingIntegrationOn) {All.ErrTolForceAcc = 0.005; All.ErrTolIntAccuracy = 0.05;}
    All.MaxNumNgbDeviation = All.DesNumNgb / 640.;
#ifdef GALSF
    All.MaxNumNgbDeviation = All.DesNumNgb / 64.;
#endif
    if(All.MaxNumNgbDeviation < 0.05) All.MaxNumNgbDeviation = 0.05;
#ifdef EOS_ELASTIC
    All.MaxNumNgbDeviation /= 20.0;
#endif
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    All.AGS_MaxNumNgbDeviation = All.AGS_DesNumNgb / 640.;
#ifdef GALSF
    All.AGS_MaxNumNgbDeviation = All.AGS_DesNumNgb / 64.;
#endif
    if(All.AGS_MaxNumNgbDeviation < 0.05) All.AGS_MaxNumNgbDeviation = 0.05;
#endif
#endif // closes DEVELOPER_MODE check //
#ifdef BH_WIND_SPAWN
    All.AGNWindID = 1913298393;       // this seems weird, but is the bitshifted version of 1234568912345 for not long IDs.
#endif

#ifdef GALSF
    All.CritOverDensity = 1000.0;
    /* this just needs to be some number >> 1, or else we get nonsense.
     In cosmological runs, star formation is not allowed below this overdensity, to prevent spurious
     star formation at very high redshifts */
#endif
#ifdef GALSF_EFFECTIVE_EQS
    All.CritPhysDensity = 0.0; /* this will be calculated by the code below */
#endif
    All.TypeOfOpeningCriterion = 1;
    /* determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion: this
     should only be changed if you -really- know what you're doing! */

#if defined(MAGNETIC) || defined(HYDRO_MESHLESS_FINITE_VOLUME) || defined(BH_WIND_SPAWN)
    if(All.CourantFac > 0.2) {All.CourantFac = 0.2;}
    /* (PFH) safety factor needed for MHD calc, because people keep using the same CFac as hydro! */
#endif

#if defined(PIC_MHD) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
    All.Grain_Internal_Density=1; All.Grain_Size_Min=1; All.Grain_Size_Max=1; All.Grain_Size_Spectrum_Powerlaw=1; /* in this case these are never used, so we treat them as dummy variables */
#endif

    /* now we're going to do a bunch of checks */
    if((All.ErrTolIntAccuracy<=0)||(All.ErrTolIntAccuracy>0.05))
    {
        if(ThisTask==0) {printf("ErrTolIntAccuracy must be >0 and <0.05 to ensure stability \n"); endrun(1);}
    }
    if((All.ErrTolTheta<=0.1)||(All.ErrTolTheta>=0.9))
    {
        if(ThisTask==0) {printf("ErrTolTheta must be >0.1 and <0.9 to ensure stability \n"); endrun(1);}
    }
    if((All.CourantFac<=0)||(All.CourantFac>0.5))
    {
        if(ThisTask==0) {printf("CourantFac must be >0 and <0.5 to ensure stability \n"); endrun(1);}
    }
    if((All.ErrTolForceAcc<=0)||(All.ErrTolForceAcc>=0.01))
    {
        if(ThisTask==0) {printf("ErrTolForceAcc must be >0 and <0.01 to ensure stability \n"); endrun(1);}
    }
    if((All.MaxRMSDisplacementFac<=0)||(All.MaxRMSDisplacementFac>0.25))
    {
        if(ThisTask==0) {printf("MaxRMSDisplacementFac must be >0 and <0.25 to ensure stability \n"); endrun(1);}
    }
#ifdef HYDRO_SPH
    if((All.ArtBulkViscConst<=0.5)||(All.ArtBulkViscConst>=2.0))
    {
        if(ThisTask==0) {printf("ArtBulkViscConst must be >0.5 and <2 to ensure stability \n"); endrun(1);}
    }
#ifdef SPHAV_ARTIFICIAL_CONDUCTIVITY
    if((All.ArtCondConstant<=0)||(All.ArtCondConstant>0.5))
    {
        if(ThisTask==0) {printf("For SPH-mode runs, ArtCondConstant must be >0 and <0.5"); endrun(1);}
    }
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
    if((All.ViscosityAMin<=0.025)||(All.ViscosityAMin>=All.ViscosityAMax)||(All.ViscosityAMin>1.0))
    {
        if(ThisTask==0) {printf("For SPH-mode runs, ViscosityAMin must be >0.025 (stability) and <MIN(1,ViscosityAMax)"); endrun(1);}
    }
    if((All.ViscosityAMax<1))
    {
        if(ThisTask==0) {printf("For SPH-mode runs, ViscosityAMax must be >1"); endrun(1);}
    }
#endif
#ifdef SPH_TP12_ARTIFICIAL_RESISTIVITY
    if((All.ArtMagDispConst<1)||(All.ArtMagDispConst>2))
    {
        if(ThisTask==0) {printf("For SPH-mode runs, ArtificialResistivityMax must be >1 and <2"); endrun(1);}
    }
#endif
#endif
#ifdef DIVBCLEANING_DEDNER
    if((All.DivBcleanParabolicSigma<0.1)||(All.DivBcleanParabolicSigma>1))
    {
        if(ThisTask==0) {printf("Divergence-Cleaning Damping Parameter DivBcleaningParabolicSigma must be >0.1 and <1"); endrun(1);}
    }
    if((All.DivBcleanHyperbolicSigma<0.5)||(All.DivBcleanHyperbolicSigma>2))
    {
        if(ThisTask==0) {printf("Divergence-Cleaning Damping Parameter DivBcleanHyperbolicSigma must be >0.5 and <2"); endrun(1);}
    }
#endif
    if((All.MaxNumNgbDeviation<=0)||(All.MaxNumNgbDeviation>0.1*All.DesNumNgb))
    {
        if(ThisTask==0) {printf("MaxNumNgbDeviation must be >0 and <0.1*DesNumNgb \n"); endrun(1);}
    }
    if(!isnan(All.DesNumNgb))
    {
        if((All.DesNumNgb<KERNEL_NMIN)||(All.DesNumNgb>KERNEL_NMAX))
        {
            if(ThisTask==0) {printf("For the kernel chosen, proper sampling and stability requires DesNumNgb must be >%d and <%d \n",KERNEL_NMIN,KERNEL_NMAX); endrun(1);}
        }
    }
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    if((All.AGS_MaxNumNgbDeviation<=0)||(All.AGS_MaxNumNgbDeviation>0.1*All.AGS_DesNumNgb))
    {
        if(ThisTask==0) {printf("AGS_MaxNumNgbDeviation must be >0 and <0.1*AGS_DesNumNgb \n"); endrun(1);}
    }
    if(!isnan(All.AGS_DesNumNgb))
    {
        if((All.AGS_DesNumNgb<KERNEL_NMIN)||(All.AGS_DesNumNgb>KERNEL_NMAX))
        {
            if(ThisTask==0) {printf("For the kernel chosen, proper sampling and stability requires AGS_DesNumNgb must be >%d and <%d \n",KERNEL_NMIN,KERNEL_NMAX); endrun(1);}
        }

    }
#endif


    for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

    if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
        if(ThisTask == 0) {printf("NumFilesWrittenInParallel MUST be a power of 2\n"); endrun(1);}
    }

    if(All.NumFilesWrittenInParallel > NTask)
    {
        if(ThisTask == 0) {printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n"); endrun(1);}
    }

#if defined(BOX_LONG_X) ||  defined(BOX_LONG_Y) || defined(BOX_LONG_Z)
#if !defined(SELFGRAVITY_OFF) && !defined(GRAVITY_NOT_PERIODIC) && (defined(BOX_PERIODIC) || defined(PMGRID))
    if(ThisTask == 0)
    {
        printf("Code was compiled with BOX_LONG_X/Y/Z and either BOX_PERIODIC or PMGRID, but not with SELFGRAVITY_OFF or GRAVITY_NOT_PERIODIC.\n");
        printf("The gravitational solver does not allow stretched-periodic boxes (cubic-box periodic or non-periodic gravity required).\n");
        endrun(1);
    }
#endif
#endif


#ifdef GR_TABULATED_COSMOLOGY_W
#ifndef GR_TABULATED_COSMOLOGY
    if(ThisTask == 0) {fprintf(stdout, "Code was compiled with GR_TABULATED_COSMOLOGY_W, but not with GR_TABULATED_COSMOLOGY; this is not allowed.\n"); endrun(1);}
#endif
#endif





#ifdef PTHREADS_NUM_THREADS
#ifdef _OPENMP
    if(ThisTask == 0) {printf("PTHREADS_NUM_THREADS is incompatible with enabling OpenMP in the compiler options \n"); endrun(1);}
#endif
#endif

#undef REAL
#undef STRING
#undef INT
#undef MAXTAGS

}


#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
int read_agetracerlist(char *fname)
{
    FILE *fd; int count,i=0; char buf[512];
    if(!(fd = fopen(fname, "r"))) {printf("can't read age tracer list in file '%s'\n", fname); return 1;}
    while(1)
    {
      if(fgets(buf, 500, fd) != buf) {break;}
      count = sscanf(buf, " %lg", &All.AgeTracerTimeBins[i]);
      if(count == 1 || count == 2)
      {
          if(i >= NUM_AGE_TRACERS+1) {PRINT_WARNING("Too many entries in age tracer list. You should increase NUM_AGE_TRACERS=%d",(int)NUM_AGE_TRACERS); endrun(314);}
          i++;
      }
    }
    if(i < NUM_AGE_TRACERS+1) {PRINT_WARNING("Not enough entries in age tracer list. Found %d entries, but we need %d\n", i, NUM_AGE_TRACERS+1); endrun(314);}
    fclose(fd);
    if(ThisTask==0) {printf("Read age tracer bin set. Found %d age tracer bin edges in age tracer list.\n", i); fflush(stdout);}
    return 0;
}
#endif


/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;
  int count, flag;
  char buf[512];

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;

  while(1)
    {
      if(fgets(buf, 500, fd) != buf)
	break;

      count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);

      if(count == 1)
	flag = 1;

      if(count == 1 || count == 2)
	{
	  if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
	    {
	      if(ThisTask == 0)
		printf("\ntoo many entries in output-list. You should increase MAXLEN_OUTPUTLIST=%d.\n",
		       (int) MAXLEN_OUTPUTLIST);
	      endrun(13);
	    }

	  All.OutputListFlag[All.OutputListLength] = flag;
	  All.OutputListLength++;
	}
    }

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i; long long ti_end;

  if(sizeof(long long) != 8)
    {if(ThisTask == 0) {printf("\nType 'long long' is not 64 bit on this platform; this will produce segfaults: need to exit.\n\n");} endrun(555);}

  if(ThisTask == 0)
    {
      printf("\n TimeMax (Time_at_End_of_Simulation) has been augmented to be larger in the parameterfile;\n");
      printf("  We need to adjust integer timeline, which perturbs all the structure of particle timesteps. Usually this is ok, but with some config flags on, your run will suddently be extremely slow (because the code cannot correctly reorder the timeline). In those cases, restarting from a snapshot is recommended.\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {if(ThisTask == 0) {printf("\n You cannot reduce TimeMax (Time_at_End_of_Simulation) in the parameterfile, in a restart [this breaks the integer timeline]. Simply stop the run when desired, instead. Quitting.\n");} endrun(556);}

  if(All.ComovingIntegrationOn) {ti_end = (long long) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);}
    else {ti_end = (long long) ((TimeMax_new - All.TimeBegin) / All.Timebase_interval);}

  while(ti_end > TIMEBASE)
  {
      All.Timebase_interval *= 2.0;
      ti_end /= 2;
      All.Ti_Current /= 2;
#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif
#ifdef TURB_DRIVING
      StTPrev /= 2;
#endif

    for(i = 0; i < NumPart; i++)
	{
        P[i].Ti_begstep /= 2;
        P[i].Ti_current /= 2;
        if(P[i].TimeBin > 0)
	    {
	      P[i].TimeBin--;
	      if(P[i].TimeBin <= 0) {printf("Attempted to restructure integer timeline but ran into an error in readjust_timebase(). The minimum timebin for particle %d has been reached -- need smaller timesteps. Exiting.\n", i); endrun(8765);}
	    }
	}
    All.Ti_nextlineofsight /= 2;
  }
  All.TimeMax = TimeMax_new;
}
