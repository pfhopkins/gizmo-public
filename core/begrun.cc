#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"


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

/*! This file was originally part of the GADGET3 code developed by
 * Volker Springel. The parser is still fundamentally the same, but the code has been modified
 * heavily by Phil Hopkins (phopkins@caltech.edu) and Mike Grudic for GIZMO. The modifications
 * mostly center on quality-of-life improvements, adding the option for default flags to exist,
 * added functionality for new modules, elimination of unnecessary variables, new extensions to
 * the original parser, implementing the DEVELOPER_MODE options, and rewriting the read order
 * to allow easier manipulation on restarts.
 */



void begrun(void)
{
  struct global_data_all_processes all;
  memset(&All, 0, sizeof(struct global_data_all_processes));
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
      printf("Size of hydro-cell structure   %d  [bytes]\n\n", (int) sizeof(struct gas_cell_data));

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

  /* Ensure the HDF5 deflate (gzip) filter is available before any IC /
   * snapshot read — registers a zlib-backed replacement only if the linked
   * HDF5 library was built without zlib (no-op otherwise). See
   * file_io/hdf5_deflate_filter.cc. */
  gizmo_register_hdf5_deflate_filter();

  mymalloc_init();

#ifdef GR_TABULATED_COSMOLOGY
#ifdef GR_TABULATED_COSMOLOGY_W
  fwa_init();
#endif
#endif

  set_units();
  set_cosmo_factors_for_current_time();
  All.Time = All.TimeBegin;

  gizmo_gpu_sync_all();  /* sync All to GPU managed copies before any EOS/cooling calls */

#ifdef COOLING
  InitCool();
#endif

#ifdef BOX_PERIODIC
  ewald_init();
#endif

    /* boxSize / boxHalf / boxSize_[XYZ] / boxHalf_[XYZ] are macros into
       All.BoxSize — no per-startup sync needed (Step 5 Phase E0). */

#ifdef BOX_SHEARING
    double L_box_towrap = All.BoxSize;
#ifdef BOX_LONG_X
    L_box_towrap *= BOX_LONG_X;
#endif
    All.Shearing_Box_Vel_Offset = BOX_SHEARING_Q * BOX_SHEARING_OMEGA_BOX_CENTER * L_box_towrap;
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




  gizmo_rng_init(&random_generator, (uint64_t)(42 + ThisTask));

#ifdef PMGRID
  if(RestartFlag != 3 && RestartFlag != 4) {long_range_init();}
#endif

#ifdef SUBFIND
  GrNr = -1;
#endif

#if defined(EOS_TABULATED) && !defined(EOS_ANEOS)
    int ierr = eos_init(All.EosTable);
    if(ierr) {printf("error initializing the eos"); endrun(1);}
#endif

#ifdef EOS_TILLOTSON
    tillotson_eos_init();
#endif

#ifdef EOS_ANEOS
    for(int aneos_k = 0; aneos_k < All.AneosNumMaterials; aneos_k++) {
        int aneos_ierr = aneos_read_table(All.AneosTableFiles[aneos_k], aneos_k);
        if(aneos_ierr) {printf("ANEOS: error loading table %d from %s\n", aneos_k, All.AneosTableFiles[aneos_k]); endrun(1);}
    }
#endif

#ifdef NUCLEAR_NETWORK
    InitNuclearNetwork();
#endif

#ifdef TURB_DRIVING
    init_turb();
#endif

#if defined(DM_SIDM)
    init_geofactor_table();
#endif

#if defined(CRFLUID_EVOLVE_SPECTRUM)
    CR_initialize_multibin_quantities(); // initialize the global variables and look-up tables //
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

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.SnapFormat = all.SnapFormat;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeOfFirstSnapshot = all.TimeOfFirstSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MinGasKernelRadiusFractional = all.MinGasKernelRadiusFractional;
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

        All.MaxKernelRadius = all.MaxKernelRadius;
        All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

        All.ErrTolForceAcc = all.ErrTolForceAcc;
        All.NumFilesPerSnapshot = all.NumFilesPerSnapshot;
        All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
        All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;
#ifdef MHD_MODIFIED_GRADIENT
        All.ActiveFractionForMGSweep = all.ActiveFractionForMGSweep;
        All.Flag_SkipMGSolve = 0; /* first timestep always runs the MG global solve */
#endif
#ifdef TWO_TEMPERATURE_PLASMA
        All.TwoTemp_InitialTeOverTgas = all.TwoTemp_InitialTeOverTgas;
        All.TwoTemp_ShockElectronFraction = all.TwoTemp_ShockElectronFraction;
#endif

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
#ifdef SINK_PARTICLES
        All.SinkEddingtonFactor = all.SinkEddingtonFactor;
        All.SeedSinkMass = all.SeedSinkMass;
        All.SinkNgbFactor = all.SinkNgbFactor;
        All.SinkMaxAccretionRadius = all.SinkMaxAccretionRadius;
        All.SinkRadiativeEfficiency = all.SinkRadiativeEfficiency;
        All.SinkFeedbackFactor = all.SinkFeedbackFactor;
#if defined(SINK_SEED_FROM_FOF) || defined(SINK_SEED_FROM_LOCALGAS)
        All.SeedSinkMassSigma = all.SeedSinkMassSigma;
        All.SeedSinkMinRedshift = all.SeedSinkMinRedshift;
#ifdef SINK_SEED_FROM_LOCALGAS
        All.SeedSinkPerUnitMass = all.SeedSinkPerUnitMass;
#endif
#endif
#ifdef SINK_ALPHADISK_ACCRETION
        All.SeedReservoirMass = all.SeedReservoirMass;
#endif
#ifdef SINK_SEED_FROM_FOF
        All.MinFoFMassForNewSeed = all.MinFoFMassForNewSeed;
#endif
#if defined(SINK_WIND_KICK) || defined(SINK_WIND_SPAWN)
        All.Sink_accreted_fraction = all.Sink_accreted_fraction;
        All.Sink_outflow_velocity = all.Sink_outflow_velocity;
#endif
#if defined(SINGLE_STAR_FB_JETS)
        All.Sink_outflow_jetlaunchvelscaling = all.Sink_outflow_jetlaunchvelscaling;
#endif
#if defined(SINK_COSMIC_RAYS)
        All.Sink_CosmicRay_Injection_Efficiency = all.Sink_CosmicRay_Injection_Efficiency;
#endif
#ifdef SINK_WIND_SPAWN
        All.Sink_outflow_temperature = all.Sink_outflow_temperature;
        All.Sink_outflow_particlemass = all.Sink_outflow_particlemass; // dangeous to change this, as it is also part of the merger criterion!
#endif
#ifdef SINK_PHOTONMOMENTUM
        All.Sink_Rad_MomentumFactor = all.Sink_Rad_MomentumFactor;
#endif
#endif // sinks
#ifdef GALSF_FB_FIRE_RT_LOCALRP
        All.RP_Local_Momentum_Renormalization = all.RP_Local_Momentum_Renormalization;
#endif
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
        All.HIIRegion_fLum_Coupled = all.HIIRegion_fLum_Coupled;
#endif
#ifdef RT_LEBRON
        All.PhotonMomentum_Coupled_Fraction = all.PhotonMomentum_Coupled_Fraction;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        All.PhotonMomentum_fUV = all.PhotonMomentum_fUV;
        All.PhotonMomentum_fOPT = all.PhotonMomentum_fOPT;
#endif
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
        All.SNe_Energy_Renormalization = all.SNe_Energy_Renormalization;
        All.StellarMassLoss_Rate_Renormalization = all.StellarMassLoss_Rate_Renormalization;
        All.StellarMassLoss_Energy_Renormalization = all.StellarMassLoss_Energy_Renormalization;
#endif
#ifdef COSMIC_RAY_FLUID
      All.CosmicRayDiffusionCoeff = all.CosmicRayDiffusionCoeff;
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
#ifdef CR_DYNAMICAL_INJECTION_IN_SNE
        All.CosmicRay_SNeFraction = all.CosmicRay_SNeFraction;
#endif

#ifdef GALSF_ISMDUSTCHEM_MODEL
        All.ISMDustChem_SNeIIDustScaling = all.ISMDustChem_SNeIIDustScaling;
        All.ISMDustChem_SNeIaDustScaling = all.ISMDustChem_SNeIaDustScaling;
        All.ISMDustChem_AGBDustScaling = all.ISMDustChem_AGBDustScaling;
        All.ISMDustChem_DustAccretionScaling = all.ISMDustChem_DustAccretionScaling;
        All.ISMDustChem_ThermalSputteringScaling = all.ISMDustChem_ThermalSputteringScaling;
        All.ISMDustChem_AccretionTcutoffScaling = all.ISMDustChem_AccretionTcutoffScaling;
        All.ISMDustChem_SNeGasClearedOfDustScaling = all.ISMDustChem_SNeGasClearedOfDustScaling;
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
        All.ISMDustChem_SNeShatteringScaling = all.ISMDustChem_SNeShatteringScaling;
        All.ISMDustChem_SNeSputteringScaling = all.ISMDustChem_SNeSputteringScaling;
        All.ISMDustChem_ShatteringScaling = all.ISMDustChem_ShatteringScaling;
        All.ISMDustChem_CoagDensityEnhancementScaling = all.ISMDustChem_CoagDensityEnhancementScaling;
        All.ISMDustChem_VCoagScaling = all.ISMDustChem_VCoagScaling;
        All.ISMDustChem_CoagulationScaling = all.ISMDustChem_CoagulationScaling;
        All.ISMDustChem_GrainVelocityScaling = all.ISMDustChem_GrainVelocityScaling;
        All.ISMDustChem_PhotodestructionScaling = all.ISMDustChem_PhotodestructionScaling;
#endif
#endif

#ifdef GR_TABULATED_COSMOLOGY
      All.DarkEnergyConstantW = all.DarkEnergyConstantW;
#endif

      All.MaxNumNgbDeviation = all.MaxNumNgbDeviation;
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
      /* Allow the tolerance over the number of neighbours to vary during the run:
        If it was initially set to a very strict value, convergence in ngb-iteration may at some point fail */
      All.AGS_MaxNumNgbDeviation = all.AGS_MaxNumNgbDeviation;
#endif

      if(All.ComovingIntegrationOn==0) {set_softenings();}

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
#ifdef NUCLEAR_NETWORK
      strcpy(All.NuclearNetworkDataFile, all.NuclearNetworkDataFile);
#endif

#if defined(EOS_TABULATED) && !defined(EOS_ANEOS)
        strcpy(All.EosTable, all.EosTable);
#endif
#ifdef EOS_ANEOS
        All.AneosNumMaterials = all.AneosNumMaterials;
        for(int aneos_k = 0; aneos_k < ANEOS_MAX_MATERIALS; aneos_k++)
            strcpy(All.AneosTableFiles[aneos_k], all.AneosTableFiles[aneos_k]);
#endif

      if(All.TimeMax != all.TimeMax) {readjust_timebase(All.TimeMax, all.TimeMax);}
    }

#ifdef GALSF_EFFECTIVE_EQS
  init_clouds();
#endif

  char contfname[DEFAULT_PATH_BUFFERSIZE_TOUSE];
  snprintf(contfname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%scont", All.OutputDir);
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
      if(P[i].Type == 0) {CellP[i].VelPred[2] = CellP[i].HydroAccel[2] = 0;}
    }
#endif
#endif

#if defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    single_star_SN_init_directions();
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
    {int i; for(i=0;i<NumPart;i++) {CellP[i].Rad_Accel = {};}}
#endif
#ifdef RADTRANSFER
#if defined(RT_EVOLVE_INTENSITIES)
    rt_init_intensity_directions();
#endif
#ifdef RT_CHEM_PHOTOION
    rt_get_sigma();
#endif
#endif
    
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
    rt_define_effective_frequencies_in_bands();
#endif
#ifdef RADTRANSFER
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

  gizmo_gpu_sync_all();  /* re-sync All to GPU after ALL initialization (init, InitCool, rt_set_simple_inits, etc.).
                            The early sync at line 100 runs before init() sets All.SolarAbundances and other fields.
                            Without this re-sync, GPU TUs see stale All_dev with SolarAbundances[0]==0,
                            causing division-by-zero in gas_dust_heating_coeff on the first kick. */

#ifdef GIZMO_NLR_ITER_HARNESS_TEST
  /* Phase 4.B.0 step 3: synthetic iterative harness opt-in.
   *
   * KNOWN HARNESS LIMITATION (codex 2026-05-10 step-3 review):
   * Iterative Mode A ghost-writeback reverse-comm CANNOT be validated
   * from this begrun-end placement: at this point the simulation has
   * not yet entered run() and the ghost-exchange/tree machinery isn't
   * fully primed, so `gizmo_request_filtered_ghost_import_fresh`
   * returns 0 ghosts regardless of where rank 0 places its actives in
   * space (probe loop over 64 rank-1 candidates verified).
   *
   * The harness's IterHarnessGhostSpec subtest detects this case via
   * its probe loop and emits VP 11 SKIP with a clear message rather
   * than a misleading PASS or FAIL. Real Mode A reverse-comm coverage
   * falls to ags_density 3d.4 (the first iterative+ghost_writeback
   * production user). All other VPs (1, 2, 3, 4, 5, 7, 8, 9, 11
   * ctx-reset-symmetry, 11 Mode B suppression+parity) ARE validated
   * by the harness here.
   *
   * Future extension to make Mode A reverse-comm testable: move this
   * invocation to core/run.cc after first-step setup, OR build a
   * uniform mini-IC for the harness. Both are non-trivial; deferred. */
  {
    extern bool gizmo_nlr_iter_harness_run_enabled(void);
    extern void run_iter_harness_tests(void);
    if (gizmo_nlr_iter_harness_run_enabled()) {
      run_iter_harness_tests();  /* never returns; endrun(0) inside */
    }
  }
#endif
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
    
#if !defined(SELFGRAVITY_OFF)
    if(ThisTask==0)
    {
        printf("\nThe chosen interaction and softening kernel function is =%d \n",(int)KERNEL_FUNCTION);
        printf("  Gravity will be exactly Newtonian beyond =%g times the given Plummer-equivalent softenings \n",(1./KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER));
        printf("\n");
    }
#endif

    double meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC); /* assumes fully-atomic otherwise */
#ifdef COOL_LOW_TEMPERATURES
    meanweight = 1. / ( HYDROGEN_MASSFRAC*0.5 + (1-HYDROGEN_MASSFRAC)/4. + 1./(16.+12.)); /* assumes fully-molecular if low-temp cooling enabled */
#endif
    All.MinEgySpec = All.MinGasTemp / (meanweight * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS);

#ifdef DISK_BETA_COOL
    All.BetaCool_u_irr = All.BetaCool_Tirr / (meanweight * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS);
#endif
#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
    All.GrainPromotion_MassThresh = All.GrainPromotion_MassThresh_cgs / UNIT_MASS_IN_CGS;
#endif
#ifdef PLANET_HEATING
    {
        double cgs_to_code_rate = UNIT_TIME_IN_CGS / UNIT_SPECEGY_IN_CGS;
        All.PlanetHeating_RadQ0  = All.PlanetHeating_RadQ0_cgs  * cgs_to_code_rate;
        All.PlanetHeating_RadTau = (All.PlanetHeating_RadTau_cgs > 0.0)
                                   ? All.PlanetHeating_RadTau_cgs / UNIT_TIME_IN_CGS : 0.0;
        All.PlanetHeating_AccQ0  = All.PlanetHeating_AccQ0_cgs  * cgs_to_code_rate;
    }
#endif


#if defined(GALSF)
  /* for historical reasons, we need to convert to "All.MaxSfrTimescale", defined as the SF timescale in code units at the critical physical
     density given above. use the dimensionless SfEffPerFreeFall (which has been read in) to calculate this. This must be done -BEFORE- calling set_units_sfr) */
#ifndef GALSF_EFFECTIVE_EQS
    All.MaxSfrTimescale = (1/All.MaxSfrTimescale) * sqrt(3.*M_PI / (32. * All.G * (All.CritPhysDensity / UNIT_DENSITY_IN_NHCGS)));
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
  char mode[2], buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
  if(RestartFlag == 0) {strcpy(mode, "w");} else {strcpy(mode, "a");}
  if(ThisTask == 0) {mkdir(All.OutputDir, 02755);}
  MPI_Barrier(MPI_COMM_WORLD);
  char prefix_char[100] = "###";

#ifdef SINK_PARTICLES /* Note: This is done by everyone [all tasks can write to these log-files], even if it might be empty */
  if(ThisTask == 0) {snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%ssink_details", All.OutputDir); mkdir(buf, 02755);}
  MPI_Barrier(MPI_COMM_WORLD);
#if defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(SINK_OUTPUT_MOREINFO)
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%ssink_details/sink_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdSinksDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif
#ifdef OUTPUT_SINK_ACCRETION_HIST
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%ssink_details/sink_swallow_%d.txt", All.OutputDir, ThisTask);
  if(!(FdSinkSwallowDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0 && ThisTask == 0) {
      fprintf(FdSinkSwallowDetails,"%s Logfiles [1 per MPI task] for sink-gas capture/swallow events. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (2,3,4) Particle ID, child ID number, ID generation number of primary (sink) \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (5) Particle mass of primary (sink) [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (6,7,8) Coordinate (x,y,z) position of primary (sink) [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (9,10,11) Particle ID, child ID number, ID generation number of secondary \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (12) Initial (pre-swallow) particle mass of secondary [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (13,14,15) Coordinate (x,y,z) position of secondary, relative to primary (Pos[sec]-Pos[prim]) [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (16,17,18) Velocity (x,y,z) difference of secondary, relative to primary (Vel[sec]-Vel[prim]) [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (19) Specific internal energy of secondary gas [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (20,21,22) Magnetic field of secondary gas [code units] \n",prefix_char);
      fprintf(FdSinkSwallowDetails,"%s   (23) Density of secondary gas [code units] \n",prefix_char);
  }
#endif
#ifdef OUTPUT_SINK_FORMATION_PROPS
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%ssink_details/sink_formation_%d.txt", All.OutputDir, ThisTask);
  if(!(FdSinkFormationDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0 && ThisTask == 0) {
      fprintf(FdSinkFormationDetails,"%s Logfiles [1 per MPI task] for sink-formation events. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (2) Particle ID becoming a sink \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (3) Mass [code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (4,5,6) x,y,z Coordinates [code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (7,8,9) x,y,z Velocity [code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (10,11,12) x,y,z Magnetic field [code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (13) Specific internal energy [code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (14) Density [physical code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (15) Effective sound Speed [code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (16) Cell linear size [physical code units] \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (17) Column density [physical code units of mass/area] integrated to infinity \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (18) Turbulent velocity variance dv*dv [physical code units] from deviatoric shear tensor \n",prefix_char);
      fprintf(FdSinkFormationDetails,"%s   (19) Minimum distance to nearest (other) sink [code units] \n",prefix_char);
  }
#endif
#ifdef SINK_OUTPUT_MOREINFO
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%ssink_details/sink_mergers_%d.txt", All.OutputDir, ThisTask);
  if(!(FdSinkMergerDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0 && ThisTask == 0) {
      fprintf(FdSinkMergerDetails,"%s Logfiles [1 per MPI task] for sink merger events. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdSinkMergerDetails,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdSinkMergerDetails,"%s   (2,3,4) Particle ID, child ID number, ID generation number of primary \n",prefix_char);
      fprintf(FdSinkMergerDetails,"%s   (5) Sink mass of primary [code units] \n",prefix_char);
      fprintf(FdSinkMergerDetails,"%s   (6,7,8) Coordinate (x,y,z) position of primary (sink) [code units] \n",prefix_char);
      fprintf(FdSinkMergerDetails,"%s   (9,10,11) Particle ID, child ID number, ID generation number of secondary \n",prefix_char);
      fprintf(FdSinkMergerDetails,"%s   (12) Sink mass of secondary [code units] \n",prefix_char);
      fprintf(FdSinkMergerDetails,"%s   (13,14,15) Coordinate (x,y,z) position of secondary [code units] \n",prefix_char);
  }
#ifdef SINK_WIND_KICK
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%ssink_details/sink_winds_%d.txt", All.OutputDir, ThisTask);
  if(!(FdSinkWindDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0 && ThisTask == 0) {
      fprintf(FdSinkWindDetails,"%s Logfiles [1 per MPI task] for sink wind-kick [SINK_WIND_KICK] events. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (2) Particle ID of kicked (gas) element \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (3) Mass of kicked element [code units] \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (4,5,6) Coordinate (x,y,z) position of kicked element [code units] \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (7,8,9) Velocity (x,y,z) of kicked element [code units] \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (10,11,12) Directional unit vector (x,y,z) of kick \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (13,14,15) Coordinate (x,y,z) position of secondary [code units] \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (16) ID of kicking sink \n",prefix_char);
      fprintf(FdSinkWindDetails,"%s   (17,18,19) Coordinate (x,y,z) position of kicking sink [code units] \n",prefix_char);
  }
#endif
#endif // output-more-info if
#endif // SINK_PARTICLES if

    if(ThisTask != 0) {return;}	/* only the root processors writes to the log files listed below */

    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "cpu.txt");
    if(!(FdCPU = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}

#ifdef OUTPUT_ADDITIONAL_RUNINFO
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "timebin.txt");
    if(!(FdTimebin = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "info.txt");
    if(!(FdInfo = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "energy.txt");
    if(!(FdEnergy = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "timings.txt");
    if(!(FdTimings = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}

    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "balance.txt");
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
    fprintf(FdBalance, "Find-timesteps = '%c' / '%c'\n", CPU_Symbol[CPU_FIND_TIMESTEPS], CPU_SymbolImbalance[CPU_FIND_TIMESTEPS]);
    fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL], CPU_SymbolImbalance[CPU_POTENTIAL]);
    fprintf(FdBalance, "PM-gravity     = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
    fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT], CPU_SymbolImbalance[CPU_SNAPSHOT]);
    fprintf(FdBalance, "Sink           = '%c' / '%c'\n", CPU_Symbol[CPU_SINKS], CPU_SymbolImbalance[CPU_SINKS]);
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
    /* Step-15-era buckets surfaced during the GPU optimization pass. */
    fprintf(FdBalance, "GPU ngb-build  = '%c' / '%c'\n", CPU_Symbol[CPU_GPU_NGB_BUILD],     CPU_SymbolImbalance[CPU_GPU_NGB_BUILD]);
    fprintf(FdBalance, "SIDX refresh   = '%c' / '%c'\n", CPU_Symbol[CPU_SIDX_REFRESH],      CPU_SymbolImbalance[CPU_SIDX_REFRESH]);
    fprintf(FdBalance, "Lazy drift hk  = '%c' / '%c'\n", CPU_Symbol[CPU_LAZY_DRIFT],        CPU_SymbolImbalance[CPU_LAZY_DRIFT]);
    fprintf(FdBalance, "GPU kernel     = '%c' / '%c'\n", CPU_Symbol[CPU_GPU_KERNEL],        CPU_SymbolImbalance[CPU_GPU_KERNEL]);
    fprintf(FdBalance, "Grav precomp   = '%c' / '%c'\n", CPU_Symbol[CPU_GRAV_PRECOMP],      CPU_SymbolImbalance[CPU_GRAV_PRECOMP]);
    fprintf(FdBalance, "Force-upd-tree = '%c' / '%c'\n", CPU_Symbol[CPU_FORCE_UPDATE_TREE], CPU_SymbolImbalance[CPU_FORCE_UPDATE_TREE]);
    fprintf(FdBalance, "Sink env       = '%c' / '%c'\n", CPU_Symbol[CPU_SINK_ENV],          CPU_SymbolImbalance[CPU_SINK_ENV]);
    fprintf(FdBalance, "Sink feed/swk  = '%c' / '%c'\n", CPU_Symbol[CPU_SINK_FEEDSWK],      CPU_SymbolImbalance[CPU_SINK_FEEDSWK]);
    fprintf(FdBalance, "\n");
#endif

#ifdef GALSF
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0) {
      fprintf(FdSfr,"%s Star formation log-file [GALSF]. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdSfr,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdSfr,"%s   (2) Expectation value of stellar mass formed this timestep [code units] \n",prefix_char);
      fprintf(FdSfr,"%s   (3) Total star formation rate [SFR] summed over gas cells [code units] \n",prefix_char);
      fprintf(FdSfr,"%s   (4) Discretized mass formed / timestep, in solar masses per year \n",prefix_char);
      fprintf(FdSfr,"%s   (5) Total mass actually formed this timestep [code units] \n",prefix_char);
  }
#endif

#ifdef GALSF_FB_FIRE_RT_LOCALRP
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "MomWinds.txt");
    if(!(FdMomWinds = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    else if(RestartFlag == 0) {
        fprintf(FdMomWinds,"%s Local Photon-momentum model log-file [GALSF_FB_FIRE_RT_LOCALRP]. See User Guide for details. Columns represent: \n",prefix_char);
        fprintf(FdMomWinds,"%s   (1) Simulation time [code units] \n",prefix_char);
        fprintf(FdMomWinds,"%s   (2) Number of cells affected by short-range radiation pressure \n",prefix_char);
        fprintf(FdMomWinds,"%s   (3) Total photon momentum (L/c x dt) for which this is being calculated (code units) in timestep \n",prefix_char);
        fprintf(FdMomWinds,"%s   (4) Total momentum actually coupled to gas [code units] \n",prefix_char);
        fprintf(FdMomWinds,"%s   (5) Average velocity of the discretized kicks assigned to gas cells from photon momentum \n",prefix_char);
        fprintf(FdMomWinds,"%s   (6) Mean infrared optical depth of incident gas being illuminated \n",prefix_char);
    }
#endif

#ifdef GALSF_FB_FIRE_RT_HIIHEATING
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "HIIheating.txt");
    if(!(FdHIIHeating = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    else if(RestartFlag == 0) {
        fprintf(FdHIIHeating,"%s Stochastic HII region model log-file [GALSF_FB_FIRE_RT_HIIHEATING]. See User Guide for details. Columns represent: \n",prefix_char);
        fprintf(FdHIIHeating,"%s   (1) Simulation time [code units] \n",prefix_char);
        fprintf(FdHIIHeating,"%s   (2) Number of sources emitting ionizing photons (very young stars active that timestep) \n",prefix_char);
        fprintf(FdHIIHeating,"%s   (3) Total ionizing photon emission in photons/s \n",prefix_char);
        fprintf(FdHIIHeating,"%s   (4) Total number of gas cells ionized (all or in part) by those stars this timestep \n",prefix_char);
        fprintf(FdHIIHeating,"%s   (5) Total gas mass (in solar) ionized this timestep \n",prefix_char);
        fprintf(FdHIIHeating,"%s   (6) Average size of the HII regions being created (distance from source to ionized gas element, in code units) \n",prefix_char);
    }
#endif

#ifdef GALSF_FB_MECHANICAL
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "MechFeedbackEvents.txt");
    if(!(FdSNeFBLogFile = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
    else if(RestartFlag == 0) {
        fprintf(FdSNeFBLogFile,"%s Mechanical feedback log-file [GALSF_FB_MECHANICAL]. See User Guide for details. Columns represent: \n",prefix_char);
        fprintf(FdSNeFBLogFile,"%s   (1) Simulation time [code units] \n",prefix_char);
        fprintf(FdSNeFBLogFile,"%s   (2) Number of active particles which could (potentially) have mechanical injection events \n",prefix_char);
        fprintf(FdSNeFBLogFile,"%s   (3) Number of particles with at least one such event this timestep \n",prefix_char);
        fprintf(FdSNeFBLogFile,"%s   (4) Total number of injection events (one particle can have multiple, and different types) \n",prefix_char);
        fprintf(FdSNeFBLogFile,"%s   (5) Integral of the rate functions over timestep: number should (on average) reflect this \n",prefix_char);
        fprintf(FdSNeFBLogFile,"%s   (6) Mean timestep of active particles [code units] \n",prefix_char);
        fprintf(FdSNeFBLogFile,"%s   (7) Mean event rates of active particles [code units] \n",prefix_char);
    }
#endif

#if defined(RT_CHEM_PHOTOION) && defined(OUTPUT_ADDITIONAL_RUNINFO)
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "rt_photoion_chem.txt");
  if(!(FdPhotoIonChemStats = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
#endif

#if defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
  FILE *FdSinkSNDetails;
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "SN_details.txt");
  if(!(FdSinkSNDetails = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0) {
      fprintf(FdSinkSNDetails,"%s Sink Supernova/Explosion log-file [SINK_PARTICLES]. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdSinkSNDetails,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdSinkSNDetails,"%s   (2) Sink ID to explode/implode \n",prefix_char);
      fprintf(FdSinkSNDetails,"%s   (3) Sink mass [code units] \n",prefix_char);
      fprintf(FdSinkSNDetails,"%s   (4,5,6) Sink x,y,z position [code units] \n",prefix_char);
      fprintf(FdSinkSNDetails,"%s   (7,8,9) Sink x,y,z velocity [code units] \n",prefix_char);
  }
  fclose(FdSinkSNDetails); /* this is special and we close it because we'll risk collisions re-opening it from each task when needed, since it's rare enough we don't need separate files for each */
#endif

#ifdef SINK_PARTICLES
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "sinks.txt");
  if(!(FdSinks = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0) {
      fprintf(FdSinks,"%s Sink particle log-file [SINK_PARTICLES]. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdSinks,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdSinks,"%s   (2) Number of sink particles in the simulation \n",prefix_char);
      fprintf(FdSinks,"%s   (3) Total mass [code units] of the sinks (using the sink central mass, not particle/reservoir masses) \n",prefix_char);
      fprintf(FdSinks,"%s   (4) Total sink accretion rate [code units] (using sink central Mdot) \n",prefix_char);
      fprintf(FdSinks,"%s   (5) Total sink accretion rate in solar masses per year \n",prefix_char);
      fprintf(FdSinks,"%s   (6) Total sink particle dynamical mass in simulation [code units] \n",prefix_char);
      fprintf(FdSinks,"%s   (7) Mean Eddington ratio of accretion onto all sinks \n",prefix_char);
  }
#endif

#if defined(TURB_DRIVING) && defined(OUTPUT_ADDITIONAL_RUNINFO)
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "turb.txt");
  if(!(FdTurb = fopen(buf, mode))) {printf("error in opening file '%s'\n", buf); endrun(1);}
  else if(RestartFlag == 0) {
      fprintf(FdTurb,"%s Turbulent driving log [TURB_DRIVING]. See User Guide for details. Columns represent: \n",prefix_char);
      fprintf(FdTurb,"%s   (1) Simulation time [code units] \n",prefix_char);
      fprintf(FdTurb,"%s   (2) Energy-weighted rms Mach number in the simulation \n",prefix_char);
      fprintf(FdTurb,"%s   (3) Mean specific total thermal plus kinetic energy [code units] \n",prefix_char);
      fprintf(FdTurb,"%s   (4) Mean specific kinetic energy injection rate per unit mass [code units] \n",prefix_char);
      fprintf(FdTurb,"%s   (5) Mean specific dissipated energy rate per unit mass [code units] \n",prefix_char);
      fprintf(FdTurb,"%s   (6) Integrated kinetic energy injection per unit mass [code units] \n",prefix_char);
      fprintf(FdTurb,"%s   (7) Integrated dissipated energy per unit mass [code units] \n",prefix_char);
  }
#endif

#if defined(GR_TABULATED_COSMOLOGY) && defined(OUTPUT_ADDITIONAL_RUNINFO)
  snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "darkenergy.txt");
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
    char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE], buf1[DEFAULT_PATH_BUFFERSIZE_TOUSE], buf2[DEFAULT_PATH_BUFFERSIZE_TOUSE], buf3[DEFAULT_PATH_BUFFERSIZE_TOUSE], tag[MAXTAGS][50], alternate_tag[MAXTAGS][50];
    int i, j, nt, id[MAXTAGS], pnum, errorFlag = 0;
    void *addr[MAXTAGS];
    double safe_memorypertask = mpi_report_comittable_memory(0,0); /* used for some checks below */
#ifdef CHIMES
    double Tdust_buf, Tmol_buf, relTol_buf, absTol_buf, expTol_buf, z_reion_buf;
#endif

    if(sizeof(long long) != 8)
    {
        if(ThisTask == 0) {printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");}
        endrun(0);
    }
    
    if(sizeof(int) != 4)
    {
        if(ThisTask == 0) {printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");}
        endrun(0);
    }
    
    if(sizeof(float) != 4)
    {
        if(ThisTask == 0) {printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");}
        endrun(0);
    }
    
    if(sizeof(double) != 8)
    {
        if(ThisTask == 0) {printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");}
        endrun(0);
    }
    

  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;
      for(j=0;j<MAXTAGS;j++) {strcpy(alternate_tag[j], "-");}

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
      strcpy(alternate_tag[nt], "Use_Automatic_Shell_Resubmission");
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
      strcpy(alternate_tag[nt], "Time_Between_Internal_Diagnostic_Statistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "TimeBegin");
      strcpy(alternate_tag[nt], "Initial_Simulation_Time");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      strcpy(alternate_tag[nt], "Final_Simulation_Time");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      strcpy(alternate_tag[nt], "ScaleFac_Between_Snapshots");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      strcpy(alternate_tag[nt], "UnitVelocity_in_centimeterspersecond");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      strcpy(alternate_tag[nt], "UnitLength_in_centimeters");
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
      strcpy(alternate_tag[nt], "TreeRebuild_ActiveFraction");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = REAL;

      strcpy(tag[nt], "LETAllocFactor");
      strcpy(alternate_tag[nt], "LET_ForeignNode_HeadroomFactor");
      addr[nt] = &All.LETAllocFactor;
      id[nt++] = REAL;

#ifdef MHD_MODIFIED_GRADIENT
      strcpy(tag[nt], "ActiveFractionForMGSweep");
      strcpy(alternate_tag[nt], "MGSweep_ActiveFraction");
      addr[nt] = &All.ActiveFractionForMGSweep;
      id[nt++] = REAL;
#endif

#ifdef TWO_TEMPERATURE_PLASMA
      strcpy(tag[nt], "TwoTemp_InitialTeOverTgas");
      strcpy(alternate_tag[nt], "TwoTemp_InitTeRatio");
      addr[nt] = &All.TwoTemp_InitialTeOverTgas;
      id[nt++] = REAL;

      strcpy(tag[nt], "TwoTemp_ShockElectronFraction");
      strcpy(alternate_tag[nt], "TwoTemp_fE_Shock");
      addr[nt] = &All.TwoTemp_ShockElectronFraction;
      id[nt++] = REAL;
#endif

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
#if !defined(PIC_MHD) || defined(GRAIN_FLUID_AND_PIC_BOTH_DEFINED)
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
        strcpy(alternate_tag[nt],"Initial_StellarAge_inICs");
        addr[nt] = &All.InitStellarAgeinGyr;
        id[nt++] = REAL;
#endif

#if defined(GALSF_ISMDUSTCHEM_MODEL)
        strcpy(tag[nt],"Initial_ISMDustChem_Depletion");
        strcpy(alternate_tag[nt],"Initial_ISMDustChem_Dust_Depletion");
        addr[nt] = &All.Initial_ISMDustChem_Depletion;
        id[nt++] = REAL;

        strcpy(tag[nt],"Initial_ISMDustChem_SiltoCarbRatio");
        strcpy(alternate_tag[nt],"Initial_ISMDustChem_Silicate_to_Carbon_Dust_Ratio");
        addr[nt] = &All.Initial_ISMDustChem_SiliconToCarbonRatio;
        id[nt++] = REAL;

        strcpy(tag[nt],"SNeIIDustScaling");
        addr[nt] = &All.ISMDustChem_SNeIIDustScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"SNeIaDustScaling");
        addr[nt] = &All.ISMDustChem_SNeIaDustScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"AGBDustScaling");
        addr[nt] = &All.ISMDustChem_AGBDustScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"DustAccretionScaling");
        addr[nt] = &All.ISMDustChem_DustAccretionScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"ThermalSputteringScaling");
        addr[nt] = &All.ISMDustChem_ThermalSputteringScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"SNeGasClearedOfDustScaling");
        addr[nt] = &All.ISMDustChem_SNeGasClearedOfDustScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"AccretionTcutoffScaling");
        addr[nt] = &All.ISMDustChem_AccretionTcutoffScaling;
        id[nt++] = REAL;        

#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)

        strcpy(tag[nt], "UnitGrainNumber");
        addr[nt] = &All.UnitGrainNumber;
        id[nt++] = REAL;

        strcpy(tag[nt], "UnitGrainLength_in_cm");
        addr[nt] = &All.UnitGrainLength_in_cm;
        id[nt++] = REAL;

        strcpy(tag[nt],"ISMDustChem_Grain_Size_Min");
        addr[nt] = &All.ISMDustChem_Grain_Size_Min;
        id[nt++] = REAL;    

        strcpy(tag[nt],"ISMDustChem_Grain_Size_Max");
        addr[nt] = &All.ISMDustChem_Grain_Size_Max;
        id[nt++] = REAL;    

        strcpy(tag[nt],"SNeShatteringScaling");
        addr[nt] = &All.ISMDustChem_SNeShatteringScaling;
        id[nt++] = REAL; 

        strcpy(tag[nt],"SNeSputteringScaling");
        addr[nt] = &All.ISMDustChem_SNeSputteringScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"GrainShatteringScaling");
        addr[nt] = &All.ISMDustChem_ShatteringScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"GrainCoagulationScaling");
        addr[nt] = &All.ISMDustChem_CoagulationScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"VCoagScaling");
        addr[nt] = &All.ISMDustChem_VCoagScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"CoagDensityEnhancementScaling");
        addr[nt] = &All.ISMDustChem_CoagDensityEnhancementScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"GrainVelocityScaling");
        addr[nt] = &All.ISMDustChem_GrainVelocityScaling;
        id[nt++] = REAL;

        strcpy(tag[nt],"PhotodestructionScaling");
        addr[nt] = &All.ISMDustChem_PhotodestructionScaling;
        id[nt++] = REAL;
#endif
#endif

#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
        strcpy(tag[nt], "SNeIIEnergyFrac");
        strcpy(alternate_tag[nt], "SNe_Energy_Renormalization");
        addr[nt] = &All.SNe_Energy_Renormalization;
        id[nt++] = REAL;

        strcpy(tag[nt],"GasReturnFraction");
        strcpy(alternate_tag[nt], "StellarMassLoss_Rate_Renormalization");
        addr[nt] = &All.StellarMassLoss_Rate_Renormalization;
        id[nt++] = REAL;

        strcpy(tag[nt],"GasReturnEnergy");
        strcpy(alternate_tag[nt], "StellarMassLoss_Energy_Renormalization");
        addr[nt] = &All.StellarMassLoss_Energy_Renormalization;
        id[nt++] = REAL;
#endif

#ifdef CR_DYNAMICAL_INJECTION_IN_SNE
        strcpy(tag[nt], "CosmicRay_SNeFraction");
        strcpy(alternate_tag[nt], "CosmicRay_SNeEnergyFraction");
        addr[nt] = &All.CosmicRay_SNeFraction;
        id[nt++] = REAL;
#endif
        
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        strcpy(tag[nt], "CosmicRay_Subgrid_Kappa_0");
        strcpy(alternate_tag[nt], "CosmicRay_Subgrid_Kappa_0");
        addr[nt] = &All.CosmicRay_Subgrid_Kappa_0;
        id[nt++] = REAL;

        strcpy(tag[nt], "CosmicRay_Subgrid_Vstream_0");
        strcpy(alternate_tag[nt], "CosmicRay_Subgrid_Vstream_0");
        addr[nt] = &All.CosmicRay_Subgrid_Vstream_0;
        id[nt++] = REAL;
#endif
        

#ifdef GALSF_FB_FIRE_RT_LOCALRP
        strcpy(tag[nt], "WindMomentumLoading");
        strcpy(alternate_tag[nt], "RP_Local_Momentum_Renormalization");
        addr[nt] = &All.RP_Local_Momentum_Renormalization;
        id[nt++] = REAL;
#endif

#ifdef GALSF_FB_FIRE_RT_HIIHEATING
        strcpy(tag[nt], "HIIRegion_fLum_Coupled");
        addr[nt] = &All.HIIRegion_fLum_Coupled;
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

#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        strcpy(tag[nt], "PhotonMomentum_fUV");
        strcpy(alternate_tag[nt], "PhotonMomentum_min_fUV_in_SED");
        addr[nt] = &All.PhotonMomentum_fUV;
        id[nt++] = REAL;

        strcpy(tag[nt], "PhotonMomentum_fOPT");
        strcpy(alternate_tag[nt], "PhotonMomentum_min_fOPT_in_SED");
        addr[nt] = &All.PhotonMomentum_fOPT;
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


        strcpy(tag[nt], "MinGasKernelRadiusFractional");
        strcpy(alternate_tag[nt], "Minimum_Gas_KernelRadius_RelativetoSoftening");
        addr[nt] = &All.MinGasKernelRadiusFractional;
        id[nt++] = REAL;

        strcpy(tag[nt], "MaxKernelRadius");
        strcpy(alternate_tag[nt], "Maximum_KernelRadius_CodeUnits");
        addr[nt] = &All.MaxKernelRadius;
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

#ifdef NUCLEAR_NETWORK
        strcpy(tag[nt], "NuclearNetworkDataFile");
        addr[nt] = All.NuclearNetworkDataFile;
        id[nt++] = STRING;

        strcpy(tag[nt], "NuclearBurningFloor_T");
        addr[nt] = &All.NuclearBurningFloor_T;
        id[nt++] = REAL;

        strcpy(tag[nt], "NuclearBurningFloor_rho");
        addr[nt] = &All.NuclearBurningFloor_rho;
        id[nt++] = REAL;

        strcpy(tag[nt], "NuclearNSE_T_threshold");
        addr[nt] = &All.NuclearNSE_T_threshold;
        id[nt++] = REAL;
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
      strcpy(alternate_tag[nt], "GravityConstant_inCodeUnits");
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

#ifdef DISK_BETA_COOL
      strcpy(tag[nt], "BetaCool_Beta");
      addr[nt] = &All.BetaCool_Beta;
      id[nt++] = REAL;

      strcpy(tag[nt], "BetaCool_Tirr");
      addr[nt] = &All.BetaCool_Tirr;
      id[nt++] = REAL;
#endif
#ifdef PLANET_HEATING
      strcpy(tag[nt], "PlanetHeating_RadQ0_cgs");
      addr[nt] = &All.PlanetHeating_RadQ0_cgs;
      id[nt++] = REAL;

      strcpy(tag[nt], "PlanetHeating_RadTau_cgs");
      addr[nt] = &All.PlanetHeating_RadTau_cgs;
      id[nt++] = REAL;

      strcpy(tag[nt], "PlanetHeating_AccQ0_cgs");
      addr[nt] = &All.PlanetHeating_AccQ0_cgs;
      id[nt++] = REAL;
#endif

#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
      strcpy(tag[nt], "GrainPromotion_MassThresh_cgs");
      addr[nt] = &All.GrainPromotion_MassThresh_cgs;
      id[nt++] = REAL;

      strcpy(tag[nt], "GrainPromotion_DustGasRatioThresh");
      addr[nt] = &All.GrainPromotion_DustGasRatioThresh;
      id[nt++] = REAL;
#endif
#ifdef GRAIN_EVOLUTION
      strcpy(tag[nt], "GrainEvolution_StickingCoeff");
      addr[nt] = &All.GrainEvolution_StickingCoeff;
      id[nt++] = REAL;

      strcpy(tag[nt], "GrainEvolution_VelThreshFrag");
      addr[nt] = &All.GrainEvolution_VelThreshFrag;
      id[nt++] = REAL;

      strcpy(tag[nt], "GrainEvolution_VelThreshShat");
      addr[nt] = &All.GrainEvolution_VelThreshShat;
      id[nt++] = REAL;

      strcpy(tag[nt], "GrainEvolution_ThermalSputteringScaling");
      addr[nt] = &All.GrainEvolution_ThermalSputteringScaling;
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
        strcpy(tag[nt], "CosmicRayDiffusionCoeff");
        addr[nt] = &All.CosmicRayDiffusionCoeff;
        id[nt++] = REAL;
#endif


#if (defined(SINK_PARTICLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
      strcpy(tag[nt], "TimeBetOnTheFlyFoF");
      addr[nt] = &All.TimeBetOnTheFlyFoF;
      id[nt++] = REAL;
#endif

#ifdef SINK_PARTICLES
        strcpy(tag[nt], "SinkAccretionFactor");
        addr[nt] = &All.SinkAccretionFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "SinkEddingtonFactor");
        addr[nt] = &All.SinkEddingtonFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "SeedSinkMass");
        addr[nt] = &All.SeedSinkMass;
        id[nt++] = REAL;

        strcpy(tag[nt], "SinkNgbFactor");
        addr[nt] = &All.SinkNgbFactor;
        id[nt++] = REAL;

        strcpy(tag[nt], "SinkMaxAccretionRadius");
        addr[nt] = &All.SinkMaxAccretionRadius;
        id[nt++] = REAL;

        strcpy(tag[nt], "SinkRadiativeEfficiency");
        addr[nt] = &All.SinkRadiativeEfficiency;
        id[nt++] = REAL;

        strcpy(tag[nt], "SinkFeedbackFactor");
        addr[nt] = &All.SinkFeedbackFactor;
        id[nt++] = REAL;

#if defined(SINK_SEED_FROM_FOF) || defined(SINK_SEED_FROM_LOCALGAS)
        strcpy(tag[nt], "SeedSinkMassSigma");
        addr[nt] = &All.SeedSinkMassSigma;
        id[nt++] = REAL;

        strcpy(tag[nt], "SeedSinkMinRedshift");
        addr[nt] = &All.SeedSinkMinRedshift;
        id[nt++] = REAL;

#ifdef SINK_SEED_FROM_LOCALGAS
        strcpy(tag[nt], "SeedSinkPerUnitMass");
        addr[nt] = &All.SeedSinkPerUnitMass;
        id[nt++] = REAL;
#endif
#endif

#ifdef SINK_ALPHADISK_ACCRETION
        strcpy(tag[nt], "SeedReservoirMass");
        addr[nt] = &All.SeedReservoirMass;
        id[nt++] = REAL;
#endif

#ifdef SINK_SEED_FROM_FOF
        strcpy(tag[nt], "MinFoFMassForNewSeed");
        addr[nt] = &All.MinFoFMassForNewSeed;
        id[nt++] = REAL;
#endif

#if defined(SINK_WIND_KICK) || defined(SINK_WIND_SPAWN)
        strcpy(tag[nt],"Sink_accreted_fraction");
        strcpy(alternate_tag[nt], "Sink_f_accretion");
        addr[nt] = &All.Sink_accreted_fraction;
        id[nt++] = REAL;

        strcpy(tag[nt],"Sink_outflow_velocity");
        strcpy(alternate_tag[nt], "Sink_v_outflow");
        addr[nt] = &All.Sink_outflow_velocity;
        id[nt++] = REAL;
#endif

#if defined(SINGLE_STAR_FB_JETS)
        strcpy(tag[nt],"Sink_outflow_jetlaunchvelscaling");
        strcpy(alternate_tag[nt], "Sink_f_launch_v");
        addr[nt] = &All.Sink_outflow_jetlaunchvelscaling;
        id[nt++] = REAL;
#endif

#if defined(SINK_COSMIC_RAYS)
        strcpy(tag[nt],"Sink_CosmicRay_Injection_Efficiency");
        strcpy(alternate_tag[nt], "Sink_CosmicRay_Injection_Efficiency");
        addr[nt] = &All.Sink_CosmicRay_Injection_Efficiency;
        id[nt++] = REAL;
#endif


#ifdef SINK_WIND_SPAWN
        strcpy(tag[nt], "Sink_outflow_temperature");
        strcpy(alternate_tag[nt], "Cell_Spawn_Temperature");
        addr[nt] = &All.Sink_outflow_temperature;
        id[nt++] = REAL;
        strcpy(tag[nt], "Sink_outflow_particlemass");
        strcpy(alternate_tag[nt], "Cell_Spawn_Mass_ratio");
        addr[nt] = &All.Sink_outflow_particlemass;
        id[nt++] = REAL;
#ifdef SINGLE_STAR_FB_WINDS
        strcpy(tag[nt], "Cell_Spawn_Mass_ratio_MS");
        strcpy(alternate_tag[nt], "Sink_outflow_particlemass_MS");
        addr[nt] = &All.Cell_Spawn_Mass_ratio_MS;
        id[nt++] = REAL;
#endif
#endif

#ifdef SINK_PHOTONMOMENTUM
        strcpy(tag[nt],"Sink_FluxMomentumFactor");
        strcpy(alternate_tag[nt], "Sink_Rad_MomentumFactor");
        addr[nt] = &All.Sink_Rad_MomentumFactor;
        id[nt++] = REAL;
#endif

#endif /* SINK_PARTICLES */


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
      strcpy(alternate_tag[nt], "B_initial_x");
      addr[nt] = &All.BiniX;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniY");
      strcpy(alternate_tag[nt], "B_initial_y");
      addr[nt] = &All.BiniY;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniZ");
      strcpy(alternate_tag[nt], "B_initial_z");
      addr[nt] = &All.BiniZ;
      id[nt++] = REAL;
#endif
#endif /* MAGNETIC */

#ifdef SINK_WIND_SPAWN_SET_BFIELD_POLTOR
      strcpy(tag[nt], "Sink_spawn_injection_radius");
      addr[nt] = &All.Sink_spawn_injectionradius;
      id[nt++] = REAL;

      strcpy(tag[nt], "Sink_spawn_poloidal_B");
      addr[nt] = &All.B_spawn_pol;
      id[nt++] = REAL;

      strcpy(tag[nt], "Sink_spawn_toroidal_B");
      addr[nt] = &All.B_spawn_tor;
      id[nt++] = REAL;
#endif
#ifdef SINK_WIND_SPAWN_SET_JET_PRECESSION
      strcpy(tag[nt], "Sink_jet_precession_degree");
      addr[nt] = &All.Sink_jet_precess_degree;
      id[nt++] = REAL;

      strcpy(tag[nt], "Sink_jet_precession_period");
      addr[nt] = &All.Sink_jet_precess_period;
      id[nt++] = REAL;
#endif

#if defined(EOS_TABULATED) && !defined(EOS_ANEOS)
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

#ifdef EOS_ANEOS
        strcpy(tag[nt], "AneosNumMaterials");
        addr[nt] = &All.AneosNumMaterials;
        id[nt++] = INT;

        for(int aneos_k = 0; aneos_k < ANEOS_MAX_MATERIALS; aneos_k++) {
            char aneos_tagname[64];
            snprintf(aneos_tagname, sizeof(aneos_tagname), "AneosTable%d", aneos_k);
            strcpy(tag[nt], aneos_tagname);
            addr[nt] = All.AneosTableFiles[aneos_k];
            id[nt++] = STRING;
        }
#endif

#ifdef EOS_ELASTIC
        strcpy(tag[nt], "Tillotson_EOS_params_mu");
        addr[nt] = &All.Tillotson_EOS_params[0][10];
        id[nt++] = REAL;

        strcpy(tag[nt], "Tillotson_EOS_params_Y0");
        addr[nt] = &All.Tillotson_EOS_params[0][11];
        id[nt++] = REAL;
#endif

#if defined(RT_CHEM_PHOTOION) && !(defined(GALSF_FB_FIRE_RT_HIIHEATING) || defined(GALSF))
        strcpy(tag[nt], "IonizingLuminosityPerSolarMass_cgs");
        addr[nt] = &All.IonizingLuminosityPerSolarMass_cgs;
        id[nt++] = REAL;

        strcpy(tag[nt], "star_Teff");
        addr[nt] = &All.star_Teff;
        id[nt++] = REAL;
#endif

#ifdef RT_ISRF_BACKGROUND
	strcpy(tag[nt], "InterstellarRadiationFieldStrength");
	strcpy(alternate_tag[nt], "ISRF");
        addr[nt] = &All.InterstellarRadiationFieldStrength;
        id[nt++] = REAL;

	strcpy(tag[nt], "Redshift_RT_Background");
	strcpy(alternate_tag[nt], "z_background");
        addr[nt] = &All.RadiationBackgroundRedshift;
        id[nt++] = REAL;
#endif
#ifdef RT_INFRARED
    strcpy(tag[nt], "InitRadiationTemp");
	strcpy(alternate_tag[nt], "InitRadTemp");
    addr[nt] = &All.InitRadiationTemp;
    id[nt++] = REAL;
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
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

        strcpy(tag[nt], "ST_SpectForm"); // driving pwr-spec: 0=Ek~k^-1; 1=sharp-peak at kc; 2=Ek~k^(-5/3); 3=Ek~k^-2
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
            snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", fname, "-usedvalues");
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
                        fprintf(stdout, "Possible warning to be aware of from file %s:   Tag '%s' was specified, but it is being ignored -- make sure this is intended!\n", fname, buf1);
                    }
                }
                fclose(fd);
                fclose(fdout);
                printf("\n");

                i = strlen(All.OutputDir);
                if(i > 0) {if(All.OutputDir[i - 1] != '/') {strcat(All.OutputDir, "/");}}

                snprintf(buf1, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", fname, "-usedvalues");
                snprintf(buf2, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s%s", All.OutputDir, "parameters-usedvalues");
                /* silence cp stderr when output dir does not yet exist; mkdir -p ensures the dest exists */
                snprintf(buf3, DEFAULT_PATH_BUFFERSIZE_TOUSE, "mkdir -p %s 2>/dev/null; cp %s %s 2>/dev/null", All.OutputDir, buf1, buf2);
                int ret; ret = system(buf3);
            }
        }
        else
        {
            printf("Parameter file %s not found.\n", fname);
            errorFlag = 1;
        }

        /* do two loops over the parameters, first to check 'high-level' parameters that others depend on for defaults in the second check, because we aren't gauranteed a certain desired order of parameter checks */
        for(i = 0; i < nt; i++)
        {
            if(*tag[i])
            {
                if(strcmp("ComovingIntegrationOn",tag[i])==0) {*((int *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume this is a NON-cosmological (static-spacetime) simulation (=%d) \n",tag[i],alternate_tag[i],All.ComovingIntegrationOn); continue;}
                if(strcmp("OutputListOn",tag[i])==0) {*((int *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume the snapshots will be output at user-specified intervals rather than adopting a list from a given file (=%d) \n",tag[i],alternate_tag[i],All.OutputListOn); continue;}
                if(strcmp("BufferSize",tag[i])==0) {*((double *)addr[i])=100; printf("Tag %s (%s) not set in parameter file: defaulting to 100MB MPI Buffer allocation (adjust if needed for your machine) (=%g) \n",tag[i],alternate_tag[i],All.BufferSize); continue;}
            }
        }
        /* ok now safe to do a loop that has statements (e.g. if's) which can depend on some of the variables above */
        for(i = 0; i < nt; i++)
        {
            if(*tag[i])
            {
                /* skip the flags we already addressed above so we dont throw an unintentional error */
                if(strcmp("ComovingIntegrationOn",tag[i])==0) {continue;}
                if(strcmp("OutputListOn",tag[i])==0) {continue;}
                if(strcmp("BufferSize",tag[i])==0) {continue;}
                /* now move on to the flags we have not addressed already */
                //if(strcmp("MaxMemSize",tag[i])==0) {*((int *)addr[i])=(int)(0.99*safe_memorypertask); printf("Tag %s (%s) not set in parameter file: We will try to assign a memory-per-MPI task according to the number of MPI tasks in total and average number of tasks per node, and the minimum available memory per node. This gives %d MB per task. Depending on your configuration, and system memory overhead, you may need to increase or decrease this.\n",tag[i],alternate_tag[i],All.MaxMemSize); continue;}
                /* below more like what is needed for safe runs on Frontera, notoriously picky about memory for the system */
                //if(strcmp("MaxMemSize",tag[i])==0) {*((int *)addr[i])=(int)(0.93*safe_memorypertask-All.BufferSize); printf("Tag %s (%s) not set in parameter file: We will try to assign a memory-per-MPI task according to the number of MPI tasks in total and average number of tasks per node, and the minimum available memory per node. This gives %d MB per task. Depending on your configuration, and system memory overhead, you may need to increase or decrease this.\n",tag[i],alternate_tag[i],All.MaxMemSize); continue;}
                if(strcmp("MaxMemSize",tag[i])==0) {*((int *)addr[i])=(int)(0.90*safe_memorypertask); printf("Tag %s (%s) not set in parameter file: We will try to assign a memory-per-MPI task according to the number of MPI tasks in total and average number of tasks per node, and the minimum available memory per node. This gives %d MB per task. Depending on your configuration, and system memory overhead, you may need to increase or decrease this.\n",tag[i],alternate_tag[i],All.MaxMemSize); continue;}
                if(strcmp("ICFormat",tag[i])==0) {*((int *)addr[i])=3; printf("Tag %s (%s) not set in parameter file: defaulting to standard hdf5 ICs format (=%d) - change this if needed for your ICs (many codes generate ICs in the old GADGET unformatted binary format, which requires value=1 here) \n",tag[i],alternate_tag[i],All.ICFormat); continue;}
                if(strcmp("NumFilesWrittenInParallel",tag[i])==0) {*((int *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to only main-task writes (=%d) \n",tag[i],alternate_tag[i],All.NumFilesWrittenInParallel); continue;}
                if(strcmp("NumFilesPerSnapshot",tag[i])==0) {*((int *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to single-file snapshots (=%d) \n",tag[i],alternate_tag[i],All.NumFilesPerSnapshot); continue;}
                if(strcmp("SnapFormat",tag[i])==0) {*((int *)addr[i])=3; printf("Tag %s (%s) not set in parameter file: defaulting to standard hdf5 snapshot format (=%d) \n",tag[i],alternate_tag[i],All.SnapFormat); continue;}
                if(strcmp("RestartFile",tag[i])==0) {strcpy((char *)addr[i],"restart"); printf("Tag %s (%s) not set in parameter file: defaulting to value = 'restart' \n",tag[i],alternate_tag[i]); continue;}
                if(strcmp("SnapshotFileBase",tag[i])==0) {strcpy((char *)addr[i],"snapshot"); printf("Tag %s (%s) not set in parameter file: defaulting to value = 'snapshot' \n",tag[i],alternate_tag[i]); continue;}
                if(All.ComovingIntegrationOn==0)
                {
                    if(strcmp("TimeBegin",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: since this is a non-cosmological simulation, defaulting to assume time starts at t = %g \n",tag[i],alternate_tag[i],All.TimeBegin); continue;}
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
                    if(strcmp("TimeOfFirstSnapshot",tag[i])==0) {*((double *)addr[i])=(All.TimeBegin); printf("Tag %s (%s) not set in parameter file: defaulting to All.TimeBegin (=%g) \n",tag[i],alternate_tag[i],All.TimeOfFirstSnapshot); continue;}
                    if(All.ComovingIntegrationOn==0)
                    {
                        if(strcmp("MaxSizeTimestep",tag[i])==0) {*((double *)addr[i])=(DMIN(1.e-3*All.TimeMax,1.e-2*All.TimeBetSnapshot)); printf("Tag %s (%s) not set in parameter file: it is often unsafe to not set a maximum timestep. We will default to assume 0.1 percent of the maximum time or 1 percent of the time between snapshots (=%g), but this may need to be set lower. \n",tag[i],alternate_tag[i],All.MaxSizeTimestep); continue;}
                    } else {
                        if(strcmp("MaxSizeTimestep",tag[i])==0) {*((double *)addr[i])=(DMIN(1.e-3*All.TimeMax,1.e-2*(All.TimeBetSnapshot-1.))); printf("Tag %s (%s) not set in parameter file: it is often unsafe to not set a maximum timestep. We will default to assume 0.1 percent of the maximum time or 1 percent of the time between snapshots (=%g), but this may need to be set lower. \n",tag[i],alternate_tag[i],All.MaxSizeTimestep); continue;}
                    }
                } else {
                    if(strcmp("TimeOfFirstSnapshot",tag[i])==0) {*((double *)addr[i])=0; continue;}
                    if(strcmp("TimeBetSnapshot",tag[i])==0) {*((double *)addr[i])=1.1; continue;}
                    if(strcmp("MaxSizeTimestep",tag[i])==0) {*((double *)addr[i])=(1.e-3*All.TimeMax); printf("Tag %s (%s) not set in parameter file: it is often unsafe to not set a maximum timestep. We will default to assume 0.1 percent of the maximum time (=%g), but this may need to be set lower. \n",tag[i],alternate_tag[i],All.MaxSizeTimestep); continue;}
                }
                if(strcmp("DesNumNgb",tag[i])==0) {*((double *)addr[i])=(0.5*(KERNEL_NMIN+KERNEL_NMAX)); printf("Tag %s (%s) not set in parameter file: you did not set a target effective neighbor number for the interaction kernel. Trying to set a reasonable guess of =%g based on the kernel specified, but PLEASE CHECK that this is intended and experiment with different values or set your own for safety. \n",tag[i],alternate_tag[i],All.DesNumNgb); continue;}
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
                if(strcmp("AGS_DesNumNgb",tag[i])==0) {*((double *)addr[i])=(0.5*(KERNEL_NMIN+KERNEL_NMAX)); printf("Tag %s (%s) not set in parameter file: you did not set a target effective neighbor number for the adaptive-gravity (non-fluid) interaction kernel. Trying to set a reasonable guess (=%g) based on the kernel specified, but PLEASE CHECK that this is intended and experiment with different values or set your own for safety. \n",tag[i],alternate_tag[i],All.AGS_DesNumNgb); continue;}
#endif
                if(strcmp("InitGasTemp",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume temperatures defined in ICs (=%g) \n",tag[i],alternate_tag[i],All.InitGasTemp); continue;}
                if(strcmp("MinGasTemp",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume no mininum (=%g) \n",tag[i],alternate_tag[i],All.MinGasTemp); continue;}
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                if(strcmp("MinGasKernelRadiusFractional",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to assume minimum matches the specified force softening (=%g) \n",tag[i],alternate_tag[i],All.MinGasKernelRadiusFractional); continue;}
#else
                if(strcmp("MinGasKernelRadiusFractional",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to assume no mininum (=%g) \n",tag[i],alternate_tag[i],All.MinGasKernelRadiusFractional); continue;}
#endif
                if(strcmp("TreeDomainUpdateFrequency",tag[i])==0) {*((double *)addr[i])=0.005; printf("Tag %s (%s) not set in parameter file: defaulting to guess that we should re-build whenever 0.5 percent of the system is active. But this should be adjusted manually for performance and accuracy in most cases (=%g) \n",tag[i],alternate_tag[i],All.TreeDomainUpdateFrequency); continue;}
                if(strcmp("LETAllocFactor",tag[i])==0) {*((double *)addr[i])=1.0; printf("Tag %s (%s) not set in parameter file: defaulting to 1.0 (LET active with 1x MaxNodes foreign headroom; increase if LET unpack overflows) (=%g) \n",tag[i],alternate_tag[i],All.LETAllocFactor); continue;}
#ifdef MHD_MODIFIED_GRADIENT
                if(strcmp("ActiveFractionForMGSweep",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to run MG global solve when any gas is active (=%g) \n",tag[i],alternate_tag[i],All.ActiveFractionForMGSweep); continue;}
#endif
#ifdef TWO_TEMPERATURE_PLASMA
                if(strcmp("TwoTemp_InitialTeOverTgas",tag[i])==0) {*((double *)addr[i])=1.0; printf("Tag %s (%s) not set in parameter file: defaulting to LTE seed (T_e = T_gas) at startup (=%g) \n",tag[i],alternate_tag[i],All.TwoTemp_InitialTeOverTgas); continue;}
                if(strcmp("TwoTemp_ShockElectronFraction",tag[i])==0) {*((double *)addr[i])=0.0; printf("Tag %s (%s) not set in parameter file: defaulting to collisionless-shock limit f_e=0 (all hydro dissipation to ions; electrons heat only via Spitzer Coulomb exchange) (=%g) \n",tag[i],alternate_tag[i],All.TwoTemp_ShockElectronFraction); continue;}
#endif
                if(strcmp("MaxKernelRadius",tag[i])==0) {*((double *)addr[i])=MAX_REAL_NUMBER; printf("Tag %s (%s) not set in parameter file: defaulting to assume no maximum (=%g) \n",tag[i],alternate_tag[i],All.MaxKernelRadius); continue;}
                if(strcmp("GravityConstantInternal",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to calculating in terms of other specified units if needed (=%g) \n",tag[i],alternate_tag[i],All.G); continue;}
                if(strcmp("MinSizeTimestep",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to minimum allowed by memory table-size (=%g) \n",tag[i],alternate_tag[i],All.MinSizeTimestep); continue;}
                if(strcmp("TimeLimitCPU",tag[i])==0) {*((double *)addr[i])=8.6e4; printf("Tag %s (%s) not set in parameter file: defaulting to 24-hours before auto-shutdown (=%g) \n",tag[i],alternate_tag[i],All.TimeLimitCPU); continue;}
                if(strcmp("PartAllocFactor",tag[i])==0) {*((double *)addr[i])=10.0; printf("Tag %s (%s) not set in parameter file: defaulting to %g (needed for ghost particle headroom in multi-rank runs) \n",tag[i],alternate_tag[i],All.PartAllocFactor); continue;}
#if !(defined(BOX_PERIODIC) || defined(BOX_SHEARING) || defined(BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE) || defined(BOX_LONG_X) || defined(BOX_LONG_Y) || defined(BOX_LONG_Z))
                if(strcmp("BoxSize",tag[i])==0) {*((double *)addr[i])=0; continue;} /* can ignore box size */
#endif
#ifdef SELFGRAVITY_OFF /* can ignore softenings here */
                if(strcmp("SofteningGas",tag[i])==0) {*((double *)addr[i])=0; continue;} 
                if(strcmp("SofteningHalo",tag[i])==0) {*((double *)addr[i])=0; continue;}
                if(strcmp("SofteningDisk",tag[i])==0) {*((double *)addr[i])=0; continue;}
                if(strcmp("SofteningBulge",tag[i])==0) {*((double *)addr[i])=0; continue;}
                if(strcmp("SofteningStars",tag[i])==0) {*((double *)addr[i])=0; continue;}
                if(strcmp("SofteningBndry",tag[i])==0) {*((double *)addr[i])=0; continue;}
#endif
                if(strcmp("CpuTimeBetRestartFile",tag[i])==0) {*((double *)addr[i])=3450.; printf("Tag %s (%s) not set in parameter file: defaulting to write restart checkpoints just under every hour (=%g) \n",tag[i],alternate_tag[i],All.CpuTimeBetRestartFile); continue;}
#if !defined(COOLING) && !defined(GALSF) && !defined(EOS_HELMHOLTZ) && !defined(EOS_ELASTIC) && !defined(EOS_TILLOTSON)
                if(strcmp("UnitLength_in_cm",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitLength_in_cm); continue;}
                if(strcmp("UnitMass_in_g",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitMass_in_g); continue;}
                if(strcmp("UnitVelocity_in_cm_per_s",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitVelocity_in_cm_per_s); continue;}
#ifdef MAGNETIC
                if(strcmp("UnitMagneticField_in_gauss",tag[i])==0) {*((double *)addr[i])=3.5449077018110318; printf("Tag %s (%s) not set in parameter file: will default to assume code units are cgs (=%g), if conversion to physical units for e.g. cooling are needed \n",tag[i],alternate_tag[i],All.UnitMagneticField_in_gauss); continue;}
#endif
#endif
#if defined(SINGLE_STAR_FB_JETS) && defined(SINGLE_STAR_STARFORGE_DEFAULTS)
                if(strcmp("Sink_outflow_jetlaunchvelscaling",tag[i])==0) {*((double *)addr[i])=0.3; printf("Tag %s (%s) not set in parameter file: will default to =%g (this is the fraction of escape velocity outflows are launched with) \n",tag[i],alternate_tag[i],All.Sink_outflow_jetlaunchvelscaling); continue;}
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
#if defined(GALSF_ISMDUSTCHEM_MODEL)
                if(strcmp("Initial_ISMDustChem_Depletion",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to zero (=%g) \n",tag[i],alternate_tag[i],All.Initial_ISMDustChem_Depletion); continue;}
                if(strcmp("Initial_ISMDustChem_SiltoCarbRatio",tag[i])==0) {*((double *)addr[i])=0.; printf("Tag %s (%s) not set in parameter file: defaulting to zero (=%g)\n",tag[i],alternate_tag[i],All.Initial_ISMDustChem_SiliconToCarbonRatio); continue;}
                if(strcmp("SNeIIDustScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_SNeIIDustScaling); continue;}
                if(strcmp("SNeIaDustScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_SNeIaDustScaling); continue;}
                if(strcmp("AGBDustScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_AGBDustScaling); continue;}
                if(strcmp("DustAccretionScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_DustAccretionScaling); continue;}
                if(strcmp("ThermalSputteringScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_ThermalSputteringScaling); continue;}
                if(strcmp("SNeGasClearedOfDustScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_SNeGasClearedOfDustScaling); continue;}
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
                if(strcmp("ISMDustChem_Grain_Size_Min",tag[i])==0) {*((double *)addr[i])=1E-7; printf("Tag %s (%s) not set in parameter file: defaulting to 1E-7 (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_Grain_Size_Min); continue;}
                if(strcmp("ISMDustChem_Grain_Size_Max",tag[i])==0) {*((double *)addr[i])=1E-4; printf("Tag %s (%s) not set in parameter file: defaulting to 1E-4 (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_Grain_Size_Max); continue;}
                if(strcmp("SNeShatteringScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_SNeShatteringScaling); continue;}
                if(strcmp("SNeSputteringScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_SNeSputteringScaling); continue;}
                if(strcmp("GrainShatteringScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_ShatteringScaling); continue;}
                if(strcmp("GrainCoagulationScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_CoagulationScaling); continue;}
                if(strcmp("VCoagScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_VCoagScaling); continue;}
                if(strcmp("CoagDensityEnhancementScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_CoagDensityEnhancementScaling); continue;}
                if(strcmp("GrainVelocityScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_GrainVelocityScaling); continue;}
                if(strcmp("PhotodestructionScaling",tag[i])==0) {*((double *)addr[i])=1.; printf("Tag %s (%s) not set in parameter file: defaulting to one (=%g)\n",tag[i],alternate_tag[i],All.ISMDustChem_PhotodestructionScaling); continue;}
#endif
#endif
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
                if(strcmp("SNeIIEnergyFrac",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to standard stellar-evolution-defaults (=%g) \n",tag[i],alternate_tag[i],All.SNe_Energy_Renormalization); continue;}
                if(strcmp("GasReturnFraction",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to standard stellar-evolution-defaults (=%g) \n",tag[i],alternate_tag[i],All.StellarMassLoss_Rate_Renormalization); continue;}
                if(strcmp("GasReturnEnergy",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to standard stellar-evolution-defaults (=%g) \n",tag[i],alternate_tag[i],All.StellarMassLoss_Energy_Renormalization); continue;}
#endif
#ifdef GALSF_FB_FIRE_RT_LOCALRP
                if(strcmp("WindMomentumLoading",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to standard stellar-evolution-defaults (=%g) \n",tag[i],alternate_tag[i],All.RP_Local_Momentum_Renormalization); continue;}
#endif
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
                if(strcmp("HIIRegion_fLum_Coupled",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to standard stellar-evolution-defaults (=%g) \n",tag[i],alternate_tag[i],All.HIIRegion_fLum_Coupled); continue;}
#endif
#ifdef RT_LEBRON
                if(strcmp("PhotonMomentum_Coupled_Fraction",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to use the explicitly-resolved absorption (=%g) \n",tag[i],alternate_tag[i],All.PhotonMomentum_Coupled_Fraction); continue;}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
                if(strcmp("PhotonMomentum_fUV",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to use the explicitly-resolved absorption (=%g) \n",tag[i],alternate_tag[i],All.PhotonMomentum_fUV); continue;}
                if(strcmp("PhotonMomentum_fOPT",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to use the explicitly-resolved absorption (=%g) \n",tag[i],alternate_tag[i],All.PhotonMomentum_fOPT); continue;}
#endif
#if defined(FIRE_MHD)
                if(strcmp("BiniX",tag[i])==0) {*((double *)addr[i])=1.e-10; printf("Tag %s (%s) not set in parameter file: defaulting =%g (assuming units of Gauss), typical for cosmological FIRE simulations starting around z=100 \n",tag[i],alternate_tag[i],All.BiniX); continue;}
                if(strcmp("BiniY",tag[i])==0) {*((double *)addr[i])=1.e-10; printf("Tag %s (%s) not set in parameter file: defaulting =%g (assuming units of Gauss), typical for cosmological FIRE simulations starting around z=100 \n",tag[i],alternate_tag[i],All.BiniY); continue;}
                if(strcmp("BiniZ",tag[i])==0) {*((double *)addr[i])=1.e-10; printf("Tag %s (%s) not set in parameter file: defaulting =%g (assuming units of Gauss), typical for cosmological FIRE simulations starting around z=100 \n",tag[i],alternate_tag[i],All.BiniZ); continue;}
#endif
#if defined(FIRE_CRS)
#if defined(COSMIC_RAY_SUBGRID_LEBRON)
                if(strcmp("CosmicRay_Subgrid_Kappa_0",tag[i])==0) {*((double *)addr[i])=120; printf("Tag %s (%s) not set in parameter file: defaulting to calibrated value from FIRE simulations, assuming units kpc/h and km/s (=%g) \n",tag[i],alternate_tag[i],All.CosmicRay_Subgrid_Kappa_0); continue;}
                if(strcmp("CosmicRay_Subgrid_Vstream_0",tag[i])==0) {*((double *)addr[i])=20; printf("Tag %s (%s) not set in parameter file: defaulting to calibrated value from FIRE simulations, assuming units kpc/h and km/s (=%g) \n",tag[i],alternate_tag[i],All.CosmicRay_Subgrid_Vstream_0); continue;}
#endif
#if defined(CRFLUID_DIFFUSION_MODEL) && (CRFLUID_DIFFUSION_MODEL == 0)
                if(strcmp("CosmicRayDiffusionCoeff",tag[i])==0) {*((double *)addr[i])=690.; printf("Tag %s (%s) not set in parameter file: defaulting to observationally-favored diffusivity ~3e29, assuming units kpc/h and km/s (=%g) \n",tag[i],alternate_tag[i],All.CosmicRayDiffusionCoeff); continue;}
#endif
                if(strcmp("CosmicRay_SNeFraction",tag[i])==0) {*((double *)addr[i])=0.1; printf("Tag %s (%s) not set in parameter file: defaulting to observationally-favored ~10 percent conversion to CRs (=%g) \n",tag[i],alternate_tag[i],All.CosmicRay_SNeFraction); continue;}
#endif
#if defined(SINGLE_STAR_STARFORGE_DEFAULTS)
                if(strcmp("SfEffPerFreeFall",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to scaling-default of unity (=%g) \n",tag[i],alternate_tag[i],All.MaxSfrTimescale); continue;}
                if(strcmp("SinkAccretionFactor",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to Hopkins and Quataert best-estimate (=%g) \n",tag[i],alternate_tag[i],All.SinkAccretionFactor); continue;}
                if(strcmp("SinkEddingtonFactor",tag[i])==0) {*((double *)addr[i])=1e10; printf("Tag %s (%s) not set in parameter file: defaulting to no Eddington-limit in accretion from disk to sink (=%g) \n",tag[i],alternate_tag[i],All.SinkEddingtonFactor); continue;}
                if(strcmp("SeedSinkMass",tag[i])==0) {*((double *)addr[i])=1e-20; printf("Tag %s (%s) not set in parameter file: defaulting to arbitrary small value which will be ignored later (=%g) \n",tag[i],alternate_tag[i],All.SeedSinkMass); continue;}
                if(strcmp("SeedReservoirMass",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to sinks beginning their existence without an active accretion disk (=%g) \n",tag[i],alternate_tag[i],All.SeedReservoirMass); continue;}
                if(strcmp("SinkNgbFactor",tag[i])==0) {*((double *)addr[i])=1.0; printf("Tag %s (%s) not set in parameter file: defaulting to no augment of sink neighbors vs gas (=%g) \n",tag[i],alternate_tag[i],All.SinkNgbFactor); continue;}
                if(strcmp("SinkMaxAccretionRadius",tag[i])==0) {*((double *)addr[i])=5.0; printf("Tag %s (%s) not set in parameter file: defaulting to some large size of order a few in code units (=%g) \n",tag[i],alternate_tag[i],All.SinkMaxAccretionRadius); continue;}
                if(strcmp("SinkFeedbackFactor",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to follow user-defined coefficients for each mechanism (=%g) \n",tag[i],alternate_tag[i],All.SinkFeedbackFactor); continue;}
                if(strcmp("SinkRadiativeEfficiency",tag[i])==0) {*((double *)addr[i])=5.e-7; printf("Tag %s (%s) not set in parameter file: defaulting to a reference radiative efficiency, but tracks ignore this (=%g) \n",tag[i],alternate_tag[i],All.SinkRadiativeEfficiency); continue;}
#if defined(SINK_WIND_SPAWN)
                if(strcmp("Sink_accreted_fraction",tag[i])==0) {*((double *)addr[i])=0.7; printf("Tag %s (%s) not set in parameter file: defaulting to assume 70 percent is accreted onto sink versus outflow (=%g) \n",tag[i],alternate_tag[i],All.Sink_accreted_fraction); continue;}
                if(strcmp("Sink_outflow_velocity",tag[i])==0) {*((double *)addr[i])=100.; printf("Tag %s (%s) not set in parameter file: defaulting to assume mechanical outflow with 100 in code units, but tracks ignore this (=%g) \n",tag[i],alternate_tag[i],All.Sink_outflow_velocity); continue;}
                if(strcmp("Sink_outflow_temperature",tag[i])==0) {*((double *)addr[i])=1.e3; printf("Tag %s (%s) not set in parameter file: defaulting to assuming ISM-type temperatures in internal spawned elements (=%g) \n",tag[i],alternate_tag[i],All.Sink_outflow_temperature); continue;}
#ifdef SINGLE_STAR_FB_WINDS
                if(strcmp("Cell_Spawn_Mass_ratio_MS",tag[i])==0) {*((double *)addr[i])=0.0; printf("Tag %s (%s) not set in parameter file, Sink_outflow_particlemass will be used instead \n",tag[i],alternate_tag[i]); continue;}
#endif
#endif
#endif
#if defined(RT_ISRF_BACKGROUND)
                if(strcmp("InterstellarRadiationFieldStrength",tag[i])==0) {*((double *)addr[i])=1.0; printf("Tag %s (%s) not set in parameter file: defaulting to assuming Solar neighborhood (Draine) background radiation field (=%g) \n",tag[i],alternate_tag[i],All.InterstellarRadiationFieldStrength); continue;}

                if(strcmp("Redshift_RT_Background",tag[i])==0) {*((double *)addr[i])=0.0; printf("Tag %s (%s) not set in parameter file: defaulting to assuming z=0 background radiation field \n",tag[i],alternate_tag[i]); continue;}
#endif
#if defined(RT_INFRARED)
                if(strcmp("InitRadiationTemp",tag[i])==0) {*((double *)addr[i])=20.; printf("Tag %s (%s) not set in parameter file: defaulting to assuming 20K radiation temperature \n",tag[i],alternate_tag[i]); continue;}
#endif
#if defined(FIRE_PHYSICS_DEFAULTS)
                if(strcmp("SfEffPerFreeFall",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to FIRE-default of unity (=%g) \n",tag[i],alternate_tag[i],All.MaxSfrTimescale); continue;}
#if (FIRE_PHYSICS_DEFAULTS >= 3)
                if(strcmp("CritPhysDensity",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: using FIRE-3 default of 1/cm3 (=%g) \n",tag[i],alternate_tag[i],All.CritPhysDensity); continue;}
#else
                if(strcmp("CritPhysDensity",tag[i])==0) {*((double *)addr[i])=100; printf("Tag %s (%s) not set in parameter file: using FIRE-1/2 default of 100/cm3 (=%g) \n",tag[i],alternate_tag[i],All.CritPhysDensity); continue;}
#endif
#endif
#if defined(FIRE_BHS)
                if(strcmp("SinkAccretionFactor",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to Hopkins and Quataert best-estimate (=%g) \n",tag[i],alternate_tag[i],All.SinkAccretionFactor); continue;}
                if(strcmp("SinkEddingtonFactor",tag[i])==0) {*((double *)addr[i])=100; printf("Tag %s (%s) not set in parameter file: defaulting no Eddington-limit in accretion from disk to sink (=%g) \n",tag[i],alternate_tag[i],All.SinkEddingtonFactor); continue;}
                if(strcmp("SeedSinkMass",tag[i])==0) {*((double *)addr[i])=0.7e-8; printf("Tag %s (%s) not set in parameter file: defaulting to upper-limit of normal stellar BHs, assuming code mass units of 1e10 Msun/h (=%g) \n",tag[i],alternate_tag[i],All.SeedSinkMass); continue;}
                if(strcmp("SinkNgbFactor",tag[i])==0) {*((double *)addr[i])=8.0; printf("Tag %s (%s) not set in parameter file: defaulting to standard augment of sink neighbors vs gas (=%g) \n",tag[i],alternate_tag[i],All.SinkNgbFactor); continue;}
                if(strcmp("SinkMaxAccretionRadius",tag[i])==0) {*((double *)addr[i])=5.0; printf("Tag %s (%s) not set in parameter file: defaulting to typical galaxy size assuming code units of kpc/h (=%g) \n",tag[i],alternate_tag[i],All.SinkMaxAccretionRadius); continue;}
                if(strcmp("SinkRadiativeEfficiency",tag[i])==0) {*((double *)addr[i])=0.1; printf("Tag %s (%s) not set in parameter file: defaulting to canonical radiative efficiency (=%g) \n",tag[i],alternate_tag[i],All.SinkRadiativeEfficiency); continue;}
                if(strcmp("SinkFeedbackFactor",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to follow user-defined coefficients for each mechanism (=%g) \n",tag[i],alternate_tag[i],All.SinkFeedbackFactor); continue;}
                if(strcmp("SeedSinkMassSigma",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to a uniform sink seed mass (=%g) \n",tag[i],alternate_tag[i],All.SeedSinkMassSigma); continue;}
                if(strcmp("SeedSinkMinRedshift",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to allow seed creation at all times (=%g) \n",tag[i],alternate_tag[i],All.SeedSinkMinRedshift); continue;}
                if(strcmp("SeedReservoirMass",tag[i])==0) {*((double *)addr[i])=0; printf("Tag %s (%s) not set in parameter file: defaulting to sinks beginning their existence without an active accretion disk (=%g) \n",tag[i],alternate_tag[i],All.SeedReservoirMass); continue;}
                if(strcmp("Sink_FluxMomentumFactor",tag[i])==0) {*((double *)addr[i])=1; printf("Tag %s (%s) not set in parameter file: defaulting to using actual AGN spectrum for radiative feedback (=%g) \n",tag[i],alternate_tag[i],All.Sink_Rad_MomentumFactor); continue;}
                if(strcmp("Sink_accreted_fraction",tag[i])==0) {*((double *)addr[i])=0.5; printf("Tag %s (%s) not set in parameter file: defaulting to assume equal sink accretion and outflow rate intrinsically (=%g) \n",tag[i],alternate_tag[i],All.Sink_accreted_fraction); continue;}
                if(strcmp("Sink_outflow_velocity",tag[i])==0) {*((double *)addr[i])=1.e4; printf("Tag %s (%s) not set in parameter file: defaulting to assume mechanical outflow with 1e4 km/s assuming km/s code units (=%g) \n",tag[i],alternate_tag[i],All.Sink_outflow_velocity); continue;}
#if defined(SINK_WIND_SPAWN)
                if(strcmp("Sink_outflow_temperature",tag[i])==0) {*((double *)addr[i])=1.e4; printf("Tag %s (%s) not set in parameter file: defaulting to assuming ISM-type temperatures in internal spawned elements (=%g) \n",tag[i],alternate_tag[i],All.Sink_outflow_temperature); continue;}
#if defined(SINK_SCALE_SPAWNINGMASS_WITH_INITIALMASS)
                if(strcmp("Sink_outflow_particlemass",tag[i])==0) {*((double *)addr[i])=0.01; printf("Tag %s (%s) not set in parameter file: defaulting to assuming spawned cells from the sink have mass a fixed fraction (=%g) of the sink initial formation mass \n",tag[i],alternate_tag[i],All.Sink_outflow_particlemass); continue;}
#endif
#endif
#if defined(SINK_COSMIC_RAYS)
                if(strcmp("Sink_CosmicRay_Injection_Efficiency",tag[i])==0) {*((double *)addr[i])=1.e-2; printf("Tag %s (%s) not set in parameter file: defaulting to assuming CR injection efficiency of ~1 percent (=%g) \n",tag[i],alternate_tag[i],All.Sink_CosmicRay_Injection_Efficiency); continue;}
#endif
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
#ifdef EOS_ANEOS
                if(strncmp(tag[i], "AneosTable", 10)==0) {strcpy((char *)addr[i], "none"); continue;} /* unused ANEOS table slots default to 'none' */
#endif
#ifdef NUCLEAR_NETWORK
                if(strcmp("NuclearNetworkDataFile",tag[i])==0) {strcpy((char *)addr[i], ""); continue;} /* empty = built-in aprox13, no external data needed */
                if(strcmp("NuclearBurningFloor_T",tag[i])==0) {*((double *)addr[i])=1.0e8; continue;} /* default: 10^8 K */
                if(strcmp("NuclearBurningFloor_rho",tag[i])==0) {*((double *)addr[i])=0; continue;} /* default: no density floor */
                if(strcmp("NuclearNSE_T_threshold",tag[i])==0) {*((double *)addr[i])=6.0e9; continue;} /* default: 6 GK */
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
#endif
    
#ifdef DIVBCLEANING_DEDNER
#ifdef MHD_CONSTRAINED_GRADIENT
    All.DivBcleanParabolicSigma = 1.0;
#else
    All.DivBcleanParabolicSigma = 0.2;
#endif
    All.DivBcleanHyperbolicSigma = 1.0;
#endif

#ifdef TURB_DIFF_DYNAMIC
    All.TurbDynamicDiffIterations = 0; /* D. Rennehan: This has NOT been tested above 0 */
#endif
#if !defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
    if(All.ComovingIntegrationOn) {All.ErrTolForceAcc = 0.005; All.ErrTolIntAccuracy = 0.05;}
#endif
    All.MaxNumNgbDeviation = All.DesNumNgb / 640.;
#ifdef GALSF
    All.MaxNumNgbDeviation = All.DesNumNgb / 64.;
#endif
#if NUMDIMS==3    
    if(All.MaxNumNgbDeviation < 0.05) All.MaxNumNgbDeviation = 0.05;
#endif    
#ifdef EOS_ELASTIC
    All.MaxNumNgbDeviation /= 20.0;
#endif
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    All.AGS_MaxNumNgbDeviation = All.AGS_DesNumNgb / 640.;
#ifdef GALSF
    All.AGS_MaxNumNgbDeviation = All.AGS_DesNumNgb / 64.;
#endif
    if(All.AGS_MaxNumNgbDeviation < 0.05) All.AGS_MaxNumNgbDeviation = 0.05;
#endif
#endif // closes DEVELOPER_MODE check //
#ifdef SINK_WIND_SPAWN
    All.SpawnedWindCellID = 1913298393;       // this seems weird, but is the bitshifted version of 1234568912345 for not long IDs.
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

#if defined(GRAVITY_ACCURATE_FEWBODY_INTEGRATION)
    if(All.ErrTolIntAccuracy > 0.01) {All.ErrTolIntAccuracy = 0.01;}
    if(All.MaxRMSDisplacementFac > 0.125) {All.MaxRMSDisplacementFac = 0.125;}
    if(All.ErrTolTheta > 0.5) {All.ErrTolTheta = 0.5;}
    if(All.MaxNumNgbDeviation > 0.05) {All.MaxNumNgbDeviation = 0.05;}
#endif
    
#if defined(MAGNETIC) || defined(HYDRO_MESHLESS_FINITE_VOLUME) || defined(SINK_WIND_SPAWN)
    if(All.CourantFac > 0.2) {All.CourantFac = 0.2;}
    /* (PFH) safety factor needed for MHD calc, because people keep using the same CFac as hydro! */
#endif

#if defined(PIC_MHD) && !defined(GRAIN_FLUID_AND_PIC_BOTH_DEFINED)
    All.Grain_Internal_Density=1; All.Grain_Size_Min=1; All.Grain_Size_Max=1; All.Grain_Size_Spectrum_Powerlaw=1; /* in this case these are never used, so we treat them as dummy variables */
#endif

    /* now we're going to do a bunch of checks */
    if(All.MaxMemSize > safe_memorypertask)
    {
        if(ThisTask==0) {printf("WARNING: MaxMemSize (Max_Memory_Per_MPI_Task_in_MB=%d) is currently set to a larger value than the maximum safe amount of memory recommended per task, given by pinging the system for allocatable memory and dividing it among the tasks per node (=%g MB). Depending on the details of your node and core configuration and memory use, this may work, but it is not safe, and can crash if too many processes on a node try to use their full memory at once.\n",All.MaxMemSize,safe_memorypertask); fflush(stdout);}
    }
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
    if(All.LETAllocFactor <= 0)
    {
        if(ThisTask==0) {printf("LETAllocFactor must be >0 in GPU builds: the legacy gravity export fallback is retired, so LET cannot be disabled. Use the default 1.0 or increase this value if LET unpack overflows.\n"); endrun(1);}
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
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
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
