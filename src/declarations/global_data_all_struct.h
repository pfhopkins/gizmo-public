/*! \file global_data_all_struct.h
 *  \brief Definition of struct global_data_all_processes ("All").
 *
 *  Intentionally includes ONLY the headers needed to compile the struct body
 *  (GIZMO_config, typedefs, math types) — no MPI, GSL, or HDF5 headers.
 *  This lets GPU-only TUs (allvars_gpu.cu, future CUDA files) include just
 *  this header to get the struct type without pulling in host-only libraries.
 *
 *  allvars.h includes this header and adds all the other global declarations.
 */

#ifndef GLOBAL_DATA_ALL_STRUCT_H
#define GLOBAL_DATA_ALL_STRUCT_H

/* GIZMO_GPU_COMPILER: defined when compiled by any GPU device compiler.
 * Placed here (earliest universal header) so it is available before allvars.h
 * and macros.h, which matters for the All_dev guards in cooling.cc/eos.cc. */
#if defined(__CUDACC__) || defined(__HIPCC__)
#ifndef GIZMO_GPU_COMPILER
#define GIZMO_GPU_COMPILER
#endif
#endif

#include "../GIZMO_config.h"
#include "precompiler_logic.h"
#include "constants.h"
#include "macros.h"
#include "typedefs.h"
#include "../math_types/vec3.h"

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

#ifdef SINK_PARTICLES
  int TotSinks;
#endif

#if defined(DM_SIDM)
    MyDouble DM_InteractionCrossSection;  /*!< self-interaction cross-section in [cm^2/g]*/
    MyDouble DM_DissipationFactor;  /*!< dimensionless parameter governing efficiency of dissipation (1=dissipative, 0=elastic) */
    MyDouble DM_KickPerCollision;  /*!< for exo-thermic DM reactions, this determines the energy gain 'per event': kick in code units (equivalent to specific energy) associated 'per event' */
    MyDouble DM_InteractionVelocityScale; /*!< scale above which the scattering becomes velocity-dependent */
#endif

#ifdef DM_HEATING
    /* Continuous gas heating from DM annihilation and/or decay. Annihilation
     * uses self-conjugate-DM convention with no factor of 1/2 (see
     * Template_Config.sh and gizmo_documentation.md). Dirac/non-self-conjugate
     * users fold the appropriate factor into DM_AnnihilationSigmaV_over_mChi. */
    double DM_AnnihilationSigmaV_over_mChi;  /*!< <sigma v> / m_chi; input cm^3 s^-1 g^-1, stored in code units after set_units() */
    double DM_AnnihilationHeatingFraction;   /*!< f_h^ann in [0,1], fraction of rest-mass energy thermalized locally per annihilation */
    double DM_DecayRate;                     /*!< Gamma; input s^-1, stored in code units after set_units() */
    double DM_DecayHeatingFraction;          /*!< f_h^dec in [0,1], fraction of m_chi*c^2 thermalized locally per decay */
#endif

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartGas;		/*!< This gives the maxmimum number of gas cells that can be stored on one processor. */
  int MaxPartAssignable;	/*!< Ceiling on the number of LOCAL particles the domain decomposition may assign to a rank.  IDENTICAL ON EVERY RANK, and that is its contract: the decomposition compares an already-reduced load against it and gates a COLLECTIVE re-split on the result, so a rank-dependent value makes ranks disagree about whether to re-split.  Distinct from MaxPart, which is this rank's PHYSICAL storage and additionally covers imported ghosts, may differ between ranks, and may grow while the run is going.  Storage is at least this large whenever the decomposition runs, which is what makes an assignment that fits this cap always fit the storage; resize_particle_storage refuses to take storage below it.  (The one window where MaxPart is deliberately smaller than the storage behind it is the restart payload read: storage is allocated at the raised capacity while this and MaxPart both still report the writer's, so the serialized tree reads back into the layout it was written with.  That window lasts for the read and no longer -- restart() puts both in force at the end of it, before any step runs, and code that grows storage or refuses a spawn mid-step must never see the smaller value.)  Set from PartAllocFactor and the balanced particle count when the ICs are read, and worked out again from both at every restart, since the factor is defined relative to a count that moves as particles are created and destroyed; uniform by construction in both places, needing no communication, and carried through restart files with the rest of this struct. */
  int MaxPartExpandable;	/*!< Hard ceiling on MaxPart for the whole run: the largest particle capacity this rank may ever be raised to.  Fixed once when the ICs are read and carried through restart files; nothing recomputes it, so a run's index layout never moves.  It is deliberately generous -- it spends a whole per-rank share of node memory on particle slots, ignoring the tree, transport buffers and the arena, which also need room -- because its job is to bound an index space cheaply, not to predict what fits; whether the memory is really there stays the allocator's answer.  It is held down by two other terms: a wide multiple of the starting capacity, so a small problem on a large machine does not size for slots it will never hold, and a cap that keeps the tree's per-particle side arrays within a small stated share of a node.  Those side arrays are sized from this rather than from the current capacity, which is what makes a particle index that outruns them impossible to construct; they cost a couple of dozen bytes per slot against a per-slot particle record that is one to two orders of magnitude larger, so the exchange is cheap wherever the particle record is fat and is bounded by the share term where it is not.  Identical on every rank, which is what lets TreeNodeIndexBase be a shared constant. */
  int TreeNodeIndexBase;	/*!< First index of the gravity tree's node range.  The index space has THREE regions, not two: particle indices occupy [0,TreeParticleSlots), the span [TreeParticleSlots,TreeNodeIndexBase) holds NO valid object and a walk that lands there is reading a malformed tree, and tree nodes / foreign nodes / pseudo-particles all live at or above TreeNodeIndexBase.  Its ONLY contract is "larger than any local particle index any rank ever uses, identical on every rank, constant for the run"; it is a namespace separator, NOT a memory size.  Set once from MaxPart when the ICs are read, carried through restart files with the rest of this struct.  NOTHING may be sized from it -- array sizes follow real counts (see MaxNodes, MaxForeignNodes, TreeParticleSlots). */
  int TreeParticleSlots;	/*!< Particle-slot count the LIVE tree's per-particle side arrays were allocated with: Father[] spans [0,TreeParticleSlots), Nextnode[]'s particle segment likewise, and Nextnode[]'s pseudo-particle segment starts at that offset.  A COUNT of real storage, distinct from TreeNodeIndexBase (where node INDICES start).  Published by force_treeallocate from its tree_particle_slots argument BEFORE any tree payload is read or walked -- that publish-before-use is the correctness contract.  The copy carried through restart files inside this struct IS relied upon, in exactly one place: a restart read allocates the reader's tree with it, because the file's tree side arrays and its pseudo-particle offset describe the WRITER's slot count, which the capacity may since have outgrown.  It therefore has the same status as MaxNodes and MaxForeignNodes -- a serialized description of a physical layout that must be reproduced, not re-derived. */
  int ICFormat;			/*!< selects different versions of IC file-format */
  int SnapFormat;		/*!< selects different versions of snapshot file-formats */
  int NumFilesPerSnapshot;	/*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/*!< maximum number of files that may be written simultaneously when writing/reading restart-files, or when writing snapshot files */
  double CommChunkSize;		/*!< Size in MB of the chunk particle data is sent between tasks in: the snapshot and initial-conditions field blocks, the generic chunked all-to-all, and the neighbour loops' round batches.  A round-count/memory trade, not a capacity: a smaller chunk means more rounds, never a failure.  Worked out by the code from what a round has to hold, unless the parameter file asks for a size. */
  long BunchSize;     	        /*!< how many targets the gravity walk can record as not covered by the locally-built tree before it stops the run */

  double PartAllocFactor;	/*!< The largest local particle imbalance the decomposition may produce, as a multiple of the balanced count: MaxPartAssignable = PartAllocFactor * ceil(TotNumPart/NTask), the balanced count rounded up so a rank can hold the larger side of an uneven split.  It is a ceiling on ASSIGNMENT, not a storage size.  Room for imported ghosts is no longer taken from it: ghosts share the P[]/CellP[] arrays, but the amount is measured per rank from the imports that actually happen and added on top, and an import that still does not fit raises the capacity when it arrives.  So this only needs to cover the imbalance the decomposition really reaches, which is why the default is small.  Raising it still inflates MaxNodes / MaxForeignNodes and the domain caps, not ghost room. */
  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times the maximum(!) number of particles.  Note: A typical local tree for N particles needs usually about ~0.65*N nodes. */
  int DomainSegmentsPerRank;	/*!< Number of separate contiguous Peano-Hilbert segments handed to each rank by the
				  domain decomposition.  Settled once when the initial conditions are read, from the
				  particle count and the number of ranks (domain_segments_per_rank_for_particles), and
				  constant for the run thereafter.  Unrelated to TopNodeAllocFactor below, which is
				  only a capacity for the top-tree arrays. */
  double TopNodeAllocFactor;	/*!< Sizes the TOP-level (domain-decomposition) tree, not the full local tree: MaxTopNodes = TopNodeAllocFactor * MaxPart + 1.  Much smaller than TreeAllocFactor (default ~0.008 vs ~0.45) because the top-tree holds one node per top-level domain cell, not one per particle; auto-ratchets on top-node overflow. */

#ifdef DM_SCALARFIELD_SCREENING
  double ScalarBeta;
  double ScalarScreeningLength;
#endif

  /* some gas/fluid arbitrary-mesh parameters */
  double DesNumNgb;		/*!< Desired number of gas cell neighbours */
#ifdef SUBFIND
  int DesLinkNgb;       /*! < Number of neighbors used for linking and density estimation in SUBFIND */
#endif

  double MaxNumNgbDeviation;	/*!< Maximum allowed deviation neighbour number */
  double ArtBulkViscConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
#ifdef HYDRO_MULTIFLUID_DM_COOLING
  double ADM_FineStructure;     /*!< alpha_ADM (SM value 7.2973525693e-3) */
  double ADM_ProtonMass;        /*!< m_p_ADM in g (SM value PROTONMASS_CGS) */
  double ADM_ElectronMass;      /*!< m_e_ADM in g (SM value ELECTRONMASS) */
  double ADM_MolecularFraction; /*!< dark-H2 fraction by mass, 0..1 */
#endif
#ifdef DISK_BETA_COOL
  double BetaCool_Beta;         /*!< beta in t_cool = beta / Omega */
  double BetaCool_Tirr;         /*!< irradiation-floor temperature [K]; 0 disables floor */
  double BetaCool_u_irr;        /*!< derived: u corresponding to BetaCool_Tirr */
#endif
#ifdef PLANET_HEATING
  double PlanetHeating_RadQ0_cgs; /*!< input: initial radiogenic specific heating rate [erg/g/s] */
  double PlanetHeating_RadTau_cgs;/*!< input: radiogenic e-folding decay time [s] */
  double PlanetHeating_AccQ0_cgs; /*!< input: background accretional heating rate [erg/g/s]; 0=off */
  double PlanetHeating_RadQ0;     /*!< derived: RadQ0_cgs in code units [code u / code t] */
  double PlanetHeating_RadTau;    /*!< derived: RadTau_cgs in code units [code t] */
  double PlanetHeating_AccQ0;     /*!< derived: AccQ0_cgs in code units [code u / code t] */
#endif
#ifdef GRAIN_EVOLUTION
  double GrainEvolution_StickingCoeff;          /*!< global sticking-coefficient multiplier for pairwise outcomes (bits 0|1|2) and condensation (bit 5). 1.0 = use species defaults from grain_collisional_outcomes.h. */
  double GrainEvolution_VelThreshFrag;          /*!< |dv| threshold for fragmentation onset [code velocity]; 0 = use species defaults. */
  double GrainEvolution_VelThreshShat;          /*!< |dv| threshold for shattering onset [code velocity]; 0 = use species defaults. */
  double GrainEvolution_ThermalSputteringScaling; /*!< global multiplier on the Nozawa+(2006) thermal-sputter erosion rate (bit 3). Mirrors All.ISMDustChem_ThermalSputteringScaling. 1.0 = nominal. */
#endif
#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
  double GrainPromotion_MassThresh;         /*!< derived: MassThresh_cgs in code units */
  double GrainPromotion_MassThresh_cgs;     /*!< user input: grain mass promotion threshold [g] */
  double GrainPromotion_DustGasRatioThresh; /*!< grain/gas density ratio threshold; 0 = disabled */
#endif
#ifdef CHIMES
  int ChimesThermEvolOn;        /*!< Flag to determine whether to evolve the temperature in CHIMES. */
#ifdef CHIMES_STELLAR_FLUXES
  double Chimes_f_esc_ion;
  double Chimes_f_esc_G0;
#endif
#endif // CHIMES

  double MinEgySpec;		/*!< the minimum allowed temperature expressed as energy per unit mass */
#ifdef SPHAV_ARTIFICIAL_CONDUCTIVITY
  double ArtCondConstant;
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
    double MeanGasParticleMass; /*!< the mean gas particle mass */
#endif
    double MinMassForParticleMerger; /*!< the minimum mass of a gas particle below which it will be merged into a neighbor */
    double MaxMassForParticleSplit; /*!< the maximum mass of a gas particle above which it will be split into a pair */

  /* some force counters  */
  long long TotNumOfForces;	/*!< counts total number of force computations  */
  long long NumForcesSinceLastTreeBuild;	/*!< count particle updates since the last gravity tree build */

  /* various cosmological factors that are only a function of the current scale factor, and in Newtonian runs are set to 1 */
  double cf_atime, cf_a2inv, cf_a3inv, cf_hubble_a, cf_hubble_a2;

  /* system of units  */
  double UnitMass_in_g,		        /*!< factor to convert internal mass unit to grams/h */
         UnitVelocity_in_cm_per_s,	/*!< factor to convert intqernal velocity unit to cm/sec */
         UnitLength_in_cm,          /*!< factor to convert internal length unit to cm/h */
         G;                         /*!< Gravity-constant in internal units */

#ifdef MAGNETIC
  double UnitMagneticField_in_gauss; /*!< factor to convert internal magnetic field (B) unit to gauss (cgs) units */
#endif

  /* Cosmology */
  double Hubble_H0_CodeUnits;		/*!< Hubble-constant (unit-ed version: 100 km/s/Mpc) in internal units */
  double OmegaMatter,		/*!< matter density in units of the critical density (at z=0) */
    OmegaLambda,		/*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,		/*!< baryon density in units of the critical density (at z=0) */
    OmegaRadiation,     /*!< radiation [including all relativistic components] density in units of the critical density (at z=0) */
    HubbleParam;		/*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute physical values for cooling physics */

  double BoxSize;		/*!< Boxsize in case periodic boundary conditions are used */
#ifdef BOX_SHEARING
  MyDouble Shearing_Box_Vel_Offset; /*!< shearing box velocity offset (in All for GPU access) */
  MyDouble Shearing_Box_Pos_Offset; /*!< shearing box position offset (in All for GPU access) */
#endif
#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z)
  short int special_boundary_condition_xyz_def_reflect[3]; /*!< per-axis reflecting-boundary code (in All for GPU access) */
  short int special_boundary_condition_xyz_def_outflow[3]; /*!< per-axis outflow-boundary code (in All for GPU access) */
#endif

  /* Code options */
  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int ResubmitOn;		/*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;	/*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion */
  int OutputListOn;		/*!< flags that output times are listed in a specified file */

  int HighestActiveTimeBin;
  int HighestOccupiedTimeBin;

  /* parameters determining output frequency */
  int SnapshotFileCount;	/*!< number of snapshot that is written next */
  double TimeBetSnapshot,	/*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,	/*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,	/*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,	/*!< cpu-time when last restart-file was written */
    TimeBetStatistics,		/*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;		/*!< simulation time when the energy statistics was computed the last time */
  integertime NumCurrentTiStep;		/*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */
  double Time,			/*!< current time of the simulation */
    TimeBegin,			/*!< time of initial conditions of the simulation */
    TimeStep,			/*!< difference between current times of previous and current timestep */
    TimeMax;			/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */
  double Timebase_interval;	/*!< factor to convert from floating point time interval to integer timeline */
  integertime Ti_Current;		/*!< current time on integer timeline */
  integertime Previous_Ti_Current;
  integertime Ti_nextoutput;		/*!< next output time on integer timeline */
  integertime Ti_lastoutput;

#ifdef PMGRID
  integertime PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif

  int    CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];
  int LevelToTimeBin[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */
  double TimeLimitCPU;
  double CPU_Sum[CPU_PARTS];    /*!< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */
  double ErrTolTheta;		/*!< Barnes-Hut tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */
  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */
  double MinSizeTimestep,	/*!< minimum allowed timestep. Normally, the simulation terminates if the timestep determined by the timestep criteria falls below this limit. */
         MaxSizeTimestep;		/*!< maximum allowed timestep */
  double MaxRMSDisplacementFac;	/*!< this determines a global timestep criterion for cosmological simulations in comoving coordinates.  To this end, the code computes the rms velocity
				   of all particles, and limits the timestep such that the rms displacement is a fraction of the mean particle separation (determined from the particle mass and the cosmological parameters). This parameter specifies this fraction. */
  int WorkingMemoryPoolSize;		/*!< Per-MPI-task size in MB of the working memory pool, the scoped last-in-first-out arena mymalloc hands blocks out of.  NOT a total memory budget: the particle arrays, the tree, device scratch and the transport buffers are allocated outside it and accounted separately.  Worked out by the code from the tenants that actually draw on it, unless the parameter file asks for a size, in which case that size is used untouched. */
  int NeighborLoopModeBThresholdSum;	/*!< optional Mode-A/B dispatch threshold on the summed active-neighbor count; -1 = unset (use each loop's Spec::modeb_threshold_sum) */
  int NeighborLoopModeBThresholdMax;	/*!< optional Mode-A/B dispatch threshold on the max-rank active-neighbor count; -1 = unset (use each loop's Spec::modeb_threshold_max) */
  int GravityHostWalkBelowActive;	/*!< below this many RANK-LOCAL active gravity candidates the gravity tree walk runs on the host instead of the device. 0 disables the count-based routing. Default 1000. */
  double TreeUpdateFullSweep_ActiveFraction;	/*!< when the device runs the dynamic tree update, it drifts every node in the tree first (which also certifies the whole tree current for the walks that follow)
					   only if the rank-local fraction of elements being updated is at least this; below it, it brings current just the chains of nodes the kick touches. 0 = always sweep every node; above 1 = chains only, with the sweep kept as the answer to a chain list that overflows. Default 2. */
  int TreeUpdateHostBelowActive;	/*!< below this many RANK-LOCAL elements due a force update at the sync point, the dynamic tree update runs on the host, which drifts nodes lazily as it
					   touches them instead of sweeping every node up front. 0 disables the count-based routing; the host still owns the rest of any time step in which a node has already been drifted lazily, since the device sweep
					   cannot refresh those nodes' mirror. Default 1000. */
  double CourantFac;		/*!< Courant factor */
#ifdef CBE_INTEGRATOR
  double CBEMassEffFloor;	/*!< CBE timestep m_eff floor fraction: m_eff = max(m_b, CBEMassEffFloor*m_cell) in the per-basis mass-depletion + moment-accel timestep criteria, so near-empty placeholder/free-slot bases cannot force an absurdly small step. Timestep-only (does not touch flux/update). Default 0.1. */
#endif
#ifdef CBE_INTEGRATOR_COLLISIONS
  double CBECollisionCrossSection;	/*!< sigma/mu_p (cross-section per unit physical particle mass) in code units, for the intra-cell CBE collision operator + collisional Riemann term. 0 disables (no-op). */
#endif

  /* how much of the system has to be active before the tree is rebuilt, and before the particles
     are redistributed between tasks. A decomposition always brings a new tree with it, so setting
     the tree fraction above the domain fraction has no effect; setting it below rebuilds the tree
     on steps that keep the existing decomposition. */
  double TreeRebuild_ActiveFraction;	/*!< rebuild the gravity tree once this fraction of the system is active */
  double DomainBuild_ActiveFraction;	/*!< redo the domain decomposition once this fraction of the system is active */
#ifdef MHD_MODIFIED_GRADIENT
  double ActiveFractionForMGSweep;  /*!< minimum active gas fraction to trigger the global MG div(B) solve; on smaller timesteps the local CG correction is used instead */
  int Flag_SkipMGSolve;             /*!< per-timestep flag: 1 = skip MG global solve this step (use CG fallback), 0 = run MG */
#endif
#ifdef TWO_TEMPERATURE_PLASMA
  double TwoTemp_InitialTeOverTgas; /*!< initial T_e / T_gas ratio used at the LTE seed in eos.cc on the first call (when u_e_cell == 0). Default 1.0 (LTE start). Set to !=1 in the param file to start in a 2-T initial state, e.g. for the 2T_relaxation regression test (electrons cold, ions hot, watch the analytic Spitzer relaxation toward T_eq). */
  double TwoTemp_ShockElectronFraction; /*!< fraction f_e of the hydro-dissipation (DtInternalEnergy) deposited into electrons; remainder goes to ions. Default 0.0 = collisionless-shock limit (Vink+15, Ghavamian+13: ions take all of the dissipation, electrons heat only via Coulomb equilibration). f_e = 1.0 = strong-coupling limit (single-fluid behavior). Mach-dependent f_e (Ghavamian-style) deferred to a later refinement; one scalar suffices for v1. */
#endif

  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length) five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars */
    double MinGasKernelRadiusFractional; /*!< minimim allowed gas kernel length relative to force softening (what you actually set) */
    double MinKernelRadius;			/*!< minimum allowed gas kernel length */
    double MaxKernelRadius;           /*!< minimum allowed gas kernel length */

  double SofteningGas,		/*!< for type 0 */
    SofteningHalo,		/*!< for type 1 */
    SofteningDisk,		/*!< for type 2 */
    SofteningBulge,		/*!< for type 3 */
    SofteningStars,		/*!< for type 4 */
    SofteningBndry;		/*!< for type 5 */

  double SofteningGasMaxPhys,	/*!< for type 0 */
    SofteningHaloMaxPhys,	/*!< for type 1 */
    SofteningDiskMaxPhys,	/*!< for type 2 */
    SofteningBulgeMaxPhys,	/*!< for type 3 */
    SofteningStarsMaxPhys,	/*!< for type 4 */
    SofteningBndryMaxPhys;	/*!< for type 5 */

  double ForceSoftening[6];	/*!< current (comoving) gravitational softening lengths for each particle type -- multiplied by a factor 1/KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER to define the maximum kernel extent - at that scale the force is Newtonian */

  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this value, * allowing the size of the snapshot files to be reduced */
  double MassTable[6];

  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    RestartFile[100], ResubmitCommand[100], OutputListFilename[100];

#ifdef COOL_GRACKLE
    char GrackleDataFile[100];
#endif
    /*! table with desired output times */
    double OutputListTimes[MAXLEN_OUTPUTLIST];
    char OutputListFlag[MAXLEN_OUTPUTLIST];
    int OutputListLength;		/*!< number of times stored in table of desired output times */

#ifdef TURB_DRIVING
    double TurbInjectedEnergy;
    double TurbDissipatedEnergy;
#if defined(TURB_DRIVING_SPECTRUMGRID)
    double TimeBetTurbSpectrum;
    double TimeNextTurbSpectrum;
    int FileNumberTurbSpectrum;
#endif
#endif

#ifdef RADTRANSFER
    integertime Radiation_Ti_begstep;
    integertime Radiation_Ti_endstep;
#endif
#ifdef TRANSPORT_SUBCYCLE
    int Transport_Subcycle_N;               /*!< number of transport subcycles this step */
    double Transport_Subcycle_dt_fraction;   /*!< = 1.0/Transport_Subcycle_N, fraction of full dt per sub-step */
#endif
#ifdef RT_EVOLVE_INTENSITIES
    double Rad_Intensity_Direction[N_RT_INTENSITY_BINS][3];
#endif

#if (defined(SINGLE_STAR_FB_SNE) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)) || defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
    double SN_Ejecta_Direction[SINGLE_STAR_FB_SNE_N_EJECTA][3];
#endif

#if defined(RT_CHEM_PHOTOION) && !(defined(GALSF_FB_FIRE_RT_HIIHEATING) || defined(GALSF))
    double IonizingLuminosityPerSolarMass_cgs;
    double star_Teff;
#endif

#ifdef RT_ISRF_BACKGROUND
    double InterstellarRadiationFieldStrength;
    double RadiationBackgroundRedshift;
#endif
#ifdef RT_INFRARED
    double InitRadiationTemp;
#endif

#ifdef RT_LEBRON
    double PhotonMomentum_Coupled_Fraction;
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    double PhotonMomentum_fUV;
    double PhotonMomentum_fOPT;
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    double Sink_Rad_MomentumFactor;
#endif

#ifdef GRAIN_FLUID
#define GRAIN_PTYPES 8 /* default to allowed particle type for grains == 3, only, but can make this a more extended list as desired */
#ifdef GRAIN_RDI_TESTPROBLEM
#if(NUMDIMS==3)
#define GRAV_DIRECTION_RDI 2
#else
#define GRAV_DIRECTION_RDI 1
#endif
    double Grain_Charge_Parameter;
    double Dust_to_Gas_Mass_Ratio;
    double Vertical_Gravity_Strength;
    double Vertical_Grain_Accel;
    double Vertical_Grain_Accel_Angle;
#ifdef BOX_SHEARING
    double Pressure_Gradient_Accel;
#endif
#ifdef RT_OPACITY_FROM_EXPLICIT_GRAINS
    double Grain_Q_at_MaxGrainSize;
#endif
#endif // GRAIN_RDI_TESTPROBLEM
    double Grain_Internal_Density;
    double Grain_Size_Min;
    double Grain_Size_Max;
    double Grain_Size_Spectrum_Powerlaw;
#endif
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS) && defined(RT_GENERIC_USER_FREQ)
    double Grain_Absorbed_Fraction_vs_Total_Extinction;
#endif

#ifdef PIC_MHD
    double PIC_Charge_to_Mass_Ratio;
#endif

#ifdef COSMIC_RAY_FLUID
    double CosmicRayDiffusionCoeff;
#if (N_CR_PARTICLE_BINS > 2)
    double CR_global_charge_in_bin[N_CR_PARTICLE_BINS];
#endif
#if defined(CRFLUID_EVOLVE_SPECTRUM)
#define N_CR_SPECTRUM_LUT 101
    double CR_global_min_rigidity_in_bin[N_CR_PARTICLE_BINS];
    double CR_global_max_rigidity_in_bin[N_CR_PARTICLE_BINS];
    double CR_global_rigidity_at_bin_center[N_CR_PARTICLE_BINS];
    int CR_species_ID_in_bin[N_CR_PARTICLE_BINS];
    double CR_global_slope_lut[N_CR_PARTICLE_BINS][N_CR_SPECTRUM_LUT];
    int CR_secondary_species_listref[N_CR_PARTICLE_SPECIES][N_CR_PARTICLE_SPECIES];
    int CR_secondary_target_bin[N_CR_PARTICLE_BINS][N_CR_PARTICLE_SPECIES];
    double CR_frag_secondary_coeff[N_CR_PARTICLE_BINS][N_CR_PARTICLE_SPECIES];
    double CR_frag_coeff[N_CR_PARTICLE_BINS];
    double CR_rad_decay_coeff[N_CR_PARTICLE_BINS];
    int CR_species_ID_active_list[N_CR_PARTICLE_SPECIES];
#endif
#endif

#ifdef GALSF		/* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
  double MaxSfrTimescale;

#ifdef GALSF_EFFECTIVE_EQS
  double EgySpecSN;
  double FactorSN;
  double EgySpecCold;
  double FactorEVP;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double FactorForSofterEQS;
#endif

#ifdef GALSF_FB_FIRE_RT_LOCALRP
  double RP_Local_Momentum_Renormalization;
#endif

#ifdef GALSF_SUBGRID_WINDS
#ifndef GALSF_SUBGRID_WIND_SCALING
#define GALSF_SUBGRID_WIND_SCALING 0 // default to constant-velocity winds //
#endif
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelMaxTimeFactor;  /* maximum free travel time in units of the Hubble time at the current simulation redshift */
  double WindFreeTravelDensFac;
#if (GALSF_SUBGRID_WIND_SCALING>0)
  double VariableWindVelFactor;  /* wind velocity in units of the halo escape velocity */
  double VariableWindSpecMomentum;  /* momentum available for wind per unit mass of stars formed, in internal velocity units */
#endif
#endif // GALSF_SUBGRID_WINDS //

#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
    double SNe_Energy_Renormalization;
    double StellarMassLoss_Rate_Renormalization;
    double StellarMassLoss_Energy_Renormalization;
#if defined(COSMIC_RAY_FLUID) || defined(COSMIC_RAY_SUBGRID_LEBRON)
#define CR_DYNAMICAL_INJECTION_IN_SNE
    double CosmicRay_SNeFraction;
#endif
#endif
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
    double HIIRegion_fLum_Coupled;
#endif

#ifdef GALSF_FB_FIRE_AGE_TRACERS
    double AgeTracerRateNormalization;              /* Determines Fraction of time to do age tracer deposition (with checks depending on time bin width for current star) */
#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
    double AgeTracerTimeBins[NUM_AGE_TRACERS+1];    /* Bin edges (left) for stellar age passive scalar tracers when using custom (uneven) bins the final value is the right edge of the final bin, hence a total size +1 the number of tracers */
    char   AgeTracerListFilename[100];              /* file name to read ages from (in Myr) as a single column */
#else
    double AgeTracerBinStart;                       /* left bin edge of first age tracers (Myr) - for log spaced bins */
    double AgeTracerBinEnd;                         /* right bin edge of last age tracer (Myr)  - for log spaced bins */
#endif
#endif

#endif // GALSF

#ifdef GALSF_LIMIT_FBTIMESTEPS_FROM_BELOW
    double Dt_Since_LastFBCalc_Gyr; // time since last feedback event occurred, needs to be set
    double Dt_Min_Between_FBCalc_Gyr; // minimum timestep to enforce between feedback calculations, for optimization
#endif

#if (defined(GALSF) && defined(METALS)) || defined(COOL_METAL_LINES_BY_SPECIES) || defined(GALSF_FB_FIRE_RT_LOCALRP) || defined(GALSF_FB_FIRE_RT_HIIHEATING) || defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_FIRE_RT_LONGRANGE) || defined(GALSF_FB_THERMAL)
#define INIT_STELLAR_METALS_AGES_DEFINED // convenience flag for later to know these variables exist
    double InitMetallicityinSolar;
    double InitStellarAgeinGyr;
#endif

#if defined(SINK_WIND_KICK) || defined(SINK_WIND_SPAWN)
    double Sink_accreted_fraction;
    double Sink_outflow_velocity;
#endif

#if defined(SINGLE_STAR_FB_JETS)
    double Sink_outflow_jetlaunchvelscaling; // scales the amount of accretion power going into jets, we eject (1-All.Sink_accreted_fraction) fraction of the accreted mass at this value times the Keplerian velocity at the protostellar radius. If set to 1 then the mass and power loading of the jets are both (1-All.Sink_accreted_fraction)
#endif

#if defined(SINK_COSMIC_RAYS)
    double Sink_CosmicRay_Injection_Efficiency;
#endif

#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double CosmicRay_Subgrid_Vstream_0;
    double CosmicRay_Subgrid_Kappa_0;
#endif


#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
    double RHD_bins_nu_min_ev[N_RT_FREQ_BINS]; /* minimum frequency of the radiation 'bin' in eV */
    double RHD_bins_nu_max_ev[N_RT_FREQ_BINS]; /* maximum frequency of the radiation 'bin' in eV */
#endif
#ifdef RT_CHEM_PHOTOION
    double rt_ion_nu_min[N_RT_FREQ_BINS];
    double rt_nu_eff_eV[N_RT_FREQ_BINS];
    double rt_ion_precalc_stellar_luminosity_fraction[N_RT_FREQ_BINS];
    double rt_ion_sigma_HI[N_RT_FREQ_BINS];
    double rt_ion_sigma_HeI[N_RT_FREQ_BINS];
    double rt_ion_sigma_HeII[N_RT_FREQ_BINS];
    double rt_ion_G_HI[N_RT_FREQ_BINS];
    double rt_ion_G_HeI[N_RT_FREQ_BINS];
    double rt_ion_G_HeII[N_RT_FREQ_BINS];
#endif

#ifdef METALS
    double SolarAbundances[NUM_METAL_SPECIES];
#ifdef COOL_METAL_LINES_BY_SPECIES
    int SpeciesTableInUse;
#endif
#endif


#if defined(GALSF_ISMDUSTCHEM_MODEL)
    double Initial_ISMDustChem_Depletion; /* initial depletion for silicate dust species if defined */
    double Initial_ISMDustChem_SiliconToCarbonRatio; /* sets rough mass ratio between silicates are carbonaceous dust for given initial depletion */
    double ISMDustChem_AtomicMassTable[NUM_ISMDUSTCHEM_ELEMENTS]; /* atomic mass for each element in metallicity field */
    double ISMDustChem_SNeSputteringShutOffTime; /* amount of time to turn off thermal sputtering after SNe event to avoid double counting dust destruction */
    int ISMDustChem_SilicateMetallicityFieldIndexTable[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES]; /* index in metallicity field for elements which make up silicate dust (O, Mg, Si, and Fe) */
    double ISMDustChem_SilicateNumberOfAtomsTable[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES]; /* number of O, Mg, Si, and Fe in one formula unit of silicate dust */
    double ISMDustChem_EffectiveSilicateDustAtomicWeight; /* atomic weight of one formula unit of silicate dust, depends on which optional module you use */
    // Scaling arguements from parameter file used to adjust each dust process
    double ISMDustChem_SNeIIDustScaling;
    double ISMDustChem_SNeIaDustScaling;
    double ISMDustChem_AGBDustScaling;
    double ISMDustChem_DustAccretionScaling;
    double ISMDustChem_ThermalSputteringScaling;
    double ISMDustChem_AccretionTcutoffScaling;
    double ISMDustChem_SNeGasClearedOfDustScaling;
    double ISMDustChem_SpeciesBulkDens[3]; /* condensed bulk density for silicates, carbonaceous, and metallic iron */
    int ISMDustChem_TrackedSpeciesIDTable[NUM_ISMDUSTCHEM_SPECIES]; /* contains unique ID numbers for each tracked dust species which correspond to their location in ISMDustChem_SpeciesFieldIndexTable. Returns -1 for untracked species  */
    int ISMDustChem_SpeciesFieldIndexTable[NUM_ISMDUSTCHEM_SPECIES_IDS]; /* index in dust species field for given dust species. Sparse table indexed by fixed species ID (0..NUM_ISMDUSTCHEM_SPECIES_IDS-1), returns the packed field slot or -1 if untracked. */
    int ISMDustChem_Sil_Index;
    int ISMDustChem_Carb_Index;
    int ISMDustChem_FreeIron_Index ;
    int ISMDustChem_ORes_Index;
    int ISMDustChem_InclIron_Index;
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
    double UnitGrainNumber; /* factor to convert internal grain number unit to number of grains */
    double UnitGrainLength_in_cm; /* factor to convert internal grain length unit to cm */
    double ISMDustChem_SNeShatteringScaling;
    double ISMDustChem_SNeSputteringScaling;
    double ISMDustChem_ShatteringScaling;
    double ISMDustChem_CoagulationScaling;
    double ISMDustChem_VCoagScaling;
    double ISMDustChem_CoagDensityEnhancementScaling;
    double ISMDustChem_GrainVelocityScaling;
    double ISMDustChem_PhotodestructionScaling; /* dispersion of grain velocities in the ISM */
    double ISMDustChem_Grain_Size_Min;
    double ISMDustChem_Grain_Size_Max;
    double ISMDustChem_GrainBinSize; /* bin width of logarithmically spaced grain sizes */
    double ISMDustChem_GrainBinEdges[NUM_ISMDUSTCHEM_SIZE_BINS+1]; /* edges of each grain size bin */
    double ISMDustChem_GrainBinCenters[NUM_ISMDUSTCHEM_SIZE_BINS]; /* centers of each grain size bin in log space */
    double ISMDustChem_C_NiNj[NUM_ISMDUSTCHEM_SIZE_BINS][NUM_ISMDUSTCHEM_SIZE_BINS]; /* pre-computed coefficients for coagulation/shattering polynomial */
    double ISMDustChem_C_Njsi[NUM_ISMDUSTCHEM_SIZE_BINS][NUM_ISMDUSTCHEM_SIZE_BINS];
    double ISMDustChem_C_Nisj[NUM_ISMDUSTCHEM_SIZE_BINS][NUM_ISMDUSTCHEM_SIZE_BINS];
    double ISMDustChem_C_sisj[NUM_ISMDUSTCHEM_SIZE_BINS][NUM_ISMDUSTCHEM_SIZE_BINS];
#endif
#endif


#ifdef GR_TABULATED_COSMOLOGY
  double DarkEnergyConstantW;	/*!< fixed w for equation of state */
#if defined(GR_TABULATED_COSMOLOGY_W) || defined(GR_TABULATED_COSMOLOGY_G) || defined(GR_TABULATED_COSMOLOGY_H)
#ifndef GR_TABULATED_COSMOLOGY_W
#define GR_TABULATED_COSMOLOGY_W
#endif
  char TabulatedCosmologyFile[100];	/*!< tabulated parameters for expansion and/or gravity */
#ifdef GR_TABULATED_COSMOLOGY_G
  double Gini;
#endif
#endif
#endif

#ifdef SPHAV_CD10_VISCOSITY_SWITCH
  double ViscosityAMin;
  double ViscosityAMax;
#endif

#ifdef TURB_DIFFUSION
  double TurbDiffusion_Coefficient;
#ifdef TURB_DIFF_DYNAMIC
  double TurbDynamicDiffFac;
  int TurbDynamicDiffIterations;
  double TurbDynamicDiffSmoothing;
  double TurbDynamicDiffMax;
#endif
#endif

#if defined(CONDUCTION)
   double ConductionCoeff;	/*!< Thermal Conductivity */
#endif

#if defined(VISCOSITY)
   double ShearViscosityCoeff;
   double BulkViscosityCoeff;
#endif

#if defined(CONDUCTION_SPITZER) || defined(VISCOSITY_BRAGINSKII)
    double ElectronFreePathFactor;	/*!< Factor to get electron mean free path */
#endif

#ifdef MAGNETIC
#ifdef MHD_B_SET_IN_PARAMS
  double BiniX, BiniY, BiniZ;	/*!< Initial values for B */
#endif
#ifdef SPH_TP12_ARTIFICIAL_RESISTIVITY
  double ArtMagDispConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial magnetic disipation */
#endif
#ifdef DIVBCLEANING_DEDNER
  double FastestWaveSpeed;
  double FastestWaveDecay;
  double DivBcleanParabolicSigma;
  double DivBcleanHyperbolicSigma;
#endif
#endif /* MAGNETIC */

#if (defined(SINK_PARTICLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
  double TimeNextOnTheFlyFoF;
  double TimeBetOnTheFlyFoF;
#endif

#ifdef SINK_PARTICLES
  double SinkAccretionFactor;	/*!< Rescale sink accretion rate normaliation */
  double SinkFeedbackFactor;	/*!< Rescale sink feedback normalization */
  double SeedSinkMass;          /*!< Seed sink particle mass */
#if defined(SINK_SEED_FROM_FOF) || defined(SINK_SEED_FROM_LOCALGAS)
  double SeedSinkMassSigma;     /*!< Standard deviation of initial sink particle masses */
  double SeedSinkMinRedshift;   /*!< Minimum redshift where sink seeds are allowed */
#ifdef SINK_SEED_FROM_LOCALGAS
  double SeedSinkPerUnitMass;   /*!< Defines probability per unit mass of seed sink forming */
#endif
#endif
#ifdef SINK_ALPHADISK_ACCRETION
  double SeedReservoirMass;         /*!< Seed alpha disk mass */
#endif
#ifdef SINK_WIND_SPAWN
  double Sink_outflow_particlemass; /*!< target mass for feedback particles to be spawned */
  double Sink_outflow_temperature;
  MyIDType SpawnedWindCellID;
#ifdef SINGLE_STAR_FB_WINDS
  double Cell_Spawn_Mass_ratio_MS;  /*!< target mass for feedback particles to be spawned for main sequence winds in STARFORGE*/
#endif
#endif
#ifdef SINK_SEED_FROM_FOF
  double MinFoFMassForNewSeed;      /*!< Halo mass required before new seed is put in */
#endif
  double SinkNgbFactor;             /*!< Factor by which the gas neighbour count should be increased/decreased */
  double SinkMaxAccretionRadius;
  double SinkEddingtonFactor;	    /*!< Factor above Eddington */
  double SinkRadiativeEfficiency;   /*!< Radiative efficiency determined by the spin value, default value is 0.1 */
#endif

#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)
  double Tillotson_EOS_params[7][18]; /*! < holds parameters for Tillotson EOS for solids; slots 12-15 (k_Weibull, m_Weibull, mu_DP, alpha_0) and slots 16-17 (P_e, P_s) are used only under EOS_DAMAGE_POROSITY */
#endif

#ifdef EOS_TABULATED
    char EosTable[100];
#endif

#ifdef EOS_ANEOS
#ifndef ANEOS_MAX_MATERIALS
#define ANEOS_MAX_MATERIALS 7
#endif
    int  AneosNumMaterials;                       /* number of ANEOS material tables to load */
    char AneosTableFiles[ANEOS_MAX_MATERIALS][256]; /* file paths for each material's SESAME table */
#endif

#ifdef RANDOMIZE_GRAVTREE_PERIODIC
    Vec3<double> RandomShift;   /*!< how far the coordinates presently sit from the frame the initial conditions were written in, kept in [0,box). Subtracted back out wherever a physical coordinate is reported. Lives here so a restart resumes in the frame it stopped in. */
#endif

#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    Vec3<double> SpecialParticle_Position_ForRefinement[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM];
    double Mass_Accreted_By_SpecialParticle[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM];
    double Mass_of_SpecialParticle[SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM];
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
  double AGS_DesNumNgb;
  double AGS_MaxNumNgbDeviation;
#endif

#ifdef DM_FUZZY
  double ScalarField_hbar_over_mass;
#endif

#ifdef TURB_DRIVING
  double TurbDriving_Global_DecayTime;
  double TurbDriving_Global_AccelerationPowerVariable;
  double TurbDriving_Global_DtTurbUpdates;
  double TurbDriving_Global_DrivingScaleKMinVar;
  double TurbDriving_Global_DrivingScaleKMaxVar;
  double TurbDriving_Global_SolenoidalFraction;
  int    TurbDriving_Global_DrivingSpectrumKey;
  int    TurbDriving_Global_DrivingRandomNumberKey;
#endif

#if defined(COOLING) && defined(COOL_GRACKLE)
    code_units GrackleUnits;
#endif

#if defined(SINK_WIND_SPAWN_SET_BFIELD_POLTOR)
  double Sink_spawn_injectionradius;
  double B_spawn_pol;
  double B_spawn_tor;
#endif
#ifdef SINK_WIND_SPAWN_SET_JET_PRECESSION
  double Sink_jet_precess_degree;
  double Sink_jet_precess_period;
#endif
#ifdef NUCLEAR_NETWORK
    char NuclearNetworkDataFile[256];
    double NuclearBurningFloor_T;
    double NuclearBurningFloor_rho;
    double NuclearNSE_T_threshold;
#endif
};

#endif /* GLOBAL_DATA_ALL_STRUCT_H */
