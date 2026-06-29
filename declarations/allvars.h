
/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "../declarations/allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified extensively
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO (many new variables,
 * structures, and different naming conventions for some old variables).
 * The declarations have also been divided into a new scheme which improves
 * read-ability and ease of use and separates all of the macros written for
 * GIZMO.
 */


#ifndef ALLVARS_H
#define ALLVARS_H

#define GIZMO_VERSION     2026  /*!< code version (should be an int corresponding to the year) */

#include <mpi.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <hdf5.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../GIZMO_config.h"

#include "precompiler_logic.h"

#include "constants.h"

#include "../eos/eos.h"

#ifdef COOL_GRACKLE
#include <grackle.h>
#endif

#include "../system/tags.h"

#include <assert.h>
#include <algorithm>
#include <type_traits>
#include <vector>

#include "macros.h"

#include "typedefs.h"
#include "../math_types/vec3.h"
#include "../math_types/mat3.h"
#include "../math_types/symmetric_tensor2.h"

#ifdef CHIMES
#include "../cooling/chimes/chimes_proto.h"
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

#ifdef BOX_PERIODIC
extern MyDouble boxSize, boxHalf;
#else
#define boxSize (All.BoxSize)
#define boxHalf (0.5*All.BoxSize)
#endif
#ifdef BOX_LONG_X
extern MyDouble boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef BOX_LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef BOX_LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif

#ifdef BOX_SHEARING
extern MyDouble Shearing_Box_Vel_Offset;
extern MyDouble Shearing_Box_Pos_Offset;
#endif

#ifdef RANDOMIZE_GRAVTREE
/* the periodic path shifts all coordinates by a random vector mod box, which is only valid
 * where every periodic dimension is a pure translation symmetry */
#if defined(BOX_SHEARING)
#error "RANDOMIZE_GRAVTREE is incompatible with BOX_SHEARING: an x-translation couples into the shear velocity offset (see do_box_wrapping) and corrupts velocities."
#endif
#if defined(BOX_PERIODIC) && (defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z))
#error "RANDOMIZE_GRAVTREE + BOX_PERIODIC is incompatible with reflecting/outflow boundaries: those dimensions are not translation-invariant."
#endif
/* Which randomization method applies is decided by whether GRAVITY is periodic, NOT by whether
 * the box is. Under GRAVITY_NOT_PERIODIC the tree uses bare |dx| (see the GRAVITY_NGB_* macros
 * in forcetree.cc), so translating coordinates mod box is NOT a symmetry: wrapping a particle
 * across the boundary moves it a full box away gravitationally, which would corrupt the forces.
 * Those runs (STARFORGE-style setups, and anything BOX_SHEARING implies) must therefore use the
 * non-periodic method -- moving and enlarging the root node -- even though the box is periodic.
 * The cost is the root-node doubling and the bit of Peano resolution it spends. */
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
#define RANDOMIZE_GRAVTREE_PERIODIC /* use the coordinate-translation method */
#endif
#endif

#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z)
extern short int special_boundary_condition_xyz_def_reflect[3];
extern short int special_boundary_condition_xyz_def_outflow[3];
#endif

template<typename T> inline void nearest_xyz(Vec3<T>& v, int sign=1) { NEAREST_XYZ(v[0], v[1], v[2], sign); }


#ifdef CHIMES
extern struct gasVariables *ChimesGasVars;
extern struct globalVariables ChimesGlobalVars;
extern char ChimesDataPath[DEFAULT_PATH_BUFFERSIZE_TOUSE];
extern char ChimesEqAbundanceTable[DEFAULT_PATH_BUFFERSIZE_TOUSE];
extern char ChimesPhotoIonTable[DEFAULT_PATH_BUFFERSIZE_TOUSE];
extern double chimes_rad_field_norm_factor;
extern double shielding_length_factor;
extern double cr_rate;
extern int ChimesEqmMode;
extern int ChimesUVBMode;
extern int ChimesInitIonState;
extern int Chimes_incl_full_output;
extern int N_chimes_full_output_freq;
#endif // CHIMES
#ifdef CHIMES_METAL_DEPLETION
#define DEPL_N_ELEM 17
struct Chimes_depletion_data_structure
{
    double SolarAbund[DEPL_N_ELEM];
    double DeplPars[DEPL_N_ELEM][3];
    double DustToGasSaturated;
    double ChimesDepletionFactors[7];
    double ChimesDustRatio;
};
extern struct Chimes_depletion_data_structure *ChimesDepletionData;
#endif // CHIMES_METAL_DEPLETION


extern std::vector<int> ActiveParticleList;
extern unsigned char *ProcessedFlag;
extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountGas[TIMEBINS];
extern int TimeBinActive[TIMEBINS];
extern int FirstInTimeBin[TIMEBINS];
extern int LastInTimeBin[TIMEBINS];
extern std::vector<int> NextInTimeBin;
extern std::vector<int> PrevInTimeBin;
#ifdef GALSF
extern double TimeBinSfr[TIMEBINS];
#endif

#ifdef SINK_PARTICLES
extern double TimeBin_Sink_mass[TIMEBINS];
extern double TimeBin_Sink_dynamicalmass[TIMEBINS];
extern double TimeBin_Sink_Mdot[TIMEBINS];
extern double TimeBin_Sink_Medd[TIMEBINS];
#endif


extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern int PTask;		/*!< note: NTask = 2^PTask */
extern double CPUThisRun;	/*!< Sums CPU time of current process */
extern int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;
extern int NumGasUpdate;	/*!< number of active gas cells on local processor in current timestep  */
extern int MaxTopNodes;	        /*!< Maximum number of nodes in the top-level tree used for domain decomposition */
extern int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from initial conditions, 1 is resuming a run from a set of restart files, while 2 marks a restart from a snapshot file. */
extern int RestartSnapNum;
extern int SelRnd;
extern int TakeLevel;
extern int *Exportflag;	        /*!< Buffer used for flagging whether a particle needs to be exported to another process */
extern int *Exportnodecount;
extern int *Exportindex;
extern int *Send_offset, *Send_count, *Recv_count, *Recv_offset;
extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;
extern double CPU_Step[CPU_PARTS];
extern char CPU_Symbol[CPU_PARTS];
extern char CPU_SymbolImbalance[CPU_PARTS];
extern char CPU_String[CPU_STRING_LEN + 1];
extern double WallclockTime;    /*!< This holds the last wallclock time measurement for timings measurements */
extern int Flag_FullStep;	/*!< Flag used to signal that the current step involves all particles */
extern size_t HighMark_run,  HighMark_domain, HighMark_gravtree, HighMark_pmperiodic,
  HighMark_pmnonperiodic,  HighMark_gasdensity, HighMark_hydro, HighMark_GasGrad;

#ifdef TURB_DRIVING
extern size_t HighMark_turbpower;
#endif
extern int TreeReconstructFlag; /*!< the standing tree is CONDEMNED: its storage/topology must be
    replaced by a rearrange+treebuild before the next walk. Consumed at the step-start ladder
    (core/run.cc) or mid-step by gravity_tree()/compute_potential(); both decomposition paths
    self-raise it to defer the actual build to the next tree consumer. Distinct from
    DomainReconstructFlag (ownership/balance decomposition owed -- strictly more expensive) and from
    TreeMomentsStaleFlag below (moments re-derivation on intact topology -- strictly cheaper).
    Serialized in restartfiles. */
extern int TreeMomentsStaleFlag; /*!< flag to refresh tree node moments without a full tree rebuild, e.g. after star formation or sink mass changes */
extern int DomainReconstructFlag; /*!< a full domain decomposition (particle ownership/balance) is owed.
    Consumed ONLY by the step-start ladder in run.cc; cleared by the decomposition itself, which
    self-raises TreeReconstructFlag for the tree half. Physics events that merely condemn the tree
    raise TreeReconstructFlag instead and cost a rebuild, not a decomposition. Serialized in
    restartfiles. */
extern int NtotSwallowedThisStep; /*!< global (Allreduce'd) count of particles swallowed in this step's sink pass; gates the cleanup rearrange in run.cc. Rank-uniform by construction. */
extern int TreeWalkValidatePending; /*!< set (rank-local) when a maintained rearrange changed the walk threading; the next tree consumer runs force_validate_tree_links() and clears it. Cleared by force_treebuild (fresh tree). */
enum TreeOpsIndex {TREEOPS_BUILD=0, TREEOPS_REFRESH, TREEOPS_DECOMP, TREEOPS_DECOMP_LIGHT, TREEOPS_DEFER, TREEOPS_ESCALATE, TREEOPS_REARRANGE, TREEOPS_N};
extern long long ForceAddElementToTree_CallsSinceBuild; /*!< rank-local insertions accepted by the standing tree since its build; the ladder SUM-reduces this and condemns the tree past a fixed fraction of TotNumPart (insertions attach at existing nodes, so tree quality degrades with their count) */
extern long long TreeOpsCount[TREEOPS_N]; /*!< running counts of whole-tree builds, moment refreshes, full/light decompositions, deferred-rebuild ladder steps, and rearranges. All are collective events, so every rank agrees; reported in cpu.txt. Diagnostic only: not serialized, resets on resume. */
extern double TreeMaintTime_SwapPointers; /*!< rank-local cumulative seconds in swap_treewalk_pointers (the per-swap O(tree) pointer surgery MAINTAIN_TREE_IN_REARRANGE pays); MAX across ranks in cpu.txt */
extern double TreeMaintTime_Rearrange; /*!< rank-local cumulative seconds in rearrange_particle_sequence; MAX across ranks in cpu.txt */
extern long long TreeMaintCount_SwapPointers; /*!< rank-local count of swap_treewalk_pointers calls */
extern int TreeAuditMomentsPending; /*!< set when a rearrange edited the list; the next moments refresh gets a FULL audit regardless of the sampling stride (TREE_INTEGRITY_AUDITS builds) */
extern int TreeDiscardImminent; /*!< held by the decomposition paths across their merge/split+rearrange prelude: the standing tree dies before the next walk, so split daughters are deliberately not inserted and the rearrange-time audit must not treat them as orphans */
extern int GlobFlag;
extern char DumpFlag;
#ifdef WAKEUP
extern int NeedToWakeupParticles;
extern int NeedToWakeupParticles_local;
#endif

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
#ifdef SINK_WIND_SPAWN
extern double Max_Unspawned_MassUnits_fromSink;
#endif

extern long long Ntype[6];	/*!< total number of particles of each type */
extern int NtypeLocal[6];	/*!< local number of particles of each type */
extern gsl_rng *random_generator;	/*!< the random number generator used */
extern int Gas_split;           /*!< current number of newly-spawned gas particles outside block */
#ifdef GALSF
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
#endif

extern double TimeOfLastTreeConstruction;	/*!< holds what it says */
extern std::vector<int> Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */
extern double *R2ngblist;
extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int *DomainStartList, *DomainEndList;
extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountGas;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern int *DomainList, DomainNumChanged;
extern peanokey *Key, *KeySorted;


#ifdef HERMITE_INTEGRATION
extern int HermiteOnlyFlag;     /*!< flag to only do Hermite integration for applicable particles (ie. stars) in the gravity routine - set =1 on the first prediction pass and =2 on the second correction pass */
#endif

#ifdef RT_CHEM_PHOTOION
extern double rt_ion_nu_min[N_RT_FREQ_BINS];
extern double rt_nu_eff_eV[N_RT_FREQ_BINS];
extern double rt_ion_precalc_stellar_luminosity_fraction[N_RT_FREQ_BINS];
extern double rt_ion_sigma_HI[N_RT_FREQ_BINS];
extern double rt_ion_sigma_HeI[N_RT_FREQ_BINS];
extern double rt_ion_sigma_HeII[N_RT_FREQ_BINS];
extern double rt_ion_G_HI[N_RT_FREQ_BINS];
extern double rt_ion_G_HeI[N_RT_FREQ_BINS];
extern double rt_ion_G_HeII[N_RT_FREQ_BINS];
#endif

#if defined(CRFLUID_EVOLVE_SPECTRUM) /* define some global variables we will need to use semi-constantly to make reference to the CR spectra */
extern double CR_global_min_rigidity_in_bin[N_CR_PARTICLE_BINS];
extern double CR_global_max_rigidity_in_bin[N_CR_PARTICLE_BINS];
extern double CR_global_rigidity_at_bin_center[N_CR_PARTICLE_BINS];
extern double CR_global_charge_in_bin[N_CR_PARTICLE_BINS];
extern int CR_species_ID_in_bin[N_CR_PARTICLE_BINS];
#define N_CR_SPECTRUM_LUT 101 /*!< number of elements per bin in the look-up-tables we will pre-compute to use for inverting the energy-number relation to determine the spectral slope */
extern double CR_global_slope_lut[N_CR_PARTICLE_BINS][N_CR_SPECTRUM_LUT]; /*!< holder for the actual look-up-tables */
extern int CR_secondary_species_listref[N_CR_PARTICLE_SPECIES][N_CR_PARTICLE_SPECIES]; /*!< list for each type of the different secondaries to which it can decay */
extern int CR_secondary_target_bin[N_CR_PARTICLE_BINS][N_CR_PARTICLE_SPECIES]; /*!< destination bin for the secondaries produced by different primaries */
extern double CR_frag_secondary_coeff[N_CR_PARTICLE_BINS][N_CR_PARTICLE_SPECIES]; /*!< coefficients for fragmentation to the given secondaries (also pre-computed for simplicity) */
extern double CR_frag_coeff[N_CR_PARTICLE_BINS]; /*!< total coefficients for fragmentation processes (pre-compute b/c cross-sections are complicated) */
extern double CR_rad_decay_coeff[N_CR_PARTICLE_BINS]; /*!< radioactive decay coefficients (pre-computed for ease, also because of dilation dependence) */
extern int CR_species_ID_active_list[N_CR_PARTICLE_SPECIES]; /*!< holds the list of species ids to loop over */
#endif


extern int MaxNodes;        /*!< maximum allowed number of internal nodes */
extern int Numnodestree;    /*!< number of (internal) nodes in each tree */

extern int *Nextnode;        /*!< gives next node in tree walk  (nodes array) */
extern int *Father;        /*!< gives parent node in tree (Prenodes array) */

extern int maxThreads;

#ifdef TURB_DRIVING // other global variables for forcing field (need to be carried through all timesteps)
extern double* StOUPhases; // random fluctuating component of the amplitudes
extern double* StAmpl; // relative amplitude for each k
extern double* StAka; // phases (real part)
extern double* StAkb; // phases (imag part)
extern double* StMode; // k vectors
extern int StNModes; // total number of modes
extern integertime StTPrev; // time of last update (to determine when next will be)
extern gsl_rng* StRng; // random number generator key
#endif


#if defined(DM_SIDM)
#define GEOFACTOR_TABLE_LENGTH 1000    /*!< length of the table used for the geometric factor spline */
extern MyDouble GeoFactorTable[GEOFACTOR_TABLE_LENGTH];
#endif
extern int NTopnodes, NTopleaves;

/* variables for input/output , usually only used on process 0 */
extern char ParameterFile[100];    /*!< file name of parameterfile used for starting the simulation */
extern FILE
#ifdef OUTPUT_ADDITIONAL_RUNINFO
*FdTimebin,    /*!< file handle for timebin.txt log-file. */
*FdInfo,       /*!< file handle for info.txt log-file. */
*FdEnergy,     /*!< file handle for energy.txt log-file. */
*FdTimings,    /*!< file handle for timings.txt log-file. */
*FdBalance,    /*!< file handle for balance.txt log-file. */
#ifdef RT_CHEM_PHOTOION
*FdPhotoIonChemStats,        /*!< file handle for photoionization chemistry log-file. */
#endif
#ifdef TURB_DRIVING
*FdTurb,       /*!< file handle for turbulent driving log-file */
#endif
#ifdef GR_TABULATED_COSMOLOGY
*FdDE,         /*!< file handle for darkenergy.txt log-file. */
#endif
#endif
*FdCPU;        /*!< file handle for cpu.txt log-file. */
#ifdef GALSF
extern FILE *FdSfr;        /*!< file handle for star formation log-file. */
#endif
#ifdef GALSF_FB_FIRE_RT_LOCALRP
extern FILE *FdMomWinds;    /*!< file handle for local photon-momentum subgrid model log-file */
#endif
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
extern FILE *FdHIIHeating;    /*!< file handle for local HII model log-file */
#endif
#ifdef GALSF_FB_MECHANICAL
extern FILE *FdSNeFBLogFile;    /*!< file handle for mechanical feedback log-file */
#endif
#ifdef SINK_PARTICLES
extern FILE *FdSinks;    /*!< file handle for sinks.txt log-file. */
#endif
#ifdef OUTPUT_SINK_ACCRETION_HIST
extern FILE *FdSinkSwallowDetails;
#endif
#ifdef OUTPUT_SINK_FORMATION_PROPS
extern FILE *FdSinkFormationDetails;
#endif
#if (defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(SINK_OUTPUT_MOREINFO)) && defined(SINK_PARTICLES)
extern FILE *FdSinksDetails;
#ifdef SINK_OUTPUT_MOREINFO
extern FILE *FdSinkMergerDetails;
#endif
#ifdef SINK_WIND_KICK
extern FILE *FdSinkWindDetails;
#endif
#endif


extern double DriftTable[DRIFT_TABLE_LENGTH]; /*! table for the cosmological drift factors */
extern double GravKickTable[DRIFT_TABLE_LENGTH]; /*! table for the cosmological kick factor for gravitational forces */
extern void *CommBuffer;    /*!< points to communication buffer, which is used at a few places */


#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
#endif

extern peanokey *DomainKeyBuf;

/* Tree variables */
extern long Nexport, Nimport;
extern int BufferCollisionFlag;
extern int BufferFullFlag;
extern int NextParticle;
extern int NextJ;
extern int TimerFlag;




extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  MyFloat GravCost;
  int Daughter;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;



#ifdef SUBFIND
extern struct Subfind_DensityOtherPropsEval_data_out
{
    MyOutputFloat M200, R200;
#ifdef SUBFIND_ADDIO_VELDISP
    MyOutputFloat V200[3], Disp200;
#endif
#ifdef SUBFIND_ADDIO_BARYONS
    MyOutputFloat gas_mass, star_mass, temp, xlum;
#endif
}
*Subfind_DensityOtherPropsEval_DataResult, *Subfind_DensityOtherPropsEval_DataOut, *Subfind_DensityOtherPropsEval_GlobalPasser;
#endif



/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
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

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartGas;		/*!< This gives the maxmimum number of gas cells that can be stored on one processor. */
  int ICFormat;			/*!< selects different versions of IC file-format */
  int SnapFormat;		/*!< selects different versions of snapshot file-formats */
  int NumFilesPerSnapshot;	/*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/*!< maximum number of files that may be written simultaneously when writing/reading restart-files, or when writing snapshot files */
  double BufferSize;		/*!< size of communication buffer in MB */
  long BunchSize;     	        /*!< number of particles fitting into the buffer in the parallel tree algorithm  */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually NOT be balanced.  Each processor allocates memory for PartAllocFactor times the average number of particles to allow for that */
  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times the maximum(!) number of particles.  Note: A typical local tree for N particles needs usually about ~0.65*N nodes. */
  double TopNodeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times the maximum(!) number of particles.  Note: A typical local tree for N particles needs usually about ~0.65*N nodes. */

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
  long long NumForcesSinceLastDomainDecomp;	/*!< force updates since the last whole-tree BUILD (reset there -- historically at the decomposition, but the cadence it drives now requests a rebuild, so age-since-build is the meaningful quantity; name kept for restartfile layout). Past TreeDomainUpdateFrequency*TotNumPart the ladder condemns the tree. */

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
#ifdef RANDOMIZE_GRAVTREE_PERIODIC
  double RandomShift[3];		/*!< current random frame translation: coordinates are stored shifted by this (mod box) and un-shifted on output */
#endif

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

  integertime Ti_nextlineofsight;
#ifdef OUTPUT_LINEOFSIGHT
  double TimeFirstLineOfSight;
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
  int MaxMemSize;
  double CourantFac;		/*!< Courant factor */

  /* frequency of tree reconstruction/domain decomposition */
  double TreeDomainUpdateFrequency;	/*!< controls frequency of domain decompositions  */
  double DomainFacWork;     /*!< adaptive load-balance weights: moving averages over the decomposition */
  double DomainFacWorkGas;  /*!<   history, kept in All so they survive a restart (resetting them shifts */
  double DomainFacLoad;     /*!<   the next domain split) */
#ifdef MHD_MODIFIED_GRADIENT
  double ActiveFractionForMGSweep;  /*!< minimum active gas fraction to trigger the global MG div(B) solve; on smaller timesteps the local CG correction is used instead */
  int Flag_SkipMGSolve;             /*!< per-timestep flag: 1 = skip MG global solve this step (use CG fallback), 0 = run MG */
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
#ifdef TURB_DRIVING_TARGET_SIGMA
    double TurbDriving_SmoothedVRMS;     /* exponentially smoothed measured RMS velocity, time constant = 2*t_corr */
    double TurbDriving_DrivingNormFactor; /* multiplicative factor on driving amplitude from proportional control loop */
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

#ifdef KETJU_REGULARIZATION
  double KetjuRegionRadius;              /*!< physical radius of KETJU regions [code units] */
  double KetjuMinBHMass;                 /*!< minimum Type 5 stage-7 (BH/relic) mass to act as always-center [code units], 0=off */
  double KetjuMinStarMass;               /*!< minimum Type 5 stage<7 (star) mass to be a chain candidate [code units], 0=off */
  double KetjuIntegrationTolerance;      /*!< GBS relative tolerance (default 1e-9) */
  char   KetjuPNTerms[20];              /*!< PN term flags: "all","none","no_spin", or bitstring */
  int    KetjuEnableBHMergerKicks;       /*!< enable GW recoil kicks for BH mergers */
  int    KetjuUseStarStarSoftening;      /*!< keep softening between stars inside chain */
  double KetjuTimestepLimitingFactor;    /*!< radius factor for timestep limiting (default 100) */
  int    KetjuMaxStepCount;              /*!< max MSTAR integration steps before bail-out (default 100000, 0=off) */
  double KetjuExpandBinariesFactor;      /*!< expand tight binaries if period < factor*dt (default 0=off, typical 0.1) */
  double KetjuMinParticlesPerTask;       /*!< minimum particles per compute task for load balancing (default 10) */
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
    int ISMDustChem_SilicateMetallicityFieldIndexTable[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES]; /* index in metallicity field for elements which make up silicate dust (O, Mg, Si, and possibly Fe) */
    double ISMDustChem_SilicateNumberOfAtomsTable[GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES]; /* number of O, Mg, Si, and possibly Fe in one formula unit of silicate dust */
    double ISMDustChem_EffectiveSilicateDustAtomicWeight; /* atomic weight of one formula unit of silicate dust, depends on which optional module you use */
    // Scaling arguements from parameter file used to adjust each dust process
    double ISMDustChem_SNeIIDustScaling;
    double ISMDustChem_SNeIaDustScaling;
    double ISMDustChem_AGBDustScaling;
    double ISMDustChem_DustAccretionScaling;
    double ISMDustChem_ThermalSputteringScaling;
    double ISMDustChem_AccretionTcutoffScaling;
    double ISMDustChem_SNeGasClearedOfDustScaling;
#if (GALSF_ISMDUSTCHEM_MODEL & 2)
    double ISMDustChem_SpeciesBulkDens[4]; /* condensed bulk density for silicates, carbonaceous, SiC, and metallic iron */
    int ISMDustChem_TrackedSpeciesIDTable[NUM_ISMDUSTCHEM_SPECIES]; /* contains unique ID numbers for each tracked dust species which correspond to their location in ISMDustChem_SpeciesFieldIndexTable. Returns -1 for untracked species  */
    int ISMDustChem_SpeciesFieldIndexTable[NUM_ISMDUSTCHEM_SPECIES_IDS]; /* index in dust species field for given dust species. Sparse table indexed by fixed species ID (0..NUM_ISMDUSTCHEM_SPECIES_IDS-1), returns the packed field slot or -1 if untracked. */
    int ISMDustChem_Sil_Index;
    int ISMDustChem_Carb_Index;
    int ISMDustChem_SiC_Index;
    int ISMDustChem_FreeIron_Index ;
    int ISMDustChem_ORes_Index;
    int ISMDustChem_InclIron_Index;
#endif
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
  double Tillotson_EOS_params[7][12]; /*! < holds parameters for Tillotson EOS for solids */
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
}
All;



#include "particle_data.h"
#include "cell_data.h"

static_assert(std::is_trivially_copyable<particle_data>::value, "particle_data must be trivially copyable for MPI");
static_assert(std::is_trivially_copyable<gas_cell_data>::value, "gas_cell_data must be trivially copyable for MPI");




/* global state of system for checking conservation, etc. */
extern struct state_of_system
{
  double Mass,
    EnergyKin,
    EnergyPot,
    EnergyInt,
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[6],
    EnergyKinComp[6],
    EnergyPotComp[6],
    EnergyIntComp[6],
    EnergyTotComp[6],
    MomentumComp[6][4],
    AngMomentumComp[6][4],
    CenterOfMassComp[6][4];
}
SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation. */

extern struct data_index
{
  int Task;
  int Index;
  int IndexGet;
}
 *DataIndexTable;		/*!< the particles to be exported are grouped by task-number. This table allows the results to be disentangled again and to be assigned to the correct particle */


extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
*DataNodeList;


extern struct gravdata_in
{
    int Type;
    Vec3<MyFloat> Pos;
    MyFloat Soft;
    MyFloat Mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    Vec3<MyFloat> Vel;
#endif
#if defined(SINK_DYNFRICTION_FROMTREE)
    MyFloat Sink_Mass;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    MyFloat AGS_zeta;
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
    MyFloat Min_Sink_OrbitalTime;   /*!<orbital time for binary */
    Vec3<MyDouble> comp_dx;        /*!< position of binary companion */
    Vec3<MyDouble> comp_dv;        /*!< velocity of binary companion */
    MyDouble comp_Mass;         /*!< mass of binary companion */
    int is_in_a_binary;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    SymmetricTensor2<MyFloat> tidal_tensorps_prevstep;
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    int SuperTimestepFlag;  /*!< 2 if allowed to super-timestep, 1 if a candidate for super-timestepping, 0 otherwise */
#endif
#ifdef KETJU_REGULARIZATION
    MyIDType KetjuChainID;  /*!< shared chain ID for KETJU members — tree walk uses this to skip same-chain Type-5 neighbors in Min_Sink_Approach_Time/Freefall_time */
#endif
    MyFloat OldAcc;
    int NodeList[NODELISTLENGTH];
}
 *GravDataIn,			/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


extern struct gravdata_out
{
    Vec3<MyDouble> Acc;
#ifdef RT_USE_TREECOL_FOR_NH
    MyDouble ColumnDensityBins[RT_USE_TREECOL_FOR_NH];
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    MyDouble TreeMass;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyDouble SubGrid_CosmicRayEnergyDensity;
#endif
#ifdef RT_OTVET
    SymmetricTensor2<MyDouble> ET[N_RT_FREQ_BINS];
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    MyDouble Rad_Flux_UV;
    MyDouble Rad_Flux_EUV;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    MyDouble Rad_E_gamma[N_RT_FREQ_BINS];
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    Vec3<MyDouble> Rad_Flux[N_RT_FREQ_BINS];
#endif
#ifdef CHIMES_STELLAR_FLUXES
    double Chimes_G0[CHIMES_LOCAL_UV_NBINS];
    double Chimes_fluxPhotIon[CHIMES_LOCAL_UV_NBINS];
#endif
#ifdef SINK_COMPTON_HEATING
    MyDouble Rad_Flux_AGN;
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    MyDouble MencInRcrit;
#endif
#ifdef EVALPOTENTIAL
    MyDouble Potential;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    SymmetricTensor2<MyFloat> tidal_tensorps;
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    MyFloat tidal_zeta;
#endif
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    Vec3<MyDouble> GravJerk;
#endif
#ifdef SINK_CALC_DISTANCES
    MyFloat Min_Distance_to_Sink;
    Vec3<MyFloat> Min_xyz_to_Sink;
#ifdef SPECIAL_POINT_MOTION
    Vec3<MyFloat> vel_of_nearest_special;
    Vec3<MyFloat> acc_of_nearest_special;
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    MyFloat weight_sum_for_special_point_smoothing;
#endif
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
    MyFloat Min_Sink_OrbitalTime; //orbital time for binary
    Vec3<MyDouble> comp_dx; //position of binary companion
    Vec3<MyDouble> comp_dv; //velocity of binary companion
    MyDouble comp_Mass; //mass of binary companion
    int is_in_a_binary; // 1 if star is in a binary, 0 otherwise
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    MyFloat Min_Sink_Freefall_time;    // minimum value of sqrt(R^3 / G(M_SINK + M_particle)) as calculated from the tree-walk
    MyFloat Min_Sink_Approach_Time; // smallest approach time t_a = |v_radial|/r
#if (SINGLE_STAR_TIMESTEPPING > 0)
    SymmetricTensor2<MyFloat> COM_tidal_tensorps; //tidal tensor evaluated at the center of mass without contribution from the companion
    MyDouble COM_GravAccel[3]; //gravitational acceleration evaluated at the center of mass without contribution from the companion
    int COM_calc_flag; //flag that tells whether this was only a rerun to get the acceleration ad the tidal tenor at the center of mass of a binary
    int SuperTimestepFlag; // 2 if allowed to super-timestep, 1 if a candidate for super-timestepping, 0 otherwise
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    MyFloat Min_Sink_FeedbackTime; // minimum time for feedback to arrive from a star
#endif    
#endif
#endif
}
 *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct potdata_out
{
  MyDouble Potential;
}
 *PotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct info_block
{
  char label[4];
  char type[8];
  int ndim;
  int is_present[6];
}
*InfoBlock;



/* this structure needs to be defined here, because routines for feedback event rates, etc, are shared among files */
extern struct addFB_evaluate_data_in_
{
    Vec3<MyDouble> Pos, Vel; MyDouble Msne, unit_mom_SNe;
    MyFloat KernelRadius, V_i, SNe_v_ejecta;
#ifdef GALSF_FB_MECHANICAL
    MyFloat Area_weighted_sum[AREA_WEIGHTED_SUM_ELEMENTS];
#endif
#ifdef METALS
    MyDouble yields[NUM_METAL_SPECIES+NUM_ADDITIONAL_PASSIVESCALAR_SPECIES_FOR_YIELDS_AND_DIFFUSION];
#endif
    int NodeList[NODELISTLENGTH];
}
*addFB_evaluate_DataIn_, *addFB_evaluate_DataGet_;





/* Header for the standard file format */
extern struct io_header
{
  int npart[6];			    /*!< number of particles of each type in this file */
  double mass[6];           /*!< mass of particles of each type. If 0, then the masses are explicitly stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;			    /*!< time of snapshot file */
  double redshift;		    /*!< redshift of snapshot file */
  int flag_sfr;			    /*!< flags whether the simulation was including star formation */
  int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];   /*!< total number of particles of each type in this snapshot. This can be different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		    /*!< flags whether cooling was included  */
  int num_files;		    /*!< number of files in multi-file snapshot */
  double BoxSize;		    /*!< box-size of simulation in case periodic boundaries were used */
  double OmegaMatter;       /*!< matter density in units of critical density */
  double OmegaLambda;		/*!< cosmological constant parameter */
  double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
  int flag_metals;		    /*!< flags whether the file contains metallicity values for gas and star particles */

  unsigned int npartTotalHighWord[6];   /*!< High word of the total number of particles of each type (needed to combine with npartTotal to allow >2^31 particles of a given type) */
  int flag_entropy_instead_u; /*!< flag here strictly for historical compatibility with unformatted binary files from GADGET-3 era formats, which expect this flag to exist. this does nothing in gizmo */
  int flag_doubleprecision; /*!< flags that snapshot contains double-precision instead of single precision */

  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */

  char fill[18];		        /*!< fills to 256 Bytes */
  char names[15][2];
}
header;				/*!< holds header for snapshot files */







enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_CHILD_ID,
  IO_GENERATION_ID,
  IO_MASS,
  IO_U,
  IO_RHO,
  IO_NE,
  IO_NH,
  IO_KERNELRADIUS,
  IO_SFR,
  IO_AGE,
  IO_GRAINSIZE,
  IO_DUST_TO_GAS, 
  IO_GRAINTYPE,
  IO_HSMS,
  IO_Z,
  IO_DUSTCHEMZMET,
  IO_DUSTCHEMSPECIESMET,
  IO_ISMDUSTCHEMMOL,
  IO_MACHNUM,
  IO_DUSTCHEMGRAINBINNUMBERS,
  IO_DUSTCHEMGRAINBINSLOPES,
  IO_DUSTCHEMGRAINBINMASS,
  IO_DUSTCHEM_COAG_MASSRATE,
  IO_DUSTCHEM_SHAT_MASSRATE,
  IO_SINKMASS,
  IO_SINKMASSALPHA,
  IO_SINK_ANGMOM,
  IO_SINKMDOT,
  IO_SINKDUSTMASSACC,
  IO_R_PROTOSTAR,
  IO_MASS_D_PROTOSTAR,
  IO_ZAMS_MASS,
  IO_STAGE_PROTOSTAR,
  IO_AGE_PROTOSTAR,
  IO_LUM_SINGLESTAR,
  IO_SINKPROGS,
  IO_SINK_DIST,
  IO_ACRB,
  IO_SINKRAD,
  IO_SINK_FORM_MASS,
  IO_POT,
  IO_ACCEL,
  IO_HYDROACCEL,
  IO_HERMITE_POS,
  IO_HERMITE_VEL,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_UNSPMASS,
  IO_UNSPJETMASS,
  IO_CRATE,
  IO_HRATE,
  IO_NHRATE,
  IO_HHRATE,
  IO_MCRATE,
  IO_DCRATE,
  IO_PHRATE,
  IO_DTENTR,
  IO_TSTP,
  IO_BFLD,
  IO_AMBIPOLAR,
  IO_OHMIC,
  IO_HALL,
  IO_IMF,
  IO_COSMICRAY_ENERGY,
  IO_COSMICRAY_KAPPA,
  IO_COSMICRAY_ALFVEN,
  IO_COSMICRAY_SLOPES,
  IO_DIVB,
  IO_ABVC,
  IO_AMDC,
  IO_PHI,
  IO_GRADPHI,
  IO_GRADRHO,
  IO_GRADVEL,
  IO_GRADMAG,
  IO_COOLRATE,
  IO_TIDALTENSORPS,
  IO_EOSTEMP,
  IO_EOSABAR,
  IO_EOSYE,
  IO_PRESSURE,
  IO_EOSCS,
  IO_EOS_STRESS_TENSOR,
  IO_CBE_MOMENTS,
  IO_EOSCOMP,
  IO_EOSPHASE,
  IO_PARTVEL,
  IO_RADGAMMA,
  IO_RAD_TEMP,
  IO_RAD_OPACITY,
  IO_DUST_TEMP,
  IO_RAD_ACCEL,
  IO_RAD_FLUX,
  IO_EDDINGTON_TENSOR,
  IO_VDIV,
  IO_VORT,
  IO_DELAYTIME,
  IO_SOFT,
  IO_AGS_HKERN,
  IO_AGS_RHO,
  IO_AGS_QPT,
  IO_AGS_PSI_RE,
  IO_AGS_PSI_IM,
  IO_AGS_ZETA,
  IO_VSTURB_DISS,
  IO_VSTURB_DRIVE,
  IO_grHI,
  IO_grHII,
  IO_grHM,
  IO_grHeI,
  IO_grHeII,
  IO_grHeIII,
  IO_grH2I,
  IO_grH2II,
  IO_grDI,
  IO_grDII,
  IO_grHDI,
  IO_OSTAR,
  IO_DTOSTAR,
  IO_TURB_DYNAMIC_COEFF,
  IO_TURB_DIFF_COEFF,
  IO_DYNERROR,
  IO_DYNERRORDEFAULT,
  IO_CHIMES_ABUNDANCES,
  IO_CHIMES_MU,
  IO_CHIMES_REDUCED,
  IO_CHIMES_NH,
  IO_CHIMES_STAR_SIGMA,
  IO_CHIMES_FLUX_G0,
  IO_CHIMES_FLUX_ION,
  IO_DENS_AROUND_STAR,
  IO_DELAY_TIME_HII,
  IO_MOLECULARFRACTION,
  IO_SHOCKMACHNUM,
  IO_KETJU_FINAL_VEL,
  IO_KETJU_SPIN,
  IO_LASTENTRY			/* This should be kept - it signals the end of the list */
};


enum siofields
{ SIO_GLEN,
  SIO_GOFF,
  SIO_MTOT,
  SIO_GPOS,
  SIO_DELTA_MSUB,
  SIO_DELTA_RSUB,
  SIO_DELTA_DISPSUB,
  SIO_DELTA_MGASSUB,
  SIO_DELTA_MSTSUB,
  SIO_DELTA_TEMPSUB,
  SIO_DELTA_LXSUB,
  SIO_NCON,
  SIO_MCON,
  SIO_BGPOS,
  SIO_BGMTOP,
  SIO_BGRTOP,
  SIO_NSUB,
  SIO_FSUB,
  SIO_SLEN,
  SIO_SOFF,
  SIO_PFOF,
  SIO_MSUB,
  SIO_SPOS,
  SIO_SVEL,
  SIO_SCM,
  SIO_SPIN,
  SIO_DSUB,
  SIO_VMAX,
  SIO_RVMAX,
  SIO_RHMS,
  SIO_MBID,
  SIO_GRNR,
  SIO_SMST,
  SIO_SLUM,
  SIO_SLATT,
  SIO_SLOBS,
  SIO_HALODUST,
  SIO_SAGE,
  SIO_SZ,
  SIO_SSFR,
  SIO_PPOS,
  SIO_PVEL,
  SIO_PTYP,
  SIO_PMAS,
  SIO_PID,

  SIO_LASTENTRY
};


/* Variables for Tree */
/* Tree Node structure (note, the ALIGN(32) directive will effectively pad the structure size to a multiple of 32 bytes) */
extern ALIGN(32) struct NODE
{
  Vec3<MyFloat> center;		/*!< geometrical center of node */
  MyFloat len;			/*!< sidelength of treenode */

  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      Vec3<MyFloat> s;		/*!< center of mass of node */
      MyFloat mass;		/*!< mass of node */
      unsigned int bitflags;	/*!< flags certain node properties */
      int sibling;		/*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;		/*!< this gives the next node in case the current node needs to be opened */
      int father;		/*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;

    double GravCost;
    integertime Ti_current;
    long N_part;   /*!< number of particles+cells in the tree node */
    MyFloat maxsoft;        /*!< hold the maximum gravitational softening of particle in the node */
#if defined(GRAVTREE_CALCULATE_GAS_MASS_IN_NODE)
  MyFloat gasmass;
#endif
#ifdef RT_USE_GRAVTREE
  MyFloat stellar_lum[N_RT_FREQ_BINS]; /*!< luminosity in the node*/
#ifdef CHIMES_STELLAR_FLUXES
  double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS];
  double chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#endif

#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<MyFloat> rt_source_lum_s; /*!< luminosity-weighted position of RT sources in the node */
#endif
#ifdef SINK_PHOTONMOMENTUM
    MyFloat sink_lum;		    /*!< luminosity of BHs in the node */
    Vec3<MyFloat> sink_lum_grad;	/*!< gradient vector for gas around sink (for angular dependence) */
#endif
    
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyFloat cr_injection;
#endif

#ifdef SINK_CALC_DISTANCES
  MyFloat sink_mass;      /*!< holds the sink mass in the node.  Used for calculating tree based dist to closest sink */
  Vec3<MyFloat> sink_pos;    /*!< holds the mass-weighted position of the the actual sink particles within the node */
#ifdef SINK_NODE_MOTION_TRACKED
    Vec3<MyFloat> sink_vel;    /*!< holds the mass-weighted avg. velocity of sink particles in the node */
#endif
#if defined(SPECIAL_POINT_MOTION)
    Vec3<MyFloat> sink_acc; /*!< holds the mass-weighted avg. acceleration of sink particles in the node */
#endif
#ifdef SINK_NODE_MOTION_TRACKED
  int N_SINK;             /*!< holds the number of sink particles in the node. Used for refinement/search criteria */
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    MyFloat MaxFeedbackVel;
#endif
#endif


#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    SymmetricTensor2<MyFloat> tidal_tensorps_prevstep;
#endif
#ifdef DM_SCALARFIELD_SCREENING
  Vec3<MyFloat> s_dm;
  MyFloat mass_dm;
#endif
}
 *Nodes_base,			/*!< points to the actual memory allocated for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] gives the first allocated node */


#ifdef SINGLE_STAR_DIRECT_GRAVITY
/* every star in the simulation, replicated on every task, for the exact star-star sum in
   gravity/star_direct_gravity.cc. Rebuilt each time gravity is evaluated; sized O(N_star). */
extern struct star_direct_data
{
    Vec3<MyDouble> Pos;
    Vec3<MyFloat> Vel;
    MyFloat Mass;
    MyFloat Soft;
    MyIDType ID;
}
 *StarDirect;
extern int N_StarDirect;
#endif


extern struct extNODE
{
  Vec3<MyDouble> dp;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<MyDouble> rt_source_lum_dp;
    Vec3<MyFloat> rt_source_lum_vs;
#endif
#ifdef DM_SCALARFIELD_SCREENING
  Vec3<MyDouble> dp_dm;
  Vec3<MyFloat> vs_dm;
#endif
#ifdef SINK_NODE_MOTION_TRACKED
  Vec3<MyDouble> sink_dp;   /*!< momentum change of the sinks in this node, applied to sink_vel on the next drift (mirrors dp/dp_dm) */
#endif
  Vec3<MyFloat> vs;
  MyFloat vmax;
  MyFloat hmax;			/*!< maximum gas kernel length in node. Only used for gas particles */
  MyFloat divVmax;
  integertime Ti_lastkicked;
  int Flag;
}
 *Extnodes, *Extnodes_base;



#endif  /* ALLVARS_H  - please do not put anything below this line */
