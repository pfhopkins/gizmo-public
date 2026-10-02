
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
#include "gizmo_rng.h"
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

/* Pull in the struct definition + extern declaration of `All` early so the
   boxSize / boxHalf / Shearing_Box_*_Offset macros below (which expand to
   All.*) parse correctly inside templates and inline functions defined in
   this header.  The full repeat of the include + extern at the bottom of
   the file is gated on `#ifndef All` and is harmless on second pass. */
#include "global_data_all_struct.h"
#ifndef All
extern struct global_data_all_processes All;
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

/* Box-size shorthands.  Always macros reading from All so the same code
   compiles for host + device with no extern globals to keep in sync.
   In GPU TUs the device compilation pass remaps `All` to that TU's
   private `AllDeviceMirror` via gpu_all_mirror.h, so these macros work
   in both contexts unchanged. All.BoxSize is always a valid field (set
   from the param file even for non-periodic runs), so boxSize/boxHalf
   are unconditionally defined. */
#define boxSize (All.BoxSize)
#define boxHalf (0.5*All.BoxSize)
#ifdef BOX_LONG_X
#define boxSize_X (All.BoxSize * ((MyDouble)(BOX_LONG_X)))
#define boxHalf_X (0.5 * All.BoxSize * ((MyDouble)(BOX_LONG_X)))
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef BOX_LONG_Y
#define boxSize_Y (All.BoxSize * ((MyDouble)(BOX_LONG_Y)))
#define boxHalf_Y (0.5 * All.BoxSize * ((MyDouble)(BOX_LONG_Y)))
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef BOX_LONG_Z
#define boxSize_Z (All.BoxSize * ((MyDouble)(BOX_LONG_Z)))
#define boxHalf_Z (0.5 * All.BoxSize * ((MyDouble)(BOX_LONG_Z)))
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif

#ifdef BOX_SHEARING
#define Shearing_Box_Vel_Offset (All.Shearing_Box_Vel_Offset)
#define Shearing_Box_Pos_Offset (All.Shearing_Box_Pos_Offset)
#endif

#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z)
/* These live in All so device code can read them: the boundary check runs inside
   the drift, which is moving to the GPU, and a free-standing global is not visible
   there. Aliased so every existing reader and the begrun initialization are
   unchanged -- same treatment as Shearing_Box_*_Offset above. Use the bare name:
   spelling out All.special_boundary_condition_xyz_def_reflect double-expands and
   will not compile, exactly as for Shearing_Box_*_Offset. */
#define special_boundary_condition_xyz_def_reflect (All.special_boundary_condition_xyz_def_reflect)
#define special_boundary_condition_xyz_def_outflow (All.special_boundary_condition_xyz_def_outflow)
#endif

template<typename T> GIZMO_GPU_FUNCTION inline void nearest_xyz(Vec3<T>& v, int sign=1) { NEAREST_XYZ(v[0], v[1], v[2], sign); }


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
extern int NumForceUpdateAtSyncPoint;	/*!< NumForceUpdate as fixed at the sync point, before the physics of the step can add to it (star formation, sink swallows). Routing decisions taken before the active list is rebuilt must use this rather than the live counter, which several later call sites mutate. */
extern long long GlobNumForceUpdate;
extern int NumGasUpdate;	/*!< number of active gas cells on local processor in current timestep  */
extern int MaxTopNodes;	        /*!< Maximum number of nodes in the top-level tree used for domain decomposition */
extern int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from initial conditions, 1 is resuming a run from a set of restart files, while 2 marks a restart from a snapshot file. */
extern int RestartSnapNum;
extern int SelRnd;
extern int TakeLevel;
extern int *Exportflag;	        /*!< per-task flag used by the gravity LET-incompleteness detector (the export round-trip is retired) */
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
extern double CPU_ChildCharged;         /*!< running total charged to child CPU_Step buckets inside enclosing timed spans */
extern double CPU_ChildCharged_at_sync; /*!< value of the above at the last chain sync point */
extern int Flag_FullStep;	/*!< Flag used to signal that the current step involves all particles */
extern size_t HighMark_run,  HighMark_domain, HighMark_gravtree, HighMark_pmperiodic,
  HighMark_pmnonperiodic,  HighMark_gasdensity, HighMark_hydro, HighMark_GasGrad;

#ifdef TURB_DRIVING
extern size_t HighMark_turbpower;
#endif
extern int TreeReconstructFlag;
extern int DomainReconstructFlag;   /*!< set when the domain decomposition itself must be redone, not merely the tree rebuilt */
extern int TreeMomentsStaleFlag; /*!< flag to refresh tree node moments without a full tree rebuild, e.g. after star formation or sink mass changes */
extern long long ForceAddElementToTree_CallsSinceBuild; /*!< diagnostic: force_add_element_to_tree calls accumulated since last full tree build.  Insertions stale the LET / pseudo-particle moments; auto-rebuild when this exceeds 1% of TotNumPart. */
extern int GlobFlag;
extern char DumpFlag;
extern int NeedToWakeupParticles;
extern int NeedToWakeupParticles_local;
extern unsigned char *WakeupDirty;   /*!< process_wake_ups superset index (UVM, MaxPart): when WakeupDirtyValid, set for every i<NumPart with P[i].wakeup!=0 (false positives allowed). P[i].wakeup stays authoritative. */
extern int WakeupDirtyValid;         /*!< 0 => WakeupDirty stale (index remap / init); next process_wake_ups rebuilds from P[]. */

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
#ifdef SINK_WIND_SPAWN
extern double Max_Unspawned_MassUnits_fromSink;
#endif

extern long long Ntype[6];	/*!< total number of particles of each type */
extern int NtypeLocal[6];	/*!< local number of particles of each type */
extern gizmo_rng_t random_generator;	/*!< the random number generator used */
extern int Gas_split;           /*!< current number of newly-spawned gas particles outside block */
#ifdef GALSF
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
#endif
#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
extern int Grains_promoted;	/*!< current number of grain particles promoted to solid body in gas block */
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
/* Top-leaf-router geometry SSOT: TopNodeNodeIndex[topnode] -> Nodes[] slot for
 * EVERY topnode (internal + leaf), populated in force_create_empty_nodes.  Lets
 * the hierarchical router read each topnode's exact center/len from
 * Nodes[TopNodeNodeIndex[no]] instead of (wrongly) deriving it from the
 * Peano-Hilbert child offset.  For leaves, TopNodeNodeIndex[leaf] ==
 * DomainNodeIndex[TopNodes[leaf].Leaf] by construction (asserted). */
extern int *TopNodeNodeIndex;
extern int *DomainList, DomainNumChanged;
extern peanokey *Key, *KeySorted;


#ifdef HERMITE_INTEGRATION
extern int HermiteOnlyFlag;     /*!< flag to only do Hermite integration for applicable particles (ie. stars) in the gravity routine - set =1 on the first prediction pass and =2 on the second correction pass */
#endif

/* RT_CHEM_PHOTOION arrays and CRFLUID_EVOLVE_SPECTRUM arrays are now in the All struct
   (global_data_all_struct.h) for GPU visibility via __managed__ All. */


extern int MaxNodes;        /*!< maximum allowed number of internal nodes */
extern int Numnodestree;    /*!< number of (internal) nodes in each tree */
extern int MaxForeignNodes;  /*!< LET: the CEILING of the foreign-node INDEX range, published by force_treeallocate -- derived there for a fresh tree, or restored verbatim from the restart file for a serialized one, whose node pointers encode the ceiling the writer used.  0 before the first tree is allocated.  Index map: foreign nodes occupy [TreeNodeIndexBase+MaxNodes, TreeNodeIndexBase+MaxNodes+MaxForeignNodes); pseudo-particles shifted to [TreeNodeIndexBase+MaxNodes+MaxForeignNodes, ...).  This is an index ceiling, NOT the amount of foreign-node storage that exists: see AllocatedForeignNodes. */
extern int AllocatedForeignNodes;  /*!< LET: how many foreign-node slots actually have STORAGE behind them (node/extnode arena extent, foreign-leaf sidecars, GPU node mirror).  Zero when the tree is allocated; raised once per tree build to the exact import this rank received, by force_tree_grow_foreign_storage.  Per-rank, unlike the ceiling.  Invariant: Numforeignnodes <= AllocatedForeignNodes <= MaxForeignNodes. */
extern int Numforeignnodes;  /*!< LET: count of foreign nodes currently installed (<= AllocatedForeignNodes); reset on each LET exchange. */

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
extern gizmo_rng_t StRng; // random number generator key
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
#if (defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(CBE_INTEGRATOR_OUTPUT_MOREINFO)) && defined(CBE_INTEGRATOR)
extern FILE *FdCbeDiagnostics;
#endif
#ifdef CBE_INTEGRATOR
/* Per-particle-type presence flag for the VlasovMoments
 * IC dataset. Set in read_ic.cc inside the `if(hdf5_dataset >= 0)` block
 * on successful H5Dread; consumed by do_cbe_initialization in sidm/
 * cbe_integrator.cc to choose between "trust loaded moments" vs
 * "synthesize cold default" per particle. Per-type (NOT a single global)
 * because VlasovMoments lives under /PartTypeX/ in the HDF5 layout — a
 * mixed cosmological IC can carry it for DM (Type=1) and not for gas
 * (Type=0) or stars (Type=4). Single-flag would corrupt absent-type
 * particles. Reset to {0,0,0,0,0,0} at the top of read_ic() each call. */
extern int CBE_Moments_LoadedFromIC_PType[6];
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
extern double DriftTable_logTimeBegin; /*! log of the scale factor the two tables above start at */
extern double DriftTable_logTimeMax;   /*! log of the scale factor they end at */
extern void *CommBuffer;    /*!< points to communication buffer, which is used at a few places */


#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
#endif

extern peanokey *DomainKeyBuf;

/* Tree variables */
extern long Nexport;   /* gravity LET-incompleteness detector count (Nimport retired with the export round-trip) */
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



/* Struct definition lives in its own minimal header so GPU-only TUs can include
   it without pulling in MPI, GSL, or HDF5.  allvars.h includes it here so
   all existing code that includes allvars.h sees the type as before. */
#include "global_data_all_struct.h"

/* In a GPU TU during the device compilation pass, gpu_all_mirror.h defines
 * `All` to that TU's private `AllDeviceMirror`. The #ifndef guard keeps
 * this re-declaration of the host extern from being macro-expanded in
 * that pass; in the host pass and in non-GPU TUs the extern is declared
 * normally. */
#ifndef All
extern struct global_data_all_processes All;
#endif

/* SSOT for whether the dm_dispersion neighbor loop is compiled AND called.
 * Resolved AFTER global_data_all_struct.h sets GALSF_SUBGRID_WIND_SCALING
 * and BEFORE particle_data.h / cell_data.h / proto.h consume it. Future
 * physics needing rho_DM / sigma_DM at gas cells: add to this OR. */
#if (defined(GALSF_SUBGRID_WINDS) && (GALSF_SUBGRID_WIND_SCALING==2)) || defined(DM_HEATING)
#define DM_DISPERSION_LOOP_ACTIVE
#endif

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
 *DataIndexTable;		/*!< records would-be-exported targets per task for the gravity LET-incompleteness detector; never shipped (the export round-trip is retired) */


extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
*DataNodeList;


/* The gravdata_in/gravdata_out export+import buffers (GravDataIn/Get/Result/Out)
 * are RETIRED: gravity runs on the target-owning rank via the GPU pre-pass + LET,
 * so there is no MPI export round-trip. The surviving export-list structs
 * (data_index/data_nodelist) now serve only as the LET-incompleteness detector. */


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
    MyDouble yields[NUM_METAL_SPECIES];
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
  IO_REFINE_FLAG,
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
  IO_DUSTCHEMGRAINBINMASS,
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
  IO_HERMITE_POS,
  IO_HERMITE_VEL,
  IO_ACCEL,
  IO_HYDROACCEL,
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
#ifdef TWO_TEMPERATURE_PLASMA
  IO_TWOTEMP_TE,
#endif
  IO_EOSABAR,
  IO_EOSYE,
#ifdef NUCLEAR_NETWORK
  IO_NUCLEAR_COMPOSITION,
#endif
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
  IO_DAMAGE_POROSITY_DAMAGE,      /* Grady-Kipp scalar D in [0,1] */
  IO_DAMAGE_POROSITY_DISTENTION,  /* Jutzi P-alpha alpha in [1,alpha_0] */
  IO_DAMAGE_POROSITY_ACTVCRACKS,  /* Grady-Kipp active-crack radius */
  IO_FLUIDTYPE,
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
#if defined(SINK_NODE_MOTION_TRACKED)
    Vec3<MyFloat> sink_vel;    /*!< holds the mass-weighted avg. velocity of sink particles in the node. Gate matches N_SINK and every sink_vel reader/writer. */
#endif
#if defined(SPECIAL_POINT_MOTION)
    Vec3<MyFloat> sink_acc; /*!< holds the mass-weighted avg. acceleration of sink particles in the node */
#endif
#if defined(SINK_NODE_MOTION_TRACKED)
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
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.TreeNodeIndexBase] gives the first allocated node */


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
  Vec3<MyDouble> sink_dp;   /*!< momentum given to this node's special-type particles since the last
                                 drift, folded into Nodes[].sink_vel there. Mirrors dp/dp_dm; the
                                 velocity itself lives in Nodes, beside sink_pos. */
#endif
  Vec3<MyFloat> vs;
  MyFloat vmax;			/*!< how fast any member can move along an axis, per unit undilated drift interval
				   (particle_motion_speed_bound); the node's box grows by TREE_NODE_WIDENING_DELTA of
				   this and the interval since it was last drifted, so a member not yet drifted
				   is still inside it */
  MyFloat hmax;			/*!< maximum gas kernel length in node. Legacy aggregate:
				   gas-only at tree build, AGS-aware after force_update_hmax
				   in ADAPTIVE_GRAVSOFT_FORALL builds. Read by legacy
				   gravity / GPU SoA / LET pack code. */
  /* Mode B per-type kernel-radius bands.  Invariant (every type t):
   *   hmax_per_type[t] = max over Type==t in subtree of
   *                       force_hmax_per_type_particle_radius(j)
   *                    = max( DMIN(All.MaxKernelRadius,
   *                                max(P[j].KernelRadius,
   *                                    P[j].AGS_KernelRadius if AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)),
   *                           P[j].ForceSoftening )
   *
   * Kernel-like radii (KernelRadius / AGS_KernelRadius) are capped at
   * All.MaxKernelRadius per legacy invariant (those fields can blow up under
   * iterative search).  ForceSoftening is NOT capped — the leaf policy may
   * admit pairs via FS reach (sink_feed's binary-merge filter), so the node
   * band must dominate FS unconditionally.
   *
   * The band is a conservative UPPER BOUND across every leaf-policy-selectable
   * radius source for that type — so Mode B SYMMETRIC tree-prune can use
   * max_{t in mask} hmax_per_type[t] regardless of which Spec's radius_policy
   * is calling.  The exact policy filter sits at the leaf predicate
   * (mode_b_neighbor_symmetric_radius -> nlr_symmetric_radius_from_fields);
   * over-opening at the node level is safe (extra candidates filter at leaf),
   * under-opening would be a correctness bug.
   *
   * Drift behavior: per-type bands inflate UPWARD only during force_drift_node
   * (and gpu_force_drift's UVM equivalent).  Downward decay is suppressed
   * because the band includes static-ish sources (ForceSoftening) that don't
   * shrink under drift; force_update_hmax() re-seeds bands per particle each
   * call, and a transiently-large band only over-opens (safe).  Scalar `hmax`
   * keeps its legacy bidirectional decay — its semantics are unchanged.
   *
   * ONEWAY walks never read this field.
   *
   * These bands are seeded rank-local at every tree build/refresh AND maintained
   * cross-rank: DomainNODE carries them on the full-build/refresh pseudodata
   * exchange (up-propagated through the top tree by gpu_topnode_moment_resum), and
   * force_update_hmax() re-exchanges them post-density alongside the scalar hmax.
   * So a remote topleaf's (and its ancestors') per-type bands are populated and
   * kept as fresh as the scalar hmax.  Consumer = the Mode B SYMMETRIC local
   * tree-walk (mesh/mode_b_local_walker.cc): rank-local bands prune the candidate
   * walk (returns only rank-local candidates, never descends pseudo/foreign
   * subtrees); the cross-rank bands bound the sender's targeted export of remote
   * topleaves (mode_b_node_symmetric_radius), replacing the broadcast path.
   * (Historically these were rank-local + broadcast-answered; the cross-rank
   * exchange + export routing were added to bound the non-gas SYMMETRIC loops --
   * correctness for those, since the gas-biased scalar hmax cannot cover their
   * reach.) */
  MyFloat hmax_per_type[6];
  MyFloat divVmax;
  integertime Ti_lastkicked;
  int Flag;
}
 *Extnodes, *Extnodes_base;



#endif  /* ALLVARS_H  - please do not put anything below this line */
