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
 * Volker Springel. The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (many new variables,
 * structures, and different naming conventions for some old variables).
 * The declarations have also been divided into a new scheme which improves
 * read-ability and ease of use and separates all of the macros written for
 * GIZMO.
 */

#include "allvars.h"
#include "../core/timestep_functions.h"   /* HermiteWalkState / DriftKickTableView, for the pass-state globals below */




/*********************************************************/
/*  Global variables                                     */
/*********************************************************/



int ThisTask;			/*!< the number of the local processor  */
int NTask;			/*!< number of processors */
int PTask;			/*!< note: NTask = 2^PTask */

double CPUThisRun;		/*!< Sums CPU time of current process */

int NumForceUpdate;		/*!< number of active particles on local processor in current timestep  */
int NumForceUpdateAtSyncPoint;	/*!< see allvars.h: sync-point snapshot of NumForceUpdate for pre-rollover routing */
long long GlobNumForceUpdate;
int NumGasUpdate;		/*!< number of active gas cells on local processor in current timestep  */

int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from initial conditions, 1 is resuming a run from a set of restart files, while 2 marks a restart from a snapshot file. */

int RestartSnapNum;
int SelRnd;

int *Exportflag;		/*!< per-task flag used by the gravity LET-incompleteness detector (the export round-trip is retired) */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount;

int TakeLevel;

std::vector<int> ActiveParticleList;
unsigned char *ProcessedFlag;

int TimeBinCount[TIMEBINS];
int TimeBinCountGas[TIMEBINS];
int TimeBinActive[TIMEBINS];

int FirstInTimeBin[TIMEBINS];
int LastInTimeBin[TIMEBINS];
std::vector<int> NextInTimeBin;
std::vector<int> PrevInTimeBin;

size_t HighMark_run, HighMark_domain, HighMark_gravtree,
  HighMark_pmperiodic, HighMark_pmnonperiodic, HighMark_gasdensity, HighMark_hydro, HighMark_GasGrad;

/* boxSize / boxHalf / boxSize_[XYZ] / boxHalf_[XYZ] / Shearing_Box_*_Offset and
   special_boundary_condition_xyz_def_{reflect,outflow} are macros into All.* —
   see declarations/allvars.h. */

#ifdef TURB_DRIVING
size_t HighMark_turbpower;
#endif

#ifdef GALSF
double TimeBinSfr[TIMEBINS];
#endif

#ifdef SINK_PARTICLES
double TimeBin_Sink_mass[TIMEBINS];
double TimeBin_Sink_dynamicalmass[TIMEBINS];
double TimeBin_Sink_Mdot[TIMEBINS];
double TimeBin_Sink_Medd[TIMEBINS];
#endif

/* RT_CHEM_PHOTOION and CRFLUID_EVOLVE_SPECTRUM arrays moved to All struct (global_data_all_struct.h) */


char DumpFlag = 1;
size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;
double CPU_Step[CPU_PARTS];
/* Characters painted into the balance.txt work/imbalance map, one pair per
 * CPU_Step bucket. write_cpu_log() lays both shares of every bucket into a
 * single flat string, so a character must identify exactly one bucket AND one
 * of the two shares.  A character may be reused by a second bucket: the map is
 * painted in bucket order and each bucket's run is contiguous, so position
 * separates two uses that sit far apart in the list.  Reusing one between
 * neighbours is what makes a column unreadable. '#' (head marker) and '-'
 * (unaccounted tail) are painted literally by the writer and are reserved.
 *
 * Only the 43 buckets that have a charge site can ever be painted (the writer
 * skips any bucket with zero time), and 43 pairs plus the two reserved
 * characters is the most that fits in printable ASCII. Buckets with no writer
 * therefore carry the sentinel '_' rather than a plausible character that a
 * reader would hunt for and never find; balance.txt names the sentinel in its
 * legend. Charging one of those buckets means giving it a real character here,
 * in BOTH arrays, checked unique against every other entry in both. */
char CPU_Symbol[CPU_PARTS] = {
    '_', '*', '_', '_', '<', '_', '_', ':', '.', '~', '|', '+', '"', '_',  '`', ',', '_', '_', '_', '&',
    '$', '_', '(', '?', ')', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '\\', '%', '{', '}', 'Z',
    '_', '_', '_', 'a', '_', '_', 'd', 'e', 'q', ']', '_', '!', 'j', 'k', 'l', 'm',  '/', '@'};
char CPU_SymbolImbalance[CPU_PARTS] = {
    '_', 't', '_', '_', 'b', '_', '_', 'r', 'h', 'B', 'n', 'C', 'o', '_',  's', 'f', '_', '_', '_', 'D', // 20 columns here
    'x', '_', 'z', 'E', 'I', 'W', 'T', 'V', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N',  'O', 'P', 'Q', 'R',
    '_', '_', '_', 'S', '_', '_', 'U', 'X', 'Y', '\'', '_', 'y', 'c', 'g', 'i', 'p',  '>', '='};
char CPU_String[CPU_STRING_LEN + 1];
double WallclockTime;		/*!< This holds the last wallclock time measurement for timings measurements */
double CPU_ChildCharged;
double CPU_ChildCharged_at_sync;
int Flag_FullStep;		/*!< Flag used to signal that the current step involves all particles */


int TreeReconstructFlag;
int DomainReconstructFlag;
int DomainExtentOutgrownLocal;
int TreeMomentsStaleFlag;
int TypePresenceMaskTrusted = 0;
long long ForceAddElementToTree_CallsSinceBuild = 0;
int NeedToWakeupParticles;      /*!< Flags used to signal that wakeups need to be processed at the beginning of the next timestep */
int NeedToWakeupParticles_local;
unsigned char *WakeupDirty = NULL;
int WakeupDirtyValid = 0;
int GlobFlag;
#ifdef HERMITE_INTEGRATION
#ifdef HERMITE_INTEGRATION
struct HermiteWalkState   HermiteWalk;
struct DriftKickTableView HermiteWalkTables;
#endif
int HermiteOnlyFlag;            /*! Flag used to indicate whether to skip non-Hermite integrated particles in the force evaluation */
#endif

int NumPart;			/*!< number of particles on the LOCAL processor */
int N_gas;			/*!< number of gas particles on the LOCAL processor  */
#ifdef SINK_WIND_SPAWN
double  Max_Unspawned_MassUnits_fromSink;
#endif

long long Ntype[6];		/*!< total number of particles of each type */
int NtypeLocal[6];		/*!< local number of particles of each type */

gizmo_rng_t random_generator;	/*!< the random number generator used */

int Gas_split;           /*!< current number of newly-spawned gas particles outside block */
#ifdef GALSF
int Stars_converted;		/*!< current number of star particles in gas particle block */
#endif
#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
int Grains_promoted;		/*!< current number of grain particles promoted to solid body in gas block */
#endif

double TimeOfLastTreeConstruction;	/*!< holds what it says */

std::vector<int> Ngblist;			/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */
double *R2ngblist;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int *DomainStartList, *DomainEndList;



double *DomainWork;
int *DomainCount;
int *DomainCountGas;
int *DomainTask;
int *DomainNodeIndex;
int *TopNodeNodeIndex;   /* [NTopnodes] topnode -> Nodes[] slot (top-leaf router geometry SSOT) */
int *DomainList, DomainNumChanged;
peanokey *Key, *KeySorted;
struct topnode_data *TopNodes;
int NTopnodes, NTopleaves;

#ifdef SUBFIND
int GrNr;
int NumPartGroup;
#endif


/* variables for input/output , usually only used on process 0 */


char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

FILE
#ifdef OUTPUT_ADDITIONAL_RUNINFO
*FdTimebin,    /*!< file handle for timebin.txt log-file. */
*FdInfo,       /*!< file handle for info.txt log-file. */
*FdEnergy,     /*!< file handle for energy.txt log-file. */
*FdTimings,    /*!< file handle for timings.txt log-file. */
*FdBalance,    /*!< file handle for balance.txt log-file. */
#ifdef RT_CHEM_PHOTOION
*FdPhotoIonChemStats,         /*!< file handle for radtransfer.txt log-file. */
#endif
#ifdef TURB_DRIVING
*FdTurb,        /*!< file handle for turb.txt log-file */
#endif
#ifdef GR_TABULATED_COSMOLOGY
*FdDE,			/*!< file handle for darkenergy.txt log-file. */
#endif
#endif
*FdCPU;        /*!< file handle for cpu.txt log-file. */

#ifdef GALSF
FILE *FdSfr;			/*!< file handle for sfr.txt log-file. */
#endif
#ifdef GALSF_FB_FIRE_RT_LOCALRP
FILE *FdMomWinds;	/*!< file handle for MomWinds.txt log-file */
#endif
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
FILE *FdHIIHeating;	/*!< file handle for HIIheating.txt log-file */
#endif
#ifdef GALSF_FB_MECHANICAL
FILE *FdSNeFBLogFile;	/*!< file handle for SNIIheating.txt log-file */
#endif

#ifdef SINK_PARTICLES
FILE *FdSinks;		/*!< file handle for sinks.txt log-file. */
#ifdef OUTPUT_SINK_ACCRETION_HIST
FILE *FdSinkSwallowDetails;
#endif
#ifdef OUTPUT_SINK_FORMATION_PROPS
FILE *FdSinkFormationDetails;
#endif
#if defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(SINK_OUTPUT_MOREINFO)
FILE *FdSinksDetails;
#ifdef SINK_OUTPUT_MOREINFO
FILE *FdSinkMergerDetails;
#ifdef SINK_WIND_KICK
FILE *FdSinkWindDetails;
#endif
#endif
#endif
#endif
#if (defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(CBE_INTEGRATOR_OUTPUT_MOREINFO)) && defined(CBE_INTEGRATOR)
FILE *FdCbeDiagnostics;
#endif
#ifdef CBE_INTEGRATOR
/* see allvars.h for the contract. Zero-initialized
 * here; reset to {0,...} at the top of read_ic() each call; set per
 * PartType inside read_ic.cc's `if(hdf5_dataset >= 0)` block on
 * successful H5Dread of the VlasovMoments dataset. */
int CBE_Moments_LoadedFromIC_PType[6] = {0, 0, 0, 0, 0, 0};
#endif











/*! table for the cosmological drift factors */
double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
double GravKickTable[DRIFT_TABLE_LENGTH];

/*! log-time bounds the two tables above were built over */
double DriftTable_logTimeBegin;
double DriftTable_logTimeMax;

void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
/* All is defined here as the plain host global. GPU builds also create
   per-TU __managed__ mirrors (`AllDeviceMirror`, declared in each GPU
   TU that includes declarations/gpu_all_mirror.h) and redirect `All` to
   the TU's local mirror during the device compilation pass only. The
   central `gizmo_gpu_sync_all()` in cooling/cooling.cc copies host All
   into every auto-registered mirror before each GPU dispatch. Host
   code reads the extern below unconditionally; device code reads the
   freshly-synced local mirror. */
struct global_data_all_processes All;

/* Canonical out-of-line host accessors for `All.*`. Useful where explicit
 * host-extern intent helps readability or where preprocessor state may
 * be ambiguous; with the device-pass-gated mirror redirect the prior
 * "host wrapper in a GPU TU must route through here" rule is no longer
 * load-bearing. */
extern "C" {
struct global_data_all_processes *gizmo_host_all_ptr(void)
{
    return &All;
}
integertime gizmo_host_ti_current(void)
{
    return All.Ti_Current;
}
}




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P,	/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



/* the following struture holds data that is stored for each gas cell in addition to the collisionless
 * variables.
 */
struct gas_cell_data *CellP,	/*!< holds gas cell data on local processor */
 *DomainGasBuf;			/*!< buffer for gas cell data in domain decomposition */

peanokey *DomainKeyBuf;

/* global state of system
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

struct data_index *DataIndexTable;	/*!< records would-be-exported targets per task for
					   the gravity LET-incompleteness detector; never
					   shipped (the export round-trip is retired) */

struct data_nodelist *DataNodeList;

/* GravDataIn/Get/Result/Out + their gravdata_in/out structs are RETIRED (the gravity
 * MPI export round-trip is gone; gravity runs on the target-owning rank via GPU+LET). */

struct addFB_evaluate_data_in_ *addFB_evaluate_DataIn_, *addFB_evaluate_DataGet_; /*!< hold partial results of feedback calls if using various feedback algorithms >*/


struct info_block *InfoBlock;

/*! Header for the standard file format.
 */
struct io_header header;	/*!< holds header for snapshot files */


#ifdef SINK_PARTICLES
int N_active_loc_Sink=0;       /*!< number of active sink particles on the LOCAL processor */
struct sink_temp_particle_data *SinkTempInfo; /*! declare this structure, we'll malloc it below */
#endif

/*
 * Variables for Tree
 * ------------------
 */

long Nexport;   /* gravity LET-incompleteness detector count (Nimport retired with the export round-trip) */
int BufferCollisionFlag;
int BufferFullFlag;
int NextParticle;
int NextJ;
int TimerFlag;

struct NODE *Nodes_base,	/*!< points to the actual memory allocated for the nodes */
*Nodes;			            /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.TreeNodeIndexBase] gives the first allocated node */
struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/*!< maximum allowed number of internal nodes */
int Numnodestree;		/*!< number of (internal) nodes in each tree */
int MaxForeignNodes = 0;        /*!< LET: ceiling of the foreign-node index range; set in force_treeallocate. */
int AllocatedForeignNodes = 0;  /*!< LET: foreign-node slots that actually have storage; raised once per tree build to this rank's exact import. */
int Numforeignnodes = 0;        /*!< LET: foreign nodes currently installed; reset on each LET exchange. */
long long RuntimeMinLETForeignNodes = 0;  /*!< adaptive lower bound on MaxForeignNodes; ratcheted up by force_treebuild on a retryable LET overflow. Not a parameter; restart-persisted. */
long long Numforeignnodes_highwater = 0;  /*!< since-start peak of Numforeignnodes actually installed (diagnostic; memory ledger foreign-used). */


int *Nextnode;			/*!< gives next node in tree walk  (nodes array) */
int *Father;			/*!< gives parent node in tree (Prenodes array) */


int maxThreads = 1;

#if defined(DM_SIDM)
#endif
