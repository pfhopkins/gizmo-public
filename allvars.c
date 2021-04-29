/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (many new variables,
 * structures, and different naming conventions for some old variables)
 */

#include "allvars.h"




/*********************************************************/
/*  Global variables                                     */
/*********************************************************/



int ThisTask;			/*!< the number of the local processor  */
int NTask;			/*!< number of processors */
int PTask;			/*!< note: NTask = 2^PTask */

double CPUThisRun;		/*!< Sums CPU time of current process */

int NumForceUpdate;		/*!< number of active particles on local processor in current timestep  */
long long GlobNumForceUpdate;
int NumSphUpdate;		/*!< number of active SPH particles on local processor in current timestep  */

int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from initial conditions, 1 is resuming a run from a set of restart files, while 2 marks a restart from a snapshot file. */

int RestartSnapNum;
int SelRnd;

int *Exportflag;		/*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount;

int TakeLevel;

int FirstActiveParticle;
int *NextActiveParticle;
unsigned char *ProcessedFlag;

int TimeBinCount[TIMEBINS];
int TimeBinCountSph[TIMEBINS];
int TimeBinActive[TIMEBINS];

int FirstInTimeBin[TIMEBINS];
int LastInTimeBin[TIMEBINS];
int *NextInTimeBin;
int *PrevInTimeBin;

size_t HighMark_run, HighMark_domain, HighMark_gravtree,
  HighMark_pmperiodic, HighMark_pmnonperiodic, HighMark_sphdensity, HighMark_sphhydro, HighMark_GasGrad;

#ifdef BOX_PERIODIC
MyDouble boxSize, boxHalf;      /* size of the box! these variables are technically redundant but used -constantly- so very helpful */
#endif
#ifdef BOX_LONG_X
MyDouble boxSize_X, boxHalf_X;
#endif
#ifdef BOX_LONG_Y
MyDouble boxSize_Y, boxHalf_Y;
#endif
#ifdef BOX_LONG_Z
MyDouble boxSize_Z, boxHalf_Z;
#endif

#ifdef BOX_SHEARING
MyDouble Shearing_Box_Vel_Offset;
MyDouble Shearing_Box_Pos_Offset;
#endif

#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z)
short int special_boundary_condition_xyz_def_reflect[3];
short int special_boundary_condition_xyz_def_outflow[3];
#endif

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
MPI_Status mpistat;
#endif


#ifdef TURB_DRIVING
size_t HighMark_turbpower;
#endif

#ifdef GALSF
double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
double TimeBin_BH_mass[TIMEBINS];
double TimeBin_BH_dynamicalmass[TIMEBINS];
double TimeBin_BH_Mdot[TIMEBINS];
double TimeBin_BH_Medd[TIMEBINS];
#endif



char DumpFlag = 1;
size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;
double CPU_Step[CPU_PARTS];
char CPU_Symbol[CPU_PARTS] = {
    '-', '*', '=', ';', '<', '[', '^', ':', '.', '~', '|', '+', '"', '/',  '`', ',', '>', '@', '#', '&',
    '$', ']', '(', '?', ')', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '\\', '%', '{', '}', 'Z',
    'Y', 'X', 'U', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',  'n', 'o'};//, 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y, 'z',
char CPU_SymbolImbalance[CPU_PARTS] = {
    'a', 't', 'u', 'v', 'b', 'w', 'd', 'r', 'h', 'm', 'n', 'l', 'o', 'p',  's', 'f', 'i', 'g', 'c', 'e', // 20 columns here
    'x', 'y', 'z', 'A', 'I', 'W', 'T', 'V', 'B', 'C', 'D', 'E', 'F', 'G', 'H',  'I', 'J', 'K', 'L', 'Q',
    'R', 'S', 'T', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',  'N', 'O'};//, 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
char CPU_String[CPU_STRING_LEN + 1];
double WallclockTime;		/*!< This holds the last wallclock time measurement for timings measurements */
int Flag_FullStep;		/*!< Flag used to signal that the current step involves all particles */


int TreeReconstructFlag;
#ifdef WAKEUP
int NeedToWakeupParticles;      /*!< Flags used to signal that wakeups need to be processed at the beginning of the next timestep */
int NeedToWakeupParticles_local;
#endif
int GlobFlag;
#ifdef HERMITE_INTEGRATION
int HermiteOnlyFlag;            /*! Flag used to indicate whether to skip non-Hermite integrated particles in the force evaluation */
#endif

int NumPart;			/*!< number of particles on the LOCAL processor */
int N_gas;			/*!< number of gas particles on the LOCAL processor  */
#ifdef SEPARATE_STELLARDOMAINDECOMP
int N_stars;
#endif
#ifdef BH_WIND_SPAWN
double  MaxUnSpanMassBH;
#endif

long long Ntype[6];		/*!< total number of particles of each type */
int NtypeLocal[6];		/*!< local number of particles of each type */

gsl_rng *random_generator;	/*!< the random number generator used */

int Gas_split;           /*!< current number of newly-spawned gas particles outside block */
#ifdef GALSF
int Stars_converted;		/*!< current number of star particles in gas particle block */
#endif

double TimeOfLastTreeConstruction;	/*!< holds what it says */

int *Ngblist;			/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */
double *R2ngblist;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int *DomainStartList, *DomainEndList;



double *DomainWork;
int *DomainCount;
int *DomainCountSph;
int *DomainTask;
int *DomainNodeIndex;
int *DomainList, DomainNumChanged;

peanokey *Key, *KeySorted;

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;

#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
double RndTable[RNDTABLE];
#endif

#ifdef SUBFIND
int GrNr;
int NumPartGroup;
#endif


/* variables for input/output , usually only used on process 0 */


char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

FILE
#ifndef IO_REDUCED_MODE
*FdTimebin,    /*!< file handle for timebin.txt log-file. */
*FdInfo,       /*!< file handle for info.txt log-file. */
*FdEnergy,     /*!< file handle for energy.txt log-file. */
*FdTimings,    /*!< file handle for timings.txt log-file. */
*FdBalance,    /*!< file handle for balance.txt log-file. */
#ifdef RT_CHEM_PHOTOION
*FdRad,         /*!< file handle for radtransfer.txt log-file. */
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
#ifdef GALSF_FB_MECHANICAL
FILE *FdSneIIHeating;	/*!< file handle for SNIIheating.txt log-file */
#endif

#ifdef BLACK_HOLES
FILE *FdBlackHoles;		/*!< file handle for blackholes.txt log-file. */
#ifdef OUTPUT_SINK_ACCRETION_HIST
FILE *FdBhSwallowDetails;
#endif
#ifdef OUTPUT_SINK_FORMATION_PROPS
FILE *FdBhFormationDetails;
#endif
#if !defined(IO_REDUCED_MODE) || defined(BH_OUTPUT_MOREINFO)
FILE *FdBlackHolesDetails;
#ifdef BH_OUTPUT_MOREINFO
FILE *FdBhMergerDetails;
#ifdef BH_WIND_KICK
FILE *FdBhWindDetails;
#endif
#endif
#endif
#endif











/*! table for the cosmological drift factors */
double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
double GravKickTable[DRIFT_TABLE_LENGTH];

void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P,	/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *SphP,	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */

peanokey *DomainKeyBuf;

/* global state of system
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

struct data_index *DataIndexTable;	/*!< the particles to be exported are grouped
					   by task-number. This table allows the
					   results to be disentangled again and to be
					   assigned to the correct particle */

struct data_nodelist *DataNodeList;

struct gravdata_in *GravDataIn,	/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


struct gravdata_out *GravDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct potdata_out *PotDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct info_block *InfoBlock;

/*! Header for the standard file format.
 */
struct io_header header;	/*!< holds header for snapshot files */


#ifdef BLACK_HOLES
int N_active_loc_BHs=0;       /*!< number of active black holes on the LOCAL processor */
struct blackhole_temp_particle_data *BlackholeTempInfo; /*! declare this structure, we'll malloc it below */
#endif

/*
 * Variables for Tree
 * ------------------
 */

long Nexport, Nimport;
int BufferFullFlag;
int NextParticle;
int NextJ;
int TimerFlag;

struct NODE *Nodes_base,	/*!< points to the actual memory allocated for the nodes */
*Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] gives the first allocated node */
struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/*!< maximum allowed number of internal nodes */
int Numnodestree;		/*!< number of (internal) nodes in each tree */


int *Nextnode;			/*!< gives next node in tree walk  (nodes array) */
int *Father;			/*!< gives parent node in tree (Prenodes array) */


#if defined(PTHREADS_NUM_THREADS)
int maxThreads = PTHREADS_NUM_THREADS;
#else
int maxThreads = 1;
#endif


#if defined(DM_SIDM)
MyDouble GeoFactorTable[GEOFACTOR_TABLE_LENGTH];
#endif
