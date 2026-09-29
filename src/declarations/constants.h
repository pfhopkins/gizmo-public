/* values of various constants to be parsed as part of global definitions */

/* Minimum number of active particles for GPU offload to be worthwhile.
 * Below this threshold, the kernel launch + memory transfer overhead exceeds
 * the GPU compute benefit, so we fall back to the CPU (OpenMP) path.
 * Applies to cooling, nuclear burning, and any future Kokkos-parallelized loops. */
#ifndef GPU_MIN_PARTICLES_FOR_OFFLOAD
#define GPU_MIN_PARTICLES_FOR_OFFLOAD 4096
#endif
/* DEBUG BISECTION: override to force CPU fallback for cooling/nuclear/turb.
   Uncomment to disable GPU kernels while keeping GPU kernels active
   (All_dev sync, __managed__ CoolTables, #pragma omp declare target, etc.) */
//#undef GPU_MIN_PARTICLES_FOR_OFFLOAD
//#define GPU_MIN_PARTICLES_FOR_OFFLOAD 999999999

/* Minimum number of particles needing a drift for the batched drift to run on the
 * device. The drift needs its own threshold rather than the shared one above: its
 * per-element work is far smaller than cooling's, while it stages the same whole
 * structs, so the two want opposite policies and the shared constant is read by
 * seven translation units.
 *
 * Set beyond any per-rank count so the drift runs on the host by default. Measured
 * on a 64-rank galaxy run, comparing clean builds with and without the batched path
 * at matched simulation time: the device path won modestly between roughly 4k and
 * 64k particles needing a drift per rank (-0.6 s and -0.7 s), and lost above that
 * (+9.5 s), as a uniform per-step increase rather than a few expensive steps. The
 * loss tracks the number of elements staged, which points at the cost of copying
 * whole particle and cell structs both ways rather than at launch overhead.
 *
 * ⚠ That measurement is for a configuration whose equation of state is cheap. It is
 * NOT a general verdict on running the drift on the device, and must not be read as
 * one. What decides the question is the work done per element against the bytes
 * staged for it: the drift loses here because an ideal-gas pressure update is trivial
 * beside copying whole particle and cell structs both ways, whereas cooling wins
 * comfortably on the same staging because its per-element solve is large. Under a
 * tabulated or solid equation of state the drift moves to the second camp -- the
 * per-particle pressure update becomes a Newton root-find over table lookups, with a
 * bisection fallback and an iteration cap in the tens -- so the same staging is
 * amortised over far more arithmetic and the balance can invert.
 *
 * Lower it to re-enable the device path once the staged set is narrowed to the fields
 * the drift actually writes, to select the band above if that is measured to be worth
 * the dispatch complexity, or per-configuration once the tabulated equations of state
 * are themselves device-callable. */
#ifndef GPU_MIN_PARTICLES_FOR_DRIFT_OFFLOAD
#define GPU_MIN_PARTICLES_FOR_DRIFT_OFFLOAD 2000000000
#endif

/* Minimum number of sources in ONE batch for a fused neighbour walk to run on the
 * device -- the upper of the two boundaries in a three-way choice, matching what
 * cooling and the other batched loops already do:
 *
 *     sources <  max(64, 4*threads)                  one host core
 *     sources >= max(64, 4*threads), < this constant host threads
 *     sources >= this constant                       the device
 *
 * The lower boundary is not a second constant; it is the neighbour-loop runner's
 * existing OpenMP work floor, which scales with the thread count. ⚠ Setting this
 * constant down to that floor would leave the threaded tier ZERO WIDTH and quietly
 * reduce the choice to two ways, so it belongs well above it.
 *
 * The walk is pointer-chasing down a dependent chain per source, so a batch of a
 * few sources leaves the device with almost nothing to overlap while still paying
 * a launch and a fence, where a host core pays neither. Per leaf visit the device
 * measures far worse than a host core at small and intermediate batch sizes and
 * reaches parity only at the largest. That deficit is a defect in this walk, to be
 * removed rather than routed around permanently, and the boundary below is set
 * where it currently pays rather than where it ought to end up.
 *
 * Decided per batch and per rank, from that batch's own source count: a rank's own
 * actives for a walk from the root, a group's received queries for a resumed one.
 * On a clustered run those two differ by three orders of magnitude inside a single
 * call, so no step-level or global property can stand in for either.
 *
 * At 1 every non-empty batch goes to the device; raising it sends batches below it
 * to the host instead.
 *
 * ⚠ THIS VALUE IS TEMPORARY AND MUST BE RE-PRICED, NOT INHERITED.
 *
 * Measured on a 128-rank cosmological zoom at matched simulation time, with the work
 * matched to 0.002% and every other cost row flat to within a second: against sending
 * every batch to the device, this saves 231 s of a 685 s wall, 236 s of it on the
 * density row alone. Replicate spread on that vehicle is ~0.3%, about 2 s, so the move
 * is a hundred times the noise; a boundary of 64 instead was measured at 133 s. On a
 * gas-only galaxy the same change is worth only 14 s, because that problem spends its
 * neighbour time in large batches where the two paths are within a few percent of each
 * other; the gain here comes from small and intermediate batches, which is where the
 * zoom spends nearly all of its.
 *
 * ⚠ What that measures is a DEFICIT IN THIS WALK'S DEVICE PATH, not a verdict that the
 * host is the right home for neighbour finding. The tile-and-BVH search on the same
 * hardware does not show it, so the work belongs on the device and the number above is
 * the size of the bill currently being paid. Sending EVERY batch to the host was also
 * measured and buys only about 4 s beyond this, so nothing material is left on the
 * table by stopping here -- and going further would be the wrong direction anyway: it
 * would concede the device for a widening share of the code, could cancel the gain
 * from moving the particle drift back onto the device, and would make the host look
 * like the better home for work -- memory placement especially -- when it is not.
 *
 * ⇒ RE-PRICE AND LOWER OR REMOVE THIS when any of the following lands: supply-type
 * pruning of the walk, bounded-occupancy terminal leaves, neighbour discovery through
 * an index built for the search rather than the gravity tree, the device drift, or the
 * next pricing of the three neighbour-search paths against each other. Each of those
 * attacks the deficit this compensates for; a threshold left in place afterwards hides
 * its own obsolescence. */
#ifndef GPU_MIN_SOURCES_FOR_WALK_OFFLOAD
#define GPU_MIN_SOURCES_FOR_WALK_OFFLOAD 4096
#endif

/* The Saitoh & Makino (2009) timestep-limiter factor: a cell is woken when a neighbour's step is
   this much shorter than its own. A run may set it in Config.sh -- their Table 1 gives f=2 a smaller
   energy error than f=4 on Sedov, at roughly twice the cost -- and otherwise it follows the slope
   limiter, which is the choice this code has always made. */
#ifndef WAKEUP
#if (SLOPE_LIMITER_TOLERANCE > 0)
#define WAKEUP   4.1            /* allows 2 timestep bins within kernel */
#else
#define WAKEUP   2.1            /* allows only 1-separated timestep bins within kernel */
#endif
#endif

#define ARENA_SHARE_OF_TASK_MEMORY          0.90 /* most of a task's memory share the working pool may take when the code has to bring its own size down to fit the machine; the rest is for the particle arrays and the tree, which are not allocated from the pool */
#define REDUC_FAC_FOR_MEMORY_IN_DOMAIN      0.98 /* used to pad memory in domain decomposition structures, should be slightly less than unity */

#define TREE_ALLOC_FACTOR_START             0.45 /* tree nodes allocated per particle at startup, before the run ratchets it upward to fit the particle distribution.  Named here because the startup memory projection has to use it before init() assigns it */

/* Cushion on the particle storage a rank holds for the coming epoch, as a fraction of the particles
 * per rank.  It covers the amount by which the next epoch's worst ghost import may exceed the one
 * just measured; on a galaxy-scale zoom that overshoot is a few percent of the per-rank count at the
 * median and around an eighth of it at the ninetieth percentile, so this covers the common case with
 * room to spare.  The rare import that exceeds it is met by growing the storage when it arrives. */
#define DOMAIN_CAPACITY_SAFETY_FRACTION     0.15

/* Releasing storage means migrating every particle array, so it is only worth doing for a saving that
 * is a real proportion of the capacity AND a real amount of memory.  The second test is expressed
 * against the particles per rank rather than as a slot count, so that it means the same thing on a
 * small problem as on a large one; a fixed count would silently forbid ever releasing anything on a
 * problem whose whole per-rank capacity is smaller than that count. */
#define DOMAIN_CAPACITY_SHRINK_FRACTION     0.10
#define DOMAIN_CAPACITY_SHRINK_MIN_FRAC     0.10

/* these are tolerances for the slope-limiters. we define them here, because the gradient constraint routine needs to be sure to use the -same- values in both the gradients and reimann solver routines */
#if MHD_CONSTRAINED_GRADIENT
#if (MHD_CONSTRAINED_GRADIENT > 1)
#define MHD_CONSTRAINED_GRADIENT_FAC_MINMAX 7.5
#define MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV 5.0
#define MHD_CONSTRAINED_GRADIENT_FAC_MED_PM 0.25
#define MHD_CONSTRAINED_GRADIENT_FAC_MAX_PM 0.25
#else
#define MHD_CONSTRAINED_GRADIENT_FAC_MINMAX 7.5
#define MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV 1.5
#define MHD_CONSTRAINED_GRADIENT_FAC_MED_PM 0.2
#define MHD_CONSTRAINED_GRADIENT_FAC_MAX_PM 0.2
#endif
#else
#define MHD_CONSTRAINED_GRADIENT_FAC_MINMAX 2.0
#define MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV 1.0
#define MHD_CONSTRAINED_GRADIENT_FAC_MED_PM 0.20
#define MHD_CONSTRAINED_GRADIENT_FAC_MAX_PM 0.125
#endif

/* Domain-decomposition granularity.  Each rank's share of the Peano-Hilbert curve is handed to it
   as All.DomainSegmentsPerRank separate contiguous segments rather than one.  More segments let the
   decomposition balance work more finely, because the assignment has more and smaller pieces to
   distribute; they cost more top-tree refinement, one more pseudo-particle all-gather each, and a
   more spatially fragmented rank territory.  The segment count is therefore chosen automatically at
   every full decomposition to hold roughly this many particles, which is where the balance gain
   stops paying for those costs. */
#define  DOMAIN_TARGET_PARTICLES_PER_SEGMENT   10000
#define  DOMAIN_MIN_SEGMENTS_PER_RANK          1
#define  DOMAIN_MAX_SEGMENTS_PER_RANK          256   /* every segment posts its own pseudo-particle
                                                        all-gather, so this bounds the message count
                                                        on a large problem spread over few ranks */
/* A small problem gets a single segment per rank, which on its own would leave the decomposition
   almost nothing to balance with.  This is the floor on how many top-tree leaves each rank is
   refined towards regardless -- a feasibility floor, not a tuned optimum; it sits below every
   measured landing point, so it binds only in the small-problem regime it exists for. */
#define  DOMAIN_MIN_TOPLEAVES_PER_RANK         16

/* Change the segment count by this factor relative to the automatic choice above.  The result is
   still held to DOMAIN_MAX_SEGMENTS_PER_RANK, so raising this has no effect once the automatic
   choice has already reached that ceiling. */
#ifndef  DOMAIN_SEGMENTS_SCALE
#define  DOMAIN_SEGMENTS_SCALE                 1
#endif
/* How finely the top tree is refined, relative to the number of domain segments it has to fill.
   Refinement splits eight ways at a time, so the leaf count achieved overshoots the number
   requested here; both are reported at each decomposition.  Must be at least 1, or the top tree
   can end up with fewer leaves than there are segments to assign them to. */
#ifndef  DOMAIN_TOPTREE_REFINEMENT_SCALE
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
#define  DOMAIN_TOPTREE_REFINEMENT_SCALE       4
#else
#define  DOMAIN_TOPTREE_REFINEMENT_SCALE       1
#endif
#endif

#ifndef  GRAVCOSTLEVELS
#define  GRAVCOSTLEVELS      20
#endif

#define  NUMBER_OF_MEASUREMENTS_TO_RECORD  6  /* this is the number of past executions of a timebin that the reported average CPU-times average over */

#define  NODELISTLENGTH      8

#define  GRAVITY_LET_DETECTOR_ENTRIES  4096  /* how many targets the gravity walk can record as not covered by the locally-built tree. Recording even one stops the run, so this only has to be large enough to report what happened, with room for every thread to record at once */

/* Maximum number of particles a tree node may hold without being subdivided further: a node at or
 * below this count becomes a terminal "bucket" leaf whose particles are visited directly, instead of
 * being split down to one particle per leaf. Trading a few extra direct particle tests for a much
 * shallower tree is standard (GADGET-4 TREE_NUM_BEFORE_NODESPLIT, ChaNGa bucketSize, SPH-EXA/SWIFT
 * leaf sizes, and our own SFC tile size), and it also bounds the depth that near-coincident
 * particles can force.
 *
 * There is no universally good value: the best choice depends on the problem, on how aggressively
 * the opening criterion descends, and on the machine. Measured across a clustered cosmological zoom
 * and a shallower-tree galaxy problem, 4 is a win on the first and neutral on the second, so it is
 * the default; larger values suit strongly clustered problems and are limited in practice by the
 * memory needed to export the tree between ranks. 1 reproduces the historical one-particle-per-leaf
 * tree exactly, and compiles the multi-particle leaf code out entirely. Override in Config.sh. */
#ifndef TREE_LEAF_BUCKET_SIZE
#define TREE_LEAF_BUCKET_SIZE 4
#endif
#if TREE_LEAF_BUCKET_SIZE < 1
/* Below one, no child count satisfies the terminal-leaf test, so a single particle would be neither
 * stored in a leaf slot nor given a node of its own: it would drop out of the tree silently. */
#error "TREE_LEAF_BUCKET_SIZE must be at least 1 (1 = one particle per leaf, the historical tree)"
#endif

/* How many targets the host gravity walk takes through the tree together. Targets adjacent in the
 * active list are usually close in space (the list follows the particle order, which follows the
 * space-filling curve), so their walks visit nearly the same nodes; a packet traverses the tree
 * once for all of them, every member judges each node by its own opening criterion, and each
 * member then evaluates the elements it accepted, so its force is exactly what its own walk would
 * produce. What is shared is the traversal, which is most of the cost where a few targets descend
 * a locally very deep tree. 1 walks every target alone, the historical walk.
 *
 * No single value suits every problem or machine. Measured on a clustered cosmological zoom at 128
 * ranks, 8 cut the host walk by 14% and was never slower than 1 on any class of step, while 64 gave
 * the same gain on the steps that matter most but cost more on the many tiny steps, where a large
 * packet does far more opening-criterion work than the node loads it shares; so 8 is the default.
 * Override in Config.sh. */
#ifndef TREE_QUERY_PACKET_SIZE
#define TREE_QUERY_PACKET_SIZE 8
#endif
#if TREE_QUERY_PACKET_SIZE < 1
#error "TREE_QUERY_PACKET_SIZE must be at least 1 (1 = every target walks the tree alone)"
#endif


#define  EPSILON_FOR_TREERND_SUBNODE_SPLITTING (1.0e-4) /* define some number << 1; particles with less than this separation will trigger randomized sub-node splitting in the tree. we set it to a global value here so that other sub-routines will know not to force particle separations below this */

#if !defined(EOS_GAMMA)
#define EOS_GAMMA (5.0/3.0) /*!< adiabatic index of simulated gas */
#endif

#ifdef GAMMA_ENFORCE_ADIABAT
#define EOS_ENFORCE_ADIABAT (GAMMA_ENFORCE_ADIABAT) /* this allows for either term to be defined, for backwards-compatibility */
#endif


#if !defined(RT_HYDROGEN_GAS_ONLY) || defined(RT_CHEM_PHOTOION_HE)
#define  HYDROGEN_MASSFRAC 0.76 /*!< mass fraction of hydrogen, relevant only for radiative cooling */
#else
#define  HYDROGEN_MASSFRAC 1.0  /*!< mass fraction of hydrogen, relevant only for radiative cooling */
#endif


#define  MAX_REAL_NUMBER  1e56
#define  MIN_REAL_NUMBER  1e-56

#if (defined(MAGNETIC) && !defined(COOLING)) || defined(EOS_ELASTIC)
#define  CONDITION_NUMBER_DANGER  1.0e7 /*!< condition number above which we will not trust matrix-based gradients */
#else
#define  CONDITION_NUMBER_DANGER  1.0e3 /*!< condition number above which we will not trust matrix-based gradients */
#endif

/* ... often used physical constants (cgs units). note many of these are defined to better precision with different units in e.g. the GSL package, these are purely for user convenience */
#define  GRAVITY_G_CGS      (6.672e-8)
#define  SOLAR_MASS_CGS     (1.989e33)
#define  SOLAR_LUM_CGS      (3.826e33)
#define  SOLAR_RADIUS_CGS   (6.957e10)
#define  BOLTZMANN_CGS      (1.38066e-16)
#define  C_LIGHT_CGS        (2.9979e10)
#define  PROTONMASS_CGS     (1.6726e-24)
#define  ELECTRONMASS_CGS   (9.10953e-28)
#define  THOMPSON_CX_CGS    (6.65245e-25)
#define  ELECTRONCHARGE_CGS (4.8032e-10)
#define  SECONDS_PER_YEAR   (3.155e7)
#define  HUBBLE_H100_CGS    (3.2407789e-18)    /* in h/sec */
#define  ELECTRONVOLT_IN_ERGS (1.60217733e-12)
#define HABING_FLUX_CGS      (1.6e-3)
#define DRAINE_FLUX_CGS      (1.7 * HABING_FLUX_CGS)


/* and a bunch of useful unit-conversion macros pre-bundled here, to help keep the 'h' terms and other correct */
#define UNIT_MASS_IN_CGS        ((All.UnitMass_in_g/All.HubbleParam))
#define UNIT_VEL_IN_CGS         ((All.UnitVelocity_in_cm_per_s))
#define UNIT_LENGTH_IN_CGS      ((All.UnitLength_in_cm/All.HubbleParam))
#define UNIT_TIME_IN_CGS        (((UNIT_LENGTH_IN_CGS)/(UNIT_VEL_IN_CGS)))
#define UNIT_ENERGY_IN_CGS      (((UNIT_MASS_IN_CGS)*(UNIT_VEL_IN_CGS)*(UNIT_VEL_IN_CGS)))
#define UNIT_PRESSURE_IN_CGS    (((UNIT_ENERGY_IN_CGS)/(UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS)))
#define UNIT_DENSITY_IN_CGS     (((UNIT_MASS_IN_CGS)/(UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS)))
#define UNIT_SPECEGY_IN_CGS     (((UNIT_PRESSURE_IN_CGS)/(UNIT_DENSITY_IN_CGS)))
#define UNIT_SURFDEN_IN_CGS     (((UNIT_DENSITY_IN_CGS)*(UNIT_LENGTH_IN_CGS)))
#define UNIT_FLUX_IN_CGS        (((UNIT_PRESSURE_IN_CGS)*(UNIT_VEL_IN_CGS)))
#define UNIT_LUM_IN_CGS         (((UNIT_ENERGY_IN_CGS)/(UNIT_TIME_IN_CGS)))
#define UNIT_B_IN_GAUSS         ((sqrt(4.*M_PI*UNIT_PRESSURE_IN_CGS)))
#define UNIT_MASS_IN_SOLAR      (((UNIT_MASS_IN_CGS)/SOLAR_MASS_CGS))
#define UNIT_DENSITY_IN_NHCGS   (((UNIT_DENSITY_IN_CGS)/PROTONMASS_CGS))
#define UNIT_TIME_IN_YR         (((UNIT_TIME_IN_CGS)/(SECONDS_PER_YEAR)))
#define UNIT_TIME_IN_MYR        (((UNIT_TIME_IN_CGS)/(1.e6*SECONDS_PER_YEAR)))
#define UNIT_TIME_IN_GYR        (((UNIT_TIME_IN_CGS)/(1.e9*SECONDS_PER_YEAR)))
#define UNIT_LENGTH_IN_SOLAR    (((UNIT_LENGTH_IN_CGS)/SOLAR_RADIUS_CGS))
#define UNIT_LENGTH_IN_AU       (((UNIT_LENGTH_IN_CGS)/1.496e13))
#define UNIT_LENGTH_IN_PC       (((UNIT_LENGTH_IN_CGS)/3.085678e18))
#define UNIT_LENGTH_IN_KPC      (((UNIT_LENGTH_IN_CGS)/3.085678e21))
#define UNIT_PRESSURE_IN_EV     (((UNIT_PRESSURE_IN_CGS)/ELECTRONVOLT_IN_ERGS))
#define UNIT_VEL_IN_KMS         (((UNIT_VEL_IN_CGS)/1.e5))
#define UNIT_LUM_IN_SOLAR       (((UNIT_LUM_IN_CGS)/SOLAR_LUM_CGS))
#define UNIT_FLUX_IN_HABING     (((UNIT_FLUX_IN_CGS)/HABING_FLUX_CGS))
#define UNIT_EGY_DENSITY_IN_HABING ((UNIT_PRESSURE_IN_CGS)/(HABING_FLUX_CGS / C_LIGHT_CGS))
#define U_TO_TEMP_UNITS         ((PROTONMASS_CGS/BOLTZMANN_CGS)*((UNIT_ENERGY_IN_CGS)/(UNIT_MASS_IN_CGS))) /* units to convert specific internal energy to temperature. needs to be multiplied by dimensionless factor=mean_molec_weight_in_amu*(gamma_eos-1) */
#define MEAN_MOLECULAR_WEIGHT_IONIZED (0.59) /* mean molecular weight of fully-ionized primordial gas, in units of the proton mass. use this where the gas is genuinely ionized by construction, e.g. inside an HII region; where a run simply solves no chemistry, MEAN_MOLECULAR_WEIGHT_DEFAULT is the value to use, and a user can set that one */
#define MEAN_MOLECULAR_WEIGHT_ATOMIC (1.28) /* mean molecular weight of fully-atomic solar metallicity gas */
#define MEAN_MOLECULAR_WEIGHT_MOLECULAR (2.3) /* mean molecular weight of fully-molecular solar metallicity gas */
#ifndef C_LIGHT_CODE
#define C_LIGHT_CODE            ((C_LIGHT_CGS/UNIT_VEL_IN_CGS)) /* pure convenience function, speed-of-light in code units */
#endif
#define C_LIGHT_CODE_REDUCED    (((RT_SPEEDOFLIGHT_REDUCTION)*(C_LIGHT_CODE))) /* reduced speed-of-light in code units, again here as a convenience function, but just returns constant */
#define H0_CGS                  ((All.HubbleParam*HUBBLE_H100_CGS)) /* actual value of H0 in cgs */
#define COSMIC_BARYON_DENSITY_CGS ((All.OmegaBaryon*(H0_CGS)*(H0_CGS)*(3./(8.*M_PI*GRAVITY_G_CGS))*All.cf_a3inv)) /* cosmic mean baryon density [scale-factor-dependent] in cgs units */



#ifdef RT_COMOVING
#define RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS (0) /* this prefactor goes in front of various terms which vanish in the comoving frame RHD equations */
#else
#define RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS ((C_LIGHT_CODE_REDUCED)/(C_LIGHT_CODE)) /* these terms in the mixed-frame equations need to be multiplied by c_reduced/c */
#endif


#ifdef GALSF_FB_FIRE_RT_HIIHEATING
#define HIIRegion_Temp (1.0e4) /* temperature (in K) of heated gas */
#endif


#if defined(COOLING) || defined(RT_INFRARED)
#define MAX_DUST_TEMP 1.0e4 // maximum dust temperature for which we expect to call opacity or dust-to-metals ratio functions
#endif

/* some flags for the field "flag_ic_info" in the file header */
#define FLAG_ZELDOVICH_ICS     1
#define FLAG_SECOND_ORDER_ICS  2
#define FLAG_EVOLVED_ZELDOVICH 3
#define FLAG_EVOLVED_2LPT      4
#define FLAG_NORMALICS_2LPT    5


#ifndef PM_ASMTH
#define PM_ASMTH (1.25) /*! PM_ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#endif
#ifndef PM_RCUT
#define PM_RCUT (4.5) /*! PM_RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range forces are evaluated in the short-range tree walk. */
#endif
#define MAXLEN_OUTPUTLIST 1201    /*!< maxmimum number of entries in output list */
#define DRIFT_TABLE_LENGTH 1000    /*!< length of the lookup table used to hold the drift and kick factors */
#define MAXITER 150

#ifndef LINKLENGTH
#define LINKLENGTH (0.2)
#endif
#ifndef FOF_GROUP_MIN_SIZE
#ifdef FOF_GROUP_MIN_LEN
#define FOF_GROUP_MIN_SIZE FOF_GROUP_MIN_LEN
#else
#define FOF_GROUP_MIN_SIZE 32
#endif
#endif
#ifndef SUBFIND_ADDIO_NUMOVERDEN
#define SUBFIND_ADDIO_NUMOVERDEN 1
#endif




#define CPU_ALL            0
#define CPU_TREEWALK1      1
#define CPU_TREEWALK2      2
#define CPU_TREEWAIT1      3
#define CPU_TREEWAIT2      4
#define CPU_TREESEND       5
#define CPU_TREERECV       6
#define CPU_TREEMISC       7
#define CPU_TREEBUILD      8
#define CPU_TREEHMAXUPDATE 9
#define CPU_DOMAIN         10
#define CPU_DENSCOMPUTE    11
/* Slot 12: holds the GRADIENTS bracket (printed as "gradients"); the name is
 * kept so slot numbers, and the meaning of columns in old logs, do not shift. */
#define CPU_DENSWAIT       12
#define CPU_DENSCOMM       13
/* Slot 14: "misc_hydro" catch-all; also receives the drift charged by the
 * shared ghost-import helpers, whose callers include non-hydro loops. */
#define CPU_DENSMISC       14
#define CPU_HYDCOMPUTE     15
#define CPU_HYDWAIT        16
#define CPU_HYDCOMM        17
#define CPU_HYDMISC        18
#define CPU_DRIFT          19
/* Slot 20: was CPU_TIMELINE (mislabeled "kicks" in cpu.txt). It's actually
 * populated only by find_timesteps() in core/timestep.cc, so renamed to
 * reflect its real meaning. No "kicks" bucket exists in legacy CPU_Step;
 * kick wall time falls into CPU_MISC by default. */
#define CPU_FIND_TIMESTEPS 20
/* Slot 21: RETIRED. The potential is computed in-tree together with the force,
 * so no separate potential phase exists and nothing charges this bucket. The
 * slot number is kept so old logs and parsers keep their column meaning. */
#define CPU_POTENTIAL      21
#define CPU_MESH           22
#define CPU_PEANO          23
#define CPU_COOLINGSFR     24
#define CPU_SNAPSHOT       25
#define CPU_FOF            26
#define CPU_SINKS     27
#define CPU_MISC           28
#define CPU_DRAGFORCE      29
#define CPU_SNIIHEATING    30
#define CPU_HIIHEATING     31
#define CPU_LOCALWIND      32
#define CPU_COOLSFRIMBAL   33
#define CPU_AGSDENSCOMPUTE 34
#define CPU_AGSDENSWAIT    35
#define CPU_AGSDENSCOMM    36
#define CPU_AGSDENSMISC    37
#define CPU_DYNDIFFMISC       38
#define CPU_DYNDIFFCOMPUTE    39
#define CPU_DYNDIFFWAIT       40
#define CPU_DYNDIFFCOMM       41
#define CPU_IMPROVDIFFMISC    42
#define CPU_IMPROVDIFFCOMPUTE 43
#define CPU_IMPROVDIFFWAIT    44
#define CPU_IMPROVDIFFCOMM    45
#define CPU_RTNONFLUXOPS  46
/* Categories carved from the original CPU_DUMMY00..10 slots. One double
 * increment per call site — the same cost as any other CPU_Step bucket.
 * Slots marked "reserved" have no writer; they are kept so slot numbers, and
 * therefore the meaning of columns in existing logs, never shift. */
#define CPU_NGB_BUILD         47  /* neighbor-list (CSR) construction, all callers */
#define CPU_SIDX_REFRESH      48  /* sync-point spatial-index upkeep (release of the all-types index) */
#define CPU_GHOSTIMPORT_SYMM  49  /* hydro-corridor SYMMETRIC ghost import (sub-row of ghost import) */
#define CPU_MHD_MG            50  /* MHD divergence-cleaning global MG solve (took a reserved, never-charged slot, so no column moved) */
#define CPU_FORCE_UPDATE_TREE 51  /* dynamic tree update on reused-tree steps (node drift + top-node kick exchange) */
#define CPU_SINK_ENV          52  /* sink_environment GPU kernel + scatter */
#define CPU_SINK_FEEDSWK      53  /* sink_feed + sink_swallow_and_kick combined */
#define CPU_PAIR_KERNEL       54  /* neighbor-loop pair-kernel launch + fence (runner only) */
#define CPU_GHOSTIMPORT       55  /* imported-ghost exchange, all callers (hydro, gravity, sinks, feedback) */
#define CPU_KICKS             56  /* first+second half-step kicks, Hermite predict/correct */
#define CPU_LOGMSG            57  /* per-step statistics + log-file writes (a barrier-sync absorber) */

#define CPU_PARTS          58  /* this gives the number of parts above (must be last) */

#define CPU_STRING_LEN 120


#if (BOX_SPATIAL_DIMENSION==1) || defined(ONEDIM)
#define NUMDIMS 1           /* define number of dimensions and volume normalization */
#define VOLUME_NORM_COEFF_FOR_NDIMS 2.0
#elif (BOX_SPATIAL_DIMENSION==2) || defined(TWODIMS)
#define NUMDIMS 2
#define VOLUME_NORM_COEFF_FOR_NDIMS M_PI
#else
#define VOLUME_NORM_COEFF_FOR_NDIMS 4.188790204786  /* 4pi/3 */
#define NUMDIMS 3
#endif


#define CUBE_EDGEFACTOR_1 0.366025403785    /* CUBE_EDGEFACTOR_1 = 0.5 * (sqrt(3)-1) */
#define CUBE_EDGEFACTOR_2 0.86602540        /* CUBE_EDGEFACTOR_2 = 0.5 * sqrt(3) */


#if !defined(MEAN_MOLECULAR_WEIGHT_DEFAULT)
#define MEAN_MOLECULAR_WEIGHT_DEFAULT (MEAN_MOLECULAR_WEIGHT_IONIZED) /* default to ionized value but allow user-setting */
#endif
