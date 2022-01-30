
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


#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "GIZMO_config.h"
/*------- Things that are always recommended (this must follow loading GIZMO_config.h!) -------*/
#define GIZMO_VERSION     2022  /*!< code version (should be an int corresponding to the year) */
#define DOUBLEPRECISION         /* using double (not floating-point) precision */
#define PEANOHILBERT            /* sort particles on a Peano-Hilbert curve (huge optimization) */
#define WALLCLOCK               /* track timing of different routines */
#define MYSORT                  /* use our custom sort (as opposed to C default, which is compiler-dependent) */
#define ALLOWEXTRAPARAMS        /* don't crash (just warn) if there are extra lines in the input parameterfile */
#define INHOMOG_GASDISTR_HINT   /* if the gas is distributed very different from collisionless particles, this can helps to avoid problems in the domain decomposition */
#ifndef OUTPUT_ADDITIONAL_RUNINFO
#define IO_REDUCED_MODE
#endif
#ifndef IO_DISABLE_HDF5
#define HAVE_HDF5               /* default to using HDF5 */
#include <hdf5.h>
#endif
#if !defined(OUTPUT_POSITIONS_IN_DOUBLE) && defined(HAVE_HDF5)
#define OUTPUT_POSITIONS_IN_DOUBLE /* recommended to always default to recording positions in double-precision: there's not really a good reason not to do this unless we need to match unformatted binary */
#endif
#if !defined(LONG_INTEGER_TIME)
#define LONG_INTEGER_TIME   /* always recommended: on modern machines the memory overhead cost of this is negligible */
#endif


#define DO_PREPROCESSOR_EXPAND_(VAL)  VAL ## 1
#define EXPAND_PREPROCESSOR_(VAL)     DO_PREPROCESSOR_EXPAND_(VAL) /* checks for a NON-ZERO value of this parameter */
#define CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(VAL) !(EXPAND_PREPROCESSOR_(VAL) == 1) /* returns True if a non-zero int value of VAL is set */

#if !defined(SLOPE_LIMITER_TOLERANCE)
#if defined(AGGRESSIVE_SLOPE_LIMITERS)
#define SLOPE_LIMITER_TOLERANCE 2
#else
#define SLOPE_LIMITER_TOLERANCE 1
#endif
#endif

#ifndef DISABLE_SPH_PARTICLE_WAKEUP
#if (SLOPE_LIMITER_TOLERANCE > 0)
#define WAKEUP   4.1            /* allows 2 timestep bins within kernel */
#else
#define WAKEUP   2.1            /* allows only 1-separated timestep bins within kernel */
#endif
#endif

/* lock the 'default' hydro mode */
#if !(defined(HYDRO_MESHLESS_FINITE_VOLUME) || defined(HYDRO_MESHLESS_FINITE_MASS) || defined(HYDRO_DENSITY_SPH) || defined(HYDRO_PRESSURE_SPH)) // default solver is not defined
#if (defined(HYDRO_FIX_MESH_MOTION) && (HYDRO_FIX_MESH_MOTION != 7)) || defined(HYDRO_REGULAR_GRID)
#define HYDRO_MESHLESS_FINITE_VOLUME /* only makes sense to use this modules with this 'backbone' of MFV here */
#else
#define HYDRO_MESHLESS_FINITE_MASS   /* otherwise default to MFM if nothing is specified */
#endif
#endif

/* define the default mesh-motion assumption, if this is not provided by the user */
#if !defined(HYDRO_FIX_MESH_MOTION)
#if defined(HYDRO_REGULAR_GRID)
#define HYDRO_FIX_MESH_MOTION 0     /* default to non-moving for regular grids */
#else
#define HYDRO_FIX_MESH_MOTION 5     /* otherwise default to smoothed motion, only relevant for MFV (MFM/SPH will always move with flow) */
#endif
#endif

/* determine whether the mesh is adaptive via splitting/merging (refinement) or 'frozen' to the initial number of elements */
#if !defined(PREVENT_PARTICLE_MERGE_SPLIT) && (HYDRO_FIX_MESH_MOTION<5)
#define PREVENT_PARTICLE_MERGE_SPLIT  /* particle merging/splitting doesn't make sense with frozen grids */
#endif

#ifdef PARTICLE_MERGE_SPLIT_EVERY_TIMESTEP
#define MAINTAIN_TREE_IN_REARRANGE
#endif


#ifdef PMGRID
#define PM_ENLARGEREGION 1.1    /* enlarges PMGRID region as the simulation evolves */
/*
#     - PM_ENLARGEREGION: The spatial region covered by the high-res zone has a fixed
#                   size during the simulation, which initially is set to the
#                   smallest region that encompasses all high-res particles. Normally, the
#                   simulation will be interrupted, if high-res particles leave this
#                   region in the course of the run. However, by setting this parameter
#                   to a value larger than one, the high-res region can be expanded.
#                   For example, setting it to 1.4 will enlarge its side-length by
#                   40% (it remains centered on the high-res particles). Hence, with
#                   such a setting, the high-res region may expand or move by a
#                   limited amount. If in addition SYNCHRONIZATION is activated, then
#                   the code will be able to continue even if high-res particles
#                   leave the initial high-res grid. In this case, the code will
#                   update the size and position of the grid that is placed onto
#                   the high-resolution region automatically. To prevent that this
#                   potentially happens every single PM step, one should nevertheless
#                   assign a value slightly larger than 1 to PM_ENLARGEREGION.
*/
#endif



#if (defined(HYDRO_DENSITY_SPH) || defined(HYDRO_PRESSURE_SPH)) && !defined(HYDRO_SPH)
#define HYDRO_SPH               /* top-level flag for SPH: must be enabled if any SPH method is used */
#endif
#ifdef HYDRO_SPH
#if !defined(SPH_DISABLE_CD10_ARTVISC) && !(defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)) // fancy viscosity switches assume positive pressures //
#define SPHAV_CD10_VISCOSITY_SWITCH 0.05   /* Enables Cullen & Dehnen 2010 'inviscid sph' (viscosity suppression outside shocks) */
#endif
#ifndef SPH_DISABLE_PM_CONDUCTIVITY
#define SPHAV_ARTIFICIAL_CONDUCTIVITY      /* Enables mixing entropy (J.Read's improved Price-Monaghan conductivity with Cullen-Dehnen switches) */
#endif
#endif

#ifdef PERIODIC
#define BOX_PERIODIC
#endif
#ifdef BND_PARTICLES
#define BOX_BND_PARTICLES
#endif
#ifdef LONG_X
#define BOX_LONG_X LONG_X
#endif
#ifdef LONG_Y
#define BOX_LONG_Y LONG_Y
#endif
#ifdef LONG_Z
#define BOX_LONG_Z LONG_Z
#endif
#ifdef REFLECT_BND_X
#define BOX_REFLECT_X
#endif
#ifdef REFLECT_BND_Y
#define BOX_REFLECT_Y
#endif
#ifdef REFLECT_BND_Z
#define BOX_REFLECT_Z
#endif
#ifdef SHEARING_BOX
#define BOX_SHEARING SHEARING_BOX
#endif
#ifdef SHEARING_BOX_Q
#define BOX_SHEARING_Q SHEARING_BOX_Q
#endif
#ifdef ANALYTIC_GRAVITY
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(ANALYTIC_GRAVITY)
#define GRAVITY_ANALYTIC ANALYTIC_GRAVITY
#else
#define GRAVITY_ANALYTIC
#endif
#endif
#ifdef NOGRAVITY
#define SELFGRAVITY_OFF
#endif


#if defined(EOS_ELASTIC)
#if !defined(DISABLE_SURFACE_VOLCORR) && !defined(HYDRO_KERNEL_SURFACE_VOLCORR)
#define HYDRO_KERNEL_SURFACE_VOLCORR
#endif
#if !defined(DISABLE_EXPLICIT_VOLUME_INTEGRATION) && !defined(HYDRO_EXPLICITLY_INTEGRATE_VOLUME)
#define HYDRO_EXPLICITLY_INTEGRATE_VOLUME
#endif
#endif



#include "eos/eos.h"


#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(DM_FUZZY)
#define AGS_FACE_CALCULATION_IS_ACTIVE
#endif


#if defined(GRAIN_COLLISIONS)
#define DM_SIDM 8 /* use the SIDM module to handle scattering of otherwise-collisionless particles against each other -- set to Particle Type=3 here */
#endif
#if defined(PIC_MHD)
#define PIC_MHD_NEW_RSOL_METHOD /* prefer new method for dealing with RSOL, should make simulations easier if done correctly */
#ifdef GRAIN_FLUID
#else
#define GRAIN_FLUID
#endif
#ifndef GRAIN_LORENTZFORCE
#define GRAIN_LORENTZFORCE
#endif
#ifndef GRAIN_BACKREACTION
#define GRAIN_BACKREACTION
#endif
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(DM_FUZZY) || defined(AGS_FACE_CALCULATION_IS_ACTIVE) || defined(DM_SIDM)
#define AGS_HSML_CALCULATION_IS_ACTIVE
#endif







#ifdef PROTECT_FROZEN_FIRE
#define GALSF_USE_SNE_ONELOOP_SCHEME // set to use the 'base' FIRE-2 SNe coupling. if commented out, will user newer version that more accurately manages the injected energy with neighbors moving to inject a specific target
#endif

#ifdef GALSF_SFR_CRITERION // flag for pure cross-compatibility [identical functionality, just ease-of-use for galaxy simulators here]
#define SINGLE_STAR_SINK_FORMATION GALSF_SFR_CRITERION
#endif

#ifdef COSMIC_RAY_FLUID
#define GAMMA_COSMICRAY(k) (4.0/3.0)
#ifndef CRFLUID_DIFFUSION_MODEL
#define CRFLUID_DIFFUSION_MODEL 0
#endif
#ifndef N_CR_PARTICLE_BINS
#define N_CR_PARTICLE_BINS 1
#endif
#if (N_CR_PARTICLE_BINS > 2) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
#endif
#endif



#if defined(COOL_GRACKLE)
#if !defined(COOLING)
#define COOLING
#endif
#include <grackle.h>
#endif

#ifdef CHIMES
#include "./cooling/chimes/chimes_proto.h"
extern struct gasVariables *ChimesGasVars;
extern struct globalVariables ChimesGlobalVars;
extern char ChimesDataPath[256];
extern char ChimesEqAbundanceTable[196];
extern char ChimesPhotoIonTable[196];
extern double chimes_rad_field_norm_factor;
extern double shielding_length_factor;
extern double cr_rate;
extern int ChimesEqmMode;
extern int ChimesUVBMode;
extern int ChimesInitIonState;
extern int Chimes_incl_full_output;
extern int N_chimes_full_output_freq;
#ifdef CHIMES_HII_REGIONS
#endif
#ifdef CHIMES_STELLAR_FLUXES
// The following defines the stellar age bins that we will use to define the UV spectra from stars used in CHIMES.
#define CHIMES_LOCAL_UV_NBINS 8
#define CHIMES_LOCAL_UV_AGE_LOW 0.0
#define CHIMES_LOCAL_UV_DELTA_AGE_LOW 0.2
#define CHIMES_LOCAL_UV_AGE_MID 1.0
#define CHIMES_LOCAL_UV_DELTA_AGE_HI 1.0
#endif
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
#endif // CHIMES




#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL) /* options for hybrid/combined FIRE+STARFORGE simulations */
#define GALSF_SFR_IMF_SAMPLING           /* use discrete sampling of 'number of O-stars' so we can handle the intermediate-mass regime in at least a simple approximate manner */
#define COOLING              /* only physical if include cooling for both sides, using same cooling functions */
#define MAGNETIC             /* enable MHD, important for systems here */
#define CONDUCTION           /* enable conduction */
#define CONDUCTION_SPITZER   /* compute proper coefficients and anisotropy for conduction */
#define VISCOSITY            /* enable viscosity */
#define VISCOSITY_BRAGINSKII /* compute proper coefficients and anisotropy for viscosity */
#define SINGLE_STAR_FB_JETS  /* enable jets from protostars */
#define SINGLE_STAR_FB_WINDS /* enable continuous mass-loss feedback - will also enable ssp mass-loss */
#define SINGLE_STAR_FB_SNE   /* enable SNe feedback - will also enable ssp mechanical feedback */
#define SINGLE_STAR_FB_RAD   /* enable RHD feedback */
#define RT_COMOVING          /* significantly more stable and accurate formulation given the structure of the problem and method we use */
#define RT_SOURCES (16+32)   /* need to allow -both- ssp-particles and single-star particles to emit */
#define RT_SPEEDOFLIGHT_REDUCTION (0.01)   /* for many problems on these scales, need much larger RSOL than default starforge values (dynamical velocities are big, without this they will severely lag behind) */
#define ADAPTIVE_TREEFORCE_UPDATE (0.0625) /* rough typical value we use for ensuring stability */
#endif // closes hybrid FIRE+STARFORGE model settings







#ifdef SINGLE_STAR_SINK_DYNAMICS
#define GALSF // top-level switch needed to enable various frameworks
#define METALS  // metals should be active for stellar return
#define BLACK_HOLES // need to have black holes active since these are our sink particles
#define BH_INTERACT_ON_GAS_TIMESTEP // BH-gas interactions (feedback and accretion) occur with frequency set by the gas timestep
#define BH_CALC_DISTANCES // calculate distance to nearest sink in gravity tree


#ifdef SINGLE_STAR_ACCRETION // figure out flags needed for the chosen sink accretion model
#define BH_SWALLOWGAS // need to swallow gas [part of sink model]
#ifndef BH_ALPHADISK_ACCRETION
#define BH_ALPHADISK_ACCRETION (2.) // all models will use a 'reservoir' of some kind to smooth out accretion rates (and represent unresolved disk)
#endif
#if (SINGLE_STAR_ACCRETION <= 8)
#define BH_GRAVACCRETION (SINGLE_STAR_ACCRETION) // use one of these pre-built accretion models
#endif
#if (SINGLE_STAR_ACCRETION == 9)
#define BH_BONDI 0 // use 'normal' Bondi-Hoyle accretion rate
#endif
#if (SINGLE_STAR_ACCRETION == 10)
#define BH_BONDI 1 // use Bondi rate ignoring local relative velocities
#endif
#if (SINGLE_STAR_ACCRETION == 11)
#define BH_GRAVCAPTURE_GAS // use gravitational capture swallow criterion for resolved gravitational capture
#endif
#if (SINGLE_STAR_ACCRETION == 12)
#define BH_GRAVCAPTURE_GAS
#define BH_GRAVCAPTURE_FIXEDSINKRADIUS // modify grav capture to Bate-style, fixed (in time) sink radius based on SF neighbor distance, plus angular momentum criterion
#endif
#endif

#if (defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_RT_HEATING) || defined(SINGLE_STAR_FB_SNE) || defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID) || defined(SINGLE_STAR_FB_LOCAL_RP)) && defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_SINK_DYNAMICS)
#define SINGLE_STAR_FB_TIMESTEPLIMIT // general flag indicating feedback is on
#endif

#if defined(SINGLE_STAR_FB_RT_HEATING) && !(defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID))
#define BH_PHOTONMOMENTUM // enable BHs within the FIRE-RT framework. 
#define RT_DISABLE_RAD_PRESSURE
#endif

#if (defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_SNE)) && !defined(GALSF_FB_MECHANICAL)
#define GALSF_FB_MECHANICAL // we will use the mechanical wind module for low mass loss rate stars (spawning leads to issues). enable regardless if either the winds or sne module is active
#define GALSF_USE_SNE_ONELOOP_SCHEME
#endif


#ifdef SINGLE_STAR_FB_SNE
#if !(CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(SINGLE_STAR_FB_SNE)) /* no numerical value is set, so set one as our 'default' */
#undef SINGLE_STAR_FB_SNE
#define SINGLE_STAR_FB_SNE 1 // fraction of the SN energy in the kinetic energy of particles vs internal
#endif
#endif

#if defined(SINGLE_STAR_FB_JETS) || ((defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_SNE)) && defined(FLAG_NOT_IN_PUBLIC_CODE))
#define BH_WIND_SPAWN (2) // leverage the BHFB model already developed within the FIRE-BHs framework. gives accurate launching of arbitrarily-structured jets.
#define MAINTAIN_TREE_IN_REARRANGE // don't rebuild the domains/tree every time a particle is spawned - salvage the existing one by redirecting pointers as needed
#endif

#if (defined(SINGLE_STAR_FB_SNE) || defined(SINGLE_STAR_FB_WINDS)) && !defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
#define SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT 2 // determines the maximum number of ejecta particles spawned per timestep, see below
#endif
#if defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
#define SINGLE_STAR_FB_SNE_N_EJECTA (4*(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)*((SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)+1)) // maximum number of ejecta cells spawned per timestep - follows from tiling rules for rays from RHD-direct-ray method
#endif

#if defined(SINGLE_STAR_FB_LOCAL_RP) // use standard angle-weighted local coupling to impart photon momentum from stars
#if !defined(BH_PHOTONMOMENTUM)
#define BH_PHOTONMOMENTUM
#endif
#if !defined(RT_DISABLE_RAD_PRESSURE)
#define RT_DISABLE_RAD_PRESSURE // we only want the local short-ranged photon momentum, since SF sims can easily get into the badly non-photon-conserving limit where LEBRON fluxes are less accurate
#endif
#endif

#if defined(COOLING) && !defined(COOL_GRACKLE) // if not using grackle modules, need to make sure appropriate cooling is enabled
#ifndef COOL_LOW_TEMPERATURES
#define COOL_LOW_TEMPERATURES // make sure low-temperature cooling is enabled!
#endif
#ifndef COOL_METAL_LINES_BY_SPECIES
#define COOL_METAL_LINES_BY_SPECIES // metal-based cooling enabled
#endif
#define OUTPUT_TEMPERATURE
#endif

#endif // SINGLE_STAR_SINK_DYNAMICS


#if (SINGLE_STAR_SINK_FORMATION & 1) || (SINGLE_STAR_SINK_FORMATION & 2048) // figure out flags needed for the chosen sink formation model [note these CAN be used even if single-star top-level flag is off, as additional SF/sink formation criteria for e.g. GALSF sims]
#if (SINGLE_STAR_SINK_FORMATION & 2048)
#define GALSF_SFR_VIRIAL_SF_CRITERION 2
#else
#define GALSF_SFR_VIRIAL_SF_CRITERION 1
#endif
#endif
#if (SINGLE_STAR_SINK_FORMATION & 16)
#ifndef SINGLE_STAR_TIMESTEPPING
#define SINGLE_STAR_TIMESTEPPING 0
#endif
#endif
#if (SINGLE_STAR_SINK_FORMATION & 32)
#define GALSF_SFR_TIDAL_HILL_CRITERION
#endif


#ifdef GRAVITY_ACCURATE_FEWBODY_INTEGRATION /* utility flag to enable a few different extra-conservative time-integration flags for gravity */
#if !defined(GRAVITY_HYBRID_OPENING_CRIT)
#define GRAVITY_HYBRID_OPENING_CRIT // use both Barnes-Hut + relative tree opening criterion
#endif
#if !defined(STOP_WHEN_BELOW_MINTIMESTEP)
#define STOP_WHEN_BELOW_MINTIMESTEP // stop when below min timestep to prevent bad timestepping
#endif
//#define RANDOMIZE_GRAVTREE /* move the top tree node around randomly so that treeforce errors are not correlated between one treebuild and another */
#define TIDAL_TIMESTEP_CRITERION // use tidal tensor timestep criterion
#endif
#ifdef HERMITE_INTEGRATION
#define COMPUTE_JERK_IN_GRAVTREE /* needs to be computed in order to do the Hermite integration */
#ifndef TIDAL_TIMESTEP_CRITERION
#define TIDAL_TIMESTEP_CRITERION // use tidal tensor timestep criterion -- otherwise won't effectively leverage the Hermite integrator timesteps
#endif
#endif

#ifdef ADAPTIVE_TREEFORCE_UPDATE // instead of going into the tree every timestep, only update gravity with a frequency set by this fraction of dynamical timescale (default for gas only)
#ifndef TIDAL_TIMESTEP_CRITERION 
#define TIDAL_TIMESTEP_CRITERION // need this to estimate the dynamical time
#endif    
#endif    


#if (SINGLE_STAR_TIMESTEPPING > 0) /* if single-star timestepping is on, need to make sure the binary-identification flag is active */
#ifndef SINGLE_STAR_FIND_BINARIES
#define SINGLE_STAR_FIND_BINARIES
#endif
#endif

#ifdef MHD_CONSTRAINED_GRADIENT
/* make sure mid-point gradient calculation for cleaning terms is enabled */
#ifndef MHD_CONSTRAINED_GRADIENT_MIDPOINT
#define MHD_CONSTRAINED_GRADIENT_MIDPOINT
#endif
#endif
/* these are tolerances for the slope-limiters. we define them here, because the gradient constraint routine needs to
    be sure to use the -same- values in both the gradients and reimann solver routines */
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


/* force 'parent' or 'top-level' flags to be enabled for the appropriate methods, if we have enabled something using those methods */


/* ----- giant block of options for RHD modules ------ */

/* options for FIRE RT method */
#if defined(RT_LEBRON)
#define RT_USE_GRAVTREE // use gravity tree for flux propagation
#define RT_USE_GRAVTREE_SAVE_RAD_ENERGY
#define RADTRANSFER // for cross-compatibility reasons, if the FIRE version is not on, need RADTRANSFER flag also enabled
#define RT_USE_GRAVTREE_SAVE_RAD_FLUX
#endif

/* check whether we want to use the implicit solver [only usable for very special cases, not recommended] */
#if defined(RT_DIFFUSION_IMPLICIT) && (defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION)) // only modules the implicit solver works with
#define RT_DIFFUSION_CG // use our implicit solver [will crash with any other modules, hence checking this before the others below]
#endif

/* options for FLD or OTVET or M1 or Ray/Rad_Intensity modules */
#if defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID)
#ifndef RADTRANSFER
#define RADTRANSFER // RADTRANSFER is ON, obviously
#endif
#define RT_SOURCE_INJECTION // need source injection enabled to define emissivity
#if !defined(RT_DIFFUSION_CG)
#define RT_SOLVER_EXPLICIT // default to explicit solutions (much more accurate/flexible)
#endif
#endif /* end of options for our general RHD methods */

/* OTVET-specific options [uses the gravity tree to calculate the Eddington tensor] */
#if defined(RT_OTVET)
#define RT_USE_GRAVTREE // use gravity tree for Eddington tensor
#ifndef RT_SEPARATELY_TRACK_LUMPOS
#define RT_SEPARATELY_TRACK_LUMPOS // and be sure to track luminosity locations
#endif
#endif /* end of otvet-specific options */

/* M1-specific options [make sure to add the flux moment */
#if defined(RT_M1)
#define RT_EVOLVE_FLUX // evolve flux moment [not just energy moment assumed by FLD/OTVET]
#endif

/* options for direct/exact Jiang et al. method for direct evolution on an intensity grid */
#if defined(RT_LOCALRAYGRID)
#define RT_EVOLVE_INTENSITIES // evolve the intensities explicitly
#define N_RT_INTENSITY_BINS (4*(RT_LOCALRAYGRID)*((RT_LOCALRAYGRID)+1)) // define number of directional bins, used throughout
#define RT_INTENSITY_BINS_DOMEGA (4.*M_PI/((double)N_RT_INTENSITY_BINS)) // normalization coefficient (for convenience defined here)
#endif

/* check if we are -explicitly- evolving the radiation energy density [0th moment], in which case we need to carry time-derivatives of the field */
#if defined(RT_SOLVER_EXPLICIT) && !defined(RT_EVOLVE_INTENSITIES) // only needed if we are -not- evolving intensities and -are- solving explicitly
#define RT_EVOLVE_ENERGY
#if !defined(RT_EVOLVE_FLUX) && !defined(RT_DISABLE_FLUXLIMITER) // evolving energy explicitly but not flux, flux-limiting is not disabled
#define RT_FLUXLIMITER // default to include flux-limiter under these conditions
#endif
#endif

/* enable radiation pressure forces unless they have been explicitly disabled */
#if defined(RADTRANSFER) && !defined(RT_DISABLE_RAD_PRESSURE) && !defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
#define RT_RAD_PRESSURE_FORCES
#endif

#ifdef RT_SOURCE_INJECTION
#if defined(GALSF) && !defined(RT_INJECT_PHOTONS_DISCRETELY)
#define RT_INJECT_PHOTONS_DISCRETELY // modules will not work correctly with differential timestepping with point sources without discrete injection
#endif
#if defined(RT_INJECT_PHOTONS_DISCRETELY) && defined(RT_RAD_PRESSURE_FORCES) && (defined(RT_ENABLE_R15_GRADIENTFIX) || defined(GALSF))
#define RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION // adds correction for un-resolved extinction which cannot generate photon momentum with M1, FLD, OTVET, etc.
#endif
#endif

/* check if we need to explicitly calculate gradients of the radiation pressure tensor for the diffusive step */
#if (defined(RT_FLUXLIMITER) || defined(RT_RAD_PRESSURE_FORCES) || defined(RT_SOLVER_EXPLICIT)) && !defined(RT_COMPGRAD_EDDINGTON_TENSOR) //&& !defined(RT_EVOLVE_FLUX) && !defined(RT_EVOLVE_INTENSITIES))
#define RT_COMPGRAD_EDDINGTON_TENSOR
#endif

/* enable appropriate chemistry flags if we are using the photoionization modules */
#if defined(RT_CHEM_PHOTOION)
#if (RT_CHEM_PHOTOION > 1)
/* enables multi-frequency radiation transport for ionizing photons. Integration variable is the ionising intensity J_nu */
#define RT_CHEM_PHOTOION_HE
#define RT_PHOTOION_MULTIFREQUENCY // if using He-ionization, default to multi-frequency RT [otherwise doesn't make sense] //
#endif
#endif

/* enable appropriate flags for X-ray sub-modules */
#if defined(RT_XRAY)
#if (RT_XRAY == 1)
#define RT_SOFT_XRAY
#endif
#if (RT_XRAY == 2)
#define RT_HARD_XRAY
#endif
#if (RT_XRAY == 3)
#define RT_SOFT_XRAY
#define RT_HARD_XRAY
#endif
#endif

/* default to speed-of-light equal to actual speed-of-light, and stars as photo-ionizing sources */
#ifndef RT_SPEEDOFLIGHT_REDUCTION
#define RT_SPEEDOFLIGHT_REDUCTION (1.0)
#endif
#ifndef RT_SOURCES
#define RT_SOURCES 1+2+4+8+16+32 // default to allowing all types to act as sources //
#endif

/* cooling must be enabled for RT cooling to function */
#if defined(RT_COOLING_PHOTOHEATING_OLDFORMAT) && !defined(COOLING)
#define COOLING
#endif

#if !defined(RT_USE_GRAVTREE) && defined(RT_SELFGRAVITY_OFF) && !defined(SELFGRAVITY_OFF)
#define SELFGRAVITY_OFF // safely define SELFGRAVITY_OFF in this case, otherwise we act like there is gravity except in the final setting of accelerations
#endif

/* ----- end block of options for RHD modules ------ */


#if defined(GALSF) || defined(BLACK_HOLES) || defined(RADTRANSFER) || defined(OUTPUT_DENS_AROUND_STAR) || defined(CHIMES)
#define DO_DENSITY_AROUND_STAR_PARTICLES
#if !defined(ALLOW_IMBALANCED_GASPARTICLELOAD)
#define ALLOW_IMBALANCED_GASPARTICLELOAD
#endif
#endif
#if defined(GALSF_SFR_VIRIAL_SF_CRITERION)
#if (GALSF_SFR_VIRIAL_SF_CRITERION >= 5)
#define GALSF_SFR_TIDAL_HILL_CRITERION
#endif
#if (GALSF_SFR_VIRIAL_SF_CRITERION >= 2)
#define GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED
#endif
#endif


#if defined(BH_SWALLOWGAS)
#define BH_FOLLOW_ACCRETED_COM
#define BH_FOLLOW_ACCRETED_MOMENTUM
#if defined(SINGLE_STAR_SINK_DYNAMICS) || defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVCAPTURE_NONGAS)
#define BH_FOLLOW_ACCRETED_ANGMOM 0 // follow accreted AM just from explicit 'swallow' operations
#else
#define BH_FOLLOW_ACCRETED_ANGMOM 1 // follow accreted AM from 'swallowed' BH particles, and from continuous/smooth properties [mdot] of kernel gas near BH
#endif
#endif


//#define ENERGY_ENTROPY_SWITCH_IS_ACTIVE
/* this is a ryu+jones type energy/entropy switch. it can help with some problems, but can also generate significant
 errors in other types of problems. in general, even for pure hydro, this isn't recommended; use it for special problems if you know what you are doing. */




#ifdef MAGNETIC
/* recommended MHD switches -- only turn these off for de-bugging */
#define DIVBCLEANING_DEDNER         /* hyperbolic/parabolic div-cleaing (Dedner 2002), with TP improvements */
/* MHD switches specific to SPH MHD */
#ifdef HYDRO_SPH
#define SPH_TP12_ARTIFICIAL_RESISTIVITY   /* turns on magnetic dissipation ('artificial resistivity'): uses tricco switch =h*|gradB|/|B| */
#endif
#endif


#if defined(TURB_DIFF_ENERGY) || defined(TURB_DIFF_VELOCITY) || defined(TURB_DIFF_MASS) || defined(TURB_DIFF_METALS)
#define TURB_DIFFUSION /* top-level switch to calculate properties needed for scalar turbulent diffusion/mixing: must enable with any specific version */
#if defined(TURB_DIFF_VELOCITY) && !defined(VISCOSITY)
#define VISCOSITY
#endif
#if defined(TURB_DIFF_ENERGY) && !defined(CONDUCTION)
#define CONDUCTION
#endif
#endif


#if defined(OUTPUT_POTENTIAL) && !defined(EVALPOTENTIAL)
#define EVALPOTENTIAL
#endif

#if defined(BLACK_HOLES) && (defined(BH_REPOSITION_ON_POTMIN) || defined(BH_SEED_FROM_FOF))
#ifndef EVALPOTENTIAL
#define EVALPOTENTIAL
#endif
#if !defined(BH_DYNFRICTION) && (BH_REPOSITION_ON_POTMIN == 2)
#define BH_DYNFRICTION 1 // use for local damping of anomalous velocities wrt background medium //
#endif
#endif


#ifdef EVALPOTENTIAL
#ifndef COMPUTE_POTENTIAL_ENERGY
#define COMPUTE_POTENTIAL_ENERGY
#endif
#endif


#if defined(COOL_MOLECFRAC)
#if (COOL_MOLECFRAC == 6) && !defined(COOL_MOLECFRAC_NONEQM)
#define COOL_MOLECFRAC_NONEQM // estimate molecular fractions for thermochemistry+cooling with explicitly-evolved non-equilibirum H2 formation+destruction with clumping and self-shielding (Hopkins+2021, in prep)
#elif (COOL_MOLECFRAC == 5) && !defined(COOL_MOLECFRAC_LOCALEQM)
#define COOL_MOLECFRAC_LOCALEQM  // estimate molecular fractions for thermochemistry+cooling from local equilibrium H2 formation+destruction with clumping and self-shielding (Hopkins+2021, in prep)
#elif (COOL_MOLECFRAC == 4) && !defined(COOL_MOLECFRAC_KMT)
#define COOL_MOLECFRAC_KMT  // estimate f_H2 from approximate large-scale expressions from Krumholz, McKee, & Tumlinson (2009ApJ...693..216K). use the simpler Kumholz, McKee, & Tumlinson 2009 sub-grid model for molecular fractions in equilibrium, which is a function modeling spherical clouds of internally uniform properties exposed to incident radiation. Depends on column density, metallicity, and incident FUV field
#elif (COOL_MOLECFRAC == 3) && !defined(COOL_MOLECFRAC_GD)
#define COOL_MOLECFRAC_GD  // estimate f_H2 from approximate large-scale expressions from Gnedin & Draine (2014ApJ...795...37G). use the sub-grid final expression calibrated to ~60pc resolution simulations with equilibrium molecular chemistry and post-processing radiative transfer from Gnedin & Draine 2014 (Eqs. 5-7)
#elif (COOL_MOLECFRAC == 2) && !defined(COOL_MOLECFRAC_KG)
#define COOL_MOLECFRAC_KG  // estimate f_H2 with Krumholz & Gnedin 2010 fitting function, assuming simple scalings of radiation field, clumping, and other factors with basic gas properties so function only of surface density and metallicity, truncated at low values (or else it gives non-sensical answers)
#elif (COOL_MOLECFRAC == 1) && !defined(COOL_MOLECFRAC_GC)
#define COOL_MOLECFRAC_GC  // if none of the above is set, default to a wildly-oversimplified scaling set by fits to the temperature below which gas at a given density becomes molecular from cloud simulations in Glover+Clark 2012
#else
#define COOL_MOLECFRAC_GC // default if no value above set
#endif
#endif


#ifdef BOX_SHEARING
/* set default compile-time flags for the shearing-box (or shearing-sheet) boundaries */
/* shearing box boundaries: 1=r-z sheet (coordinates [0,1,2] = [r,z,phi]), 2=r-phi sheet [r,phi,z], 3=[r-phi-z] box */
#if (BOX_SHEARING==1)
#define BOX_SHEARING_PHI_COORDINATE 2
#else
#define BOX_SHEARING_PHI_COORDINATE 1
#endif
/* if the r-z or r-phi sheet is set, the code must be compiled in 2D mode */
#if (BOX_SHEARING==1) || (BOX_SHEARING==2)
#ifndef BOX_SPATIAL_DIMENSION
#define BOX_SPATIAL_DIMENSION 2
#endif
#endif
/* box must be periodic in this approximation */
#ifndef BOX_PERIODIC
#define BOX_PERIODIC
#endif
/* if not set, default to q=3/2 (q==-dlnOmega/dlnr, used for boundary and velocity corrections) */
#ifndef BOX_SHEARING_Q
#define BOX_SHEARING_Q (3.0/2.0)
#endif
/* set omega - usually we will default to always using time coordinates such that Omega = 1 at the box center */
#define BOX_SHEARING_OMEGA_BOX_CENTER 1.0
/* need analytic gravity on so we can add the appropriate source terms to the EOM */
#ifndef GRAVITY_ANALYTIC
#define GRAVITY_ANALYTIC
#endif
/* if self-gravity is on, we need to make sure the gravitational forces are not periodic. this is going to cause some errors at the x/y 'edges',
    but for now at least, the periodic gravity routines (particularly the FFT's involved) require a regular periodic map, they cannot handle the
    non-standard map that the shearing box represents. */
#ifndef GRAVITY_NOT_PERIODIC
#define GRAVITY_NOT_PERIODIC
#endif
#endif // BOX_SHEARING



#if defined(GRAVITY_ANALYTIC)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(GRAVITY_ANALYTIC)
#if (GRAVITY_ANALYTIC > 0)
#define GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE /* ok, analytic gravity is defined with a numerical value > 0, indicating we should use this flag */
#ifndef BH_CALC_DISTANCES
#define BH_CALC_DISTANCES
#endif
#endif
#endif
#endif

#if defined(EOS_SUBSTELLAR_ISM) || defined(COOL_MOLECFRAC_NONEQM)
#define EOS_GAMMA_VARIABLE
#endif

#if defined(EOS_GAMMA_VARIABLE)
#define GAMMA(i) (gamma_eos(i)) /*! use an actual function! */
#ifndef EOS_GENERAL
#define EOS_GENERAL /*! needs to be on for this to work */
#endif
#else
#define GAMMA(i) (EOS_GAMMA) /*! default to this being a universal constant */
#endif
#define GAMMA_DEFAULT (EOS_GAMMA)

#if defined(CONDUCTION) || defined(EOS_GENERAL)
#define DOGRAD_INTERNAL_ENERGY 1
#endif

#if defined(EOS_GENERAL)
#define DOGRAD_SOUNDSPEED 1
#endif






/*------- Things that are always recommended -------*/


#ifdef MPISENDRECV_CHECKSUM
#define MPI_Sendrecv MPI_Check_Sendrecv
#endif

#ifdef MPISENDRECV_SIZELIMIT
#define MPI_Sendrecv MPI_Sizelimited_Sendrecv
#endif

#include "tags.h"
#include <assert.h>


#ifdef MYSORT
#define MYSORT_DATAINDEX mysort_dataindex
#else // MYSORT
#define MYSORT_DATAINDEX qsort
#endif

#ifndef DISABLE_MEMORY_MANAGER // compiler specific data alignment hints: use only with memory manager as malloc'd memory is not sufficiently aligned
// (experimenting right now with removing this, as many compilers internal AVX optimizations appear to be doing marginally better, and can resolve crashes on some compilers)
#if defined(__xlC__) // XLC compiler
#define ALIGN(n) __attribute__((__aligned__(n)))
#elif defined(__GNUC__) // GNU compiler
#define ALIGN(n) __attribute__((__aligned__(n)))
#elif defined(__INTEL_COMPILER) // Intel Compiler
#define ALIGN(n) __declspec(align(n))
#endif
#endif
#ifndef ALIGN // Unknown Compiler or using default malloc
#define ALIGN(n)
#endif

#define ASSIGN_ADD(x,y,mode) (mode == 0 ? (x=y) : (x+=y))


#ifndef  GALSF_GENERATIONS
#define  GALSF_GENERATIONS     1	/*!< Number of star particles that may be created per gas particle */
#endif

#ifdef LONG_INTEGER_TIME
typedef  long long integertime;
static MPI_Datatype MPI_TYPE_TIME = MPI_LONG_LONG;
#define  TIMEBINS        39
#else
typedef  int integertime;
static MPI_Datatype MPI_TYPE_TIME = MPI_INT;
#define  TIMEBINS        29
#endif
#define  TIMEBASE        (((integertime) 1)<<TIMEBINS)  /* The simulated timespan is mapped onto the integer interval [0,TIMESPAN], where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds to 2^29 */
#define UNIT_INTEGERTIME_IN_PHYSICAL ((All.Timebase_interval/All.cf_hubble_a))
#define GET_INTEGERTIME_FROM_TIMEBIN(bin) ((bin ? (((integertime) 1) << bin) : 0))
#define GET_PHYSICAL_TIMESTEP_FROM_TIMEBIN(bin) ((GET_INTEGERTIME_FROM_TIMEBIN(bin) * UNIT_INTEGERTIME_IN_PHYSICAL))
#ifndef WAKEUP
#define GET_PARTICLE_INTEGERTIME(i) ((GET_INTEGERTIME_FROM_TIMEBIN(P[i].TimeBin)))
#else
#define GET_PARTICLE_INTEGERTIME(i) ((P[i].dt_step))
#endif
#define GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) ((GET_PARTICLE_INTEGERTIME(i) * UNIT_INTEGERTIME_IN_PHYSICAL))


#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
#define AGS_OUTPUTGRAVSOFT 1  /*! output softening to snapshots */
//#define AGS_OUTPUTZETA 1 /*! output correction zeta term to snapshots */
#endif

#ifdef NUCLEAR_NETWORK
#include "nuclear/nuclear_network.h"
#endif

#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
#define RT_BIN0 (-1)

#ifndef RT_CHEM_PHOTOION
#define RT_FREQ_BIN_H0 (RT_BIN0+0)
#else
#define RT_FREQ_BIN_H0 (RT_BIN0+1)
#endif

#ifndef RT_PHOTOION_MULTIFREQUENCY
#define RT_FREQ_BIN_He0 (RT_FREQ_BIN_H0+0)
#define RT_FREQ_BIN_He1 (RT_FREQ_BIN_He0+0)
#define RT_FREQ_BIN_He2 (RT_FREQ_BIN_He1+0)
#else
#define RT_FREQ_BIN_He0 (RT_FREQ_BIN_H0+1)
#define RT_FREQ_BIN_He1 (RT_FREQ_BIN_He0+1)
#define RT_FREQ_BIN_He2 (RT_FREQ_BIN_He1+1)
#endif

#define RT_FREQ_BIN_FIRE_UV (RT_FREQ_BIN_He2+0)
#define RT_FREQ_BIN_FIRE_OPT (RT_FREQ_BIN_FIRE_UV+0)
#define RT_FREQ_BIN_FIRE_IR (RT_FREQ_BIN_FIRE_OPT+0)

#ifndef RT_SOFT_XRAY
#define RT_FREQ_BIN_SOFT_XRAY (RT_FREQ_BIN_FIRE_IR+0)
#else
#define RT_FREQ_BIN_SOFT_XRAY (RT_FREQ_BIN_FIRE_IR+1)
#endif

#ifndef RT_HARD_XRAY
#define RT_FREQ_BIN_HARD_XRAY (RT_FREQ_BIN_SOFT_XRAY+0)
#else
#define RT_FREQ_BIN_HARD_XRAY (RT_FREQ_BIN_SOFT_XRAY+1)
#endif

#ifndef RT_PHOTOELECTRIC
#define RT_FREQ_BIN_PHOTOELECTRIC (RT_FREQ_BIN_HARD_XRAY+0)
#else
#define RT_FREQ_BIN_PHOTOELECTRIC (RT_FREQ_BIN_HARD_XRAY+1)
#endif

#ifndef RT_LYMAN_WERNER
#define RT_FREQ_BIN_LYMAN_WERNER (RT_FREQ_BIN_PHOTOELECTRIC+0)
#else
#define RT_FREQ_BIN_LYMAN_WERNER (RT_FREQ_BIN_PHOTOELECTRIC+1)
#endif

#ifndef RT_NUV
#define RT_FREQ_BIN_NUV (RT_FREQ_BIN_LYMAN_WERNER+0)
#else
#define RT_FREQ_BIN_NUV (RT_FREQ_BIN_LYMAN_WERNER+1)
#endif

#ifndef RT_OPTICAL_NIR
#define RT_FREQ_BIN_OPTICAL_NIR (RT_FREQ_BIN_NUV+0)
#else
#define RT_FREQ_BIN_OPTICAL_NIR (RT_FREQ_BIN_NUV+1)
#endif

#ifndef RT_FREEFREE
#define RT_FREQ_BIN_FREEFREE (RT_FREQ_BIN_OPTICAL_NIR+0)
#else
#define RT_FREQ_BIN_FREEFREE (RT_FREQ_BIN_OPTICAL_NIR+1)
#endif


#ifndef RT_GENERIC_USER_FREQ
#define RT_FREQ_BIN_GENERIC_USER_FREQ (RT_FREQ_BIN_FREEFREE+0)
#else
#define RT_FREQ_BIN_GENERIC_USER_FREQ (RT_FREQ_BIN_FREEFREE+1)
#endif



/* be sure to add all new wavebands to these lists, or else we will run into problems */
/* ALSO, the IR bin here should be the last bin: add additional bins ABOVE this line */
#ifndef RT_INFRARED
#define RT_FREQ_BIN_INFRARED (RT_FREQ_BIN_GENERIC_USER_FREQ+0)
#else
#define RT_FREQ_BIN_INFRARED (RT_FREQ_BIN_GENERIC_USER_FREQ+1)
#endif

#define N_RT_FREQ_BINS (RT_FREQ_BIN_INFRARED+1)

#endif // #if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)


#ifndef  MULTIPLEDOMAINS
#define  MULTIPLEDOMAINS     8
#endif

#ifndef  TOPNODEFACTOR
#define  TOPNODEFACTOR       4.0
#endif

#ifndef  GRAVCOSTLEVELS
#define  GRAVCOSTLEVELS      6
#endif

#define  NUMBER_OF_MEASUREMENTS_TO_RECORD  6  /* this is the number of past executions of a timebin that the reported average CPU-times average over */

#define  NODELISTLENGTH      8


#define EPSILON_FOR_TREERND_SUBNODE_SPLITTING (1.0e-4) /* define some number << 1; particles with less than this separation will trigger randomized sub-node splitting in the tree.
                                                            we set it to a global value here so that other sub-routines will know not to force particle separations below this */

#ifdef GALSF_SFR_IMF_VARIATION
#define N_IMF_FORMPROPS  13  /*!< formation properties of star particles to record for output */
#endif


typedef unsigned long long peanokey;


#define  BITS_PER_DIMENSION 21	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))

#define  BITS_PER_DIMENSION_SAVE_KEYS 10
#define  PEANOCELLS_SAVE_KEYS (((peanokey)1)<<(3*BITS_PER_DIMENSION_SAVE_KEYS))




#ifndef DISABLE_MEMORY_MANAGER
#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)
#else
#define  mymalloc(x, y)            malloc(y)
#define  mymalloc_movable(x, y, z) malloc(z)

#define  myrealloc(x, y)           realloc(x, y)
#define  myrealloc_movable(x, y)   realloc(x, y)

#define  myfree(x)                 free(x)
#define  myfree_movable(x)         free(x)

#define  report_memory_usage(x, y) printf("Memory manager disabled.\n")
#endif

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

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#if (defined(MAGNETIC) && !defined(COOLING)) || defined(EOS_ELASTIC)
#define  CONDITION_NUMBER_DANGER  1.0e7 /*!< condition number above which we will not trust matrix-based gradients */
#else
#define  CONDITION_NUMBER_DANGER  1.0e3 /*!< condition number above which we will not trust matrix-based gradients */
#endif

#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
#define  RNDTABLE 16384 /*!< this is arbitrary, but some power of 2 makes much easier */
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
#define  HUBBLE_H100_CGS    (3.2407789e-18)	/* in h/sec */
#define  ELECTRONVOLT_IN_ERGS (1.60217733e-12)

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
#define UNIT_FLUX_IN_HABING     (((UNIT_FLUX_IN_CGS)/1.6e-3))

#define U_TO_TEMP_UNITS         ((PROTONMASS_CGS/BOLTZMANN_CGS)*((UNIT_ENERGY_IN_CGS)/(UNIT_MASS_IN_CGS))) /* units to convert specific internal energy to temperature. needs to be multiplied by dimensionless factor=mean_molec_weight_in_amu*(gamma_eos-1) */
#ifndef C_LIGHT_CODE
#define C_LIGHT_CODE            ((C_LIGHT_CGS/UNIT_VEL_IN_CGS)) /* pure convenience function, speed-of-light in code units */
#endif
#define C_LIGHT_CODE_REDUCED    (((RT_SPEEDOFLIGHT_REDUCTION)*(C_LIGHT_CODE))) /* reduced speed-of-light in code units, again here as a convenience function */
#define H0_CGS                  ((All.HubbleParam*HUBBLE_H100_CGS)) /* actual value of H0 in cgs */
#define COSMIC_BARYON_DENSITY_CGS ((All.OmegaBaryon*(H0_CGS)*(H0_CGS)*(3./(8.*M_PI*GRAVITY_G_CGS))*All.cf_a3inv)) /* cosmic mean baryon density [scale-factor-dependent] in cgs units */



#ifdef RT_COMOVING
#define RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS (0) /* this prefactor goes in front of various terms which vanish in the comoving frame RHD equations */
#else
#define RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS ((C_LIGHT_CODE_REDUCED)/(C_LIGHT_CODE)) /* these terms in the mixed-frame equations need to be multiplied by c_reduced/c */
#endif




#ifdef COOL_MOLECFRAC_NONEQM
#ifndef OUTPUT_MOLECULAR_FRACTION
#define OUTPUT_MOLECULAR_FRACTION
#endif
#endif


#ifdef METALS
#ifdef GALSF_FB_FIRE_RPROCESS
#define NUM_RPROCESS_SPECIES (GALSF_FB_FIRE_RPROCESS)
#else
#define NUM_RPROCESS_SPECIES 0
#endif

#ifdef GALSF_FB_FIRE_AGE_TRACERS
#define NUM_AGE_TRACERS (GALSF_FB_FIRE_AGE_TRACERS)
#else
#define NUM_AGE_TRACERS 0
#endif

#ifdef STARFORGE_FEEDBACK_TRACERS
#define NUM_STARFORGE_FEEDBACK_TRACERS (STARFORGE_FEEDBACK_TRACERS)
#else
#define NUM_STARFORGE_FEEDBACK_TRACERS 0
#endif

#ifdef COOL_METAL_LINES_BY_SPECIES
#define NUM_LIVE_SPECIES_FOR_COOLTABLES 10
#else
#define NUM_LIVE_SPECIES_FOR_COOLTABLES 0
#endif

#define NUM_METAL_SPECIES (1+NUM_LIVE_SPECIES_FOR_COOLTABLES+NUM_RPROCESS_SPECIES+NUM_AGE_TRACERS+NUM_STARFORGE_FEEDBACK_TRACERS)
#endif // METALS //



#if defined(CRFLUID_M1)
#define CRFLUID_REDUCED_C_CODE(k) (CRFLUID_M1) // single-bin -- compiles to simply replace this macro with the M1 value, trivially
#endif // M1 cosmic rays
#if defined(CRFLUID_ALT_RSOL_FORM) && defined(CRFLUID_M1)
#define CosmicRayFluid_RSOL_Corrfac(k) (((CRFLUID_REDUCED_C_CODE(k))/(C_LIGHT_CODE))) // this needs to be defined after the code SOL for obvious reasons
#else
#define CosmicRayFluid_RSOL_Corrfac(k) (1.0) // this is always unity, macro is trivial
#endif


#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
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
#define MAXLEN_OUTPUTLIST 1201	/*!< maxmimum number of entries in output list */
#define DRIFT_TABLE_LENGTH 1000	/*!< length of the lookup table used to hold the drift and kick factors */
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

#ifndef GDE_TYPES
#define GDE_TYPES 2
#endif

#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


#ifndef DOUBLEPRECISION     /* default is single-precision */
typedef float  MyFloat;
typedef float  MyDouble;
#else                        /* everything double-precision */
typedef double  MyFloat;
typedef double  MyDouble;
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif
#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif

#ifdef OUTPUT_POSITIONS_IN_DOUBLE
typedef double MyOutputPosFloat;
#else
typedef MyOutputFloat MyOutputPosFloat;
#endif
#ifdef INPUT_POSITIONS_IN_DOUBLE
typedef double MyInputPosFloat;
#else
typedef MyInputFloat MyInputPosFloat;
#endif


struct unbind_data
{
  int index;
};


#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif

#define FLT(x) (x)
typedef MyFloat MyLongDouble;
typedef MyDouble MyBigFloat;

#define GDE_ABS(x) (fabs(x))
#define GDE_SQRT(x) (sqrt(x))
#define GDE_LOG(x) (log(x))


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
#define CPU_DENSWAIT       12
#define CPU_DENSCOMM       13
#define CPU_DENSMISC       14
#define CPU_HYDCOMPUTE     15
#define CPU_HYDWAIT        16
#define CPU_HYDCOMM        17
#define CPU_HYDMISC        18
#define CPU_DRIFT          19
#define CPU_TIMELINE       20
#define CPU_POTENTIAL      21
#define CPU_MESH           22
#define CPU_PEANO          23
#define CPU_COOLINGSFR     24
#define CPU_SNAPSHOT       25
#define CPU_FOF            26
#define CPU_BLACKHOLES     27
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
#define CPU_DUMMY00       47
#define CPU_DUMMY01       48
#define CPU_DUMMY02       49
#define CPU_DUMMY03       50
#define CPU_DUMMY04       51
#define CPU_DUMMY05       52
#define CPU_DUMMY06       53
#define CPU_DUMMY07       54
#define CPU_DUMMY08       55
#define CPU_DUMMY09       56
#define CPU_DUMMY10       57

#define CPU_PARTS          58  /* this gives the number of parts above (must be last) */

#define CPU_STRING_LEN 120

#if (BOX_SPATIAL_DIMENSION==1) || defined(ONEDIM)
#define NUMDIMS 1           /* define number of dimensions and volume normalization */
#define NORM_COEFF 2.0
#elif (BOX_SPATIAL_DIMENSION==2) || defined(TWODIMS)
#define NUMDIMS 2
#define NORM_COEFF M_PI
#else
#define NORM_COEFF 4.188790204786  /* 4pi/3 */
#define NUMDIMS 3
#endif


#define PPP P
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
#define PPPZ P
#else
#define PPPZ SphP
#endif

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

#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z)
#define BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE 1 /* flag to let the code know to use everything below */
extern short int special_boundary_condition_xyz_def_reflect[3];
extern short int special_boundary_condition_xyz_def_outflow[3];
#define BOX_VALUE_FOR_NOTHING_SPECIAL_BOUNDARY_ 20 /* define a dummy value we won't have the user set for reference below */
#endif



/****************************************************************************************************************************/
/* Here we define the box-wrapping macros NEAREST_XYZ and NGB_PERIODIC_BOX_LONG_X,NGB_PERIODIC_BOX_LONG_Y,NGB_PERIODIC_BOX_LONG_Z.
 *   The inputs to these functions are (dx_position, dy_position, dz_position, sign), where
 *     'sign' = -1 if dx_position = x_test_point - x_reference (reference = particle from which we are doing a calculation),
 *     'sign' = +1 if dx_position = x_reference - x_test_point
 *
 *   For non-periodic cases these functions are trivial (do nothing, or just take absolute values).
 *
 *   For standard periodic cases it will wrap in each dimension, allowing for a different box length in X/Y/Z.
 *      here the "sign" term is irrelevant. Also NGB_PERIODIC_BOX_LONG_X, NGB_PERIODIC_BOX_LONG_Y, NGB_PERIODIC_BOX_LONG_Z will each
 *      compile to only use the x,y, or z information, but all four inputs are required for the sake of completeness
 *      and consistency.
 *
 *   The reason for the added complexity is for shearing boxes. In this case, the Y(phi)-coordinate for particles being
 *      wrapped in the X(r)-direction must be modified by a time-dependent term. It also matters for the sign of that
 *      term "which side" of the box we are wrapping across (i.e. does the 'virtual particle' -- the test point which is
 *      not the particle for which we are currently calculating forces, etc -- lie on the '-x' side or the '+x' side)
 *      (note after all that: if very careful, sign -cancels- within the respective convention, for the type of wrapping below)
 */
/****************************************************************************************************************************/

#if defined(BOX_PERIODIC) && !(defined(BOX_REFLECT_X) || defined(BOX_OUTFLOW_X)) // x-axis is periodic
#define TMP_WRAP_X_S(x,y,z,sign) (x=((x)>boxHalf_X)?((x)-boxSize_X):(((x)<-boxHalf_X)?((x)+boxSize_X):(x))) /* normal (signed) periodic wrap */
#define NGB_PERIODIC_BOX_LONG_X(x,y,z,sign) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp) /* absolute value of normal periodic wrap */
#else // x-axis is non-periodic
#define TMP_WRAP_X_S(x,y,z,sign) /* this is an empty macro: nothing will happen to the variables input here */
#define NGB_PERIODIC_BOX_LONG_X(x,y,z,sign) (fabs(x)) /* simple absolute value */
#endif

#if defined(BOX_PERIODIC) && !(defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_Z)) // z-axis is periodic
#define TMP_WRAP_Z_S(x,y,z,sign) (z=((z)>boxHalf_Z)?((z)-boxSize_Z):(((z)<-boxHalf_Z)?((z)+boxSize_Z):(z))) /* normal (signed) periodic wrap */
#define NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign) (xtmp=fabs(z),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp) /* absolute value of normal periodic wrap */
#else // z-axis is non-periodic
#define TMP_WRAP_Z_S(x,y,z,sign) /* this is an empty macro: nothing will happen to the variables input here */
#define NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign) (fabs(z)) /* simple absolute value */
#endif

#if defined(BOX_PERIODIC) && !(defined(BOX_REFLECT_Y) || defined(BOX_OUTFLOW_Y)) // y-axis is periodic
#if (BOX_SHEARING > 1) // Shearing Periodic Box:: in this case, we have a shearing box with the '1' coordinate being phi, so there is a periodic extra wrap

#define TMP_WRAP_Y_S(x,y,z,sign) (\
y += Shearing_Box_Pos_Offset * (((x)>boxHalf_X)?(1):(((x)<-boxHalf_X)?(-1):(0))),\
y = ((y)>boxSize_Y)?((y)-boxSize_Y):(((y)<-boxSize_Y)?((y)+boxSize_Y):(y)),\
y=((y)>boxHalf_Y)?((y)-boxSize_Y):(((y)<-boxHalf_Y)?((y)+boxSize_Y):(y))) /* shear-periodic wrap in y, accounting for the position offset needed for azimuthal wrap off the radial axis */

#define NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) (\
xtmp = y + Shearing_Box_Pos_Offset * (((x)>boxHalf_X)?(1):(((x)<-boxHalf_X)?(-1):(0))),\
xtmp = fabs(((xtmp)>boxSize_Y)?((xtmp)-boxSize_Y):(((xtmp)<-boxSize_Y)?((xtmp)+boxSize_Y):(xtmp))),\
(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp) /* shear periodic wrap in y, accounting for the position offset needed for azimuthal wrap off the radial axis: absolute value here */

#else // 'normal' periodic y-axis, nothing special
#define TMP_WRAP_Y_S(x,y,z,sign) (y=((y)>boxHalf_Y)?((y)-boxSize_Y):(((y)<-boxHalf_Y)?((y)+boxSize_Y):(y))) /* normal (signed) periodic wrap */
#define NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) (xtmp=fabs(y),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp) /* absolute value of normal periodic wrap */
#endif
#else // y-axis is non-periodic
#define TMP_WRAP_Y_S(x,y,z,sign) /* this is an empty macro: nothing will happen to the variables input here */
#define NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) (fabs(y)) /* simple absolute value */
#endif

#define NEAREST_XYZ(x,y,z,sign) {\
TMP_WRAP_Y_S(x,y,z,sign);\
TMP_WRAP_X_S(x,y,z,sign);\
TMP_WRAP_Z_S(x,y,z,sign);} /* note the ORDER MATTERS here for shearing boxes: Y-wrap must precede x/z wrap to allow correct re-assignment. collect the box-wrapping terms into one function here */



#if 0 /* below is the old code block for this, replaced with the more flexible structures above, retained for de-bugging for now */

#ifdef BOX_PERIODIC
#define NGB_PERIODIC_BOX_LONG_X(x,y,z,sign) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp) // normal periodic wrap //
#define NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign) (xtmp=fabs(z),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp) // normal periodic wrap //

#if (BOX_SHEARING > 1)
/* Shearing Periodic Box::
    in this case, we have a shearing box with the '1' coordinate being phi, so there is a periodic extra wrap */
#define NEAREST_XYZ(x,y,z,sign) (\
y += Shearing_Box_Pos_Offset * (((x)>boxHalf_X)?(1):(((x)<-boxHalf_X)?(-1):(0))),\
y = ((y)>boxSize_Y)?((y)-boxSize_Y):(((y)<-boxSize_Y)?((y)+boxSize_Y):(y)),\
y=((y)>boxHalf_Y)?((y)-boxSize_Y):(((y)<-boxHalf_Y)?((y)+boxSize_Y):(y)),\
x=((x)>boxHalf_X)?((x)-boxSize_X):(((x)<-boxHalf_X)?((x)+boxSize_X):(x)),\
z=((z)>boxHalf_Z)?((z)-boxSize_Z):(((z)<-boxHalf_Z)?((z)+boxSize_Z):(z)))

#define NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) (\
xtmp = y + Shearing_Box_Pos_Offset * (((x)>boxHalf_X)?(1):(((x)<-boxHalf_X)?(-1):(0))),\
xtmp = fabs(((xtmp)>boxSize_Y)?((xtmp)-boxSize_Y):(((xtmp)<-boxSize_Y)?((xtmp)+boxSize_Y):(xtmp))),\
(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp)

#else
/* Standard Periodic Box::
    this box-wraps all three (x,y,z) separation variables when taking position differences */
#define NEAREST_XYZ(x,y,z,sign) (\
x=((x)>boxHalf_X)?((x)-boxSize_X):(((x)<-boxHalf_X)?((x)+boxSize_X):(x)),\
y=((y)>boxHalf_Y)?((y)-boxSize_Y):(((y)<-boxHalf_Y)?((y)+boxSize_Y):(y)),\
z=((z)>boxHalf_Z)?((z)-boxSize_Z):(((z)<-boxHalf_Z)?((z)+boxSize_Z):(z)))
#define NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) (xtmp=fabs(y),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp) // normal periodic wrap //

#endif

#else
/* Non-periodic box:: */
#define NEAREST_XYZ(x,y,z,sign) /* this is an empty macro: nothing will happen to the variables input here */
#define NGB_PERIODIC_BOX_LONG_X(x,y,z,sign) (fabs(x))
#define NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) (fabs(y))
#define NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign) (fabs(z))
#endif

#endif // 0


/* this function, like the NEAREST and NGB_PERIODIC functions above, does -velocity wrapping- for periodic boundary
    conditions. this is currently only relevant for shearing boxes, where the box ends in the '0' axis direction have
    systematically different (shear-periodic instead of periodic) velocities associated, so the box needs to be able to
    know how to wrap them. this takes the vector of positions of particle "i" pos_i (the particle "seeing" particle j),
    particle j position pos_j, the velocity difference vector dv_ij=v_i-v_j. last  dv_sign_flipped = 1 if dv_ij=v_i-v_j,
    but dv_sign_flipped=-1 if dv_ij=v_j-v_i (flipped from normal order) */
#ifdef BOX_SHEARING
#define NGB_SHEARBOX_BOUNDARY_VELCORR_(pos_i,pos_j,dv_ij,dv_sign_flipped) (dv_ij[BOX_SHEARING_PHI_COORDINATE] += dv_sign_flipped*Shearing_Box_Vel_Offset * ((pos_i[0]-pos_j[0]>boxHalf_X)?(1):((pos_i[0]-pos_j[0]<-boxHalf_X)?(-1):(0))))
#else
#define NGB_SHEARBOX_BOUNDARY_VELCORR_(pos_i,pos_j,dv_ij,dv_sign_flipped)
#endif


#define FACT1 0.366025403785	/* FACT1 = 0.5 * (sqrt(3)-1) */
#define FACT2 0.86602540        /* FACT2 = 0.5 * sqrt(3) */



/*****************************************************************/
/*  Utility functions used for printing status, warning, endruns */
/*****************************************************************/

#define terminate(x) {char termbuf[2000]; sprintf(termbuf, "TERMINATE issued on task=%d, function '%s()', file '%s', line %d: '%s'\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, x); fflush(stdout); printf("%s", termbuf); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, 1); exit(0);}
#define endrun(x) {if(x==0) {MPI_Finalize(); exit(0);} else {char termbuf[2000]; sprintf(termbuf, "ENDRUN issued on task=%d, function '%s()', file '%s', line %d: error level %d\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, x); fflush(stdout); printf("%s", termbuf); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, x); exit(0);}}
#define PRINT_WARNING(...) {char termbuf1[1000], termbuf2[1000]; sprintf(termbuf1, "WARNING issued on task=%d, function %s(), file %s, line %d", ThisTask, __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); fflush(stdout); printf("%s: %s\n", termbuf1, termbuf2); fflush(stdout);}
#ifdef IO_REDUCED_MODE
#define PRINT_STATUS(...) {if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) {if(ThisTask==0) {fflush(stdout); printf( __VA_ARGS__ ); printf("\n"); fflush(stdout);}}}
#else
#define PRINT_STATUS(...) {if(ThisTask==0) {fflush(stdout); printf( __VA_ARGS__ ); printf("\n"); fflush(stdout);}}
#endif

#define MACRO_NAME_CONCATENATE(A, B) MACRO_NAME_CONCATENATE_(A, B)
#define MACRO_NAME_CONCATENATE_(A, B) A##B


/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

extern int FirstActiveParticle;
extern int *NextActiveParticle;
extern unsigned char *ProcessedFlag;
extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountSph[TIMEBINS];
extern int TimeBinActive[TIMEBINS];
extern int FirstInTimeBin[TIMEBINS];
extern int LastInTimeBin[TIMEBINS];
extern int *NextInTimeBin;
extern int *PrevInTimeBin;
#ifdef GALSF
extern double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
#define BH_COUNTPROGS
/* carries a counter for each BH that gives the total number of seeds that merged into it */
#define BH_ENFORCE_EDDINGTON_LIMIT
/* put a hard limit on the maximum accretion rate (set BlackHoleEddingtonFactor>>1 to allow super-eddington) */
extern double TimeBin_BH_mass[TIMEBINS];
extern double TimeBin_BH_dynamicalmass[TIMEBINS];
extern double TimeBin_BH_Mdot[TIMEBINS];
extern double TimeBin_BH_Medd[TIMEBINS];
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(RT_BH_ANGLEWEIGHT_PHOTON_INJECTION)
#define BH_CALC_LOCAL_ANGLEWEIGHTS
#endif
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVACCRETION) || defined(BH_GRAVCAPTURE_NONGAS) || defined(BH_CALC_LOCAL_ANGLEWEIGHTS) || defined(BH_DYNFRICTION)
#define BH_NEIGHBOR_BITFLAG 63 /* allow all particle types in the BH search: 63=2^0+2^1+2^2+2^3+2^4+2^5 */
#else
#define BH_NEIGHBOR_BITFLAG 33 /* only search for particles of types 0 and 5 (gas and black holes) around a primary BH particle */
#endif
#endif

#ifdef HERMITE_INTEGRATION
extern int HermiteOnlyFlag;     /*!< flag to only do Hermite integration for applicable particles (ie. stars) in the gravity routine - set =1 on the first prediction pass and =2 on the second correction pass */
#endif
extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern int PTask;		/*!< note: NTask = 2^PTask */
extern double CPUThisRun;	/*!< Sums CPU time of current process */
extern int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;
extern int NumSphUpdate;	/*!< number of active SPH particles on local processor in current timestep  */
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
  HighMark_pmnonperiodic,  HighMark_sphdensity, HighMark_sphhydro, HighMark_GasGrad;

#ifdef TURB_DRIVING
extern size_t HighMark_turbpower;
#endif
extern int TreeReconstructFlag;
extern int GlobFlag;
extern char DumpFlag;
#ifdef WAKEUP
extern int NeedToWakeupParticles;
extern int NeedToWakeupParticles_local;
#endif

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
#ifdef SEPARATE_STELLARDOMAINDECOMP
extern int N_stars;
#endif
#ifdef BH_WIND_SPAWN
extern double MaxUnSpanMassBH;
#endif

extern long long Ntype[6];	/*!< total number of particles of each type */
extern int NtypeLocal[6];	/*!< local number of particles of each type */
extern gsl_rng *random_generator;	/*!< the random number generator used */
extern int Gas_split;           /*!< current number of newly-spawned gas particles outside block */
#ifdef GALSF
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
#endif

extern double TimeOfLastTreeConstruction;	/*!< holds what it says */
extern int *Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */
extern double *R2ngblist;
extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int *DomainStartList, *DomainEndList;
extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern int *DomainList, DomainNumChanged;
extern peanokey *Key, *KeySorted;

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

extern int NTopnodes, NTopleaves;
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
extern double RndTable[RNDTABLE];
#endif

#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
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


/* variables for input/output , usually only used on process 0 */

extern char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */
extern FILE
#ifndef IO_REDUCED_MODE
 *FdTimebin,    /*!< file handle for timebin.txt log-file. */
 *FdInfo,       /*!< file handle for info.txt log-file. */
 *FdEnergy,     /*!< file handle for energy.txt log-file. */
 *FdTimings,    /*!< file handle for timings.txt log-file. */
 *FdBalance,    /*!< file handle for balance.txt log-file. */
#ifdef RT_CHEM_PHOTOION
 *FdRad,		/*!< file handle for radtransfer.txt log-file. */
#endif
#ifdef TURB_DRIVING
 *FdTurb,       /*!< file handle for turb.txt log-file */
#endif
#ifdef GR_TABULATED_COSMOLOGY
 *FdDE,         /*!< file handle for darkenergy.txt log-file. */
#endif
#endif
 *FdCPU;        /*!< file handle for cpu.txt log-file. */
#ifdef GALSF
extern FILE *FdSfr;		/*!< file handle for sfr.txt log-file. */
#endif
#ifdef GALSF_FB_MECHANICAL
extern FILE *FdSneIIHeating;	/*!< file handle for SNIIheating.txt log-file */
#endif
#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;	/*!< file handle for blackholes.txt log-file. */
#ifdef OUTPUT_SINK_ACCRETION_HIST
extern FILE *FdBhSwallowDetails;
#endif
#ifdef OUTPUT_SINK_FORMATION_PROPS
extern FILE *FdBhFormationDetails;
#endif
#if !defined(IO_REDUCED_MODE) || defined(BH_OUTPUT_MOREINFO)
extern FILE *FdBlackHolesDetails;
#ifdef BH_OUTPUT_MOREINFO
extern FILE *FdBhMergerDetails;
#ifdef BH_WIND_KICK
extern FILE *FdBhWindDetails;
#endif
#endif
#endif
#endif


#if defined(COOLING) && defined(GALSF_EFFECTIVE_EQS)
#ifndef COOLING_OPERATOR_SPLIT
#define COOLING_OPERATOR_SPLIT /*!< the Springel-Hernquist EOS depends explicitly on the cooling time in a way that requires de-coupled hydro cooling */
#endif
#endif

extern double DriftTable[DRIFT_TABLE_LENGTH]; /*! table for the cosmological drift factors */
extern double GravKickTable[DRIFT_TABLE_LENGTH]; /*! table for the cosmological kick factor for gravitational forces */
extern void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

#ifdef BLACK_HOLES
  int TotBHs;
#endif

#if defined(DM_SIDM)
    MyDouble DM_InteractionCrossSection;  /*!< self-interaction cross-section in [cm^2/g]*/
    MyDouble DM_DissipationFactor;  /*!< dimensionless parameter governing efficiency of dissipation (1=dissipative, 0=elastic) */
    MyDouble DM_KickPerCollision;  /*!< for exo-thermic DM reactions, this determines the energy gain 'per event': kick in code units (equivalent to specific energy) associated 'per event' */
    MyDouble DM_InteractionVelocityScale; /*!< scale above which the scattering becomes velocity-dependent */
#endif

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one processor. */
  int ICFormat;			/*!< selects different versions of IC file-format */
  int SnapFormat;		/*!< selects different versions of snapshot file-formats */
  int NumFilesPerSnapshot;	/*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/*!< maximum number of files that may be written simultaneously when writing/reading restart-files, or when writing snapshot files */
  double BufferSize;		/*!< size of communication buffer in MB */
  long BunchSize;     	        /*!< number of particles fitting into the buffer in the parallel tree algorithm  */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

#ifdef DM_SCALARFIELD_SCREENING
  double ScalarBeta;
  double ScalarScreeningLength;
#endif

  /* some SPH parameters */
  double DesNumNgb;		/*!< Desired number of SPH neighbours */
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

#ifdef PM_HIRES_REGION_CLIPDM
    double MassOfClippedDMParticles; /*!< the mass of high-res DM particles which the low-res particles will target if they enter the highres region */
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
    double MeanGasParticleMass; /*!< the mean gas particle mass */
#endif
    double MinMassForParticleMerger; /*!< the minimum mass of a gas particle below which it will be merged into a neighbor */
    double MaxMassForParticleSplit; /*!< the maximum mass of a gas particle above which it will be split into a pair */

  /* some force counters  */
  long long TotNumOfForces;	/*!< counts total number of force computations  */
  long long NumForcesSinceLastDomainDecomp;	/*!< count particle updates since last domain decomposition */

  /* some variable for dynamic work-load adjustment based on CPU measurements */
  double cf_atime, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a, cf_hubble_a2;   /* various cosmological factors that are only a function of the current scale factor, and in Newtonian runs are set to 1 */

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
  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */
  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */
  double MinSizeTimestep,	/*!< minimum allowed timestep. Normally, the simulation terminates if the timestep determined by the timestep criteria falls below this limit. */
    MaxSizeTimestep;		/*!< maximum allowed timestep */
  double MaxRMSDisplacementFac;	/*!< this determines a global timestep criterion for cosmological simulations in comoving coordinates.  To this end, the code computes the rms velocity
				   of all particles, and limits the timestep such that the rms displacement is a fraction of the mean particle separation (determined from the
				   particle mass and the cosmological parameters). This parameter specifies this fraction. */
  int MaxMemSize;
  double CourantFac;		/*!< SPH-Courant factor */

  /* frequency of tree reconstruction/domain decomposition */
  double TreeDomainUpdateFrequency;	/*!< controls frequency of domain decompositions  */
#ifdef TURB_DRIVING
    double TurbInjectedEnergy;
    double TurbDissipatedEnergy;
#if defined(TURB_DRIVING_SPECTRUMGRID)
    double TimeBetTurbSpectrum;
    double TimeNextTurbSpectrum;
    int FileNumberTurbSpectrum;
#endif
#endif

  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length)
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars */
    double MinGasHsmlFractional; /*!< minimim allowed gas kernel length relative to force softening (what you actually set) */
    double MinHsml;			/*!< minimum allowed gas kernel length */
    double MaxHsml;           /*!< minimum allowed gas kernel length */

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

  double SofteningTable[6];	/*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];	/*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */

  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this value, * allowing the size of the snapshot files to be reduced */
  double MassTable[6];

  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    RestartFile[100], ResubmitCommand[100], OutputListFilename[100];
    /* EnergyFile[100], CpuFile[100], InfoFile[100], TimingsFile[100], TimebinFile[100], */
#ifdef COOL_GRACKLE
    char GrackleDataFile[100];
#endif
    /*! table with desired output times */
    double OutputListTimes[MAXLEN_OUTPUTLIST];
    char OutputListFlag[MAXLEN_OUTPUTLIST];
    int OutputListLength;		/*!< number of times stored in table of desired output times */

#ifdef RADTRANSFER
    integertime Radiation_Ti_begstep;
    integertime Radiation_Ti_endstep;
#endif
#ifdef RT_EVOLVE_INTENSITIES
    double Rad_Intensity_Direction[N_RT_INTENSITY_BINS][3];
#endif


#if defined(RT_CHEM_PHOTOION) && !(defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF))
    double IonizingLuminosityPerSolarMass_cgs;
    double star_Teff;
#endif

#ifdef RT_LEBRON
    double PhotonMomentum_Coupled_Fraction;
#endif
#ifdef BH_PHOTONMOMENTUM
    double BH_Rad_MomentumFactor;
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

#ifdef GDE_DISTORTIONTENSOR
  /* present day velocity dispersion of DM particle in cm/s (e.g. Neutralino = 0.03 cm/s) */
  double DM_velocity_dispersion;
  double TidalCorrection;
#ifdef GDE_LEAN
  double GDEInitStreamDensity;
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

#if (defined(GALSF) && defined(METALS)) || defined(COOL_METAL_LINES_BY_SPECIES) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF_FB_MECHANICAL) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF_FB_THERMAL)
#define INIT_STELLAR_METALS_AGES_DEFINED // convenience flag for later to know these variables exist
    double InitMetallicityinSolar;
    double InitStellarAgeinGyr;
#endif

#if defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK) || defined(BH_WIND_SPAWN)
    double BAL_f_accretion;
    double BAL_v_outflow;
#endif

#if defined(SINGLE_STAR_FB_JETS)
    double BAL_f_launch_v; // scales the amount of accretion power going into jets, we eject (1-All.BAL_f_accretion) fraction of the accreted mass at this value times the Keplerian velocity at the protostellar radius. If set to 1 then the mass and power loading of the jets are both (1-All.BAL_f_accretion)
#endif

#if defined(BH_COSMIC_RAYS)
    double BH_CosmicRay_Injection_Efficiency;
#endif
    

#ifdef METALS
    double SolarAbundances[NUM_METAL_SPECIES];
#ifdef COOL_METAL_LINES_BY_SPECIES
    int SpeciesTableInUse;
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

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
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

#if (defined(BLACK_HOLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
  double TimeNextOnTheFlyFoF;
  double TimeBetOnTheFlyFoF;
#endif

#ifdef BLACK_HOLES
  double BlackHoleAccretionFactor;	/*!< Fraction of BH bondi accretion rate */
  double BlackHoleFeedbackFactor;	/*!< Fraction of the black luminosity feed into thermal feedback */
  double SeedBlackHoleMass;         /*!< Seed black hole mass */
#if defined(BH_SEED_FROM_FOF) || defined(BH_SEED_FROM_LOCALGAS)
  double SeedBlackHoleMassSigma;    /*!< Standard deviation of init black hole masses */
  double SeedBlackHoleMinRedshift;  /*!< Minimum redshift where BH seeds are allowed */
#ifdef BH_SEED_FROM_LOCALGAS
  double SeedBlackHolePerUnitMass;  /*!< Defines probability per unit mass of seed BH forming */
#endif
#endif
#ifdef BH_ALPHADISK_ACCRETION
  double SeedAlphaDiskMass;         /*!< Seed alpha disk mass */
#endif
#ifdef BH_WIND_SPAWN
  double BAL_wind_particle_mass;        /*!< target mass for feedback particles to be spawned */
  double BAL_internal_temperature;
  MyIDType AGNWindID;
#endif
#ifdef BH_SEED_FROM_FOF
  double MinFoFMassForNewSeed;      /*!< Halo mass required before new seed is put in */
#endif
  double BlackHoleNgbFactor;        /*!< Factor by which the SPH neighbour should be increased/decreased */
  double BlackHoleMaxAccretionRadius;
  double BlackHoleEddingtonFactor;	/*!< Factor above Eddington */
  double BlackHoleRadiativeEfficiency;  /**< Radiative efficiency determined by the spin value, default value is 0.1 */
#endif

#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)
  double Tillotson_EOS_params[7][12]; /*! < holds parameters for Tillotson EOS for solids */
#endif

#ifdef EOS_TABULATED
    char EosTable[100];
#endif

#ifdef NUCLEAR_NETWORK
  char EosSpecies[100];
  char NetworkRates[100];
  char NetworkPartFunc[100];
  char NetworkMasses[100];
  char NetworkWeakrates[100];
  struct network_data nd;
  struct network_workspace nw;
  double NetworkTempThreshold;
#endif

#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
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

#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR)
  double BH_spawn_rinj;
  double B_spawn_pol;
  double B_spawn_tor;
#endif
#ifdef BH_WIND_SPAWN_SET_JET_PRECESSION
  double BH_jet_precess_degree;
  double BH_jet_precess_period;
#endif
#ifdef BH_DEBUG_FIX_MDOT_MBH
  double BH_fb_duty_cycle;
  double BH_fb_period;
#endif

}
All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern ALIGN(32) struct particle_data
{
    short int Type;                 /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
    short int TimeBin;
    MyIDType ID;                    /*! < unique ID of particle (assigned at beginning of the simulation) */
    MyIDType ID_child_number;       /*! < child number for particles 'split' from main (retain ID, get new child number) */
#ifndef BH_WIND_SPAWN
    int ID_generation;              /*! < generation (need to track for particle-splitting to ensure each 'child' gets a unique child number */
#else
    MyIDType ID_generation;
#endif

    integertime Ti_begstep;         /*!< marks start of current timestep of particle on integer timeline */
    integertime Ti_current;         /*!< current time of the particle */

    ALIGN(32) MyDouble Pos[3];      /*!< particle position at its current time */
    MyDouble Mass;                  /*!< particle mass */

    MyDouble Vel[3];                /*!< particle velocity at its current time */
    MyDouble dp[3];
    MyFloat Particle_DivVel;        /*!< velocity divergence of neighbors (for predict step) */

    MyDouble GravAccel[3];          /*!< particle acceleration due to gravity */
#ifdef PMGRID
    MyFloat GravPM[3];		/*!< particle acceleration due to long-range PM gravity force */
#endif
    MyFloat OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening criterion */
#ifdef HERMITE_INTEGRATION
    MyFloat Hermite_OldAcc[3];
    MyFloat OldPos[3];
    MyFloat OldVel[3];
    MyFloat OldJerk[3];
    short int AccretedThisTimestep;     /*!< flag to decide whether to stick with the KDK step for stability reasons, e.g. when actively accreting */
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    MyFloat TreeMass;  /*!< Mass seen by the particle as it sums up the gravitational force from the tree - should be equal to total mass, a useful debug diagnostic  */
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL)
    MyFloat Potential;		/*!< gravitational potential */
#if defined(EVALPOTENTIAL) && defined(PMGRID)
    MyFloat PM_Potential;
#endif
#endif
#if defined(GALSF_SFR_TIDAL_HILL_CRITERION) || defined(TIDAL_TIMESTEP_CRITERION) || defined(GDE_DISTORTIONTENSOR) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(OUTPUT_TIDAL_TENSOR) || (defined(SINGLE_STAR_TIMESTEPPING) && (SINGLE_STAR_TIMESTEPPING > 0))
#define COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    double tidal_tensorps[3][3];                        /*!< tidal tensor (=second derivatives of grav. potential) */
#ifdef PMGRID
    double tidal_tensorpsPM[3][3];                /*!< for TreePM simulations, long range tidal field */
#endif
#endif
    
#ifdef ADAPTIVE_TREEFORCE_UPDATE
    MyFloat time_since_last_treeforce;
    MyFloat tdyn_step_for_treeforce;
#ifndef COMPUTE_JERK_IN_GRAVTREE
#define COMPUTE_JERK_IN_GRAVTREE    
#endif    
#endif
    
#ifdef COMPUTE_JERK_IN_GRAVTREE
    double GravJerk[3];
#endif
#ifdef GDE_DISTORTIONTENSOR
    MyBigFloat distortion_tensorps[6][6];               /*!< phase space distortion tensor */
    MyBigFloat last_determinant;                        /*!< last real space distortion tensor determinant */
    MyBigFloat stream_density;                          /*!< physical stream density that is going to be integrated */
    float caustic_counter;                              /*!< caustic counter */
#ifndef GDE_LEAN
    MyBigFloat annihilation;                            /*!< integrated annihilation rate */
    MyBigFloat analytic_annihilation;                   /*!< analytically integrated annihilation rate */
    MyBigFloat rho_normed_cutoff_current;               /*!< current and last normed_cutoff density in rho_max/rho_init * sqrt(sigma) */
    MyBigFloat rho_normed_cutoff_last;
    MyBigFloat s_1_current, s_2_current, s_3_current;   /*! < current and last stretching factor */
    MyBigFloat s_1_last, s_2_last, s_3_last;
    MyBigFloat second_deriv_current;                    /*! < current and last second derivative */
    MyBigFloat second_deriv_last;
    double V_matrix[3][3];                              /*!< initial orientation of CDM sheet the particle is embedded in */
    float init_density;                                 /*!< initial stream density */
    float analytic_caustics;                            /*!< number of caustics that were integrated analytically */
    float a0;
#endif
#ifdef OUTPUT_GDE_LASTCAUSTIC
    MyFloat lc_Time;                                  /*!< time of caustic passage */
    MyFloat lc_Pos[3];                                /*!< position of caustic */
    MyFloat lc_Vel[3];                                /*!< particle velocity when passing through caustic */
    MyFloat lc_rho_normed_cutoff;                     /*!< normed_cutoff density at caustic */
    MyFloat lc_Dir_x[3];                              /*!< principal axis frame of smear out */
    MyFloat lc_Dir_y[3];
    MyFloat lc_Dir_z[3];
    MyFloat lc_smear_x;                               /*!< smear out length */
    MyFloat lc_smear_y;
    MyFloat lc_smear_z;
#endif
#endif // GDE_DISTORTIONTENSOR //

#ifdef GALSF
    MyFloat StellarAge;		/*!< formation time of star particle */
#endif
#ifdef METALS
    MyFloat Metallicity[NUM_METAL_SPECIES]; /*!< metallicity (species-by-species) of gas or star particle */
#endif
#ifdef GALSF_SFR_IMF_VARIATION
    MyFloat IMF_Mturnover; /*!< IMF turnover mass [in solar] (or any other parameter which conveniently describes the IMF) */
    MyFloat IMF_FormProps[N_IMF_FORMPROPS]; /*!< formation properties of star particles to record for output */
#endif
#ifdef GALSF_SFR_IMF_SAMPLING
    MyFloat IMF_NumMassiveStars; /*!< number of massive stars to associate with this star particle (for feedback) */
#endif

    MyFloat Hsml;                   /*!< search radius around particle for neighbors/interactions */
    MyFloat NumNgb;                 /*!< neighbor number around particle */
    MyFloat DhsmlNgbFactor;        /*!< correction factor needed for varying kernel lengths */
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    MyFloat DensAroundStar;         /*!< gas density in the neighborhood of the collisionless particle (evaluated from neighbors) */
    MyFloat GradRho[3];             /*!< gas density gradient evaluated simply from the neighboring particles, for collisionless centers */
#endif
#ifdef RT_USE_TREECOL_FOR_NH
    MyFloat ColumnDensityBins[RT_USE_TREECOL_FOR_NH];     /*!< angular bins for column density */
    MyFloat SigmaEff;              /*!< effective column density -log(avg(exp(-sigma))) averaged over column density bins from the gravity tree (does not include the self-contribution) */
#endif
#if defined(RT_SOURCE_INJECTION)
    MyFloat KernelSum_Around_RT_Source; /*!< kernel summation around sources for radiation injection (save so can be different from 'density') */
#endif

#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
    MyFloat SNe_ThisTimeStep; /* flag that indicated number of SNe for the particle in the timestep */
#endif
#ifdef GALSF_FB_MECHANICAL
#define AREA_WEIGHTED_SUM_ELEMENTS 11 /* number of weights needed for full momentum-and-energy conserving system */
    MyFloat Area_weighted_sum[AREA_WEIGHTED_SUM_ELEMENTS]; /* normalized weights for particles in kernel weighted by area, not mass */
#endif

#if defined(GRAIN_FLUID)
    MyFloat Grain_Size;
    MyFloat Gas_Density;
    MyFloat Gas_InternalEnergy;
    MyFloat Gas_Velocity[3];
    MyFloat Grain_AccelTimeMin;
#if defined(GRAIN_BACKREACTION)
    MyFloat Grain_DeltaMomentum[3];
#endif
#if defined(GRAIN_LORENTZFORCE)
    MyFloat Gas_B[3];
#endif
#endif
#if defined(PIC_MHD)
    short int MHD_PIC_SubType;
#endif

#if defined(BLACK_HOLES)
    MyIDType SwallowID;
    int IndexMapToTempStruc;   /*!< allows for mapping to BlackholeTempInfo struc */
#ifdef BH_WIND_SPAWN
    MyFloat unspawned_wind_mass;    /*!< tabulates the wind mass which has not yet been spawned */
#endif
#ifdef BH_COUNTPROGS
    int BH_CountProgs;
#endif
    MyFloat BH_Mass;
#if defined(BH_GRAVCAPTURE_FIXEDSINKRADIUS)
    MyFloat SinkRadius;
#endif
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
    MyFloat dt_since_last_gas_search; /* keep track of time since the sink's last neighbor search and gas interaction (for feedback/accretion) */
    short int do_gas_search_this_timestep; /* flag for deciding whether to do gas stuff for a given timestep */
#endif
#ifdef GRAIN_FLUID
    MyFloat BH_Dust_Mass;
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
    MyFloat SwallowTime; /* freefall time of a particle onto a sink particle  */
#endif
#if defined(SINGLE_STAR_TIMESTEPPING)
    MyFloat BH_SurroundingGasVel; /* Relative speed of sink to surrounding gas  */
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
    int BH_Ngb_Flag; /* Whether or not the gas live's in a sink's hydro stencil */
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat BH_Mass_AlphaDisk;
#endif
#ifdef BH_WAKEUP_GAS /* force all gas within the interaction radius of a sink to timestep at the same rate */
    int LowestBHTimeBin;
#endif
#ifdef BH_FOLLOW_ACCRETED_ANGMOM
    MyFloat BH_Specific_AngMom[3];
#endif
#ifdef BH_RETURN_BFLUX
    MyDouble B[3];
#endif
#ifdef JET_DIRECTION_FROM_KERNEL_AND_SINK
    MyFloat Mgas_in_Kernel;
    MyFloat Jgas_in_Kernel[3];
#endif
    MyFloat BH_Mdot;
    int BH_TimeBinGasNeighbor;
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
    MyFloat BH_dr_to_NearestGasNeighbor;
#endif
#ifdef BH_REPOSITION_ON_POTMIN
    MyFloat BH_MinPotPos[3];
    MyFloat BH_MinPot;
#endif
#endif  /* if defined(BLACK_HOLES) */
#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    MyFloat MencInRcrit;
#endif

#ifdef BH_CALC_DISTANCES
    MyFloat min_dist_to_bh;
    MyFloat min_xyz_to_bh[3];
#if defined(SINGLE_STAR_FIND_BINARIES) || (SINGLE_STAR_TIMESTEPPING > 0)
    MyDouble min_bh_t_orbital; //orbital time for binary
    MyDouble comp_dx[3]; //position offset of binary companion - this will be evolved in the Kepler solution while we use the Pos attribute to track the binary COM
    MyDouble comp_dv[3]; //velocity offset of binary companion - this will be evolved in the Kepler solution while we use the Vel attribute to track the binary COM velocity
    MyDouble comp_Mass; //mass of binary companion
    int is_in_a_binary; // flag whether star is in a binary or not
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    MyFloat min_bh_freefall_time;
    MyFloat min_bh_approach_time;
#if (SINGLE_STAR_TIMESTEPPING > 0)
    int SuperTimestepFlag; // >=2 if allowed to super-timestep (increases with each drift/kick), 1 if a candidate for super-timestepping, 0 otherwise
    MyDouble COM_dt_tidal; //timescale from tidal tensor evaluated at the center of mass without contribution from the companion
    MyDouble COM_GravAccel[3]; //gravitational acceleration evaluated at the center of mass without contribution from the companion
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    MyFloat MaxFeedbackVel; // maximum signal velocity of any feedback mechanism emanating from the star
    MyFloat min_bh_fb_time;  // minimum time for feedback to arrive from a star
#endif    
#endif
#endif
   


#if defined(DM_SIDM)
    double dtime_sidm; /*!< timestep used if self-interaction probabilities greater than 0.2 are found */
    long unsigned int NInteractions; /*!< Total number of interactions */
#endif

#if defined(SUBFIND)
    int GrNr;
    int SubNr;
    int DM_NumNgb;
    unsigned short targettask, origintask2;
    int origintask, submark, origindex;
    MyFloat DM_Hsml;
    union
    {
        MyFloat DM_Density;
        MyFloat DM_Potential;
    } u;
    union
    {
        MyFloat DM_VelDisp;
        MyFloat DM_BindingEnergy;
    } v;
#ifdef FOF_DENSITY_SPLIT_TYPES
    union
    {
        MyFloat int_energy;
        MyFloat density_sum;
    } w;
#endif
#endif

    float GravCost[GRAVCOSTLEVELS];   /*!< weight factor used for balancing the work-load */

#ifdef WAKEUP
    integertime dt_step;
#endif

#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    MyDouble AGS_Hsml;          /*!< smoothing length (for gravitational forces) */
    MyFloat AGS_zeta;           /*!< factor in the correction term */
    MyDouble AGS_vsig;          /*!< signal velocity of particle approach, to properly time-step */
#if defined(WAKEUP)
    short int wakeup;                     /*!< flag to wake up particle */
#endif
#endif

#ifdef DM_FUZZY
    MyFloat AGS_Density;                /*!< density calculated corresponding to AGS routine (over interacting DM neighbors) */
    MyFloat AGS_Gradients_Density[3];   /*!< density gradient calculated corresponding to AGS routine (over interacting DM neighbors) */
    MyFloat AGS_Gradients2_Density[3][3];   /*!< density gradient calculated corresponding to AGS routine (over interacting DM neighbors) */
    MyFloat AGS_Numerical_QuantumPotential; /*!< additional potential terms 'generated' by un-resolved compression [numerical diffusivity] */
    MyFloat AGS_Dt_Numerical_QuantumPotential; /*!< time derivative of the above */
#if (DM_FUZZY > 0)
    MyFloat AGS_Psi_Re;
    MyFloat AGS_Psi_Re_Pred;
    MyFloat AGS_Dt_Psi_Re;
    MyFloat AGS_Gradients_Psi_Re[3];
    MyFloat AGS_Gradients2_Psi_Re[3][3];
    MyFloat AGS_Psi_Im;
    MyFloat AGS_Psi_Im_Pred;
    MyFloat AGS_Dt_Psi_Im;
    MyFloat AGS_Gradients_Psi_Im[3];
    MyFloat AGS_Gradients2_Psi_Im[3][3];
    MyFloat AGS_Dt_Psi_Mass;
#endif
#endif
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    MyLongDouble NV_T[3][3];                                           /*!< holds the tensor used for gradient estimation */
#endif
}
 *P,				/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */


#ifndef GDE_LEAN
#define GDE_TIMEBEGIN(i) (P[i].a0)
#define GDE_VMATRIX(i, a, b) (P[i].V_matrix[a][b])
#define GDE_INITDENSITY(i) (P[i].init_density)
#else
#define GDE_TIMEBEGIN(i) (All.TimeBegin)
#define GDE_VMATRIX(i, a, b) (0.0)
#define GDE_INITDENSITY(i) (All.GDEInitStreamDensity)
#endif

#if defined(BLACK_HOLES)
#define BPP(i) P[(i)]
#endif


#if defined(HYDRO_TENSOR_FACE_CORRECTIONS_NGBITER)
#define HYDRO_TENSOR_FACE_CORRECTIONS_NUMBER_MOMWTS 15
#else
#define HYDRO_TENSOR_FACE_CORRECTIONS_NUMBER_MOMWTS 9
#endif


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
    /* the PRIMITIVE and CONSERVED hydro variables used in STATE reconstruction */
    MyDouble Density;               /*!< current baryonic mass density of particle */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    MyDouble MassTrue;              /*!< true particle mass ('mass' now is -predicted- mass */
    MyDouble dMass;                 /*!< change in particle masses from hydro step (conserved variable) */
    MyDouble DtMass;                /*!< rate-of-change of particle masses (for drifting) */
    MyDouble GravWorkTerm[3];       /*!< correction term needed for hydro mass flux in gravity */
    MyDouble ParticleVel[3];        /*!< actual velocity of the mesh-generating points */
#endif

    MyDouble Pressure;              /*!< current pressure */
    MyDouble InternalEnergy;        /*!< specific internal energy [internal thermal energy per unit mass] of cell */
    MyDouble InternalEnergyPred;    /*!< predicted value of the specific internal energy at the current time */
    //MyDouble dInternalEnergy;     /*!< change in specific internal energy from hydro step */ //manifest-indiv-timestep-debug//
    MyDouble DtInternalEnergy;      /*!< rate of change of specific internal energy */

    MyDouble VelPred[3];            /*!< predicted SPH particle velocity at the current time */
    //MyDouble dMomentum[3];        /*!< change in momentum from hydro step (conserved variable) */ //manifest-indiv-timestep-debug//
    MyDouble HydroAccel[3];         /*!< acceleration due to hydrodynamical force (for drifting) */

#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
    MyDouble Density_ExplicitInt;   /*!< explicitly integrated volume/density variable to be used if integrating the SPH-like form of the continuity directly */
#endif
    
#ifdef HYDRO_VOLUME_CORRECTIONS
    MyDouble Volume_0;              /*!< 0th-order cell volume for mesh-free (MFM/MFV-type) reconstruction at 0th-order volume quadrature */
    MyDouble Volume_1;              /*!< 1st-order cell volume for mesh-free (MFM/MFV-type) reconstruction at 1st-order volume quadrature */
#endif

#ifdef MAGNETIC
    MyDouble Face_Area[3];          /*!< vector sum of effective areas of 'faces'; this is used to check closure for meshless methods */
    MyDouble BPred[3];              /*!< current magnetic field strength */
    MyDouble B[3];                  /*!< actual B (conserved variable used for integration; can be B*V for flux schemes) */
    MyDouble DtB[3];                /*!< time derivative of B-field (of -conserved- B-field) */
    MyFloat divB;                   /*!< storage for the 'effective' divB used in div-cleaning procedure */
#ifdef DIVBCLEANING_DEDNER
    MyDouble DtB_PhiCorr[3];        /*!< correction forces for mid-face update to phi-field */
    MyDouble PhiPred;               /*!< current value of Phi */
    MyDouble Phi;                   /*!< scalar field for Dedner divergence cleaning */
    MyDouble DtPhi;                 /*!< time derivative of Phi-field */
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
    int FlagForConstrainedGradients;/*!< flag indicating whether the B-field gradient is a 'standard' one or the constrained-divB version */
#endif
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
    MyFloat Balpha;                 /*!< effective resistivity coefficient */
#endif
#endif /* MAGNETIC */

#if defined(KERNEL_CRK_FACES)
    MyFloat Tensor_CRK_Face_Corrections[16]; /*!< tensor set for face-area correction terms for the CRK formulation of SPH or MFM/V areas */
#endif
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS)
    MyFloat Tensor_MFM_Face_Corrections[9]; /*!< alternative tensor face corrections for linear consistency */
#endif

#ifdef COSMIC_RAY_FLUID
    MyFloat CosmicRayEnergy[N_CR_PARTICLE_BINS];        /*!< total energy of cosmic ray fluid (the conserved variable) */
    MyFloat CosmicRayEnergyPred[N_CR_PARTICLE_BINS];    /*!< total energy of cosmic ray fluid (the conserved variable) */
    MyFloat DtCosmicRayEnergy[N_CR_PARTICLE_BINS];      /*!< time derivative of cosmic ray energy */
    MyFloat CosmicRayDiffusionCoeff[N_CR_PARTICLE_BINS];/*!< diffusion coefficient kappa for cosmic ray fluid */
    MyFloat Face_DivVel_ForAdOps;                                 /*!< face-centered definition of the velocity divergence, needed to carefully handle adiabatic terms when Pcr >> Pgas */
#ifdef CRFLUID_M1
    MyFloat CosmicRayFlux[N_CR_PARTICLE_BINS][3];       /*!< CR flux vector [explicitly evolved] - conserved-variable */
    MyFloat CosmicRayFluxPred[N_CR_PARTICLE_BINS][3];   /*!< CR flux vector [explicitly evolved] - conserved-variable */
#endif
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
    MyFloat CosmicRayAlfvenEnergy[N_CR_PARTICLE_BINS][2];       /*!< forward and backward-traveling Alfven wave-packet energies */
    MyFloat CosmicRayAlfvenEnergyPred[N_CR_PARTICLE_BINS][2];   /*!< drifted forward and backward-traveling Alfven wave-packet energies */
    MyFloat DtCosmicRayAlfvenEnergy[N_CR_PARTICLE_BINS][2];     /*!< time derivative fof forward and backward-traveling Alfven wave-packet energies */
#endif
#endif

#ifdef SUPER_TIMESTEP_DIFFUSION
    MyDouble Super_Timestep_Dt_Explicit; /*!< records the explicit step being used to scale the sub-steps for the super-stepping */
    int Super_Timestep_j; /*!< records which sub-step if the super-stepping cycle the particle is in [needed for adaptive steps] */
#endif

#if (SINGLE_STAR_SINK_FORMATION & 4)
    MyFloat Density_Relative_Maximum_in_Kernel; /*!< hold density_max-density_i, for particle i, so we know if its a local maximum */
#endif

    /* matrix of the primitive variable gradients: rho, P, vx, vy, vz, B, phi */
    struct
    {
        MyDouble Density[3];
        MyDouble Pressure[3];
        MyDouble Velocity[3][3];
#ifdef MAGNETIC
        MyDouble B[3][3];
#ifdef DIVBCLEANING_DEDNER
        MyDouble Phi[3];
#endif
#endif
#ifdef DOGRAD_SOUNDSPEED
        MyDouble SoundSpeed[3];
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
        MyDouble InternalEnergy[3];
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        MyDouble Metallicity[NUM_METAL_SPECIES][3];
#endif
#ifdef COSMIC_RAY_FLUID
        MyDouble CosmicRayPressure[N_CR_PARTICLE_BINS][3];
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
        MyDouble Rad_E_gamma_ET[N_RT_FREQ_BINS][3];
#endif
    } Gradients;
    MyLongDouble NV_T[3][3];        /*!< holds the tensor used for gradient estimation */
    MyLongDouble ConditionNumber;   /*!< condition number of the gradient matrix: needed to ensure stability */
    MyDouble FaceClosureError;      /*!< dimensionless measure of face closure */
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
    MyDouble MaxKineticEnergyNgb;   /*!< maximum kinetic energy (with respect to neighbors): use for entropy 'switch' */
#endif

#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    MyFloat Dyield[NUM_METAL_SPECIES];
#endif

#ifdef HYDRO_SPH
    MyDouble DhsmlHydroSumFactor;   /* for 'traditional' SPH, we need the SPH hydro-element volume estimator */
#endif

#ifdef HYDRO_PRESSURE_SPH
    MyDouble EgyWtDensity;          /*!< 'effective' rho to use in hydro equations */
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(AGS_HSML_CALCULATION_IS_ACTIVE)
    MyFloat AGS_zeta;               /*!< correction term for adaptive gravitational softening lengths */
#endif

    MyFloat MaxSignalVel;           /*!< maximum signal velocity (needed for time-stepping) */


#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) 
    MyDouble IniDen;
    MyDouble IniB[3];
#endif     

#ifdef CHIMES_STELLAR_FLUXES
    double Chimes_G0[CHIMES_LOCAL_UV_NBINS];            /*!< 6-13.6 eV flux, in Habing units */
    double Chimes_fluxPhotIon[CHIMES_LOCAL_UV_NBINS];   /*!< ionising flux (>13.6 eV), in cm^-2 s^-1 */
#ifdef CHIMES_HII_REGIONS
    double Chimes_G0_HII[CHIMES_LOCAL_UV_NBINS];
    double Chimes_fluxPhotIon_HII[CHIMES_LOCAL_UV_NBINS];
#endif
#endif
#ifdef CHIMES_TURB_DIFF_IONS
    double ChimesNIons[CHIMES_TOTSIZE];
#endif
#ifdef BH_COMPTON_HEATING
    MyFloat Rad_Flux_AGN;             /*!< local AGN flux */
#endif


#if defined(TURB_DRIVING) || defined(OUTPUT_VORTICITY)
   MyFloat Vorticity[3];
   MyFloat SmoothedVel[3];
#endif

#if defined(BH_THERMALFEEDBACK)
    MyDouble Injected_BH_Energy;
#endif

#ifdef COOLING
#if !defined(COOLING_OPERATOR_SPLIT)
    int CoolingIsOperatorSplitThisTimestep; /* flag to tell us if cooling is operator split or not on a given timestep */
#endif
#ifndef CHIMES
    MyFloat Ne;  /*!< electron fraction, expressed as local electron number
		    density normalized to the hydrogen number density. Gives
		    indirectly ionization state and mean molecular weight. */
#endif
#endif
#ifdef GALSF
  MyFloat Sfr;                      /*!< particle star formation rate */
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED)
  MyFloat AlphaVirial_SF_TimeSmoothed;  /*!< dimensionless number > 0.5 if self-gravitating for smoothed virial criterion */
#endif
#endif
#ifdef GALSF_SUBGRID_WINDS
  MyFloat DelayTime;                /*!< remaining maximum decoupling time of wind particle */
#if (GALSF_SUBGRID_WIND_SCALING==1)
  MyFloat HostHaloMass;             /*!< host halo mass estimator for wind launching velocity */
#endif
#if (GALSF_SUBGRID_WIND_SCALING==2)
  MyFloat HsmlDM;                   /*!< smoothing length to find neighboring dark matter particles */
  MyDouble NumNgbDM;                /*!< number of neighbor dark matter particles */
  MyDouble DM_Vx, DM_Vy, DM_Vz, DM_VelDisp; /*!< surrounding DM velocity and velocity dispersion */
#endif
#endif

#ifdef GALSF_FB_TURNOFF_COOLING
  MyFloat DelayTimeCoolingSNe;      /*!< flag indicating cooling is suppressed b/c heated by SNe */
#endif

#ifdef TURB_DRIVING
  MyDouble DuDt_diss;               /*!< quantities specific to turbulent driving routines */
  MyDouble DuDt_drive;
  MyDouble EgyDiss;
  MyDouble EgyDrive;
  MyDouble TurbAccel[3];
#endif

#ifdef TURB_DIFFUSION
  MyFloat TD_DiffCoeff;             /*!< effective diffusion coefficient for sub-grid turbulent diffusion */
#ifdef TURB_DIFF_DYNAMIC
  MyDouble h_turb;
  MyDouble MagShear;
  MyFloat TD_DynDiffCoeff;          /*!< improved Smag. coefficient (squared) for sub-grid turb. diff. - D. Rennehan */
#endif
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
  MyFloat NV_DivVel;                /*!< quantities specific to the Cullen & Dehnen viscosity switch */
  MyFloat NV_dt_DivVel;
  MyFloat NV_A[3][3];
  MyFloat NV_D[3][3];
  MyFloat NV_trSSt;
  MyFloat alpha;
#endif

#ifdef HYDRO_SPH
  MyFloat alpha_limiter;                /*!< artificial viscosity limiter (Balsara-like) */
#endif

#ifdef CONDUCTION
    MyFloat Kappa_Conduction;                   /*!< conduction coefficient */
#endif

#if defined(OUTPUT_MOLECULAR_FRACTION) || defined(COOL_MOLECFRAC_NONEQM)
    MyFloat MolecularMassFraction;              /*!< holder for molecular mass fraction for sims where we evaluate it on-the-fly and wish to save it [different from detailed chemistry modules] */
#if defined(COOL_MOLECFRAC_NONEQM)
    MyFloat MolecularMassFraction_perNeutralH;  /*! molecular mass fraction -of-the-neutral-gas-, which we retain as a separate variable since we have a hybrid model here using implicit updates for the ionization fraction */
#endif
#endif

#ifdef MHD_NON_IDEAL
    MyFloat Eta_MHD_OhmicResistivity_Coeff;     /*!< Ohmic resistivity coefficient [physical units of L^2/t] */
    MyFloat Eta_MHD_HallEffect_Coeff;           /*!< Hall effect coefficient [physical units of L^2/t] */
    MyFloat Eta_MHD_AmbiPolarDiffusion_Coeff;   /*!< Hall effect coefficient [physical units of L^2/t] */
#endif


#if defined(VISCOSITY)
    MyFloat Eta_ShearViscosity;         /*!< shear viscosity coefficient */
    MyFloat Zeta_BulkViscosity;         /*!< bulk viscosity coefficient */
#endif


#if defined(RADTRANSFER)
    MyFloat ET[N_RT_FREQ_BINS][6];          /*!< eddington tensor - symmetric -> only 6 elements needed: this is dimensionless by our definition */
    MyFloat Rad_Je[N_RT_FREQ_BINS];         /*!< emissivity (includes sources like stars, as well as gas): units=Rad_E_gamma/time  */
    MyFloat Rad_E_gamma[N_RT_FREQ_BINS];    /*!< photon energy (integral of dRad_E_gamma/dvol*dVol) associated with particle [for simple frequency bins, equivalent to photon number] */
    MyFloat Rad_Kappa[N_RT_FREQ_BINS];      /*!< opacity [physical units ~ length^2 / mass]  */
#ifdef RT_FLUXLIMITER
    MyFloat Rad_Flux_Limiter[N_RT_FREQ_BINS]; /*!< dimensionless flux-limiter (0<lambda<1) */
#endif
#ifdef RT_EVOLVE_INTENSITIES
    MyFloat Rad_Intensity[N_RT_FREQ_BINS][N_RT_INTENSITY_BINS]; /*!< intensity values along different directions, for each frequency */
    MyFloat Rad_Intensity_Pred[N_RT_FREQ_BINS][N_RT_INTENSITY_BINS]; /*!< predicted [drifted] values of intensities */
    MyFloat Dt_Rad_Intensity[N_RT_FREQ_BINS][N_RT_INTENSITY_BINS]; /*!< time derivative of intensities */
#endif
#ifdef RT_EVOLVE_FLUX
    MyFloat Rad_Flux[N_RT_FREQ_BINS][3];    /*!< photon energy flux density (energy/time/area), for methods which track this explicitly (e.g. M1) */
    MyFloat Rad_Flux_Pred[N_RT_FREQ_BINS][3];/*!< predicted photon energy flux density for drift operations (needed for adaptive timestepping) */
    MyFloat Dt_Rad_Flux[N_RT_FREQ_BINS][3]; /*!< time derivative of photon energy flux density */
#else
#define Rad_Flux_Pred Rad_Flux
#endif
#ifdef RT_EVOLVE_ENERGY
    MyFloat Rad_E_gamma_Pred[N_RT_FREQ_BINS]; /*!< predicted Rad_E_gamma for drift operations (needed for adaptive timestepping) */
    MyFloat Dt_Rad_E_gamma[N_RT_FREQ_BINS]; /*!< time derivative of photon number in particle (used only with explicit solvers) */
#else
#define Rad_E_gamma_Pred Rad_E_gamma        /*! define a useful shortcut for use throughout code so we don't have to worry about Pred-vs-true difference */
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
    MyFloat Rad_Accel[3];
#endif
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
    MyDouble Interpolated_Opacity[N_RT_FREQ_BINS]; /* opacity values interpolated to gas positions */
#endif
#ifdef RT_INFRARED
    MyFloat Radiation_Temperature; /* IR radiation field temperature (evolved variable ^4 power, for convenience) */
    MyFloat Dt_Rad_E_gamma_T_weighted_IR; /* IR radiation temperature-weighted time derivative of photon energy (evolved variable ^4 power, for convenience) */
    MyFloat Dust_Temperature; /* Dust temperature (evolved variable ^4 power, for convenience) */
    MyFloat LambdaDust; /* Dust cooling rate */
#endif
#ifdef RT_CHEM_PHOTOION
    MyFloat HI;                  /* HI fraction */
    MyFloat HII;                 /* HII fraction */
#ifndef COOLING
    MyFloat Ne;               /* electron fraction */
#endif
#ifdef RT_CHEM_PHOTOION_HE
    MyFloat HeI;                 /* HeI fraction */
    MyFloat HeII;                 /* HeII fraction */
    MyFloat HeIII;                 /* HeIII fraction */
#endif
#endif // end of chem-photoion
#endif // end of radtransfer
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) && !defined(RADTRANSFER)
    MyFloat Rad_E_gamma[N_RT_FREQ_BINS];
#define Rad_E_gamma_Pred Rad_E_gamma
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX) && !defined(RT_EVOLVE_FLUX)
    MyFloat Rad_Flux[N_RT_FREQ_BINS][3];
#define Rad_Flux_Pred Rad_Flux
#endif
    


#ifdef EOS_GENERAL
    MyFloat SoundSpeed;                   /* Sound speed */
#ifdef EOS_CARRIES_TEMPERATURE
    MyFloat Temperature;                  /* Temperature */
#endif
#ifdef EOS_CARRIES_YE
    MyFloat Ye;                           /* Electron fraction */
#endif
#ifdef EOS_CARRIES_ABAR
    MyFloat Abar;                         /* Average atomic weight (in atomic mass units) */
#endif
#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)
    int CompositionType;                  /* define the composition of the material */
#endif
#ifdef EOS_ELASTIC
    MyDouble Elastic_Stress_Tensor[3][3]; /* deviatoric stress tensor */
    MyDouble Elastic_Stress_Tensor_Pred[3][3];
    MyDouble Dt_Elastic_Stress_Tensor[3][3];
#endif
#endif

#ifdef NUCLEAR_NETWORK
    MyDouble Temperature;
    MyDouble xnuc[EOS_NSPECIES], dxnuc[EOS_NSPECIES];
#endif

#if defined(WAKEUP) && !defined(AGS_HSML_CALCULATION_IS_ACTIVE)
    short int wakeup;                     /*!< flag to wake up particle */
#endif

#if defined(OUTPUT_COOLRATE_DETAIL) && defined(COOLING)
    MyFloat CoolingRate;
    MyFloat HeatingRate;
    MyFloat NetHeatingRateQ;
    MyFloat HydroHeatingRate;
    MyFloat MetalCoolingRate;
#endif

#if defined(COOLING) && defined(COOL_GRACKLE)
#if (COOL_GRACKLE_CHEMISTRY >= 1)
    gr_float grHI;
    gr_float grHII;
    gr_float grHM;
    gr_float grHeI;
    gr_float grHeII;
    gr_float grHeIII;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 2)
    gr_float grH2I;
    gr_float grH2II;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 3)
    gr_float grDI;
    gr_float grDII;
    gr_float grHDI;
#endif
#endif

#ifdef TURB_DIFF_DYNAMIC
  MyDouble VelShear_bar[3][3];
  MyDouble MagShear_bar;
  MyDouble Velocity_bar[3];
  MyDouble Velocity_hat[3];
  MyFloat FilterWidth_bar;
  MyFloat MaxDistance_for_grad;
  MyDouble Norm_hat;
  MyDouble Dynamic_numerator;
  MyDouble Dynamic_denominator;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
  MyDouble TD_DynDiffCoeff_error;
  MyDouble TD_DynDiffCoeff_error_default;
#endif
#endif

}
  *SphP,				/*!< holds SPH particle data on local processor */
  *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


extern peanokey *DomainKeyBuf;

/* global state of system
*/
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


/* Various structures for communication during the gravity computation.
 */

extern struct data_index
{
  int Task;
  int Index;
  int IndexGet;
}
 *DataIndexTable;		/*!< the particles to be exported are grouped
                                  by task-number. This table allows the
                                  results to be disentangled again and to be
                                  assigned to the correct particle */

extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
*DataNodeList;

extern struct gravdata_in
{
    MyFloat Pos[3];
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING)
    MyFloat Mass;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(BH_DYNFRICTION_FROMTREE)
    MyFloat Vel[3];
#endif
    int Type;
#if defined(BH_DYNFRICTION_FROMTREE)
    MyFloat BH_Mass;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    MyFloat Soft;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    MyFloat AGS_zeta;
#endif
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
    MyFloat min_bh_t_orbital;   /*!<orbital time for binary */
    MyDouble comp_dx[3];        /*!< position of binary companion */
    MyDouble comp_dv[3];        /*!< velocity of binary companion */
    MyDouble comp_Mass;         /*!< mass of binary companion */
    int is_in_a_binary;
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    int SuperTimestepFlag;  /*!< 2 if allowed to super-timestep, 1 if a candidate for super-timestepping, 0 otherwise */
#endif
    MyFloat OldAcc;
    int NodeList[NODELISTLENGTH];
}
 *GravDataIn,			/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


extern struct gravdata_out
{
    MyLongDouble Acc[3];
#ifdef RT_USE_TREECOL_FOR_NH
    MyDouble ColumnDensityBins[RT_USE_TREECOL_FOR_NH];
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    MyLongDouble TreeMass;
#endif
#ifdef RT_OTVET
    MyLongDouble ET[N_RT_FREQ_BINS][6];
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    MyDouble Rad_E_gamma[N_RT_FREQ_BINS];
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    MyDouble Rad_Flux[N_RT_FREQ_BINS][3];
#endif
#ifdef CHIMES_STELLAR_FLUXES
    double Chimes_G0[CHIMES_LOCAL_UV_NBINS];
    double Chimes_fluxPhotIon[CHIMES_LOCAL_UV_NBINS];
#endif
#ifdef BH_COMPTON_HEATING
    MyLongDouble Rad_Flux_AGN;
#endif
#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    MyLongDouble MencInRcrit;
#endif
#ifdef EVALPOTENTIAL
    MyLongDouble Potential;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    MyLongDouble tidal_tensorps[3][3];
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    MyLongDouble GravJerk[3];
#endif
#ifdef BH_CALC_DISTANCES
    MyFloat min_dist_to_bh;
    MyFloat min_xyz_to_bh[3];
#ifdef SINGLE_STAR_FIND_BINARIES
    MyFloat min_bh_t_orbital; //orbital time for binary
    MyDouble comp_dx[3]; //position of binary companion
    MyDouble comp_dv[3]; //velocity of binary companion
    MyDouble comp_Mass; //mass of binary companion
    int is_in_a_binary; // 1 if star is in a binary, 0 otherwise
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    MyFloat min_bh_freefall_time;    // minimum value of sqrt(R^3 / G(M_BH + M_particle)) as calculated from the tree-walk
    MyFloat min_bh_approach_time; // smallest approach time t_a = |v_radial|/r
#if (SINGLE_STAR_TIMESTEPPING > 0)
    MyLongDouble COM_tidal_tensorps[3][3]; //tidal tensor evaluated at the center of mass without contribution from the companion
    MyDouble COM_GravAccel[3]; //gravitational acceleration evaluated at the center of mass without contribution from the companion
    int COM_calc_flag; //flag that tells whether this was only a rerun to get the acceleration ad the tidal tenor at the center of mass of a binary
    int SuperTimestepFlag; // 2 if allowed to super-timestep, 1 if a candidate for super-timestepping, 0 otherwise
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    MyFloat min_bh_fb_time; // minimum time for feedback to arrive from a star
#endif    
#endif
#endif
}
 *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct potdata_out
{
  MyLongDouble Potential;
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






/*! Header for the standard file format */
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

  char fill[18];		/*!< fills to 256 Bytes */
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
  IO_HSML,
  IO_SFR,
  IO_AGE,
  IO_GRAINSIZE,
  IO_GRAINTYPE,
  IO_HSMS,
  IO_Z,
  IO_BHMASS,
  IO_BHMASSALPHA,
  IO_BH_ANGMOM,
  IO_BHMDOT,
  IO_BHDUSTMASS,
  IO_R_PROTOSTAR,
  IO_MASS_D_PROTOSTAR,
  IO_ZAMS_MASS,
  IO_STAGE_PROTOSTAR,
  IO_AGE_PROTOSTAR,
  IO_LUM_SINGLESTAR,
  IO_BHPROGS,
  IO_BH_DIST,
  IO_ACRB,
  IO_SINKRAD,
  IO_POT,
  IO_ACCEL,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_IDEN,
  IO_INIB,
  IO_UNSPMASS,
  IO_CRATE,
  IO_HRATE,
  IO_NHRATE,
  IO_HHRATE,
  IO_MCRATE,
  IO_DTENTR,
  IO_TSTP,
  IO_BFLD,
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
  IO_COOLRATE,
  IO_TIDALTENSORPS,
  IO_GDE_DISTORTIONTENSOR,
  IO_FLOW_DETERMINANT,
  IO_PHASE_SPACE_DETERMINANT,
  IO_ANNIHILATION_RADIATION,
  IO_STREAM_DENSITY,
  IO_EOSTEMP,
  IO_EOSABAR,
  IO_EOSYE,
  IO_PRESSURE,
  IO_EOSCS,
  IO_EOS_STRESS_TENSOR,
  IO_CBE_MOMENTS,
  IO_EOSCOMP,
  IO_PARTVEL,
  IO_RADGAMMA,
  IO_RAD_ACCEL,
  IO_RAD_FLUX,
  IO_EDDINGTON_TENSOR,
  IO_LAST_CAUSTIC,
  IO_SHEET_ORIENTATION,
  IO_INIT_DENSITY,
  IO_CAUSTIC_COUNTER,
  IO_VDIV,
  IO_VORT,
  IO_DELAYTIME,
  IO_AGS_SOFT,
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
  SIO_DUST,
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

/*
 * Variables for Tree
 * ------------------
 */

extern long Nexport, Nimport;
extern int BufferFullFlag;
extern int NextParticle;
extern int NextJ;
extern int TimerFlag;

// note, the ALIGN(32) directive will effectively pad the structure size
// to a multiple of 32 bytes
extern ALIGN(32) struct NODE
{
  MyFloat center[3];		/*!< geometrical center of node */
  MyFloat len;			/*!< sidelength of treenode */

  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      MyFloat s[3];		/*!< center of mass of node */
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
#ifdef RT_USE_TREECOL_FOR_NH
  MyFloat gasmass;
#endif
#ifdef RT_USE_GRAVTREE
  MyFloat stellar_lum[N_RT_FREQ_BINS]; /*!< luminosity in the node*/
#ifdef CHIMES_STELLAR_FLUXES
  double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS];
  double chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#endif

#ifdef BH_PHOTONMOMENTUM
    MyFloat bh_lum;		    /*!< luminosity of BHs in the node */
    MyFloat bh_lum_grad[3];	/*!< gradient vector for gas around BH (for angular dependence) */
#endif
    

#ifdef BH_CALC_DISTANCES
  MyFloat bh_mass;      /*!< holds the BH mass in the node.  Used for calculating tree based dist to closest bh */
  MyFloat bh_pos[3];    /*!< holds the mass-weighted position of the the actual black holes within the node */
#if defined(SINGLE_STAR_TIMESTEPPING)
  MyFloat bh_vel[3];    /*!< holds the mass-weighted avg. velocity of black holes in the node */
  int N_BH;             /*!< holds the number of BH particles in the node. Used for refinement/search criteria */
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
  MyFloat MaxFeedbackVel;
#endif
#endif    
#endif

#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MyFloat rt_source_lum_s[3];     /*!< center of luminosity for sources in the node*/
#endif

  MyFloat maxsoft;		/*!< hold the maximum gravitational softening of particle in the node */

#ifdef DM_SCALARFIELD_SCREENING
  MyFloat s_dm[3];
  MyFloat mass_dm;
#endif
}
 *Nodes_base,			/*!< points to the actual memory allocated for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
				   gives the first allocated node */


extern struct extNODE
{
  MyLongDouble dp[3];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MyLongDouble rt_source_lum_dp[3];
    MyFloat rt_source_lum_vs[3];
#endif
#ifdef DM_SCALARFIELD_SCREENING
  MyLongDouble dp_dm[3];
  MyFloat vs_dm[3];
#endif
  MyFloat vs[3];
  MyFloat vmax;
  MyFloat hmax;			/*!< maximum gas kernel length in node. Only used for gas particles */
  MyFloat divVmax;
  integertime Ti_lastkicked;
  int Flag;
}
 *Extnodes, *Extnodes_base;


extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


extern int *Nextnode;		/*!< gives next node in tree walk  (nodes array) */
extern int *Father;		/*!< gives parent node in tree (Prenodes array) */

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

#endif  /* ALLVARS_H  - please do not put anything below this line */
