/*! \file precompiler_logic.h
 *  \brief logic statements that need to be near the top of variable definitions.
 */

/*------- Things that are always recommended (this must follow loading GIZMO_config.h!) -------*/
/* also many logical options to force 'parent' or 'top-level' flags to be enabled for the appropriate methods, if we have enabled something using those methods */

/* macro definition we will use throughout below, make sure this is defined up-top here */
#ifdef NOGRAVITY
#define SELFGRAVITY_OFF
#endif
#define DO_PREPROCESSOR_EXPAND_(VAL)  VAL ## 1
#define EXPAND_PREPROCESSOR_(VAL)     DO_PREPROCESSOR_EXPAND_(VAL) /* checks for a NON-ZERO value of this parameter */
#define CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(VAL) !(EXPAND_PREPROCESSOR_(VAL) == 1) /* returns True if a non-zero int value of VAL is set */

/* set default slope limiters */
#if !defined(SLOPE_LIMITER_TOLERANCE)
#define SLOPE_LIMITER_TOLERANCE 1
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

#if defined(PARTICLE_MERGE_SPLIT_EVERY_TIMESTEP) && !defined(DISABLE_MAINTAIN_TREE_IN_REARRANGE)
#define MAINTAIN_TREE_IN_REARRANGE
#endif

#if (defined(HYDRO_DENSITY_SPH) || defined(HYDRO_PRESSURE_SPH)) && !defined(HYDRO_SPH)
#define HYDRO_SPH               /* top-level flag for SPH: must be enabled if any SPH method is used */
#endif
#ifdef HYDRO_SPH
#if !defined(SPH_DISABLE_CD10_ARTVISC) && !(defined(EOS_TILLOTSON) || defined(EOS_ELASTIC) || defined(EOS_ANEOS)) // fancy viscosity switches assume positive pressures //
#define SPHAV_CD10_VISCOSITY_SWITCH        /* Enables Cullen & Dehnen 2010 'inviscid sph' (viscosity suppression outside shocks) */
#endif
#ifndef SPH_DISABLE_PM_CONDUCTIVITY
#define SPHAV_ARTIFICIAL_CONDUCTIVITY      /* Enables mixing entropy (J.Read's improved Price-Monaghan conductivity with Cullen-Dehnen switches) */
#endif
#endif

#if defined(DILATION_FOR_STELLAR_KINEMATICS_ONLY) && !defined(USE_TIMESTEP_DILATION_FOR_ZOOMS)
#define USE_TIMESTEP_DILATION_FOR_ZOOMS
#endif
#if defined(SPECIAL_POINT_WEIGHTED_MOTION) && !defined(SPECIAL_POINT_MOTION)
#define SPECIAL_POINT_MOTION (SPECIAL_POINT_WEIGHTED_MOTION) /* doesn't seem to be working ??? */
#endif
#if defined(SPECIAL_POINT_MOTION) && !defined(SINK_CALC_DISTANCES)
#define SINK_CALC_DISTANCES
#endif

#if defined(EOS_ELASTIC)
#if !defined(DISABLE_SURFACE_VOLCORR) && !defined(HYDRO_KERNEL_SURFACE_VOLCORR)
#define HYDRO_KERNEL_SURFACE_VOLCORR
#endif
#if !defined(DISABLE_EXPLICIT_VOLUME_INTEGRATION) && !defined(HYDRO_EXPLICITLY_INTEGRATE_VOLUME)
#define HYDRO_EXPLICITLY_INTEGRATE_VOLUME
#endif
#endif

#if (defined(EOS_TILLOTSON) || defined(EOS_ELASTIC) || defined(EOS_HELMHOLTZ) || defined(EOS_ANEOS) || defined(COSMIC_RAY_FLUID) || defined(RT_RADPRESSURE_IN_HYDRO) || defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE) || defined(EOS_GMC_BAROTROPIC) || defined(COSMIC_RAY_SUBGRID_LEBRON)) && !defined(EOS_GENERAL)
#define EOS_GENERAL
#endif

#if defined(CBE_INTEGRATOR) || defined(DM_FUZZY)
#define AGS_FACE_CALCULATION_IS_ACTIVE
#endif

#if defined(CBE_INTEGRATOR)
#define CBE_INTEGRATOR_NBASIS CBE_INTEGRATOR
#ifdef CBE_INTEGRATOR_SECONDMOMENT
#if (BOX_SPATIAL_DIMENSION==1) || defined(ONEDIM)
#define CBE_INTEGRATOR_NMOMENTS 3  /* [0-norm,1-mom,1-second] */
#elif (BOX_SPATIAL_DIMENSION==2) || defined(TWODIMS)
#define CBE_INTEGRATOR_NMOMENTS 6  /* [0-norm,2-mom,3-symm-tensor] */
#else
#define CBE_INTEGRATOR_NMOMENTS 10  /* 10 non-trivial moments [0th/norm, 3-mom, 6-symm-tensor] */
#endif
#else
#define CBE_INTEGRATOR_NMOMENTS ((BOX_SPATIAL_DIMENSION)+1)
#endif
#endif

#if defined(GRAIN_COLLISIONS)
#define DM_SIDM 8 /* use the SIDM module to handle scattering of otherwise-collisionless particles against each other -- set to Particle Type=3 here */
#endif
#if defined(PIC_MHD)
#define PIC_MHD_NEW_RSOL_METHOD /* prefer new method for dealing with RSOL, should make simulations easier if done correctly */
#ifdef GRAIN_FLUID
#define GRAIN_FLUID_AND_PIC_BOTH_DEFINED /* keyword used later to know what we need to read in */
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

#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(CBE_INTEGRATOR) || defined(DM_FUZZY) || defined(AGS_FACE_CALCULATION_IS_ACTIVE) || defined(DM_SIDM)
#define AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORALL)
#define ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING /* comment out to revert to behavior of taking 'greater' softening in pairwise kernel interactions with adaptive softenings enabled. really only needed currently for this particular AGS model given how it computes zeta terms (could be made optional with one more loop for those as well) */
#endif

#if defined(COOL_GRACKLE) && !defined(COOLING)
#define COOLING
#endif



/* set 'default' options for FIRE module packages as a whole */
#ifdef FIRE_PHYSICS_DEFAULTS

#define COOLING                             /*! top-level switch for cooling */
#define COOL_LOW_TEMPERATURES               /*! include low-temperature (<1e4 K) cooling */
#define COOL_METAL_LINES_BY_SPECIES         /*! include high-temperature metal-line cooling, species-by-species */
#define GALSF                               /*! top-level switch for star formation */
#define METALS                              /*! follow metals as passive scalars, use in cooling, etc */
#define TURB_DIFF_METALS                    /*! explicit sub-grid diffusivity for metals/passive scalars */
#define TURB_DIFF_METALS_LOWORDER           /*! memory-saving custom mod */
#define STOP_WHEN_BELOW_MINTIMESTEP         /*! this is general good practice */
#if !defined(MULTIPLEDOMAINS)
#define MULTIPLEDOMAINS 32                  /*! slightly closer to our usual default, but users should feel free to adjust */
#endif
#if !defined(DOMAIN_LIGHTWEIGHT_REPARTITION)
#define DOMAIN_LIGHTWEIGHT_REPARTITION      /*! enable lightweight domain repartition that reuses the top tree between full decompositions */
#endif
#define OPENMP_TREE_UPDATE

#define GALSF_FB_MECHANICAL                 /*! top-level switch for mechanical feedback modules */
#define GALSF_FB_FIRE_STELLAREVOLUTION (FIRE_PHYSICS_DEFAULTS) /*! turns on default FIRE processes+lookup tables including gas return, SNe, R-process, etc. this carries a number matching the defaults set you choose */
#if !(defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID))
/* Sub-grid FIRE RT modules: only enabled when no explicit RT solver is active. When a real RT
   method is on, photoionization / local radiation pressure / long-range starlight propagation
   are handled by the explicit solver -- enabling both produces duplicate symbols and
   double-counted physics. */
#define GALSF_FB_FIRE_RT_HIIHEATING         /*! gas within HII regions around young stars is photo-heated to 10^4 K - local stromgren approximation */
#define GALSF_FB_FIRE_RT_LOCALRP            /*! turn on local radiation pressure coupling to gas - account for local multiple-scattering and isotropic local absorption */
#define GALSF_FB_FIRE_RT_LONGRANGE          /*! continuous acceleration from starlight (uses luminosity tree) to propagate FIRE RT */
#endif
#define GALSF_FB_FIRE_AGE_TRACERS 16        /*! tracks a set of passive scalars corresponding to stellar ages for chemical evolution model postprocessing */

#if !(defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL))
#define ADAPTIVE_GRAVSOFT_FORGAS            /*! default choice is adaptive force softening for gas, but not stars [since ambiguously defined] */
#endif

#if (FIRE_PHYSICS_DEFAULTS <= 2)
#define GALSF_USE_SNE_ONELOOP_SCHEME // set to use the 'base' FIRE-2 SNe coupling. if commented out, will user newer version that more accurately manages the injected energy with neighbors moving to inject a specific target
#define GALSF_SFR_CRITERION (0+256) // 0=density threshold, 1=virial criterion (strict), 2=convergent flow, 4=local extremum, 8=no sink in kernel, 16=not falling into sink, 32=hill (tidal) criterion, 64=Jeans criterion, 128=converging flow along all principle axes, 256=self-shielding/molecular, 512=multi-free-fall (smooth dependence on virial), 1024='catch' for un-resolvable densities
#define GALSF_SFR_VIRIAL_SCALING (-1) // old threshold modification with exception when densities 100x above threshold density
#endif

#if (FIRE_PHYSICS_DEFAULTS == 3)
#define COOL_MOLECFRAC_NONEQM
#define OUTPUT_MOLECULAR_FRACTION
#define OUTPUT_COOLRATE
#if !(defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID))
#define RT_USE_GRAVTREE_SAVE_RAD_FLUX  /* LEBRON-style gravtree flux save; mutex with explicit RT solvers (same reasoning as the GALSF_FB_FIRE_RT_* / SINK_*_HEATING blocks) */
#endif
#define OUTPUT_POTENTIAL
#define GALSF_SFR_CRITERION (0+1+2+64) // 0=density threshold, 1=virial criterion (strict), 2=convergent flow, 4=local extremum, 8=no sink in kernel, 16=not falling into sink, 32=hill (tidal) criterion, 64=Jeans criterion, 128=converging flow along all principle axes, 256=self-shielding/molecular, 512=multi-free-fall (smooth dependence on virial), 1024='catch' for un-resolvable densities
#define ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT (0.1/UNIT_LENGTH_IN_KPC)
#define GALSF_SFR_IMF_SAMPLING /* use the IMF-sampling discrete number of O-star scheme, no penalty at low mass-res */
#define FIRE_SNE_ENERGY_METAL_DEPENDENCE_EXPERIMENT (1) /* ramp the SNe rate and massive stellar feedback fraction of total mass (essentially L/M) at low metallicities, leaves no dwarf stars below [Z/H]<-7 or so ramping down to -5 */
#endif // defaults = 3

#if defined(FIRE_MHD)
#define MAGNETIC            /* top-level flag */
#define MHD_B_SET_IN_PARAMS /* B-field must be set in ICs */
#define CONDUCTION          /* enable conduction */
#define CONDUCTION_SPITZER  /* compute proper coefficients and anisotropy for conduction */
#define VISCOSITY           /* enable viscosity */
#define VISCOSITY_BRAGINSKII /* compute proper coefficients and anisotropy for viscosity */
#endif // FIRE_MHD

#if defined(FIRE_CRS)
#if !defined(CRFLUID_SPEEDOFLIGHT_REDUCTION)
#define CRFLUID_SPEEDOFLIGHT_REDUCTION (1.0e8/C_LIGHT_CGS) /*! maximum CR transport speed: 1000 km/s safe for our default diffusivities in constant/variable-kappa model */
#endif
#if (FIRE_CRS == -2)
#define COSMIC_RAY_SUBGRID_LEBRON   /*! this simply uses the sub-grid model */
#else
#define COSMIC_RAY_FLUID            /*! use 'explicit' CR integration in one of the code formulations */
#endif // closes whether to do cr fluid or subgrid
#if !defined(CRFLUID_EVOLVE_SPECTRUM) /*! check to enable flags for which CR spectrum or single-bin to evolve */
#if (FIRE_CRS >= 0) || (FIRE_CRS <= 1)
#define CRFLUID_EVOLVE_SPECTRUM 1   /*! evolve proton + electron spectrum */
#endif
#if (FIRE_CRS >= 2)
#define CRFLUID_EVOLVE_SPECTRUM 2   /*! evolve full set of 10 different CR species */
#endif
#endif // closes which spectrum check
#if !defined(CRFLUID_DIFFUSION_MODEL)
#if (FIRE_CRS <= 0)
#define CRFLUID_DIFFUSION_MODEL 0   /*! constant diffusivity (set by params file) */
#else
#define CRFLUID_DIFFUSION_MODEL 8   /*! best-guess for an empirical model which varies as a function of CR gradient scale length from ISM-to-CGM, but not with real plasma properties. lots of options for that. */
//#define CRFLUID_DIFFUSION_MODEL 6   /*! best-guess for variable-kappa model, combining updated SC+ET */
//#define CRFLUID_SET_SC_MODEL (7)    /*! set mode for SC model using best-estimate of fQLT and fCAS, and best model for extrinsic driving of CRs */
//#define CRFLUID_SET_ET_MODEL (-1)   /*! set mode for ET model using best-estimate of fturb from Alfven-wave scattering */
#endif
#endif // closes which diffusion model[s] to enable
#endif // FIRE_CRS

#if defined(FIRE_BHS)
#define SINK_PARTICLES                 /* top-level flag */
#define SINK_SEED_FROM_LOCALGAS       /* seed BHs locally in SF-ing gas */
#define SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA /* use the total surface-density criterion, not just gas */
#define SINK_CALC_DISTANCES           /* use this for various checks, particularly in seeding */
#if (defined(GALSF_SFR_IMF_SAMPLING) || (FIRE_PHYSICS_DEFAULTS > 2)) && !defined(SINK_REPOSITION_ON_POTMIN) && !defined(SINK_DYNFRICTION_FROMTREE)
#define SINK_DYNFRICTION_FROMTREE /* use the dynamical friction model instead of pinning/forcing BHs to potential minimum */
#endif
#if !defined(SINK_REPOSITION_ON_POTMIN) && !defined(SINK_DYNFRICTION_FROMTREE)
#define SINK_REPOSITION_ON_POTMIN   /* anchor BHs to centers smoothly */
#endif
#if (FIRE_PHYSICS_DEFAULTS > 2)
#define SINK_SCALE_SPAWNINGMASS_WITH_INITIALMASS /* purely a convention-choice when doing spawning, to use fraction of original BH mass -- this one more useful if using multi-resolution (hyper-refinement) techniques */
#endif
#define SINK_SWALLOWGAS               /* allow BHs to accrete in principle */
#if !defined(SINK_GRAVACCRETION)
#define SINK_GRAVACCRETION 1          /* accrete following our standard gravitational torques model */
#define SINK_GRAVACCRETION_STELLARFBCORR          /* account for additional acceleration-dependent retention from stellar FB in Mdot */
#endif
#if !defined(SINK_ALPHADISK_ACCRETION)
#define SINK_ALPHADISK_ACCRETION (1.e10) /* smooth out accretion + allow super-eddington capture with alpha-disk model */
#endif
#if !(defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID))
/* AGN-source counterparts of the FIRE sub-grid stellar RT modules. They share the
   LEBRON/gravity-tree numerical infrastructure (RT_USE_GRAVTREE,
   valid_gas_particle_for_rt accumulators, etc.) and are therefore mutually exclusive
   with explicit RT solvers, exactly as their stellar GALSF_FB_FIRE_RT_* counterparts. */
#define SINK_PHOTONMOMENTUM           /* allow AGN radiation pressure */
#define SINK_COMPTON_HEATING          /* allow Compton heating from AGN spectrum */
#define SINK_HII_HEATING              /* allow photo-ionization heating from AGN spectrum */
#endif
#define SINK_FB_COLLIMATED            /* BHFB directed along collimated axis following BH ang. mom */
#define SINK_WIND_SPAWN (2)           /* spawn module: N=min num spawned/step */
#if !defined(DISABLE_MAINTAIN_TREE_IN_REARRANGE)
#define MAINTAIN_TREE_IN_REARRANGE    /* avoid constant domain decompositions in bottom timebin each time a spawn occurs */
#endif
#if defined(COSMIC_RAY_FLUID) || defined(COSMIC_RAY_SUBGRID_LEBRON)
#define SINK_COSMIC_RAYS              /* allow CR injection from AGN */
#endif
#endif // FIRE_BHS

#if defined(PMGRID)
#if !defined(PM_PLACEHIGHRESREGION)
#if (FIRE_PHYSICS_DEFAULTS == 3)
#define PM_PLACEHIGHRESREGION 2 /* 2 -- if we want to have potentially resolved hyper-velocity stars/bhs/cells, need to avoid them restructuring the PMGRID completely */
#else
#if defined(SINK_PARTICLES)
#define PM_PLACEHIGHRESREGION 51 /* 1+2+16+32 */
#else
#define PM_PLACEHIGHRESREGION 19 /* 1+2+16 */
#endif
#endif
#endif
#endif // PMGRID check

#endif // FIRE_PHYSICS_DEFAULTS clauses


/* SINK_COMPTON_HEATING reads the same angle-weighted sink luminosity tree data as
 * SINK_PHOTONMOMENTUM, so it needs that infrastructure compiled in; the momentum kick itself can
 * still be disabled at runtime (PhotonMomentum_Coupled_Fraction). Kept outside the
 * FIRE_PHYSICS_DEFAULTS block so it applies to any Config that sets SINK_COMPTON_HEATING. */
#if defined(SINK_COMPTON_HEATING) && !defined(SINK_PHOTONMOMENTUM)
#define SINK_PHOTONMOMENTUM
#endif


#ifdef GALSF_SFR_CRITERION // flag for pure cross-compatibility [identical functionality, just ease-of-use for galaxy simulators here]
#define SINGLE_STAR_SINK_FORMATION GALSF_SFR_CRITERION
#endif


/* set default modules for different cosmic ray transport settings/packages */
#ifdef COSMIC_RAY_FLUID
#if !defined(CRFLUID_SPEEDOFLIGHT_REDUCTION)
#define CRFLUID_SPEEDOFLIGHT_REDUCTION (1.0)
#endif
#if !defined(CRFLUID_EVOLVE_SPECTRUM)
#define GAMMA_COSMICRAY(k) (4.0/3.0)
#endif
#ifndef CRFLUID_DIFFUSION_MODEL
#define CRFLUID_DIFFUSION_MODEL 0
#endif
#if defined(CRFLUID_EVOLVE_SPECTRUM)
#define GAMMA_COSMICRAY(k) ((4.0+gamma_eos_of_crs_in_bin(k))/3.0)
#if (CRFLUID_EVOLVE_SPECTRUM == 2)
#define N_CR_PARTICLE_BINS 70       /*<! set default bin number here -- needs to match hard-coded list in function 'CR_spectrum_define_bins', for now> */
#define N_CR_PARTICLE_SPECIES 8     /*<! total number of CR species to be evolved. must be set here because of references below*/
#else
#define N_CR_PARTICLE_BINS 19       /*<! set default bin number here -- needs to match hard-coded list in function 'CR_spectrum_define_bins', for now> */
#define N_CR_PARTICLE_SPECIES 2     /*<! total number of CR species to be evolved. must be set here because of references below*/
#endif
#endif
#ifndef N_CR_PARTICLE_BINS
#define N_CR_PARTICLE_BINS 1
#endif
#if (N_CR_PARTICLE_BINS > 2) && !defined(CRFLUID_EVOLVE_SPECTRUM)
#define CRFLUID_EVOLVE_SPECTRUM
#endif
#endif


#ifdef COSMIC_RAY_SUBGRID_LEBRON
#define N_CR_PARTICLE_BINS 1
#endif





/* set default options for hybrid single star and ssp model packages */
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL_DEFAULTS) /* options for hybrid/combined FIRE+STARFORGE simulations */
#define SINGLE_STAR_STARFORGE_DEFAULTS   /* parent flag enabling the STARFORGE modules, which themselves enable the FIRE modules */
#define SINGLE_STAR_AND_SSP_HYBRID_MODEL (SINGLE_STAR_AND_SSP_HYBRID_MODEL_DEFAULTS) /* do single-star routines below this mass resolution in solar, FIRE-like above */
#define GALSF_SFR_IMF_SAMPLING           /* use discrete sampling of 'number of O-stars' so we can handle the intermediate-mass regime in at least a simple approximate manner */
#define COOLING              /* only physical if include cooling for both sides, using same cooling functions */
#define MAGNETIC             /* enable MHD, important for systems here */
#define CONDUCTION           /* enable conduction */
#define CONDUCTION_SPITZER   /* compute proper coefficients and anisotropy for conduction */
#define VISCOSITY            /* enable viscosity */
#define VISCOSITY_BRAGINSKII /* compute proper coefficients and anisotropy for viscosity */
#define SINGLE_STAR_FB_JETS  /* enable jets from protostars */
#define SINGLE_STAR_FB_WINDS 0 /* enable continuous mass-loss feedback - will also enable ssp mass-loss */
#define SINGLE_STAR_FB_SNE   /* enable SNe feedback - will also enable ssp mechanical feedback */
#define SINGLE_STAR_FB_RAD   /* enable RHD feedback */
#define RT_COMOVING          /* significantly more stable and accurate formulation given the structure of the problem and method we use */
#define RT_SOURCES (16+32)   /* need to allow -both- ssp-particles and single-star particles to emit */
#if !defined(RT_SPEEDOFLIGHT_REDUCTION)
#define RT_SPEEDOFLIGHT_REDUCTION (0.1)   /* for many problems on these scales, need much larger RSOL than default starforge values (dynamical velocities are big, without this they will severely lag behind) */
#endif
#define ADAPTIVE_TREEFORCE_UPDATE (0.0625) /* rough typical value we use for ensuring stability */
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
#define OUTPUT_ACCELERATION
#define OUTPUT_HYDROACCELERATION
#define OUTPUT_MOLECULAR_FRACTION
#define OUTPUT_TEMPERATURE
#define OUTPUT_ADDITIONAL_RUNINFO
#endif
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
#define PARTICLE_EXCISION
#define OUTPUT_GRADIENT_RHO
#define OUTPUT_GRADIENT_VEL
#define OUTPUT_RT_RAD_FLUX
#define OUTPUT_RT_RAD_OPACITY
#define RT_RAD_PRESSURE_OUTPUT
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES
#define GRAVITY_ANALYTIC
#endif
#endif
#endif // closes hybrid FIRE+STARFORGE model settings




/* set default options for STARFORGE module packages */
#ifdef SINGLE_STAR_STARFORGE_DEFAULTS /* bunch of options -NOT- strictly required here, but this is a temporary convenience block */
#define INPUT_POSITIONS_IN_DOUBLE
#define OUTPUT_POTENTIAL
#ifndef SINGLE_STAR_AND_SSP_HYBRID_MODEL
#define IO_GRADUAL_SNAPSHOT_RESTART
#endif
#ifndef IO_SINKS_ONLY_SNAPSHOT_FREQUENCY
#define IO_SINKS_ONLY_SNAPSHOT_FREQUENCY 0 /* Determines the number of snapshots with reduced data (stars only) per full snapshots (gas+stars), e.g., setting it to 2 means 2/3 of the snapshots will be reduced, 1/3 will have full data. Setting it to 0 disables it.  */
#endif
#define SINGLE_STAR_SINK_DYNAMICS
#if !defined(PMGRID) && !defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) && !defined(USE_TIMESTEP_DILATION_FOR_ZOOMS) && !defined(DISABLE_HERMITE_INTEGRATION)
#define HERMITE_INTEGRATION 32 // bitflag for which particles to do 4th-order Hermite integration
#endif
#define ADAPTIVE_GRAVSOFT_FORGAS
#define GRAVITY_ACCURATE_FEWBODY_INTEGRATION
#define SINGLE_STAR_TIMESTEPPING 0
#define SINGLE_STAR_ACCRETION 12
#define SINGLE_STAR_SINK_FORMATION (0+1+2+4+8+16+32+64+2048) // 0=density threshold, 1=virial criterion, 2=convergent flow, 4=local extremum, 8=no sink in kernel, 16=not falling into sink, 32=hill (tidal) criterion, 64=Jeans criterion, 128=converging flow along all principle axes, 256=self-shielding/molecular, 512=multi-free-fall (smooth dependence on virial), 1024=numerical escape if too dense, 2048=virial is time-averaged
#ifndef SINGLE_STAR_DIRECT_GRAVITY_RADIUS
#define SINGLE_STAR_DIRECT_GRAVITY_RADIUS (1000.) // distance inside of which star-star gravitational interactions are calculated exactly, in AU
#endif
#ifndef ADAPTIVE_TREEFORCE_UPDATE
#define ADAPTIVE_TREEFORCE_UPDATE (0.0625) // optimization
#endif
#ifndef SINK_MDOT_RESERVOIR_CAP
#define SINK_MDOT_RESERVOIR_CAP 30 /* throttle the reservoir drain of marginally-resolved sinks; see set_sink_mdot in sinks/sink.cc */
#endif
#if !defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM) && !defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT)
#ifndef IO_SUPPRESS_TIMEBIN_STDOUT
#define IO_SUPPRESS_TIMEBIN_STDOUT 16 // only prints outputs to log file if the highest active timebin index is within n of the highest timebin (dt_bin=2^(-N)*dt_bin,max)
#endif
#endif
#define OUTPUT_SINK_ACCRETION_HIST // save accretion histories
#define OUTPUT_SINK_FORMATION_PROPS // save at-formation properties of sink particles
#define INPUT_READ_SINKPROPS
#if ( defined(STARFORGE_GMC_TURBINIT) || defined(STARFORGE_FILAMENT_TURBINIT) ) // these flags should be given numerical values equal to the desired virial parameter
#define TURB_DRIVING
#define GRAVITY_ANALYTIC
#define SELFGRAVITY_OFF
#endif
#define SINK_ALPHADISK_ACCRETION (1.0e6)
#ifdef GRAIN_FLUID
#define SINK_GRAVCAPTURE_NONGAS
#endif
#ifdef MAGNETIC
#define MHD_CONSTRAINED_GRADIENT 1
#endif
#if ( defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_RT_HEATING) || defined(SINGLE_STAR_FB_SNE) || defined(SINGLE_STAR_FB_RAD) || defined(SINGLE_STAR_FB_LOCAL_RP))
#define SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION 2 //we are using the protostellar evolution model from ORION
#endif
#if ( defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_SNE) ) //enable diffusion for metals and enable tracers for different feedback channels
#define STARFORGE_FEEDBACK_TRACERS 3 // 0 for jets, 1 for winds, 2 for SNe
#define TURB_DIFF_METALS
#define TURB_DIFF_METALS_LOWORDER
#endif
#ifdef SINGLE_STAR_FB_RAD
#define SINGLE_STAR_RT_DEFAULTS
#define RT_OPTICAL_NIR
#define RT_NUV
#define RT_PHOTOELECTRIC
#ifndef RT_CHEM_PHOTOION
#define RT_CHEM_PHOTOION 1
#endif
#define RT_INFRARED
#if !defined(RT_ISRF_BACKGROUND) && !defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
#define RT_ISRF_BACKGROUND
#endif
#endif
#if (defined(COOLING) && !defined(COOL_LOWTEMP_THIN_ONLY) && !defined(RT_INFRARED))
#define RT_USE_TREECOL_FOR_NH 6 /* This gives a better approximation for column density than the usual scale-length estimator, but is overkill for typical 1e-3msun-resolving simulations that only marginally resolve the opacity limit. Enable for high (<1e-5msun) resolution sims */
#endif
#ifdef COOLING
#define SIMPLE_STEADYSTATE_CHEMISTRY
#define COOL_MOLECFRAC_NONEQM
#define EOS_SUBSTELLAR_ISM
#define OUTPUT_MOLECULAR_FRACTION
#if defined(MAGNETIC) && !defined(CONDUCTION) && !defined(VISCOSITY) // if we have cooling and magnetic fields, enable conduction + viscosity
#define CONDUCTION           /* enable conduction */
#define CONDUCTION_SPITZER   /* compute proper coefficients and anisotropy for conduction */
#define VISCOSITY            /* enable viscosity */
#define VISCOSITY_BRAGINSKII /* compute proper coefficients and anisotropy for viscosity */
#define DIFFUSION_OPTIMIZERS
#endif // MAGNETIC
#if !defined(RT_ISRF_BACKGROUND) && !defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
#define RT_ISRF_BACKGROUND  // Draine 1978 ISRF for photoelectric heating (appropriate for solar circle, must be re-scaled for different environments)
#endif
#endif // COOLING
#if defined(SINGLE_STAR_FB_WINDS) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION) || defined(COOLING)
#define GALSF_FB_FIRE_STELLAREVOLUTION 3 // enable multi-loop feedback from such sources [this is specific to the DG-MG implementations here, not for public use right now!]. for now set to =2, which should force the code version to match previous iterations, as compared to the newer implementations.
#endif
#if defined(RT_ISRF_BACKGROUND)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(RT_ISRF_BACKGROUND)
#if (RT_ISRF_BACKGROUND <= 0)
#undef RT_ISRF_BACKGROUND /* use the negative or zero value above as a key to specifically -undefine- this variable, otherwise it will cause problems below by calling 0 to reset quantities it should not */
#endif
#endif
#endif
#endif // closes SINGLE_STAR_STARFORGE_DEFAULTS settings

#ifdef SINGLE_STAR_RT_DEFAULTS
#define RT_M1
#define RT_COMOVING
#define OUTPUT_RT_RAD_FLUX
#ifndef IO_SUPPRESS_OUTPUT_EDDINGTON_TENSOR
#define IO_SUPPRESS_OUTPUT_EDDINGTON_TENSOR
#endif
#ifndef RT_SOURCES
#define RT_SOURCES 32
#endif
#ifndef RT_SPEEDOFLIGHT_REDUCTION
#define RT_SPEEDOFLIGHT_REDUCTION (3.0e-4)
#endif
#define RT_REPROCESS_INJECTED_PHOTONS
#define RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION
#define RT_REINJECT_ACCRETED_PHOTONS
#define RT_INJECT_PHOTONS_DISCRETELY
#endif // closes SINGLE_STAR_RT_DEFAULTS settings


#ifdef SINGLE_STAR_SINK_DYNAMICS
#define GALSF // top-level switch needed to enable various frameworks
#define METALS  // metals should be active for stellar return
#define SINK_PARTICLES // need to have sink particles active since these are our sink particles
#define SINK_INTERACT_ON_GAS_TIMESTEP // BH-gas interactions (feedback and accretion) occur with frequency set by the gas timestep
#define SINK_CALC_DISTANCES // calculate distance to nearest sink in gravity tree

#ifdef SINGLE_STAR_ACCRETION // figure out flags needed for the chosen sink accretion model
#define SINK_SWALLOWGAS // need to swallow gas [part of sink model]
#ifndef SINK_ALPHADISK_ACCRETION
#define SINK_ALPHADISK_ACCRETION (2.) // all models will use a 'reservoir' of some kind to smooth out accretion rates (and represent unresolved disk)
#endif
#if (SINGLE_STAR_ACCRETION <= 10)
#define SINK_GRAVACCRETION (SINGLE_STAR_ACCRETION) // use one of these pre-built accretion models
#endif
#if (SINGLE_STAR_ACCRETION == 11)
#define SINK_GRAVCAPTURE_GAS // use gravitational capture swallow criterion for resolved gravitational capture
#endif
#if (SINGLE_STAR_ACCRETION == 12)
#define SINK_GRAVCAPTURE_GAS
#define SINK_GRAVCAPTURE_FIXEDSINKRADIUS // modify grav capture to Bate-style, fixed (in time) sink radius based on SF neighbor distance, plus angular momentum criterion
#endif
#endif

#if (defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_RT_HEATING) || defined(SINGLE_STAR_FB_SNE) || defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID) || defined(SINGLE_STAR_FB_LOCAL_RP)) && defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_SINK_DYNAMICS)
#define SINGLE_STAR_FB_TIMESTEPLIMIT // general flag indicating feedback is on
#endif

#if defined(SINGLE_STAR_FB_RT_HEATING) && !(defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION) || defined(RT_M1) || defined(RT_LOCALRAYGRID))
#define GALSF_FB_FIRE_RT_LONGRANGE  // turn on FIRE RT approximation: no Type-4 particles so don't worry about its approximations
#define SINK_PHOTONMOMENTUM // enable BHs within the FIRE-RT framework.
#define RT_DISABLE_RAD_PRESSURE
#endif

#if (defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_SNE)) && !defined(GALSF_FB_MECHANICAL)
#define GALSF_FB_MECHANICAL // we will use the mechanical wind module for low mass loss rate stars (spawning leads to issues). enable regardless if either the winds or sne module is active
#define GALSF_USE_SNE_ONELOOP_SCHEME
#endif

#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION) && !defined(GALSF_FB_FIRE_STELLAREVOLUTION)
#define GALSF_FB_FIRE_STELLAREVOLUTION 3 /* this particular set of modules depends on the fire modules explicitly now, using them for yields and other infrastructure */
#endif

#if defined(SINGLE_STAR_FB_JETS) || ((defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_SNE)) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION))
#define SINK_WIND_SPAWN (2) // leverage the BHFB model already developed within the FIRE-BHs framework. gives accurate launching of arbitrarily-structured jets.
#if !defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM) && !defined(DISABLE_MAINTAIN_TREE_IN_REARRANGE)
#define MAINTAIN_TREE_IN_REARRANGE // don't rebuild the domains/tree every time a particle is spawned - salvage the existing one by redirecting pointers as needed
#endif
#endif

/* TREE_INTEGRITY_AUDITS must be set explicitly in Config.sh; nothing implies it. It enables fatal
   deep tree audits (exactly-once walk, Father/threading cross-check, moment re-derivation) at
   build/refresh/mutation points. It is deliberately not tied to DEVELOPER_MODE, which many test
   configs set only to expose accuracy parameters. Without it only the always-on checks run
   (force_validate_tree_links and the treebuild bookkeeping reduce). */

#if defined(HERMITE_INTEGRATION) && defined(SINK_PARTICLES) && !defined(MAINTAIN_TREE_IN_REARRANGE)
/* Hermite is the one consumer that walks the STANDING tree while deliberately skipping the rebuild
   it was asked for (the HermiteOnlyFlag gate in gravtree.cc), so any Hermite run whose particle list
   is rearranged mid-step must have the full per-swap maintenance -- the parent-only carry is not
   enough for a tree that keeps being walked. Spawning rearranges, and so does every swallow: the sink
   pass's cleanup rearrange (run.cc) eliminates the victims even with no spawning at all. Without the
   maintenance the Hermite corrector refreshes moments on, and takes sink forces from, a condemned
   tree whose threading still points at the vacated slots. This deliberately OVERRIDES the
   nuclear-zoom exclusion above. (cf. gizmo-cpp a06e0073) */
#ifdef DISABLE_MAINTAIN_TREE_IN_REARRANGE
#warning "DISABLE_MAINTAIN_TREE_IN_REARRANGE overridden: HERMITE_INTEGRATION + SINK_PARTICLES requires tree maintenance (this config cannot serve as the no-MAINTAIN A/B arm)"
#endif
#define MAINTAIN_TREE_IN_REARRANGE
#endif

#if defined(SINGLE_STAR_FB_LOCAL_RP) // use standard angle-weighted local coupling to impart photon momentum from stars
#if !defined(SINK_PHOTONMOMENTUM)
#define SINK_PHOTONMOMENTUM
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


#if (SINGLE_STAR_SINK_FORMATION & 16)
#ifndef SINGLE_STAR_TIMESTEPPING
#define SINGLE_STAR_TIMESTEPPING 0
#endif
#endif
#if (SINGLE_STAR_SINK_FORMATION & 32)
#define GALSF_SFR_TIDAL_HILL_CRITERION
#endif
#if (SINGLE_STAR_SINK_FORMATION & 512)
#define GALSF_SFR_VIRIAL_SCALING 2
#endif
#if (SINGLE_STAR_SINK_FORMATION & 1)
#ifndef GALSF_SFR_VIRIAL_SCALING
#define GALSF_SFR_VIRIAL_SCALING 0
#endif
#endif
#if (SINGLE_STAR_SINK_FORMATION & 2048)
#define GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED
#endif

#ifdef GRAVITY_ACCURATE_FEWBODY_INTEGRATION /* utility flag to enable a few different extra-conservative time-integration flags for gravity */
#if !defined(GRAVITY_HYBRID_OPENING_CRIT)
#define GRAVITY_HYBRID_OPENING_CRIT // use both Barnes-Hut + relative tree opening criterion
#endif
#if !defined(STOP_WHEN_BELOW_MINTIMESTEP)
#define STOP_WHEN_BELOW_MINTIMESTEP // stop when below min timestep to prevent bad timestepping
#endif
#define TIDAL_TIMESTEP_CRITERION // use tidal tensor timestep criterion
#endif
#ifdef HERMITE_INTEGRATION
#define COMPUTE_JERK_IN_GRAVTREE /* needs to be computed in order to do the Hermite integration */
#ifndef HERMITE_CORRECTOR_ITERATIONS
#define HERMITE_CORRECTOR_ITERATIONS 1 /* P(EC)^n: n>1 re-evaluates gravity at the corrected state and re-applies the corrector (Kokubo, Yoshinaga & Makino 1998) */
#endif
#ifndef TIDAL_TIMESTEP_CRITERION
#define TIDAL_TIMESTEP_CRITERION // use tidal tensor timestep criterion -- otherwise won't effectively leverage the Hermite integrator timesteps
#endif
#if defined(HERMITE_SYMMETRIC_TIMESTEP) && !defined(FORCE_EQUAL_TIMESTEPS)
/* inert without HERMITE_INTEGRATION or on equal steps. The code tests this name rather than relying on
   an #undef, which the re-inclusion of GIZMO_config.h (via eos.h) after this file would undo */
#define HERMITE_SYMMETRIC_TIMESTEP_ACTIVE
#endif
#endif



#ifdef RT_USE_TREECOL_FOR_NH
#if !defined(GRAVTREE_CALCULATE_GAS_MASS_IN_NODE)
#define GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
#endif
#endif


#if (defined(SINGLE_STAR_FB_SNE) || defined(SINGLE_STAR_FB_WINDS)) && !defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
#define SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT 2 // determines the maximum number of ejecta particles spawned per timestep, see below
#endif
#if defined(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)
#define SINGLE_STAR_FB_SNE_N_EJECTA (4*(SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)*((SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT)+1)) // maximum number of ejecta cells spawned per timestep - follows from tiling rules for rays from RHD-direct-ray method
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


#ifdef SINGLE_STAR_DIRECT_GRAVITY /* exact O(N_star^2) star-star gravity, replacing the tree for those pairs */
#if !defined(SINK_CALC_DISTANCES)
#define SINK_CALC_DISTANCES /* the tree's node-level sink_mass is what tells the walk which nodes to drop for star targets */
#endif
#if defined(SINGLE_STAR_FIND_BINARIES)
/* SINGLE_STAR_FIND_BINARIES identifies a binary companion during the tree walk, from the very
   star-star node interactions this flag removes, and the direct sum does not reconstruct it.
   Silently losing binary detection would change the integration, so refuse the combination. */
#error "SINGLE_STAR_DIRECT_GRAVITY does not yet supply the companion search that SINGLE_STAR_FIND_BINARIES (and hence SINGLE_STAR_TIMESTEPPING > 0) needs; use one or the other."
#endif
#endif


#ifdef MHD_CONSTRAINED_GRADIENT
/* make sure mid-point gradient calculation for cleaning terms is enabled */
#ifndef MHD_CONSTRAINED_GRADIENT_MIDPOINT
#define MHD_CONSTRAINED_GRADIENT_MIDPOINT
#endif
#endif

#ifdef MHD_MODIFIED_GRADIENT
#ifndef MAGNETIC
#define MAGNETIC
#endif
/* MG builds on top of the CG infrastructure (FaceCrossX, FaceDotB, Face_Area),
   but replaces the iterative correction with a global sparse-matrix solve.
   If user set CG explicitly alongside MG, warn them. Otherwise, auto-enable CG=1. */
#ifdef MHD_CONSTRAINED_GRADIENT
#undef MHD_CONSTRAINED_GRADIENT
#endif
#define MHD_CONSTRAINED_GRADIENT 1
#endif


#if defined(SPECIAL_POINT_MOTION)
#if !defined(SINK_CALC_DISTANCES)
#define SINK_CALC_DISTANCES /* make sure the distance tracking is actually enabled */
#endif
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(SPECIAL_POINT_MOTION)
#define SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES (SPECIAL_POINT_MOTION) /* set the special particle type to be used for tracking */
#endif
#endif
#if defined(SINK_CALC_DISTANCES) && !defined(SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
#define SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES (5) /* default to type = 5 for this module */
#endif

#if defined(GRAVITY_ANALYTIC)
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(GRAVITY_ANALYTIC)
#if (GRAVITY_ANALYTIC > 0)
#define GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE /* ok, analytic gravity is defined with a numerical value > 0, indicating we should use this flag */
#ifndef SINK_CALC_DISTANCES
#define SINK_CALC_DISTANCES
#endif
#endif
#endif
#endif


/* drift and kick the sink node moments between tree rebuilds, as the main moments already are. Requires
   sink_vel to exist to drift with, hence the condition below -- which has to be tested here, after every
   site that can turn SINK_CALC_DISTANCES on. SINGLE_STAR_DIRECT_GRAVITY is in the list because it
   subtracts the sink monopole from mixed nodes and so needs sink_pos on the same clock as u.d.s; it can
   be set without SINGLE_STAR_TIMESTEPPING, and where it is, the sink_vel declaration and the moment sums
   that fill it widen to match this macro. */
#if defined(SINK_CALC_DISTANCES) && (defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES) || defined(SPECIAL_POINT_MOTION) || defined(SINGLE_STAR_DIRECT_GRAVITY))
#define SINK_NODE_MOTION_TRACKED
#endif


#if defined(COOL_MOLECFRAC)
#if (COOL_MOLECFRAC == 6) && !defined(COOL_MOLECFRAC_NONEQM)
#define COOL_MOLECFRAC_NONEQM // estimate molecular fractions for thermochemistry+cooling with explicitly-evolved non-equilibirum H2 formation+destruction with clumping and self-shielding (Hopkins et al arXiv:2203.00040)
#endif
#endif


#if defined(EOS_SUBSTELLAR_ISM) || defined(COOL_MOLECFRAC_NONEQM)
#define EOS_GAMMA_VARIABLE
#endif


#if defined(EOS_GAMMA_VARIABLE)
#ifndef EOS_GENERAL
#define EOS_GENERAL /*! needs to be on for this to work */
#endif
#endif
#define GAMMA_DEFAULT (EOS_GAMMA)


#if defined(CONDUCTION) || defined(EOS_GENERAL) || defined(TURB_DIFF_ENERGY)
#define DOGRAD_INTERNAL_ENERGY 1
#endif


#if defined(EOS_GENERAL)
#define DOGRAD_SOUNDSPEED 1
#endif


#ifndef  GALSF_GENERATIONS
#define  GALSF_GENERATIONS     1    /*!< Number of star particles that may be created per gas particle */
#endif


#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
#define OUTPUT_SOFTENING  /*! output softening to snapshots */
#endif


#ifdef GALSF_SFR_IMF_VARIATION
#define N_IMF_FORMPROPS  13  /*!< formation properties of star particles to record for output */
#endif


#define CRFLUID_REDUCED_C_CODE(k) (return_CRbin_M1speed(k)) // allow for bin-to-bin variations in RSOL


#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif


#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif


#if defined(ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)
#define OUTPUT_TIDAL_TENSOR
#define OUTPUT_SOFTENING
#endif


#ifdef RT_INFRARED
#define COOL_LOWTEMP_THIN_ONLY // don't want to double-count trapping of radiation if we're doing it self-consistently
#endif



/* ----- giant block of options for RHD modules ------ */

/* options for FIRE RT method */
#if defined(GALSF_FB_FIRE_RT_LONGRANGE)
#define RT_LEBRON // this flag requires lebron for rhd
#endif
#if defined(RT_LEBRON)
#define RT_USE_GRAVTREE // use gravity tree for flux propagation
#define RT_USE_GRAVTREE_SAVE_RAD_ENERGY
#if !defined(GALSF_FB_FIRE_RT_LONGRANGE)
#define RADTRANSFER // for cross-compatibility reasons, if the FIRE version is not on, need RADTRANSFER flag also enabled
#define RT_USE_GRAVTREE_SAVE_RAD_FLUX
#endif
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
#if !defined(RT_EVOLVE_FLUX) // evolving energy explicitly but not flux, flux-limiting is not disabled
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

#if !defined(RT_USE_GRAVTREE) && defined(RT_SELFGRAVITY_OFF) && !defined(SELFGRAVITY_OFF)
#define SELFGRAVITY_OFF // safely define SELFGRAVITY_OFF in this case, otherwise we act like there is gravity except in the final setting of accelerations
#endif

/* turn on outputs appropriately */
#ifdef RADTRANSFER
#if !defined(OUTPUT_EDDINGTON_TENSOR) && !defined(IO_SUPPRESS_OUTPUT_EDDINGTON_TENSOR)
#define OUTPUT_EDDINGTON_TENSOR
#endif
#endif


/* set up the 'bins' for different frequencies for the RHD -- there's probably a more elegant way to do this with inline functions
    instead of precompiler logic, but for now this is the solution */
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

#ifdef RT_CHEM_PHOTOION
#define RT_BAND_IS_IONIZING(k) ((k==RT_FREQ_BIN_H0) || (k==RT_FREQ_BIN_He0) || (k==RT_FREQ_BIN_He1) || (k==RT_FREQ_BIN_He2))
#endif

#ifndef GALSF_FB_FIRE_RT_LONGRANGE
#define RT_FREQ_BIN_FIRE_UV (RT_FREQ_BIN_He2+0)
#define RT_FREQ_BIN_FIRE_OPT (RT_FREQ_BIN_FIRE_UV+0)
#define RT_FREQ_BIN_FIRE_IR (RT_FREQ_BIN_FIRE_OPT+0)
#else
#define RT_FREQ_BIN_FIRE_UV (RT_FREQ_BIN_He2+1)
#define RT_FREQ_BIN_FIRE_OPT (RT_FREQ_BIN_FIRE_UV+1)
#define RT_FREQ_BIN_FIRE_IR (RT_FREQ_BIN_FIRE_OPT+1)
#endif

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


/* ----- end block of options for RHD modules ------ */



#if defined(GALSF) || defined(SINK_PARTICLES) || defined(RADTRANSFER) || defined(OUTPUT_DENS_AROUND_NONGAS) || defined(CHIMES) || defined(RT_REPROCESS_INJECTED_PHOTONS)
#define DO_DENSITY_AROUND_NONGAS_PARTICLES
#if !defined(ALLOW_IMBALANCED_GASPARTICLELOAD)
#define ALLOW_IMBALANCED_GASPARTICLELOAD
#endif
#endif

#if defined(SINK_SWALLOWGAS)
#define SINK_FOLLOW_ACCRETED_COM
#define SINK_FOLLOW_ACCRETED_MOMENTUM
#if defined(SINGLE_STAR_SINK_DYNAMICS) || defined(SINK_GRAVCAPTURE_GAS)
#define SINK_FOLLOW_ACCRETED_ANGMOM 0 // follow accreted AM just from explicit 'swallow' operations
#else
#define SINK_FOLLOW_ACCRETED_ANGMOM 1 // follow accreted AM from 'swallowed' BH particles, and from continuous/smooth properties [mdot] of kernel gas near BH
#endif
#endif

#ifdef MAGNETIC /* recommended MHD switches -- only turn these off for de-bugging */
#define DIVBCLEANING_DEDNER         /* hyperbolic/parabolic div-cleaing (Dedner 2002), with TP improvements */
#ifdef HYDRO_SPH /* MHD switches specific to SPH MHD */
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

#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS) || defined(GALSF_ISMDUSTCHEM_MODEL) || defined(RT_INFRARED)
#define OUTPUT_DUST_TO_GAS_RATIO // helpful if these special modules are on to see this output and save it for use in analysis
#endif

#if defined(OUTPUT_POTENTIAL) && !defined(EVALPOTENTIAL)
#define EVALPOTENTIAL
#endif

#if defined(SINK_PARTICLES) && (defined(SINK_REPOSITION_ON_POTMIN) || defined(SINK_SEED_FROM_FOF))
#ifndef EVALPOTENTIAL
#define EVALPOTENTIAL
#endif
#endif

#ifdef EVALPOTENTIAL
#ifndef COMPUTE_POTENTIAL_ENERGY
#define COMPUTE_POTENTIAL_ENERGY
#endif
#endif

#ifdef COOL_MOLECFRAC_NONEQM
#ifndef OUTPUT_MOLECULAR_FRACTION
#define OUTPUT_MOLECULAR_FRACTION
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



/* block for metals and other passive scalars, should stay in this order. like the RHD, probably a more elegant way to do this with functions, but designed here to use compiler logic instead */
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


#define NUM_ADDITIONAL_PASSIVESCALAR_SPECIES_FOR_YIELDS_AND_DIFFUSION 0 /* placeholder for arbitrary number of additional species to be used for different operations like yields, etc. */


#if defined(GALSF_ISMDUSTCHEM_MODEL) /* define some global and other useful variables for dust chemistry modules which also utilize the metals info above */
#ifndef GALSF_ISMDUSTCHEM_SILICATE_COMPOSITION
#define GALSF_ISMDUSTCHEM_SILICATE_COMPOSITION (1) /* default: olivine-pyroxene mix; the module reads this as a plain integer, so it must always be defined when GALSF_ISMDUSTCHEM_MODEL is active */
#endif
#if defined(COOLING)
#define GALSF_ISMDUSTCHEM_HIGHTEMPDUSTCOOLING // optional, can turn off
#endif
#ifndef OUTPUT_MOLECULAR_FRACTION
#define OUTPUT_MOLECULAR_FRACTION // useful since we track dense molecular fraction
#endif

#define NUM_ISMDUSTCHEM_ELEMENTS (1+NUM_LIVE_SPECIES_FOR_COOLTABLES) // number of metal species evolved for dust
#define NUM_ISMDUSTCHEM_SOURCES (4) // Sources of dust creation/growth 0=gas-dust accretion, 1=SNe Ia, 2=SNe II, 3=AGB outflows
#define NUM_ISMDUSTCHEM_SPECIES (2*(GALSF_ISMDUSTCHEM_MODEL & 1) + (GALSF_ISMDUSTCHEM_MODEL & 2)/2 + (GALSF_ISMDUSTCHEM_MODEL & 4)/4 + (GALSF_ISMDUSTCHEM_MODEL & 8)/8) /* number of dust species tracked, depends on the model */
#define NUM_ISMDUSTCHEM_SPECIES_IDS (5) /* number of possible species IDs (0=silicate,1=carbon,2=free iron,3=O reservoir,4=iron inclusions); the sparse SpeciesFieldIndexTable is indexed by these fixed IDs, not the packed tracked-species count */
#if !(GALSF_ISMDUSTCHEM_MODEL & 1)
#error "GALSF_ISMDUSTCHEM_MODEL requires bit 1 (silicates+carbon are the base species always tracked); e.g. use (1), (1+2), ..."
#endif
#if (GALSF_ISMDUSTCHEM_MODEL & 8) && !(GALSF_ISMDUSTCHEM_MODEL & 2)
#error "GALSF_ISMDUSTCHEM_MODEL bit 8 (iron inclusions in silicates) requires bit 2 (metallic iron)"
#endif
#define GALSF_ISMDUSTCHEM_VAR_ELEM_IN_SILICATES 4 /* O, Mg, Si, and Fe needed to make silicates */
#undef NUM_ADDITIONAL_PASSIVESCALAR_SPECIES_FOR_YIELDS_AND_DIFFUSION
#ifdef GALSF_ISMDUSTCHEM_GRAINSIZEEVO
#ifndef OUTPUT_MACH_NUMBER // mach number is used for subresolved clumping in grain size evo model
#define OUTPUT_MACH_NUMBER
#endif
#define NUM_ISMDUSTCHEM_SIZE_BINS (GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
#define UNIT_GRAIN_NUMBER      (All.UnitGrainNumber)
#define UNIT_GRAIN_LENGTH      (All.UnitGrainLength_in_cm)
#define NUM_ADDITIONAL_PASSIVESCALAR_SPECIES_FOR_YIELDS_AND_DIFFUSION (NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES+NUM_ISMDUSTCHEM_SPECIES+(2*NUM_ISMDUSTCHEM_SPECIES*NUM_ISMDUSTCHEM_SIZE_BINS))
#else
#define NUM_ADDITIONAL_PASSIVESCALAR_SPECIES_FOR_YIELDS_AND_DIFFUSION (NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES+NUM_ISMDUSTCHEM_SPECIES)
#endif
#endif // end of GALSF_ISMDUSTCHEM_MODEL block

/* end of metals block */




/* ??? this should be done with an inline function inside the particle data, not a compiler directive */
#define SinkParticle_GravityKernelRadius (All.ForceSoftening[5])


#ifdef SINK_PARTICLES
#define SINK_COUNTPROGS /* carries a counter for each BH that gives the total number of seeds that merged into it */
#define SINK_ENFORCE_EDDINGTON_LIMIT /* put a hard limit on the maximum accretion rate (set SinkEddingtonFactor>>1 to allow super-eddington) */
#if defined(SINK_PHOTONMOMENTUM) || defined(RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION)
#define SINK_CALC_LOCAL_ANGLEWEIGHTS
#endif
#if defined(SINK_GRAVCAPTURE_GAS) || defined(SINK_GRAVACCRETION) || defined(SINK_GRAVCAPTURE_NONGAS) || defined(SINK_CALC_LOCAL_ANGLEWEIGHTS) || defined(SINK_REPOSITION_ON_POTMIN)
#define SINK_NEIGHBOR_BITFLAG 63 /* allow all particle types in the BH search: 63=2^0+2^1+2^2+2^3+2^4+2^5 */
#else
#define SINK_NEIGHBOR_BITFLAG 33 /* only search for particles of types 0 and 5 (gas and sinks) around a primary BH particle */
#endif
#endif


#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_X) || defined(BOX_OUTFLOW_Y) || defined(BOX_OUTFLOW_Z)
#define BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE 1 /* flag to let the code know to use everything below */
#define BOX_VALUE_FOR_NOTHING_SPECIAL_BOUNDARY_ 20 /* define a dummy value we won't have the user set for reference below */
#endif


#if defined(CHIMES_HII_REGIONS) && !defined(GALSF_FB_FIRE_RT_HIIHEATING)
#define GALSF_FB_FIRE_RT_HIIHEATING // must be on, this module uses the same code
#endif
#ifdef CHIMES_STELLAR_FLUXES /* The following defines the stellar age bins that we will use to define the UV spectra from stars used in CHIMES. */
#define CHIMES_LOCAL_UV_NBINS 8
#define CHIMES_LOCAL_UV_AGE_LOW 0.0
#define CHIMES_LOCAL_UV_DELTA_AGE_LOW 0.2
#define CHIMES_LOCAL_UV_AGE_MID 1.0
#define CHIMES_LOCAL_UV_DELTA_AGE_HI 1.0
#endif


#if defined(COOLING) && defined(GALSF_EFFECTIVE_EQS)
#ifndef COOLING_OPERATOR_SPLIT
#define COOLING_OPERATOR_SPLIT /*!< the Springel-Hernquist EOS depends explicitly on the cooling time in a way that requires de-coupled hydro cooling */
#endif
#endif


#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
#if !(CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)) // allow to be set to integer value to represent a >1 number of special zoom sites
#undef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
#define SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM (1)
#endif
#endif

#if defined(SINK_WIND_SPAWN) && !defined(WAKEUP_REVERSE_KICK_ON_DEMOTION) && !defined(WAKEUP_TRUNCATE_STEP_ON_DEMOTION)
#define WAKEUP_TRUNCATE_STEP_ON_DEMOTION // the WAKEUP kick reversal injects energy and only fires with spawned cells present, so truncate the step instead whenever we are spawning. must come after every SINK_WIND_SPAWN derivation above. WAKEUP_REVERSE_KICK_ON_DEMOTION opts back out
#endif

#if defined(SINK_WIND_SPAWN) && !defined(MERGE_SPLIT_DISCARD_ENERGY) && !defined(MERGE_SPLIT_CONSERVE_ENERGY)
#define MERGE_SPLIT_CONSERVE_ENERGY // merging sets the velocity momentum-conservingly, so the pair's COM-frame kinetic energy leaves the budget; discarding it costs ~23% of the injected wind energy in the adiabatic wind_singlestar test, so return it as heat whenever we are spawning. pairs with the truncation above: the two errors have opposite sign and were cancelling, which is why the discard went unnoticed. MERGE_SPLIT_DISCARD_ENERGY opts back out
#endif

/* The two above fix the ENERGY BOOKKEEPING of a merge; the two below govern WHICH pairs may merge at
   all, which is what preserves the structure a sustained outflow builds. Keyed on the jet/wind modules
   rather than on SINK_WIND_SPAWN, because the concern is a long-lived contact discontinuity between
   outflow and ambient gas: SINGLE_STAR_FB_SNE also spawns cells but is a single impulsive event with
   no such interface to maintain. With the bookkeeping fixed but these off, the adiabatic
   wind_singlestar bubble still misses Weaver's coefficient by 2.4% and the radiative one by 5.7%,
   losing ~11 points of retained energy and expanding closer to momentum- than energy-driven. */
#if (defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS)) && !defined(SINK_SPAWN_MERGE_ANY_NEIGHBOR) && !defined(SINK_SPAWN_NO_MERGE) && !defined(SINK_SPAWN_MERGE_WHEN_AMBIENT)
#define SINK_SPAWN_MERGE_WHEN_AMBIENT // retire a spawned cell only once its own state has equilibrated with its kernel, kinematically and thermally, rather than as soon as it is subsonic wrt whichever neighbour is picked as target. A contact discontinuity is in pressure equilibrium with no velocity jump, so outflow material sitting at the contact passes every kinematic test while remaining a distinct phase. SINK_SPAWN_MERGE_ANY_NEIGHBOR opts back out; SINK_SPAWN_NO_MERGE supersedes it
#endif

#if (defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS)) && !defined(MERGE_SPLIT_ALLOW_KINETIC_DISSIPATION) && !defined(MERGE_SPLIT_LIMIT_KINETIC_DISSIPATION)
#define MERGE_SPLIT_LIMIT_KINETIC_DISSIPATION // the shock-side counterpart: never merge a pair whose merge would thermalize the majority of the energy involved, which cooling then radiates away, making the merge an artificial unresolved shock. MERGE_SPLIT_ALLOW_KINETIC_DISSIPATION opts back out
#endif


