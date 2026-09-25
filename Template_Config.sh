#!/bin/bash            # this line only there to enable syntax highlighting in this file

####################################################################################################
#  Enable/Disable compile-time options as needed: this is where you determine how the code will act
#  From the list below, please activate/deactivate the
#       options that apply to your run. If you modify any of these options,
#       make sure that you recompile the whole code by typing "make clean; make".
#
#  Consult the User Guide before enabling any option. Some modules are proprietary -- access to them
#    must be granted separately by the code authors (just having the code does NOT grant permission).
#    Even public modules have citations which must be included if the module is used for published work,
#    these are all given in the User Guide.
#
# This file was originally part of the GADGET3 code developed by Volker Springel. It has been modified
#   substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO (to add new modules, change
#   naming conventions, restructure, add documention, and match GIZMO conventions)
#
####################################################################################################



####################################################################################################
# --------------------------------------- Boundary Conditions & Dimensions
####################################################################################################
#BOX_SPATIAL_DIMENSION=3    # sets number of spatial dimensions evolved (default=3). Switch for 1D/2D test problems: if =1, code only follows the x-line (all y=z=0), if =2, only xy-plane (all z=0). requires SELFGRAVITY_OFF
#BOX_PERIODIC               # Use this if periodic/finite boundaries are needed (otherwise an infinite box [no boundary] is assumed)
#BOX_BND_PARTICLES          # particles with ID=0 are forced in place (their accelerations are set =0): use for special boundary conditions where these particles represent fixed "walls"
#BOX_SHEARING=1             # shearing box boundaries: 1=r-z sheet (r,z,phi coordinates), 2=r-phi sheet (r,phi,z), 3=r-phi-z box, 4=as 3, with vertical gravity
#BOX_SHEARING_Q=(3./2.)     # shearing box q=-dlnOmega/dlnr; will default to 3/2 (Keplerian) if not set
#BOX_LONG_X=140             # modify box dimensions (non-square finite box): multiply X (not compatible with periodic gravity: if BOX_PERIODIC or PMGRID is active, make sure SELFGRAVITY_OFF or GRAVITY_NOT_PERIODIC is on)
#BOX_LONG_Y=1               # modify box dimensions (non-square finite box): multiply Y
#BOX_LONG_Z=1               # modify box dimensions (non-square finite box): multiply Z
#BOX_REFLECT_X=0            # make the x-boundary reflecting (assumes a box 0<x<BoxSize_X, where BoxSize_X=BoxSize*BOX_LONG_X, if BOX_LONG_X is set); if no value set or =0, both x-boundaries reflect, if =-1, only lower-x (x=0) boundary reflects, if =+1, only upper-x (x=BoxSize) boundary reflects
#BOX_REFLECT_Y              # make the y-boundary reflecting (assumes a box 0<y<BoxSize_Y); if no value set or =0, both y-boundaries reflect, if =-1, only lower-y (y=0) boundary reflects, if =+1, only upper-y (y=BoxSize) boundary reflects
#BOX_REFLECT_Z              # make the z-boundary reflecting (assumes a box 0<z<BoxSize_Z); if no value set or =0, both z-boundaries reflect, if =-1, only lower-z (z=0) boundary reflects, if =+1, only upper-z (z=BoxSize) boundary reflects
#BOX_OUTFLOW_X=0            # make the x-boundary outflowing (assumes a box 0<x<BoxSize_X, where BoxSize_X=BoxSize*BOX_LONG_X, if BOX_LONG_X is set); if no value set or =0, both x-boundaries outflow, if =-1, only lower-x (x=0) boundary outflows, if =+1, only upper-x (x=BoxSize) boundary outflows
#BOX_OUTFLOW_Y              # make the y-boundary outflowing (rules follow BOX_OUTFLOW_X, for the y-axis here). note that outflow boundaries are usually not needed, with Lagrangian methods, but may be useful in special cases.
#BOX_OUTFLOW_Z              # make the z-boundary outflowing (rules follow BOX_OUTFLOW_X, for the z-axis here)
####################################################################################################



####################################################################################################
# --------------------------------------- Hydro solver method
####################################################################################################
# --------------------------------------- Finite-volume Godunov methods (choose one, or SPH)
#HYDRO_MESHLESS_FINITE_MASS     # solve hydro using the mesh-free Lagrangian (fixed-mass) finite-volume Godunov method
#HYDRO_MESHLESS_FINITE_VOLUME   # solve hydro using the mesh-free (quasi-Lagrangian) finite-volume Godunov method (control mesh motion with HYDRO_FIX_MESH_MOTION)
#HYDRO_REGULAR_GRID             # solve hydro equations on a regular (recti-linear) Cartesian mesh (grid) with a finite-volume Godunov method
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- Options to explicitly control the mesh motion (for use with the MFV or grid solvers): only set for non-standard behavior
#HYDRO_FIX_MESH_MOTION=0        # mesh with arbitrarily-defined mesh-generating velocities: (0=non-moving, 1=fixed-v [set in ICs] cartesian, 2=fixed-v [ICs] cylindrical, 3=fixed-v [ICs] spherical, 4=analytic function, 5=smoothed-Lagrangian, 6=glass-generating, 7=fully-Lagrangian)
#HYDRO_GENERATE_TARGET_MESH     # use for IC generation (can be used with -any- hydro method: MFM/MFV/SPH/grid): this allows you to specify in the functions 'return_user_desired_target_density' and 'return_user_desired_target_pressure' (in eos.c) the desired initial density/pressure profile, and the code will try to evolve towards this.
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- SPH methods (enable one of these flags to use SPH):
#HYDRO_PRESSURE_SPH             # solve hydro using SPH with the 'pressure-sph' formulation ('P-SPH')
#HYDRO_DENSITY_SPH              # solve hydro using SPH with the 'density-sph' formulation (GADGET-2 & GASOLINE SPH)
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- Kernel Options
#KERNEL_FUNCTION=3              # Choose the kernel function (2=quadratic peak, 3=cubic spline [default], 4=quartic spline, 5=quintic spline, 6=Wendland C2, 7=Wendland C4, 8=2-part quadratic, 9=Wendland C6)
#KERNEL_CRK_FACES               # Use the consistent reproducing kernel [higher-order tensor corrections to kernel above, compared to our usual matrix formalism] from Frontiere, Raskin, and Owen to define the faces in MFM/MFV methods. can give more accurate closure, potentially improved accuracy in MHD problems. remains experimental for now.
####################################################################################################



####################################################################################################
# --------------------------------------- Additional Fluid Physics
####################################################################################################
## ----------------------------------------------------------------------------------------------------
# --------------------------------------- Gas (or Material) Equations-of-State [some EOS options for specific regimes, like galaxy or star formation simulations, are also described in the blocks below for those sections]
#EOS_GAMMA=(5.0/3.0)            # Polytropic Index of Gas (for an ideal gas law): if not set and no other (more complex) EOS set, defaults to GAMMA=5/3
#EOS_HELMHOLTZ                  # Use Timmes & Swesty 2000 EOS (for e.g. stellar or degenerate equations of state); if additional tables needed, download at http://www.tapir.caltech.edu/~phopkins/public/helm_table.dat (or the GitHub site)
#EOS_TILLOTSON                  # Use Tillotson (1962) EOS (for solid/liquid+vapor bodies, impacts); custom EOS params can be specified or pre-computed materials used. see User Guide and Deng et al., arXiv:1711.04589
#EOS_ANEOS                      # Use ANEOS/SESAME tabulated EOS (for planetary impacts/collisions); specify table file paths in parameterfile (AneosTable0, AneosTable1, etc.)
#EOS_ELASTIC                    # treat fluid as elastic or plastic (or visco-elastic) material, obeying Hooke's law with full stress terms and von Mises yield model. custom EOS params can be specified or pre-computed materials used.
## -----------------------------------------------------------------------------------------------------
# --------------------------------- Magneto-Hydrodynamics
# ---------------------------------  these modules are public, but if used, the user should also cite the MHD-specific GIZMO methods paper
# ---------------------------------  (Hopkins 2015: 'Accurate, Meshless Methods for Magneto-Hydrodynamics') as well as the standard GIZMO paper
#MAGNETIC                       # top-level switch for MHD, regardless of which Hydro solver is used
#MHD_B_SET_IN_PARAMS            # set initial fields (Bx,By,Bz) in parameter file
#MHD_NON_IDEAL                  # enable non-ideal MHD terms: Ohmic resistivity, Hall effect, and ambipolar diffusion (solved explicitly); Users should cite Hopkins 2017, MNRAS, 466, 3387, in addition to the MHD paper
#MHD_CONSTRAINED_GRADIENT=1     # use CG method (in addition to cleaning, optional!) to maintain low divB: set this value to control how aggressive the div-reduction is:
                                # 0=minimal (safest), 1=intermediate (recommended), 2=aggressive (less stable), 3+=very aggressive (less stable+more expensive). [Please cite Hopkins, MNRAS, 2016, 462, 576]
#MHD_MODIFIED_GRADIENT          # use MG method (Tu, Wang, Gao & Tang 2026) to correct B-field gradients via global sparse-matrix solve for exact divB=0 (to machine precision). applied after gradient calculation + slope limiting, before hydro force loop. mutually exclusive with MHD_CONSTRAINED_GRADIENT. can also enable MHD_MODIFIED_GRADIENT_USE_PARDISO to use MKL PARDISO direct sparse solver instead of Hypre AMG for the MG matrix solve. requires Intel MKL. gathers matrix to rank 0 (best for small-medium problems).
#MHD_NON_IDEAL_CORRECTIONTERMS  # enable approximate corrections for anomalous resistivity and Epstein-like (drift/slip-dependent) cross sections in non-ideal MHD coefficients. Please cite Hopkins et al., https://arxiv.org/abs/2405.06026, where the scalings here are derived and presented
## ----------------------------------------------------------------------------------------------------
# -------------------------------------- Conduction
# ----------------------------------------- [Please cite and read the methods paper Hopkins 2017, MNRAS, 466, 3387]
#CONDUCTION                     # Thermal conduction solved *explicitly*: isotropic if MAGNETIC off, otherwise anisotropic
#CONDUCTION_SPITZER             # Spitzer conductivity accounting for saturation: otherwise conduction coefficient is constant  [cite Su et al., 2017, MNRAS, 471, 144, in addition to the conduction methods paper above].  Requires COOLING to calculate local thermal state of gas.
## ----------------------------------------------------------------------------------------------------
# -------------------------------------- Viscosity
# ----------------------------------------- [Please cite and read the methods paper Hopkins 2017, MNRAS, 466, 3387]
#VISCOSITY                      # Navier-stokes equations solved *explicitly*: isotropic coefficients if MAGNETIC off, otherwise anisotropic
#VISCOSITY_BRAGINSKII           # Braginskii viscosity tensor for ideal MHD [cite Su et al., 2017, MNRAS, 471, 144, in addition to the viscosity methods paper above]. Requires COOLING to calculate local thermal state of gas.
## ----------------------------------------------------------------------------------------------------
# -------------------------------------- Smagorinsky Turbulent Eddy Diffusion Model
# --------------------------------------- Users of these modules should cite Hopkins et al. 2017 (arXiv:1702.06148) and Colbrook et al. (arXiv:1610.06590)
#TURB_DIFF_METALS               # turbulent diffusion of metals (passive scalars); requires METALS
#TURB_DIFF_ENERGY               # turbulent diffusion of internal energy (conduction with effective turbulent coefficients)
#TURB_DIFF_VELOCITY             # turbulent diffusion of momentum (viscosity with effective turbulent coefficients)
#TURB_DIFF_DYNAMIC              # replace Smagorinsky-style eddy diffusion with the 'dynamic localized Smagorinsky' model from Rennehan et al. (arXiv:1807.11509 and 2104.07673): cite those papers for all methods. more accurate but more complex and expensive.
## ----------------------------------------------------------------------------------------------------
# --------------------------------------- Aerodynamic Particles
# ----------------------------- This is developed by P. Hopkins, who requests that you inform him of planned projects with these modules
# ------------------------------  because he is supervising several students using them as well, and there are some components still in active development.
# ------------------------------  Users should cite: Hopkins & Lee 2016, MNRAS, 456, 4174, and Lee, Hopkins, & Squire 2017, MNRAS, 469, 3532, for the numerical methods (plus other papers cited or listed below, for each of the appropriate modules as described below or in the User Guide)
#GRAIN_FLUID                    # aerodynamically-coupled grains (particle type 3 are grains); default is Epstein drag. Cite papers above.
#GRAIN_EPSTEIN_STOKES=1         # uses the cross section for molecular hydrogen (times this number) to calculate Epstein-Stokes drag; need to set GrainType=1 (will use calculate which applies and use appropriate value); if used with GRAIN_LORENTZFORCE and GrainType=2, will also compute Coulomb drag. Cite Hopkins et al., 2020, MNRAS, 496, 2123
#GRAIN_BACKREACTION             # account for momentum of grains pushing back on gas (from drag terms); users should cite Moseley et al., 2018, arXiv:1810.08214.
#GRAIN_LORENTZFORCE             # charged grains feel Lorentz forces (requires MAGNETIC); if used with GRAIN_EPSTEIN_STOKES flag, will also compute Coulomb drag (grain charges self-consistently computed from gas properties). Need to set GrainType=2. Please cite Seligman et al., 2019, MNRAS 485 3991
#GRAIN_COLLISIONS               # model collisions between grains (super-particles; so this is stochastic). Default = hard-sphere scattering, with options for inelastic or velocity-dependent terms. Approved users please cite papers above and Rocha et al., MNRAS 2013, 430, 81
####################################################################################################



####################################################################################################
# ------------------------------------- Driven turbulence (for turbulence tests, large-eddy sims)
# ------------------------------- users of these routines should cite Bauer & Springel 2012, MNRAS, 423, 3102. Thanks to A. Bauer for providing the core algorithms
####################################################################################################
#TURB_DRIVING                   # turns on turbulent driving/stirring. see begrun for parameters that must be set
#TURB_DRIVING_SPECTRUMGRID=128  # activates on-the-fly calculation of the turbulent velocity, vorticity, and smoothed-velocity power spectra, evaluated on a grid of linear-size TURB_DRIVING_SPECTRUMGRID elements. Requires BOX_PERIODIC
#TURB_DRIVING_TARGET_SIGMA=0    # enables proportional control loop adjusting driving amplitude to converge measured gas RMS velocity to TurbDrive_ApproxRMSVturb, smoothed over 2 coherence times. If value > 0, only gas cells within that distance (code units) of box center contribute to the measured velocity dispersion; value=0 uses all gas cells.
####################################################################################################



####################################################################################################
## ------------------------ Gravity & Cosmological Integration Options ---------------------------------
####################################################################################################
# --------------------------------------- TreePM Options (recommended for cosmological sims)
#PMGRID=512                     # adds Particle-Mesh grid for faster (but less accurate) long-range gravitational forces: value sets resolution (e.g. a PMGRID^3 grid will overlay the box, as the 'top level' grid)
#PM_PLACEHIGHRESREGION=1+2+16   # adds a second-level (nested) PM grid before the tree: value denotes particle types (via bit-mask) to place high-res PMGRID around. Requires PMGRID.
## -----------------------------------------------------------------------------------------------------
# ---------------------------------------- Adaptive Grav. Softening (including Lagrangian conservation terms!)
#ADAPTIVE_GRAVSOFT_FORGAS       # allows variable softening length for gas particles (scaled with local inter-element separation), so gravity traces same density field seen by hydro
#ADAPTIVE_GRAVSOFT_FORALL=1+2   # enable adaptive gravitational softening lengths for designated particle types (ADAPTIVE_GRAVSOFT_FORGAS should be disabled). the softening is set to the distance
                                # enclosing a neighbor number set in the parameter file. flag value = bitflag like PM_PLACEHIGHRESREGION, which determines which particle types are adaptive (others use fixed softening). cite Hopkins et al., arXiv:1702.06148
#ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION=2 # enable adaptive gravitational softening lengths for designated particle types (value = bitflag like ADAPTIVE_GRAVSOFT_FORALL), but setting the softening to the tidal-tensor based value instead of nearest-neighbor values, following Hopkins et al. (arXiv:2212.06851). this is cross-compatible with ADAPTIVE_GRAVSOFT_FORGAS, or even ADAPTIVE_GRAVSOFT_FORALL so long as the same particle type does not appear in both.
## -----------------------------------------------------------------------------------------------------
#SELFGRAVITY_OFF                # turn off self-gravity (compatible with GRAVITY_ANALYTIC)
#GRAVITY_NOT_PERIODIC           # self-gravity is not periodic, even though the rest of the box is periodic
## -----------------------------------------------------------------------------------------------------
#GRAVITY_ANALYTIC               # specific analytic gravitational force to use instead of or with self-gravity. If set to a numerical value
                                #  > 0 (e.g. =1), then SINK_CALC_DISTANCES will be enabled, and it will use the nearest BH particle as the center for analytic gravity computations
                                #  (edit "gravity/analytic_gravity.h" to actually assign the analytic gravitational forces). 'ANALYTIC_GRAVITY' gives same functionality
## ----------------------------------------------------------------------------------------------------
# -------------------------------------- Self-Interacting DM (Rocha et al. 2012) and Scalar-field DM and Fuzzy DM; users please cite Rocha et al., MNRAS 2013, 430, 81 and Robles et al, 2017 (arXiv:1706.07514)
#DM_SIDM=2                      # self-interacting particle types (specify the particle types which are self-interacting DM with a bit mask, as for PM_PLACEHIGHRESREGION above (see description); previous "DMDISK_INTERACTIONS" is identical to setting DM_SIDM=2+4  [cite Rocha et al., MNRAS 2013, 430, 81 and Robles et al, 2017 (arXiv:1706.07514)]
#DM_SCALARFIELD_SCREENING       # gravity is mediated by a long-range scalar field, with dynamical screening (primarily alternative DE models) [cite Rocha et al., MNRAS 2013, 430, 81 and Robles et al, 2017 (arXiv:1706.07514)]
#DM_FUZZY=0                     # DM particles (Type=1) are described by Bose-Einstein Condensate: within gravity kernel (adaptive), solves quantum pressure tensor for non-linear terms arising from Schroedinger equation for a given particle mass. Still some testing recommended, while this is public we encourage discussion with PFH over different use cases and applications of the method. The value here is: 0=fully-conservative Madelung method. 1=mass-conserving direct SPE integration. 2=direct SPE (non mass-conserving). Cite Hopkins et al., 2019MNRAS.489.2367H if used.
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- Pure-Tree Options for Direct N-body of small-N groups (recommended for hard binaries, etc)
#GRAVITY_ACCURATE_FEWBODY_INTEGRATION # enables a suite: GRAVITY_HYBRID_OPENING_CRIT, TIDAL_TIMESTEP_CRITERION, to more accurately follow few-body point-like dynamics in the tree. currently compatible only with pure-tree gravity.
## ----------------------------------------------------------------------------------------------------
# -------------------------------------- arbitrary time-dependent dark energy equations-of-state, expansion histories, or gravitational constants
#GR_TABULATED_COSMOLOGY         # enable reading tabulated cosmological/gravitational parameters (top-level switch)
#GR_TABULATED_COSMOLOGY_W       # read pre-tabulated dark energy equation-of-state w(z)
#GR_TABULATED_COSMOLOGY_H       # read pre-tabulated hubble function (expansion history) H(z)
#GR_TABULATED_COSMOLOGY_G       # read pre-tabulated gravitational constant G(z) [also rescales H(z) appropriately]
## ----------------------------------------------------------------------------------------------------
#EOS_TRUELOVE_PRESSURE          # adds artificial pressure floor force Jeans length above resolution scale (means you can get the wrong answer, but things will look smooth).  cite Robertson & Kravtsov 2008, ApJ, 680, 1083
####################################################################################################



####################################################################################################
# --------------------------------------- On the fly FOF groupfinder
# ----------------- This is originally developed as part of GADGET-3 by V. Springel
# ----------------- Users of any of these modules should cite Springel et al., MNRAS, 2001, 328, 726 for the numerical methods.
####################################################################################################
## ----------------------------------------------------------------------------------------------------
# ------------------------------------- Friends-of-friends on-the-fly finder options (source in fof.c)
# -----------------------------------------------------------------------------------------------------
#FOF                                # top-level switch: enable FoF searching on-the-fly and outputs (set parameter LINKLENGTH=x to control LinkingLength; default=0.2)
#FOF_PRIMARY_LINK_TYPES=2           # bitflag: sum of 2^type for the primary type used to define initial FOF groups (use a common type to ensure 'start' in reasonable locations)
#FOF_SECONDARY_LINK_TYPES=1+16+32   # bitflag: sum of 2^type for the seconary types which can be linked to nearest primaries (will be 'seen' when calculating group properties)
#FOF_DENSITY_SPLIT_TYPES=1+2+16+32  # bitflag: sum of 2^type for which the densities should be calculated seperately (i.e. if 1+2+16+32, fof densities are separately calculated for types 0,1,4,5, and shared for types 2,3)
#FOF_GROUP_MIN_SIZE=32              # minimum number of identified members required to qualify as a 'group': default is 32
## ----------------------------------------------------------------------------------------------------
# -------------------------------------  Subhalo on-the-fly finder options (uses "subfind" source code).
## ----------------------------------------------------------------------------------------------------
#SUBFIND                            # top-level switch to enable substructure-finding with the SubFind algorithm
#SUBFIND_ADDIO_NUMOVERDEN=1         # for M200,R200-type properties, compute values within in this number of different overdensities (default=1=)
#SUBFIND_ADDIO_VELDISP              # add the mass-weighted 1D velocity dispersions to properties computed in parent group[s], within the chosen overdensities
#SUBFIND_ADDIO_BARYONS              # add gas mass, mass-weighted temperature, and x-ray luminosity (assuming ionized primoridal gas), and stellar masses, to properties computed in parent group[s], within the chosen overdensities
## ----------------------------------------------------------------------------------------------------
#SUBFIND_REMOVE_GAS_STRUCTURES      # delete (do not save) any structures which are entirely gas (or have fewer than target number of elements which are non-gas, with the rest in gas)
#SUBFIND_SAVE_PARTICLEDATA          # save all particle positions,velocity,type,mass in subhalo file (in addition to IDs: this is highly redundant with snapshots, so makes subhalo info more like a snapshot)
####################################################################################################



####################################################################################################
# ----------------- Galaxy formation & Galactic Star formation
####################################################################################################
## ---------------------------------------------------------------------------------------------------
#GALSF                           # top-level switch for galactic star formation model: enables SF, stellar ages, generations, etc. [cite Springel+Hernquist 2003, MNRAS, 339, 289]
## ----------------------------------------------------------------------------------------------------
# --- star formation law/particle spawning (additional options: otherwise all star particles will reflect IMF-averaged populations and form strictly based on a density criterion) ---- #
## ----------------------------------------------------------------------------------------------------
#GALSF_SFR_CRITERION=(0+1+2)     # mix-and-match SF criteria with a bitflag: 0=density threshold, 1=virial criterion, 2=convergent flow, 4=local extremum, 8=no sink in kernel, 16=not falling into sink, 32=hill (tidal) criterion, 64=Jeans criterion, 128=converging flow along all principle axes, 256=self-shielding/molecular, 512=multi-free-fall (smooth dependence on virial), 1024=adds a 'catch' which weakens some kinematic criteria when forces become strongly non-Newtonian (when approach minimum force-softening), 2048=uses time-averaged virial criterion
#GALSF_SFR_VIRIAL_SCALING=2      # instead of a threshold, implements a semi-continuous SF efficiency as a function of alpha_vir. set 0=step function above 1 to 0; 1=Padoan 2012 prescription; 2=multi-free-fall model, as in e.g. Federrath+Klessen 2012/2013 ApJ 761,156; 763,51 (similar to that implemented in e.g. Kretschmer+Teyssier 2020), based on the analytic models in Hopkins MNRAS 2013, 430 1653, with correct virial parameter
#GALSF_SFR_IMF_VARIATION         # determines the stellar IMF for each particle from the Guszejnov/Hopkins/Hennebelle/Chabrier/Padoan theory. Cite Guszejnov, Hopkins, & Ma 2017, MNRAS, 472, 2107
#GALSF_SFR_IMF_SAMPLING          # discretely sample the IMF: simplified model with quantized number of massive stars. Cite Kung-Yi Su, Hopkins, et al., Hayward, et al., 2017, "Discrete Effects in Stellar Feedback: Individual Supernovae, Hypernovae, and IMF Sampling in Dwarf Galaxies". 
#GALSF_GENERATIONS=1             # the number of star particles a gas particle may spawn (defaults to 1, set otherwise if desired)
## ----------------------------------------------------------------------------------------------------------------------------
# ---- sub-grid models (for large-volume simulations or modest/low resolution galaxy simulations) -----------------------------
# -------- the SUBGRID_WINDS models are variations of the Springel & Hernquist 2005 sub-grid models for the ISM, star formation, and winds.
# -------- Volker has granted permissions for their use, provided users properly cite the sources for the relevant models and scalings (described below)
#GALSF_EFFECTIVE_EQS            # Springel-Hernquist 'effective equation of state' model for the ISM and star formation [cite Springel & Hernquist, MNRAS, 2003, 339, 289]
#GALSF_SUBGRID_WINDS            # sub-grid winds ('kicks' as in Oppenheimer+Dave,Springel+Hernquist,Boothe+Schaye,etc): enable this top-level switch for basic functionality [cite Springel & Hernquist, MNRAS, 2003, 339, 289]
#GALSF_SUBGRID_WIND_SCALING=0   # set wind velocity scaling: 0 (default)=constant v [and mass-loading]; 1=velocity scales with halo mass (cite Oppenheimer & Dave, 2006, MNRAS, 373, 1265), requires FOF modules; 2=scale with local DM dispersion as Vogelsberger 13 (cite Zhu & Li, ApJ, 2016, 831, 52)
#GALSF_WINDS_ORIENTATION=0      # directs wind orientation [0=isotropic/random, 1=polar, 2=along density gradient]
#GALSF_FB_TURNOFF_COOLING       # turn off cooling for SNe-heated particles (as Stinson+ 2006 GASOLINE model, cite it); requires GALSF_FB_THERMAL
## ----------------------------------------------------------------------------------------------------------------------------
# ---- explicit thermal/kinetic stellar models: i.e. models which track individual 'events' (SNe, stellar mass loss, etc) and inject energy/mass/metals/momentum directly from star particles into neighboring gas
# -------- these modules explicitly evolve individual stars+stellar populations. Event rates (SNe rates, mass-loss rates) and associated yields, etc, are all specified in 'stellar_evolution.c'. the code will then handle the actual injection and events.
# -------- users are encouraged to explore their own stellar evolution models and include various types of feedback (e.g. SNe, stellar mass-loss, NS mergers, etc)
#GALSF_FB_MECHANICAL            # explicit algorithm including thermal+kinetic/momentum terms from Hopkins+ 2018 (MNRAS, 477, 1578): manifestly conservative+isotropic, and accounts properly for un-resolved PdV work+cooling during blastwave expansion. cite Hopkins et al. 2018, MNRAS, 477, 1578, and Hopkins+ 2014 (MNRAS 445, 581)
#GALSF_FB_THERMAL               # simple 'pure thermal energy dump' feedback: mass, metals, and thermal energy are injected locally in simple kernel-weighted fashion around young stars. tends to severely over-cool owing to lack of mechanical/kinetic treatment at finite resolution (better algorithm is mechanical)
#GALSF_FB_FIRE_AGE_TRACERS=16   # model for arbitrary tracers of different age-bins of stellar yields (number here = number of log-spaced bins), which can be re-convolved in post-processing. developed by A. Emerick, paper in prep by A. Wetzel, meantime cite arXiv:2203.00040
## ----------------------------------------------------------------------------------------------------
# ----- FIRE simulation modules for mechanical+radiative FB with full evolution+yield tracks (Hopkins et al. 2014, Hopkins et al., 2017a, arXiv:1702.06148 and 2203.00040) ------ ##
# -------- Use of these modules as part of the public code is now allowed with appropriate citations to the specific methods papers above. These should be referred to as using the methods "from the FIRE public code (as in citations)", not as FIRE collaboration papers or FIRE simulations. FIRE simulations/papers follow FIRE collaboration guidelines, and you should reach out to members of the collaboration if you wish to write FIRE papers or access still in-development (non-public) FIRE codes/outputs/etc.
#FIRE_PHYSICS_DEFAULTS=3        # enable standard set of FIRE physics packages (number=fire version). note use policy above. convenience flags for MHD (FIRE_MHD), BHs (FIRE_BHS), and CRs (FIRE_CRS=X, where X=-2 is toy-model sub-grid LEBRON-approximation CRs; -1=single-bin CRs (defaults to constant scattering rate); 0=full-spectrum p+e, constant (power-law) scattering rate; 1=full-spectrum p+e, variable scattering rate 2509.07104; 2=full-spectrum 10-species treatment)
#FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT # super-lagrangian refinement based on jeans mass or other criteria. this is a generic flag to be used for high-resolution massive-galaxy simulations using hyper-refinement to achieve 'standard' high-resolution FIRE quality in galaxies, without trillion-particle loads
#GALSF_FB_FIRE_RT_HIIHEATING    # enable the generic stochastic HII-region model when HII regions are not extremely well-resolved, ionize and heat to local equilibrium temperature in immediate vicinity of O-stars probabilistically so they can be semi-resolved in a statistical sense. cite FIRE methods papers. do not use if enabling FIRE defaults flags OR any explicit RHD modules of your own (FIRE modules will work with any other RHD modules/bands where they are enabled)
#GALSF_FB_FIRE_RT_LONGRANGE     # enables LEBRON radiative transport with continuous acceleration from starlight (uses luminosity tree) to propagate approximate RT, also enables simple uv/euv/nuv/optical/nir/fir band set, for photo-ionization, photo-electric, and dust radiation heating. do not use if enabling FIRE defaults flags OR any explicit RHD modules of your own (FIRE modules will work with any other RHD modules/bands where they are enabled)
#GALSF_FB_FIRE_RT_LOCALRP       # enables local radiation pressure coupling to gas - account for local multiple-scattering and isotropic local absorption, assuming unresolved local-source approximation, designed to be used when those are not treated explicitly (i.e. if and only if an RHD module like LEBRON is used, as opposed to M1 or some other RHD module. do not use if enabling FIRE defaults flags OR any explicit RHD modules of your own (FIRE modules will work with any other RHD modules/bands where they are enabled)
############################################################################################################################



############################################################################################################################
## ----------------------------------------------------------------------------------------------------
# --------------- Star+Planet+Compact Object Formation (Sink Particle + Explicit/Keplerian N-Body Dynamics)
# -------------------- (unlike GALSF options, these sinks are individual accretors, not populations). Much in common with sink particle modules below.
# -------------------- Most of the 'core' modules here are now public. Some specific stellar evolution tracks and modifications to the public modules made for the STARFORGE project remain in development, and permissions from the authors (Mike Grudic) is required for their use if they are only in the development code and -not- in the public code: please contact Mike Grudic or Stella Offner or Claude-Andre Faucher-Giguere or Phil Hopkins if you wish to run simulations as part of the STARFORGE project/suite/collaboration
## ----------------------------------------------------------------------------------------------------
#SINGLE_STAR_SINK_DYNAMICS      # top-level switch to enable any other modules in this section
## ----------------------------------------------------------------------------------------------------
# ----- time integration, regularization, and explicit small-N-body dynamical treatments (for e.g. hard binaries, etc)
## ----------------------------------------------------------------------------------------------------
#SINGLE_STAR_TIMESTEPPING=1     # use additional timestep criteria to ensure resolved binaries/multiples dont dissolve in close encounters. 0=most conservative. 1=super-timestep hard binaries by operator-splitting the binary orbit. 2=more aggressive super-timestep. cite Grudic et al., arXiv:2010.11254, for the methods here.
#HERMITE_INTEGRATION=32         # Instead of the usual 2nd order DKD Leapfrog timestep, do 4th order Hermite integration for particles matching the bitflag. Allows longer timesteps and higher accuracy collisional dynamics. cite Grudic et al., arXiv:2010.11254, for the methods here.
#HERMITE_CORRECTOR_ITERATIONS=2 # P(EC)^n Hermite: after the corrector, re-evaluate gravity at the corrected state and re-apply the corrector, n passes in all (default 1 = standard PEC). n=2 removes the secular energy drift of a fixed-step orbit (~1000x on circular pairs) for one extra Hermite force pass per step, but on its own makes eccentric orbits on adaptive block steps WORSE (test/binary |dE/E| 8.9e-5 -> 1.9e-3). Off by default: its cost in production runs is unmeasured
## ----------------------------------------------------------------------------------------------------
# ----- sink creation and accretion/growth/merger modules
## ----------------------------------------------------------------------------------------------------
#SINGLE_STAR_SINK_FORMATION=(0+1+2+4+8+16+32+64) # form new sinks on the fly, criteria from bitflag: 0=density threshold, 1=virial criterion, 2=convergent flow, 4=local extremum, 8=no sink in kernel, 16=not falling into sink, 32=hill (tidal) criterion, 64=Jeans criterion, 128=converging flow along all principle axes, 256=self-shielding/molecular, 512=multi-free-fall (smooth dependence on virial). cite Grudic et al., arXiv:2010.11254, for the methods here.
#SINGLE_STAR_ACCRETION=7        # sink accretion [details in BH info below]: 0-10: use SINK_GRAVACCRETION=X, 11: SINK_GRAVCAPTURE_GAS, 12: SINK_GRAVCAPTURE_GAS modified with Bate-style FIXEDSINKRADIUS.  cite Grudic et al., arXiv:2010.11254, for the methods here.
## ----------------------------------------------------------------------------------------------------
# ----- star (+planet) formation-specific modules (feedback, jets, radiation, protostellar evolution, etc)
## ----------------------------------------------------------------------------------------------------
#SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION=1 # sinks are assumed to be proto-stars and follow protostellar evolution tracks as they accrete to evolve radii+luminosities, determines proto-stellar feedback properties. 1=simple model [PFH], 2=fancy model [DG+MG]. Cite Grudic et al., arXiv:2010.11254, for the methods here.
#SINGLE_STAR_FB_JETS            # kinematic jets from sinks: outflow rate+velocity set by Sink_accreted_fraction+Sink_outflow_velocity. for now cite Angles-Alcazar et al., 2017, MNRAS, 464, 2840 (for algorithm, developed for sink particle jets), though now using SPAWN algorithm developed by KY Su. cite Su et al, arXiv:2102.02206, and Grudic et al., arXiv:2010.11254, for the methods here.
#SINGLE_STAR_FB_WINDS=0+1+2     # enable continuous main-sequence mechanical feedback from single stellar sources accounting for OB/AGB/WR winds. with STARFORGE parent flag[s] enabled, this will following STARFORGE methods (Grudic+ arXiv:2010.11254). Otherwise, this will follow a simpler Castor, Abbot, & Klein scaling, for type=4 particles representing single stars, using the standard GALSF_FB_MECHNICAL algorithms in code, for which you should cite Hopkins et al. 2018MNRAS.477.1578H. With STARFORGE model enabled, the numerical value encodes a bitflag toggling use of the Vink 2001 mass loss prescription (&1) and Eddington-limit floor by Sabhahit et al. arXiv:2205.09125 (&2)
#SINGLE_STAR_WIND_MODE=2        # override the automatic wind mode selection (1=spawn particles, 2=mechanical FB injection). if not set, the code chooses based on the expected number of spawned particles per wind crossing time
#SINGLE_STAR_FB_SNE             # enable supernovae from single stellar sources at end of main-sequence lifetime. with STARFORGE flag[s] enabled this will use particle spawning in shells following STARFORGE methods (Grudic+ arXiv:2010.11254). Otherwise, this will act uniformly and in a single timestep at the end of the stellar main-sequence lifetime for type=4 particles representing single stars, using the standard GALSF_FB_MECHNICAL algorithms in code, for which you should cite Hopkins et al. 2018MNRAS.477.1578H
#SINGLE_STAR_FB_RAD             # enable radiative feedback from stars, hooking into the standard radiation hydrodynamics algorithms. you need to determine how the effective temperatures of the stars scale, which will be used to determine their input fluxes into the different explicitly-evolved bands.
#EOS_SUBSTELLAR_ISM             # allows for the local equation of state polytropic index to vary between 7/5 and 5/3 and outside this range following the detailed fit from Vaidya et al. A&A 580, A110 (2015) for n_H ~ 10^7, which accounts for collisional dissociation at 2000K and ionization at 10^4K, and take the fmol-weighted average with 5./3 at the end to interpolate between atomic/not self-shielding and molecular/self-shielding. Gamma should technically really come from calculating the species number-weighted specific heats, but fmol is very approximate so this should be OK. See the code and uncomment the noted lines in EOS.c if you want to use the exact version from Vaidya+15, which rolls the heat of ionization into the EOS. cite Grudic+ arXiv:2010.11254
#EOS_GMC_BAROTROPIC             # Barotropic EOS calibratied to Masunaga & Inutsuka 2000; useful for test problems in small-scale star formation such as cloud collapse, jet launching. See Federrath et al. 2014ApJ...790..128F. Can also set to a numerical value =1 to instead use EOS used in Bate Bonnell & Bromm 2003
## ----------------------------------------------------------------------------------------------------
# ----- optional and de-bugging modules (intended for specific behaviors)
## ----------------------------------------------------------------------------------------------------
#SINK_RETURN_ANGMOM_TO_GAS      # BH/sink particles return accreted angular momentum to surrounding gas (following Hubber+13) to represent AM transfer (loss in accreting material)
#SINGLE_STAR_FIND_BINARIES      # manually enable identification of close binaries (normally enabled automatically if actually used for e.g. hermite timestepping). cite Grudic et al. arXiv:2010.11254
#SINGLE_STAR_FB_LOCAL_RP        # approximate local radiation pressure from single-star sources, using the same LEBRON-type approximation as in FIRE - useage follows the FIRE collaboration policies
#SINGLE_STAR_FB_RT_HEATING      # proto-stellar heating: luminosity determined by SinkRadiativeEfficiency (typical ~5e-7). This particular module used without radiation-hydrodynamics uses FIRE modules, so permissions follow those. But by enabling explicit radiation-hydrodynamics, this is not needed, and the user can treat full radiative feedback in the public code.
#SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT=2  # determines the maximum number of ejecta particles spawned per timestep in the supernova shell approximation for spawning. only needs to be modified for testing purposes.
############################################################################################################################



####################################################################################################
# ---------------- sink particles (Sink-Particles with Accretion and Feedback: Black Holes, Stars, Planets, etc.)
####################################################################################################
#SINK_PARTICLES                   # top-level switch to enable sink particle modules
## ----------------------------------------------------------------------------------------------------
# ----- seeding / BH-particle spawning
## ----------------------------------------------------------------------------------------------------
#SINK_SEED_FROM_FOF=0             # use FOF on-the-fly to seed BHs in massive FOF halos; =0 uses DM groups, =1 uses stellar [type=4] groups; requires FOF with linking type including relevant particles (cite Angles-Alcazar et al., MNRAS, 2017, arXiv:1707.03832)
#SINK_SEED_FROM_LOCALGAS          # BHs seeded on-the-fly from dense, low-metallicity gas (no FOF), like star formation; criteria define-able in modular fashion in sfr_eff.c (function return_probability_of_this_forming_sink_from_seed_model). cite Grudic et al. (arXiv:1612.05635) and Lamberts et al. (MNRAS, 2016, 463, L31). Requires GALSF and METALS
#SINK_INCREASE_DYNAMIC_MASS=100   # increase the particle dynamical mass by this factor at the time of BH seeding
## ----------------------------------------------------------------------------------------------------
# ----- dynamics (when BH mass is not >> other particle masses, it will artificially get kicked and not sink via dynamical friction; these 're-anchor' the BH for low-res sims)
## ----------------------------------------------------------------------------------------------------
#SINK_DYNFRICTION_FROMTREE        # compute dynamical friction forces on BH following the discrete DF estimator in Linhao Ma et al., arXiv:2101.02727 and arXiv:2208.12275. This is a more flexible, general, and less noisy and more accurate version of the traditional Chandrasekhar dynamical friction formula. Cite L Ma et al. 2021 and 2022 if used, and contact author L. Ma for applications as testing is still ongoing.
#SINK_REPOSITION_ON_POTMIN        # reposition sink particle on potential minimum. moves smoothly with damped velocity to most-bound particle in kernel. Cite Wellons et al. arXiv:2203.06201
## ----------------------------------------------------------------------------------------------------
# ----- accretion models (modules for gas or other particle accretion)
## ----------------------------------------------------------------------------------------------------
#SINK_SWALLOWGAS                  # 'top-level switch' for accretion (should always be enabled if accretion is on). enables BH to actually eliminate gas particles and take their mass.
#SINK_ALPHADISK_ACCRETION=(10)    # gas accreted goes into a 'virtual' alpha-disk (mass reservoir), which then accretes onto the BH at the viscous rate (determining luminosity, etc). cite GIZMO methods. should be set to a value, which limits the maximum mass of the reservoir to that multiple of the central sink mass
#SINK_MDOT_RESERVOIR_CAP=30       # with SINK_ALPHADISK_ACCRETION: while a sink holds fewer than SINK_MDOT_RESERVOIR_CAP_NRESOLVED (default 100) cells' worth of mass, drain at most M_sink/N of the reservoir per t_acc (floored at one gas cell). ON with N=30 under SINGLE_STAR_STARFORGE_DEFAULTS; set only to change N
#SINK_SUBGRIDBHVARIABILITY        # model variability below resolved dynamical time for BH (convolve accretion rate with a uniform power spectrum of fluctuations on timescales below the minimum resolved dynamical time). cite Hopkins & Quataert 2011, MNRAS, 415, 1027. Requires GALSF.
#SINK_GRAVCAPTURE_NONGAS          # accretion determined only by resolved gravitational capture by the BH, for non-gas particles (can be enabled with other accretion models for gas). cite Hopkins et al., 2016, MNRAS, 458, 816
## ----
#SINK_GRAVCAPTURE_GAS             # accretion determined only by resolved gravitational capture by the BH (for gas particles). cite Hopkins et al., 2016, MNRAS, 458, 816
#SINK_GRAVACCRETION=1             # family of gravitational/torque/angular-momentum-driven accretion models from Hopkins & Quataert (2011): cite Hopkins & Quataert 2011, MNRAS, 415, 1027 and Angles-Alcazar et al. 2017, MNRAS, 464, 2840. see `notes_blackholes` for details:
#                                 # [=0] evaluate at density kernel radius, [=1] evaluate at fixed physical radius, [=2] fixed efficiency per FF time at physical radius, [=3] gravito-turbulent scaling, [=4] fixed per FF at BH radius of influence, [=5] hybrid scaling (switch to Bondi if circularization radius small),
#                                 # [=6] modified bondi-hoyle/fixed accretion in sonic point for rho~r^-1 profile, [=7] shu+pressure+turbulence solution for isothermal sphere (self-similar isothermal sphere solution with these terms), [=8] hubber+13 estimator of local inflow (limited by 'external alpha-disk' and 'internal bondi' estimates)
#                                 # [=9] pure Bondi-Hoyle (ignore cooling/angular momentum), [=10] pure Bondi (dont use gas velocity with sound speed), [=11] variable-alpha tweak (Booth & Schaye 2009; requires GALSF)
#SINK_GRAVACCRETION_STELLARFBCORR # account for additional acceleration-dependent retention from stellar FB in Mdot. cite Hopkins et al., arXiv:2103.10444, for both the analytic derivation of these scalings and the numerical methods/implementation.
## ----------------------------------------------------------------------------------------------------
# ----- feedback models/options
## ----------------------------------------------------------------------------------------------------
#SINK_FB_COLLIMATED               # BH feedback is narrowly collimated along the axis defined by the angular momentum accreted thus far in the simulation. Cite Su et al. arXiv:2102.02206
#SINK_THERMALFEEDBACK             # thermal (pure thermal energy injection around BH particle, proportional to BH accretion rate). constant fraction of luminosity coupled in kernel around BH. cite Springel, Di Matteo, and Hernquist, 2005, MNRAS, 361, 776
#SINK_WIND_KICK=1                 # mechanical (wind from accretion disk/BH with specified mass/momentum/energy-loading relative to accretion rate). gas in kernel given stochastic 'kicks' at fixed velocity. (>0=isotropic, <0=collimated, absolute value sets momentum-loading in L/c units). cite Angles-Alcazar et al., 2017, MNRAS, 464, 2840
#SINK_WIND_SPAWN=2                # mechanical (wind from accretion disk/BH with specified mass/momentum/energy-loading relative to accretion rate). spawn virtual 'wind' particles to carry BH winds out. value=min number spawned per spawn-step. Cite Torrey et al 2020MNRAS.497.5292T, Su et al., arXiv:2102.02206, and Grudic et al. arXiv:2010.11254 for use and numerical methods and tests
#SINK_COSMIC_RAYS                 # cosmic ray: explicitly inject and transport CRs from BH. set injection energy efficiency. injected alongside mechanical energy (params file sets ratios of energy in different mechanisms). these currently build on the architecture of the SINK_WIND modules, one of those must be enabled, along with the usual cosmic-ray physics set of modules for CR transport. the same restrictions apply to this as CR modules. developed by P. Hopkins
# --- radiative: [LEBRON] these currently are built on the architecture of the FIRE stellar FB modules, and require some of those be active. their use therefore follows FIRE policies (see details above). however, if explicit-radiation-hydrodynamics is enabled, users can achieve this functionality entirely in the public code, with appropriate hooks in the cooling functions
#SINK_COMPTON_HEATING             # enable Compton heating/cooling from BHs in cooling function (needs SINK_PHOTONMOMENTUM). cite Hopkins et al., 2016, MNRAS, 458, 816
#SINK_HII_HEATING                 # photo-ionization feedback from BH (needs GALSF_FB_FIRE_RT_HIIHEATING). cite Hopkins et al., arXiv:1702.06148
#SINK_PHOTONMOMENTUM              # continuous long-range IR radiation pressure acceleration from BH (needs GALSF_FB_FIRE_RT_LONGRANGE). cite Hopkins et al., arXiv:1702.06148
## ----------------------------------------------------------------------------------------------------
# ----- output options
## ----------------------------------------------------------------------------------------------------
#SINK_OUTPUT_MOREINFO             # output additional info to "sink_details" on timestep-level, following Angles-Alcazar et al. 2017, MNRAS 472, 109 (use caution: files can get very large if many BHs exist)
#SINK_CALC_DISTANCES              # calculate distances for all particles to closest BH for, e.g., refinement, external potentials, etc. cite Garrison-Kimmel et al., MNRAS, 2017, 471, 1709
####################################################################################################



####################################################################################################
# ---- Radiative Cooling & Thermo-Chemistry
# ------ Modules designed to follow radiative cooling in optically thin/thick limits, with ionized/atomic/molecular gas-phase chemistry.
# ------  These are generally designed to be applicable at densities << 1e-6 g/cm^3, or nH << 10^18 atoms/cm^3 -- i.e. densities from
# ------  low-density inter-galactic medium through proto-planetary/stellar disks, but not planetary or stellar interiors (for those, other modules are more appropriate).
# ------  Proper citations are below and in User Guide; all users should cite Hopkins et al. 2017 (arXiv:1702.06148), where Appendix B details the cooling physics
####################################################################################################
## ----------------------------------------------------------------------------------------------------
#COOLING                        # top-level switch to enable radiative cooling and heating. if nothing else enabled, uses Hopkins et al. arXiv:1702.06148 cooling physics. if GALSF, also external UV background read from file "TREECOOL" (included in the cooling folder; be sure to cite its source as well, given in the TREECOOL file)
#METALS                         # top-level switch to enable tracking metallicities / different heavy elements (with multiple species optional) for gas and stars [must be included in ICs or injected via dynamical feedback; needed for some routines]
## ----------------------------------------------------------------------------------------------------
# ---- additional cooling physics options within the default COOLING (Hopkins et al. 2017) module
## ----------------------------------------------------------------------------------------------------
#COOL_METAL_LINES_BY_SPECIES    # use full multi-species-dependent cooling tables ( http://www.tapir.caltech.edu/~phopkins/public/spcool_tables.tgz, or the GitHub site); requires METALS on; cite Wiersma et al. 2009 (MNRAS, 393, 99) in addition to Hopkins et al. 2017 (arXiv:1702.06148)
#COOL_LOW_TEMPERATURES          # allow fine-structure and molecular cooling to ~10 K; account for optical thickness and line-trapping effects with proper opacities [requires METALS]. attempts to interpolate between optically-thin and optically-thick cooling limits even if explicit rad-hydro not enabled. Cite Hopkins et al. arXiv:1702.06148
#COOL_MOLECFRAC=6               # track molecular H2 fractions for use in COOL_LOW_TEMPERATURES and thermochemistry using different estimators: (1) simplest, fit to density+temperature from Glover+Clark 2012; (2) Krumholz+Gnedin 2010 fit vs. column+metallicity; (3) Gnedin+Draine 2014 fit vs column+metallicity+MW radiation field; (4) Krumholz, McKee, & Tumlinson 2009 local equilibrium cloud model vs column, metallicity, incident FUV; (5) explicit local equilibrium H2 fraction explicitly tracking rates, metals, clumping, shielding, UV [cite Hopkins et al. 2023MNRAS.519.3154H]; (6) explicit non-equilibrium integration of rates in level 5 [cite Hopkins et al. 2023MNRAS.519.3154H]
## ----------------------------------------------------------------------------------------------------
# ---- GRACKLE: alternative chemical network using external libraries for solving thermochemistry+cooling. These treat molecular hydrogen, in particular, in more detail than our default networks, and are more accurate for 'primordial' (e.g. 1st-star) gas. But they have less-accurate treatment of
# ----            effects such as dust-gas coupling and radiative feedback (Compton and photo-electric and local ionization heating) and high-optical-depth effects, so are usually less accurate for low-redshift, metal-rich star formation or planet formation simulations.
## ----------------------------------------------------------------------------------------------------
#COOL_GRACKLE                   # enable Grackle: cooling+chemistry package (requires COOLING above; https://grackle.readthedocs.org/en/latest ); see Grackle code for their required citations
#COOL_GRACKLE_CHEMISTRY=1       # choose Grackle cooling chemistry: (0)=tabular, (1)=Atomic, (2)=(1)+H2+H2I+H2II, (3)=(2)+DI+DII+HD. Modules with dust and/or metal-line cooling require METALS also
#COOL_GRACKLE_APIVERSION=1      # set the version of the grackle api: =1 (default) is compatible with versions of grackle below 2.2. After 2.2 significant changes to the grackle api were made which require different input formats, which require setting this to =2 or larger. note newest grackle apis may not yet be compatible with the hooks here!
## ----------------------------------------------------------------------------------------------------
# ---- CHIMES: alternative non-equilibrium chemical (ion+atomic+molecular) network, developed by Alex Richings. The core methods are laid out in 2014MNRAS.440.3349R, 2014MNRAS.442.2780R. These should be cited in any paper that uses the modules below.
# ----   Per permission from Alex Richings, the CHIMES modules are now public. However recall that Alex Richings is the lead developer of CHIMES, please contact Alex or Joop Schaye, or Ben Oppenheimer to obtain the relevant permissions to port beyond GIZMO or questions about CHIMES
# ----   This implementation of CHIMES has additional hooks to use the various gizmo radiation fields if desired. The modules solve a large molecular and ion network, so can trace predictive chemistry for species in dense ISM gas in much greater detail than the other modules above (at additional CPU cost)
## ----------------------------------------------------------------------------------------------------
#CHIMES                         # top-level switch to enable CHIMES. Requires COOLING above. Also, requires COOL_METAL_LINES_BY_SPECIES to include metals.
#CHIMES_SOBOLEV_SHIELDING       # enables local self-shielding for different species, using a Sobolev-like length scale
#CHIMES_HII_REGIONS             # disables shielding withing HII region (requires FIRE modules for radiation transport/coupling: uses GALSF_FB_FIRE_RT_HIIHEATING, and permissions follow those modules)
#CHIMES_STELLAR_FLUXES          # couple UV fluxes from the luminosity tree to CHIMES (requires FIRE modules for radiation transport/coupling: use permissions follow those modules)
#CHIMES_TURB_DIFF_IONS          # turbulent diffusions of CHIMES abundances. Requires TURB_DIFF_METALS and TURB_DIFF_METALS_LOWORDER (see modules for metal diffusion above: use/citation policy follows those)
#CHIMES_METAL_DEPLETION         # uses density-dependent metal depletion factors (Jenkins 2009, De Cia et al. 2016) to obtain gas-phase abundances for chemical network
## ------------ CHIMES de-bugging and special behaviors ------------------------------------------------------------------------
#CHIMES_HYDROGEN_ONLY           # hydrogen-only. This is ignored if METALS are also set.
#CHIMES_REDUCED_OUTPUT          # full CHIMES abundance array only output in some snapshots
#CHIMES_NH_OUTPUT               # write out column densities of gas particles to snapshots
#CHIMES_INITIALISE_IN_EQM       # initialise CHIMES abundances in equilibrium at the start of the simulation
## ----------------------------------------------------------------------------------------------------
# ---- JACO: code-generated microphysics solver. Uses the jaco Python package (https://github.com/mikegrudic/jaco)
# ----   to generate a backward-Euler Newton-Raphson solver with analytic Jacobian for a user-specified cooling/chemistry model.
# ----   Set JACO=model_name where model_name is a jaco model module (e.g. wind_comparison, starforge). The build system runs
# ----   the jaco codegen to produce microphysics_func_jac.cc, jaco_eos.cc, and associated headers in cooling/.
# ----   When enabled, JACO replaces the default COOLING solver for thermochemistry (the call is made from do_the_cooling_for_particle).
# ----   Requires the jaco Python package to be installed (pip install jaco). Requires COOLING to be enabled.
## ----------------------------------------------------------------------------------------------------
#JACO=model_name            # enable the jaco  microphysics solver with the specified model. requires COOLING.
####################################################################################################



############################################################################################################################
# -------------------------------------- Radiative Transfer & Radiation Hydrodynamics:
# -------------------------------------------- modules developed by PFH with David Khatami, Mike Grudic, and Nathan Butcher (special  thanks to Alessandro Lupi)
# --------------------------------------------  these are now public, but if used, cite the appropriate paper[s] for their methods/implementation in GIZMO
############################################################################################################################
# -------------------- methods for calculating photon propagation (one, and only one, of these MUST be on for RT). whatever method is used, you must cite the appropriate methods paper.
#RT_FLUXLIMITEDDIFFUSION                # RT solved using moments-based 0th-order flux-limited diffusion approximation (constant, always-isotropic Eddington tensor). cite Hopkins & Grudic, 2018, arXiv:1803.07573
#RT_M1                                  # RT solved using moments-based 1st-order M1 approximation (solve fluxes and tensors with M1 closure; gives better shadowing; currently only compatible with explicit diffusion solver). cite Hopkins & Grudic, 2018, arXiv:1803.07573
#RT_OTVET                               # RT solved using moments-based 0th-order OTVET approximation (optically thin Eddington tensor, but interpolated to thick when appropriate). cite Hopkins & Grudic, 2018, arXiv:1803.07573
#RT_LOCALRAYGRID=1                      # RT solved using exact method of Jiang et al. (each cell carries a mesh in phase space of the intensity directions, rays directly solved over the 6+1D direction-space-frequency-time mesh [value=number of polar angles per octant: N_rays=4*value*(value+1)]. this is still in development, DO NOT USE without contacting PFH
#RT_LEBRON                              # RT solved using ray-based LEBRON approximation (locally-extincted background radiation in optically-thin networks; default in the FIRE simulations). cite Hopkins et al. 2012, MNRAS, 421, 3488 and Hopkins et al. 2018, MNRAS, 480, 800 [former developed methods and presented tests, latter details all algorithmic aspects explicitly]
# -------------------- solvers (numerical) --------------------------------------------------------
#RT_SPEEDOFLIGHT_REDUCTION=1            # set to a number <1 to use the 'reduced speed of light' approximation for photon propagation (C_eff=C_true*RT_SPEEDOFLIGHT_REDUCTION)
#RT_COMOVING                            # solve RHD equations formulated in the comoving frame, as compared to the default mixed-frame formulation; see Mihalas+Mihalas 84
#RT_DIFFUSION_IMPLICIT                  # solve the diffusion part of the RT equations (if needed) implicitly with Conjugate Gradient iteration (Petkova+Springel): less accurate and only works with some methods, but allows larger timesteps [otherwise more accurate explicit used]
# -------------------- physics: wavelengths+coupled RT-chemistry networks (if any of these is used, cite Hopkins et al. 2018, MNRAS, 480, 800) -----------------------------------
#RT_SOURCES=1+16+32                     # source types for radiation given by bitflag (1=2^0=gas,16=2^4=new stars,32=2^5=BH)
#RT_XRAY=3                              # x-rays: 1=soft (0.5-2 keV), 2=hard (>2 keV), 3=soft+hard; used for Compton-heating
#RT_CHEM_PHOTOION=2                     # ionizing photons: 1=H-only [single-band], 2=H+He [four-band]
#RT_LYMAN_WERNER                        # lyman-werner [narrow H2 dissociating] band
#RT_PHOTOELECTRIC                       # far-uv (8-13.6eV): track photo-electric heating photons + their dust interactions
#RT_NUV                                 # near-UV: 1550-3600 Angstrom (where direct stellar emission dominates)
#RT_OPTICAL_NIR                         # optical+near-ir: 3600 Angstrom-3 micron (where direct stellar emission dominates)
#RT_FREEFREE                            # scattering from Thompson, absorption+emission from free-free, appropriate for fully-ionized plasma
#RT_INFRARED                            # infrared: photons absorbed in other bands are down-graded to IR: IR radiation + dust + gas temperatures evolved independently. Requires METALS and COOLING.
#RT_GENERIC_USER_FREQ                   # example of an easily-customizable, grey or narrow band: modify this to add your own custom wavebands easily!
#RT_OPACITY_FROM_EXPLICIT_GRAINS        # calculate opacities back-and-forth from explicitly-resolved grain populations. Cite Hopkins et al., arXiv:2107.04608, if used.
# -------------------- radiation pressure options -------------------------------------------------
#RT_DISABLE_RAD_PRESSURE                # turn off radiation pressure forces (included by default)
#RT_RAD_PRESSURE_OUTPUT                 # print radiation pressure to file (requires some extra variables to save it)
#RT_ENABLE_R15_GRADIENTFIX              # for moments [FLD/OTVET/M1]: enable the Rosdahl+ 2015 approximate 'fix' (off by default) for gradients under-estimating flux when under-resolved by replacing it with E_nu*c
## ----------------------------------------------------------------------------------------------------
# ----------- alternative, test-problem, or special behavior options
## ----------------------------------------------------------------------------------------------------
#RT_SELFGRAVITY_OFF                     # turn off gravity: if using an RT method that needs the gravity tree (FIRE, OTVET), use this -instead- of SELFGRAVITY_OFF to safely turn off gravitational forces
#RT_USE_TREECOL_FOR_NH=6                # uses the TreeCol method to estimate effective optical depth using non-local information from the gravity tree; cite Clark, Glover & Klessen 2012 MNRAS 420 754. Value specifies the number of angular bins on the sky for ray-tracing column density.
#RT_INJECT_PHOTONS_DISCRETELY           # do photon injection in discrete packets, instead of sharing a continuous source function. works better with adaptive timestepping (default with GALSF)
#RT_USE_GRAVTREE_SAVE_RAD_FLUX          # save radiative fluxes incident on each cell if using RHD methods that propagate fluxes through the gravity tree when these wouldn't be saved by default
#RT_REPROCESS_INJECTED_PHOTONS          # re-process photon energy while doing the discrete injection operation conserving photon energy, put only the un-absorbed component of the current band into that band, putting the rest in its "donation" bin (ionizing->optical, all others->IR). This would happen anyway during the routine for resolved absorption, but this may more realistically handle situations where e.g. your dust destruction front is at totally unresolved scales and you don't want to spuriously ionize stuff on larger scales. Assume isotropic re-radiation, so inject only energy for the donated bin and not net flux/momentum. follows STARFORGE methods (Grudic+ arXiv:2010.11254) - cite this
#RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION   # uses a solid-angle as opposed to simple kernel weight (requires extra passes) for depositing radiation from sinks/BHs when the direct deposition is used. also ensures the sink uses a 2-way search to ensure overlapping diffuse gas gets radiation. cite Grudic+ arXiv:2010.11254
#RT_ISRF_BACKGROUND=1                   # include Draine 1978 ISRF for photoelectric heating (appropriate for solar circle, must be re-scaled for different environments); rescaled by a constant normalization given by this constant, if defined
## ----------------------------------------------------------------------------------------------------
# ----------- Transport subcycling: allow RT and/or CR transport to take multiple smaller steps per hydro step  [HIGHLY EXPERIMENTAL]
## ----------------------------------------------------------------------------------------------------
#TRANSPORT_SUBCYCLE=100                 # max number of transport (RT+CR) subcycles per hydro step. decouples the transport CFL from the hydro timestep.
#TRANSPORT_SUBCYCLE_COOLING             # also subcycle cooling within the transport subcycle loop (otherwise cooling runs once per hydro step
####################################################################################################



####################################################################################################
## --------------------------------------------------------------------------------------------------
# --------- Cosmic Rays & Relativistic Particles: MHD-PIC and Cosmic Ray-MHD simulations
# ---------  This is developed by P. Hopkins, with major contributions from TK Chan for the CR-fluid and S Ji for the MHD-PIC implementations. The fluid modules solve the CR transport in the two moment limit [second-order expansion of the collisionless boltzmann eqn] as in Hopkins et al. arXiv:2103.10443,2002.06211,2202.05283 (cite that paper for the modern version, but also Chan et al. 2019MNRAS.488.3716C, Hopkins et al. arXiv:2002.06211); can set here the reduced speed of light (maximum free-streaming speed) in code units. requires MAGNETIC for proper behavior for fully-anisotropic equations
# ---------  Note this formulation with appropriate physical terms is valid in both the frequent-scattering and free-streaming limits, it does NOT require CRs be in the 'fluid-like' or frequent-scattering limit, only that resolved scales are large compared to CR gyro radii. For resolved scales smaller than CR gyro radii, use the PIC methods below.
## --------------------------------------------------------------------------------------------------
#PIC_MHD                          #  hybrid MHD-PIC simulations for relativistic particles / cosmic rays (particle type=3). need to set 'subtype'. cite Ji, Squire, & Hopkins, arXiv:2112.00752
#PIC_SPEEDOFLIGHT_REDUCTION=1     #  factor to reduce the speed-of-light for mhd-pic simulations (relative to true value of c). requires PIC_MHD. cite Ji & Hopkins, arXiv:2111.14704
#GRAIN_FLUID_AND_PIC_BOTH_DEFINED #  this tells the code that both GRAIN_FLUID (dust) and MHD_PIC (for e.g. cosmic rays or other applications) are simultaneously active. cite Ji, Squire, & Hopkins, arXiv:2112.00752
## -----
#COSMIC_RAY_FLUID                 #  top-level switch to evolve the distribution function (continuum limit) of a population of CRs. includes losses/gains, coupling to gas, streaming, diffusion. this uses the multi-moment expansion valid generally in free-streaming and fluid-like limits so long as CR gyro radii are small, as derived and implemented in Hopkins et al. arXiv:2103.10443,2002.06211,2202.05283 (cite these). all operators will be anisotropic as it should unless MHD is turned off
#CRFLUID_DIFFUSION_MODEL=0        #  determine how coefficients for CR transport scale. 0=spatial/temporal constant diffusivity (power law in rigidity), -1=no diffusion (but stream at vAlfven), values >=1 correspond to different literature scalings for the coefficients (see user guide). cite Hopkins et al. arXiv:2002.06211 for all derivations, models here (and see refs therein for sources for some of the models re-derived and implemented here)
#CRFLUID_EVOLVE_SPECTRUM=2        #  follow a spectrally-resolved CR population of (=1:e-/p, =2:e-,e+,p,anti-p,B,C/N/O,stable (7-9)Be, unstable 10Be) from ~MeV-TeV (and other species if extended network is enabled), including injection and adiabatic+hadronic/catastrophic/pionic/fragmentation+inverse compton+ionization+coulomb+bremstrahhlung+gyroresonant/streaming+synchrotron+annihilation+radioactive losses and (optionally) re-acceleration. cite Hopkins et al. 2022MNRAS.516.3470H and arXiv:2202.05283 for implementation+electrons+loss/gain terms and Girichidis+ 2020MNRAS.491..993G for the fundamental spectral bin-to-bin method
#CRFLUID_EVOLVE_SCATTERINGWAVES   #  follows Zweibel+13,17 and Thomas+Pfrommer 18 to explicitly evolve gyro-resonant wave packets which define the (gyro-averaged) CR scattering rates; requires MAGNETIC and COOLING for detailed MHD and ionization+thermal states. cite Hopkins et al. arXiv:2002.06211 for numerical implementation into various solvers here
#CRFLUID_SPEEDOFLIGHT_REDUCTION=1 #  reduce the speed of light for CRs specifically (can be done separately from e.g. radiation), following the self-consistent formulation in Hopkins et al. arXiv:2103.10443, which should be cited for these methods
## -----
#COSMIC_RAY_SUBGRID_LEBRON        #  uses a simplified sub-grid LEBRON-type model for CR transport for cheap approximation in galaxy simulations. Cite Hopkins et al. arXiv:2211.05811 for any use (see there for methods details)
## -----
#CRFLUID_ALT_RSOL_FORM            #  enable the alternative reduced-speed-of-light formulation where 1/reduced_c appears in front of all D/Dt terms, as opposed to only in the flux equation. converges more slowly but accurately. will get made default, we think, once de-bugged
####################################################################################################



####################################################################################################
# --------------------------------------- Multi-Threading and Parallelization options
####################################################################################################
#OPENMP                         # top-level switch for explicit OpenMP implementation (can turn on here, or enable in Makefile for your machine)
#MULTIPLEDOMAINS=16             # Multi-Domain option for the top-tree level (alters load-balancing)
#DOMAIN_TIMEBINS=0              # Domain timebin cost weighting: 0=frequency-weighted costs, 1=full per-timebin balancing (Gadget-4 scheme). Omit for unweighted.
####################################################################################################



####################################################################################################
# --------------------------------------- Input/Output options
####################################################################################################
#OUTPUT_ADDITIONAL_RUNINFO      # enables extended simulation output data (can slow down machines significantly in massively-parallel runs)
#OUTPUT_IN_DOUBLEPRECISION      # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION       # input files assumed to be in double precision (otherwise float is assumed)
#INPUT_POSITIONS_IN_DOUBLE      # as above, but specific to the ICs file
#OUTPUT_POTENTIAL               # forces code to compute+output potentials in snapshots
#OUTPUT_TIDAL_TENSOR            # writes tidal tensor (computed in gravity) to snapshots
#OUTPUT_ACCELERATION            # output physical acceleration of each particle in snapshots
#OUTPUT_HYDROACCELERATION       # output the 'hydrodynamic' (includes -all- stress tensor terms) acceleration. if enabled with 'OUTPUT_ACCELERATION', that will output the gravitational acceleration, so the sum of the two is the total
#IO_HERMITE_SYNC                # output HermiteSyncCoordinates/HermiteSyncVelocities: a mutually consistent (r,v) pair at the snapshot time for Hermite-integrated particles. The ordinary Coordinates/Velocities are NOT such a pair (positions are drifted to the output time, velocities are left at the last kick), so vis-viva, orbital elements and kinetic energy are all evaluated on a state that never existed. Predicted, not corrected: 4th order in position, 3rd in velocity
#OUTPUT_CHANGEOFENERGY          # outputs rate-of-change of internal energy of gas particles in snapshots
#OUTPUT_VORTICITY               # outputs the vorticity vector
#OUTPUT_GRADIENT_RHO            # outputs the gradients of the gas density field
#OUTPUT_GRADIENT_VEL            # outputs the full velocity gradient tensor field for the gas
#OUTPUT_BFIELD_DIVCLEAN_INFO    # outputs the phi, phi-gradient, and numerical div-B fields used for de-bugging MHD simulations
#OUTPUT_TIMESTEP                # outputs timesteps for each particle
#OUTPUT_SOFTENING               # outputs force softening for each particle
#OUTPUT_COOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_COOLRATE_DETAIL         # outputs cooling rate term by term [saves all individually to snapshot]
#OUTPUT_LINEOFSIGHT				# enables on-the-fly output of Ly-alpha absorption spectra. requires METALS and COOLING.
#OUTPUT_LINEOFSIGHT_SPECTRUM    # computes power spectrum of these (requires additional code integration)
#OUTPUT_LINEOFSIGHT_PARTICLES   # computes power spectrum of these (requires additional code integration)
#OUTPUT_POWERSPEC               # compute and output cosmological power spectra. requires BOX_PERIODIC and PMGRID.
#OUTPUT_RECOMPUTE_POTENTIAL     # update potential every output even it EVALPOTENTIAL is set
#OUTPUT_DENS_AROUND_NONGAS      # output gas density in neighborhood of stars [collisionless particle types], not just gas
#OUTPUT_DELAY_TIME_HII          # output DelayTimeHII. Requires GALSF_FB_FIRE_RT_HIIHEATING (and corresponding flags/permissions set)
#OUTPUT_MOLECULAR_FRACTION      # output the code-estimated molecular mass fraction [needs COOLING], for e.g. approximate molecular fraction estimators (as opposed to detailed chemistry modules, which already output this)
#OUTPUT_TEMPERATURE             # output the in-code gas temperature
#OUTPUT_SINK_ACCRETION_HIST     # save full accretion histories of sink (BH/star/etc) particles
#OUTPUT_SINK_FORMATION_PROPS    # save at-formation properties of sink particles
#OUTPUT_SINK_DISTANCES          # saves the distance to the nearest sink, if SINK_CALC_DISTANCES is enabled, to snapshots
#OUTPUT_RT_RAD_FLUX             # save flux vector for radiation methods that explictly evolve the flux (e.g. M1)
#OUTPUT_RT_RAD_OPACITY          # save opacities for the different bands for explicit radiation-hydro methods
#OUTPUT_UNSPAWNED_SINKMASS      # save the unspawned mass variable used for sink cell-spawning modules
#INPUT_READ_KERNELRADIUS        # force reading rkern from IC file (instead of re-computing them; in general this is redundant but useful if special guesses needed)
#INPUT_READ_SINKPROPS           # force reading sink properties including sink radius, zams mass, luminosity, age, etc, from ICs file if it includes sink particles and the IC is designed for use with the single-star modules
#INPUT_READ_TEMPERATURE         # force reading the gas Temperature block from gadget-binary ICs (ICFormat 1/2), which are skipped by default because a missing block cannot be detected there and the read would fail with end-of-file. Only enable if your binary ICs genuinely contain one. HDF5 ICs are unaffected: they always read it when present
#OUTPUT_TWOPOINT_ENABLED        # allows user to calculate mass 2-point function by enabling and setting restartflag=5
#IO_COMPRESS_HDF5     		    # write HDF5 in compressed form (will slow down snapshot I/O and may cause issues on old machines, but reduce snapshots 2x)
#IO_SUPPRESS_TIMEBIN_STDOUT=10  # only prints timebin-list to log file if highest active timebin index is within N (value set) of the highest timebin (dt_bin=2^(-N)*dt_bin,max)
#IO_SUBFIND_READFOF_FROMIC      # try read already existing FOF files associated with a run instead of recomputing them: not de-bugged
#OUTPUT_TURB_DIFF_DYNAMIC_ERROR # save error terms from localized dynamic Smagorinsky model to snapshots
#IO_MOLECFRAC_NOT_IN_ICFILE     # special flag needed if using certain molecular modules with restart flag=2 where molecular data was not in that snapshot, to tell code not to read it
#IO_REDUNDANT_BACKUP_RESTARTFILE_FREQUENCY=3  # keep an extra set of backup files that are IO_REDUNDANT_BACKUP_RESTARTFILE_FREQUENCY number of restarts old (allows for soft restarts from an older position)
#IO_GRADUAL_SNAPSHOT_RESTART    # when restarting from a snapshot (flag=2) start every element on the shortest possible timestep - can reduce certain transient behaviors from the restart procedure
#IO_SINKS_ONLY_SNAPSHOT_FREQUENCY=0 # determines the number of snapshots with reduced data (sinks only) per full snapshots (gas+sinks+other), e.g., setting this to 2 means 2/3 of the snapshots will be reduced, 1/3 will have full data. Setting this to 0 disables it. developed by DG.
####################################################################################################



####################################################################################################
# -------------------------------------------- De-Bugging & special (usually test-problem only) behaviors
####################################################################################################
# --------------------
# ----- General De-Bugging and Special Behaviors
#DEVELOPER_MODE                    # allows you to modify various numerical parameters (courant factor, etc) at run-time
#FORCE_EQUAL_TIMESTEPS             # force the code to use a single universal timestep (can change in time, but all particles advance together). chosen as minimum of any particle that step.
#STOP_WHEN_BELOW_MINTIMESTEP       # forces code to quit when stepsize wants to go below MinSizeTimestep specified in the parameterfile
# --------------------
# ----- Hydrodynamics (and MHD)
#FREEZE_HYDRO                      # zeros all fluxes from RP and doesn't let particles move (for testing additional physics layers)
#EOS_ENFORCE_ADIABAT=(1.0)         # if set, this forces gas to lie -exactly- along the adiabat P=EOS_ENFORCE_ADIABAT*(rho^GAMMA)
#HYDRO_REPLACE_RIEMANN_KT          # replaces the hydro Riemann solver (HLLC) with a Kurganov-Tadmor flux derived in Panuelos, Wadsley, and Kevlahan, 2019. works with MFM/MFV/fixed-grid methods [-without- MHD active, but other modules are fine]. more diffusive, but smoother, and more stable convergence results
#SLOPE_LIMITER_TOLERANCE=1         # sets the slope-limiters used. higher=more aggressive (less diffusive, but less stable). 1=default. 0=conservative. use on problems where sharp density contrasts in poor particle arrangement may cause errors. 2=use the original GIZMO paper (more aggressive) slope-limiters. more accurate for smooth problems, but these can introduce numerical instability in problems with poorly-resolved large noise or density contrasts (e.g. multi-phase, self-gravitating flows)
#ENERGY_ENTROPY_SWITCH_IS_ACTIVE   # enable energy-entropy switch as described in GIZMO methods paper. This can greatly improve performance on some problems where the the flow is very cold and highly super-sonic. it can cause problems in multi-phase flows with strong cooling, though, and is not compatible with non-barytropic equations of state
#FORCE_ENTROPIC_EOS_BELOW=(0.01)   # set (manually) the alternative energy-entropy switch which is enabled by default in MFM/MFV: if relative velocities are below this threshold, it uses the entropic EOS
#HYDRO_KERNEL_SURFACE_VOLCORR      # attempt to correct SPH/MFM/MFV cell volumes for free-surface effects, using the estimated boundary correction for the Wendland C2 kernel (works with others but most accurate for this) based on asymmetry of neighbors within kernel, as calibrated in Reinhardt & Stadel 2017 (arXiv:1701.08296), see e.g. their Fig 3
#DISABLE_SURFACE_VOLCORR           # disables HYDRO_KERNEL_SURFACE_VOLCORR if it would be set by default (e.g. if EOS_ELASTIC is enabled)
#HYDRO_PARTITION_UNITY_IMPROVE_FD  # apply the first-derivative (FD) partition-of-unity volume correction from Massaro Acha, Alonso Asensio & Dalla Vecchia 2026, which corrects cell volume estimates using the spatial gradient of the kernel support size. This is a low-cost correction (no extra neighbor loop) that improves the cubature rule accuracy by a factor ~3-4x. Only relevant for meshless (MFM/MFV) methods. Complementary to HYDRO_VOLUME_CORRECTIONS (which uses a separate neighbor loop for a different quadrature correction)
#HYDRO_EXPLICITLY_INTEGRATE_VOLUME # explicitly integrate the kernel continuity equation for cell volumes (giving e.g. densities), as in e.g. Monaghan 2000, but with a term that relaxes the integrated cell volume back to the explicitly evaluated kernel calculation on a timescale ~10 t_cross where t_cross ~ MAX(H_kernel , L_grad) / MIN(cs_eff) where L_grad is the density gradient scale length and cs_eff the minimum sound/torsion/tension wave speed. This module ONLY makes sense for strictly fixed-mass (SPH/MFM) methods
#DISABLE_EXPLICIT_VOLUME_INTEGRATION # disables HYDRO_EXPLICITLY_INTEGRATE_VOLUME if it would be set by default (e.g. if EOS_ELASTIC is enabled)
#SPH_DISABLE_CD10_ARTVISC          # for SPH only: Disable Cullen & Dehnen 2010 'inviscid sph' (viscosity suppression outside shocks); just use Balsara switch
#SPH_DISABLE_PM_CONDUCTIVITY       # for SPH only: Disable mixing entropy (J.Read's improved Price-Monaghan conductivity with Cullen-Dehnen switches)
#MHD_ALTERNATIVE_LEAPFROG_SCHEME   # use alternative leapfrog where magnetic fields are treated like potential/positions (per Federico Stasyszyn's suggestion): still testing
# --------------------
# ----- Cooling and Additional Fluid Physics
#COOLING_OPERATOR_SPLIT            # do the hydro heating/cooling in operator-split fashion from chemical/radiative. slightly more accurate when tcool >> tdyn, but much noisier when tcool << tdyn
#COOL_LOWTEMP_THIN_ONLY            # in the COOL_LOW_TEMPERATURES module, neglect the suppression of cooling at very high surface densities due to the opacity limit (disables limiter in Eqs B29-B30, Hopkins et al arXiv:1702.06148)
#SUPER_TIMESTEP_DIFFUSION          # use super-timestepping to accelerate integration of diffusion operators [for testing or if there are stability concerns]
#TURB_DRIVING_UPDATE_FORCE_ON_TURBUPDATE # if this is enabled, we only update as frequently as the driving phases are recomputed, as set by TurbDrive_TimeBetweenTurbUpdates. Only enable as an optimization if the cost of evaluating the turbulent force is large. To avoid large errors, TurbDrive_TimeBetweenTurbUpdates must be set by-hand to be << lambda_min / V where V is the typical turbulent velocity and lambda_min is the smallest driven wavelength.
# --------------------
# ----- Gravity (force/timestep/potential/gravity-tree related options)
#EVALPOTENTIAL                     # computes gravitational potential (even if not otherwise needed)
#GRAVITY_HYBRID_OPENING_CRIT       # use -both- Barnes-Hut + relative angle opening criterion for the gravity tree (normally choose one or the other)
#TIDAL_TIMESTEP_CRITERION          # replace standard acceleration-based timestep criterion with one based on the tidal tensor norm, which is more accurate and adaptive (testing, but may be promoted to default code)
#ADAPTIVE_TREEFORCE_UPDATE=0.06    # use the tidal timescale to estimate how often gravity needs to be updated, updating a gas cell's gravity no more often than ADAPTIVE_TREEFORCE_UPDATE * dt_tidal, the factor N_f in Grudic 2020 arxiv:2010.13792 (cite this). Smaller is more accurate, larger is faster, should be tuned for your problem if used.
#RANDOMIZE_GRAVTREE                # move the top tree node around randomly so that treeforce errors are not correlated between one treebuild and another. cite Grudic+ arXiv:2010.11254
#GRAVITY_SPHERICAL_SYMMETRY=0      # modifies the tree gravity solver to give the solution assuming spherical symmetry about the origin (if BOX_PERIODIC is not enabled) or the box center. Useful for IC generation and test problems. Numerical value specifies a minimum softening length. (cite Lane et al., arXiv:2110.14816)
#SINGLE_STAR_DIRECT_GRAVITY_RADIUS=1000. # enforce direct gravity summation for star-star gravity interactions *and* tree searches (e.g. for timestepping, binarity checks) within this radius *in AU*
#SINGLE_STAR_DIRECT_GRAVITY        # exact star-star (type 5) gravity at ALL separations, by brute-force summation over a table of every star broadcast to every task, with the tree dropping star mass for star targets so nothing is counted twice. Unlike the _RADIUS option above there is no cutoff, so the tree never approximates one star's pull on another. Costs O(N_star^2) per gravity step and stores O(N_star) per task, so it is for runs with few stars; incompatible with SINGLE_STAR_TIMESTEPPING>0 (whose binary search reads the star-star tree interactions this removes)
#ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT=(1) # impose a hard upper limit (arbitrarily) to the softening kernel size for adaptive gravitational softening for gas cells
#ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING # use the 'average the forces' implementation of force symmtrization in the gravity solver, instead of the default 'calculate forces for the larger softening', when two particles or cells with different force softening values interact inside each others kernels
# --------------------
# ----- Particle IDs
#TEST_FOR_IDUNIQUENESS             # explicitly check if particles have unique id numbers (only use for special behaviors)
#ASSIGN_NEW_IDS                    # assign IDs on startup instead of reading from ICs
#NO_CHILD_IDS_IN_ICS               # IC file does not have child IDs: do not read them (used for compatibility with snapshot restarts from old versions of the code)
# --------------------
# ----- Particle Merging/Splitting/Deletion/Boundaries
#MAINTAIN_TREE_IN_REARRANGE        # don't rebuild the domains/tree every time a particle is spawned - salvage the existing one by redirecting pointers as needed. cite Grudic+ arXiv:2010.11254
#PREVENT_PARTICLE_MERGE_SPLIT      # don't allow gas particle splitting/merging operations
#PREVENT_PARTICLE_MERGE            # don't allow gas MERGING but still allow splitting (cf. PREVENT_PARTICLE_MERGE_SPLIT, which disables both). Gives a diffusion-free, finest-resolution reference for spawned outflows; cf. SINK_SPAWN_NO_MERGE, which protects only spawned cells.

#PARTICLE_EXCISION                 # enable dynamical excision (remove particles within some radius)
#MERGESPLIT_HARDCODE_MAX_MASS=(1.0e-6)   # manually set maximum mass for particle merge-split operations (in code units): useful for snapshot restarts and other special circumstances
#MERGESPLIT_HARDCODE_MIN_MASS=(1.0e-7)   # manually set minimum mass for particle merge-split operations (in code units): useful for snapshot restarts and other special circumstances
#PARTICLE_MERGE_SPLIT_EVERY_TIMESTEP     # force merge/split operations to occur every timestep, instead of only on domain decomposition steps
#MHD_CONSERVE_B_ON_REFINEMENT      # redefine B after density step after a refinement/de-refinement operation so as to exactly conserve "B" between the steps, as compared to conserving the code "VB" between the steps [the default conserved integrated quantity]
#ENERGY_BUDGET_DIAGNOSTIC           # energy diagnostics, no behaviour change. Reports (a) the total gas energy (thermal+kinetic) at every full synchronization, i.e. when All.HighestActiveTimeBin == All.HighestOccupiedTimeBin -- the only point at which Vel and InternalEnergy share a clock and the sum is a genuine conserved quantity, so the 'Energy (synced)' line is the real conservation test and in an adiabatic run should track the injected energy exactly; and (b) cumulative wakeup tallies from timestep.cc -- the residual left by the kick reversal, and processed wakeups split by trigger (hydro / merge / split / newly-spawned) and demotion depth. Do NOT substitute snapshots for the synced total: io.cc writes Velocities from P[].Vel (kick-time) but InternalEnergy from InternalEnergyPred (drift-time), so a snapshot sum is mixed-clock and not conserved. Costs one MPI_Allreduce per synchronization. The per-phase breakdown that once accompanied this (kick1/drift/grav/hydro/mergesplit/kick2/spawn/srcterms, one Allreduce per mark) was hoisted out to scripts/energy_budget_phases.patch -- apply that to localize an error to a specific stage of the step.
#  MERGE_SPLIT_ENERGY_DIAGNOSTIC was hoisted out to scripts/merge_energy_diagnostic.patch -- apply that to get it back; setting the flag alone does nothing. It reports the cumulative energy budget of merge operations: dE (change in the true InternalEnergy+KE), dE_pred (same for the predicted variables the hydro actually integrates), KE_com (what a momentum-conserving merge must dissipate), resid (the correction merge_particles_ij computes), how often the InternalEnergyPred floor fires, and the TimeBin/Ti_current relationship of merged pairs. The patch header records what it established, including the dead ends worth not repeating.
#WAKEUP_TRUNCATE_STEP_ON_DEMOTION   # AUTO-ENABLED whenever the spawn module (SINK_WIND_SPAWN, itself implied by SINGLE_STAR_FB_JETS/WINDS/SNE) is on; set it explicitly only to force it without spawning, and see WAKEUP_REVERSE_KICK_ON_DEMOTION to opt out. When the WAKEUP timestep limiter demotes a particle mid-step, truncate its current step to end at the system time instead of reversing the kick it already received. GIZMO's default calls do_the_kick with reversed time bounds to un-do part of the applied kick, but that re-derives the increment from state which has changed since, so it does not cancel and instead injects energy -- of order 20% of the injected energy in the adiabatic wind_singlestar test and a few percent of E_SN in SN_singlestar, in both cases from essentially all hydro-triggered wakeups. Saitoh & Makino (2009), whose limiter this is, explicitly avoid "the backward integration of the entire system" (their eq 3) and instead set the new time consistent with the system time with the step as the difference -- which is what this does. Eliminates the error (reversal dE -> 0 to roundoff) and is inert in sedov/noh/interactblast/evrard, where the reversal never fires. CAVEAT: dt_step is then not a power-of-two bin length, so it disagrees with TimeBin; code deriving a step from TimeBin rather than integertime_step() will over-report it.
#WAKEUP_REVERSE_KICK_ON_DEMOTION    # opt back out of the WAKEUP_TRUNCATE_STEP_ON_DEMOTION auto-enable above and restore the legacy kick reversal, which is not energy-conserving. For reproducing older spawning runs, bisecting against them, and regression-testing both paths -- not recommended for production.
#MERGE_SPLIT_CONSERVE_ENERGY        # make particle merging energy-conserving. A merge sets the new velocity momentum-conservingly, so the pair's COM-frame kinetic energy (1/2)*mu*|dv|^2 leaves the kinetic budget; this returns it as heat. The residual has always been computed in merge_particles_ij but the line applying it had been commented out since the earliest revision in this repo (marked '//???'), so historically every merge silently destroyed that energy -- ~23% of the injected energy in the adiabatic wind_singlestar test, where E_kin+E_th=L*t is exact. AUTO-ENABLED with SINK_WIND_SPAWN (see precompiler_logic.h) alongside WAKEUP_TRUNCATE_STEP_ON_DEMOTION, since the two errors have opposite sign and fixing either alone unmasks the other; MERGE_SPLIT_DISCARD_ENERGY opts back out. With both active that test conserves to a few percent; the remainder is merge-associated but is not a bookkeeping error.
#MERGE_SPLIT_DISCARD_ENERGY         # opt back out of MERGE_SPLIT_CONSERVE_ENERGY where it is auto-enabled, restoring the legacy behaviour of discarding the COM-frame kinetic energy at every merge. For reproducing older results only -- it reintroduces a ~23% energy sink in spawning runs.
#MERGE_SPLIT_LIMIT_KINETIC_DISSIPATION   # never merge a gas pair if that would convert the majority of the energy involved from kinetic into heat. Merging thermalizes KE_com=(1/2)*mu*|dv|^2, which cooling then radiates away, making such a merger an artificial unresolved shock; vetoing it lets the solver follow the shock instead. Compared against the thermal+magnetic energy of the lighter cell (the one destroyed), which keeps the test mass-ratio-neutral. Threshold: MERGE_SPLIT_MAX_KINETIC_DISSIPATION_FRACTION.
#MERGE_SPLIT_MAX_KINETIC_DISSIPATION_FRACTION=(0.5)  # max allowed KE_com/(KE_com+E_light) under MERGE_SPLIT_LIMIT_KINETIC_DISSIPATION. 0.5 = "majority is kinetic"; lower is more conservative (~0.34 = Mach 1 on the lighter cell).
#MERGE_SPLIT_ALLOW_KINETIC_DISSIPATION   # opt back out of the MERGE_SPLIT_LIMIT_KINETIC_DISSIPATION auto-enable, permitting merges that thermalize the majority of the energy involved. For reproducing older results only.
# --------------------
# ----- Radiation-Hydrodynamics Special Options for Test Problems + Disabled or Other Special Features
#RT_DISABLE_UV_BACKGROUND          # disable extenal UV background in cooling functions (to isolate pure effects of local RT, or if simulating the background directly)
#RT_SEPARATELY_TRACK_LUMPOS        # keep luminosity vs. mass positions separate in tree. not compatible with Tree-PM mode, but it can be slightly more accurate and useful for testing in tree-only mode with LEBRON or OTVET algorithms.
#RT_HYDROGEN_GAS_ONLY              # sets hydrogen fraction to 1.0 (used for certain idealized chemistry calculations)
#RT_TIMESTEP_LIMIT_RECOMBINATION   # limit timesteps to the explicit recombination time when transporting ionizing photons. note our chemistry solvers are all implicit and can handle larger timesteps, but no gaurantee of transport accuracy for much larger steps since opacities depend on ionization states.
#RT_ENHANCED_NUMERICAL_DIFFUSION   # option which increases numerical diffusion, to get smoother solutions (akin to using HLL instead of HLLC+E fluxes), if desired; akin to slopelimiters~0 model
#RT_COMPGRAD_EDDINGTON_TENSOR      # forces computation of eddington tensor even when not needed by the code
#RT_REINJECT_ACCRETED_PHOTONS      # when sink particles are used, photons lost when a gas cell is accreted are reinjected into the lowest-energy frequency bin on the following photon injection from that sink
# --------------------
# ----- Sink particle/sink particle special options
#SINK_WIND_SPAWN_SET_BFIELD_POLTOR  # set poloridal and toroidal magnetic field for spawn particles (should work for all particle spawning). Cite Su et al., arXiv:2102.02206, for methods.
#SINK_WIND_SPAWN_SET_JET_PRECESSION # manually set precession in parameter file (does not work for cosmological simulations).  Cite Su et al., arXiv:2102.02206, for methods.
#JET_DIRECTION_FIXED_Z          # TESTING ONLY: pin the collimated spawn axis to global +/-z, skipping SINK_FOLLOW_ACCRETED_ANGMOM reorientation and SINK_WIND_SPAWN_SET_JET_PRECESSION, so launch geometry can be validated against a known axis. Not physical for production.
#SINK_SPAWN_NO_MERGE            # never merge spawned jet/wind cells, in either direction, holding the spawned population fixed (e.g. to compare merged vs unmerged outflows). Diagnostic: it removes retirement entirely, so cell counts grow without bound. cf. SINK_SPAWN_MERGE_WHEN_AMBIENT.
#SINK_SPAWN_MERGE_WHEN_AMBIENT   # retire spawned jet/wind cells only once they have joined the ambient medium, rather than as soon as their relative velocity dips below one neighbor's sound speed. Requires (a) the target to be a normal cell, never another spawned cell (restoring the restriction from 9ce8be04, disabled for STARFORGE by 22c755df), and (b) the cell to be subsonic wrt the mass-weighted mean velocity of non-spawned gas in its kernel; no ambient gas in kernel = still in the outflow. Pair with MERGE_SPLIT_LIMIT_KINETIC_DISSIPATION for the energetic half.
#SINK_SPAWN_MERGE_MAX_THERMAL_CONTRAST=(1.0)  # under SINK_SPAWN_MERGE_WHEN_AMBIENT, max |ln(u_cell/u_ambient)| in specific internal energy for a spawned cell to count as having joined the ISM, u_ambient being the mass-weighted geometric mean over non-spawned gas in its kernel. Subsonic alone is not sufficient: a contact discontinuity is in pressure equilibrium with no velocity jump, so outflow material at the contact passes every kinematic test while still being a distinct phase -- in the adiabatic wind_singlestar test every otherwise-eligible cell sits at |dln u| > 7. Deliberately u and not an entropic function: A=P/rho^gamma is a state function only for a single-gamma ideal gas, so with gamma varying by composition and temperature (EOS_GAMMA_VARIABLE, COOL_MOLECFRAC_NONEQM) equal A does not mean equal entropy and the test would be ill-founded, not merely imprecise. u needs no EOS assumption and is what the merge actually averages. Lower retires later. The flag reports the |dln u| histogram for calibration on a problem where wind genuinely mixes; that regime is empty in wind_singlestar, whose exact solution has no mixing, so any value up to ~7 behaves identically there.
#SINK_SPAWN_MERGE_ANY_NEIGHBOR   # opt back out of the SINK_SPAWN_MERGE_WHEN_AMBIENT auto-enable, restoring the legacy behaviour of retiring a spawned cell as soon as it is subsonic wrt whichever single neighbour is picked as its merge target. For reproducing older jet/wind runs; it lets outflow material sitting at the contact discontinuity be absorbed into ambient gas.
#SINK_SCALE_SPAWNINGMASS_WITH_INITIALMASS # rescale the cell spawning mass criterion to scale with the initial sink mass for any sink and sink feedback model (instead of setting the spawning mass to a fixed universal constant in code units).
#SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA # modifies SINK_SEED_FROM_LOCALGAS to require that the local acceleration scale (including all mass/gravity sources) exceeds the critical value defined in arXiv:2103.10444. cite that paper for implementation
#SINK_INTERACT_ON_GAS_TIMESTEP      # force sink particles to be active in the timestep hierarchy on a timestep no larger than the minimum timestep of any gas cell inside the sink particle interaction/accretion/neighbor kernel
#SINK_GRAVCAPTURE_FIXEDSINKRADIUS   # uses a fixed sink radius for gravitational capture/accretion onto sink particles, taken to be the force softening kernel radius of the sink particle
#SINGLE_STAR_FB_TIMESTEPLIMIT       # limit the timesteps in gas cells potentially interacting with sink particles in the SINGLE_STAR modules, so that they cannot take timesteps large compared to quantities like the time for jets/winds to cross from the sink to the gas cell
# --------------------
# ----- Cosmic ray special options
#CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART=1  # allows restart from a snapshot (flag=2) where single-bin CR model was used, for runs with CR spectra: the spectra are populated with the energy of the single-bin snapshot and fixed initial spectral shapes/ratios
#CRFLUID_ALT_DISABLE_LOSSES         # turn off CR heating/cooling interactions with gas (catastrophic losses, hadronic interactions, etc; only adiabatic/D:GradU work terms remain)
#CRFLUID_ALT_REACCEL_ONLY_DIFFUSIVE # replaces the correct form of the re-acceleration terms calculated from the focused CR transport equation with the more ad-hoc diffusive re-acceleration assumption that CR scattering/diffusion is dominated by an undamped, perfectly-symmetric, isotropic extrinsic turbulent cascade (not likely valid below ~TeV), following e.g. Drury+Strong 2017A&A...597A.117D; requires CRFLUID_EVOLVE_SPECTRUM.
#CRFLUID_ALT_VARIABLE_RSOL          # allows a variable (CR energy-dependent) reduced speed of light to be used for CRs, which is set in the function return_CRbin_M1speed defined by the user. cite Hopkins et al. 2021, arXiv:2103.10443
# --------------------
# ----- FIRE and STARFORGE sub-module special options
#GALSF_FB_FIRE_RPROCESS             # optional module to inject some tracers representing R-process or S-process or any other species into the system, with an arbitrary yield or rate model. designed to be customized
#FIRE_SNE_ENERGY_METAL_DEPENDENCE_EXPERIMENT=0 # experiment with modified supernova energies as a function of metallicity - freely modify this module as desired. used for numerical experiments only.
#SINGLE_STAR_AND_SSP_HYBRID_MODEL=1 # cells with mass less than this (in solar) are treated with the single-stellar evolution models, larger mass with ssp models. needs user to specify a refinement criterion, and a criterion for when one module or another will be used. still in testing.
#SINGLE_STAR_AND_SSP_HYBRID_MODEL_DEFAULTS=1 # uses default settings for SINGLE_STAR_AND_SSP_HYBRID_MODEL=SINGLE_STAR_AND_SSP_HYBRID_MODEL_DEFAULTS (plus lots of other flag settings) from zoom-in experiments around special particles
#STARFORGE_GMC_TURBINIT=1           # special flag for custom ICs behavior in star formation simulations for turbulent clouds. adds an analytic uniform sphere harmonic potential + r^-3 halo outside to confine stirred turbulent gas, during the 'stirring' phase. cite Lane et al., 2022MNRAS.510.4767L
#STARFORGE_FILAMENT_TURBINIT        # special flag for custom ICs behavior in star formation simulations for turbulent clouds. adds an analytic potential of an finitite cylinder with a Plummer density profile, truncated at the ends of the cylinder, to confine stirred turbulent gas, during the 'stirring' phase. cite Lane et al., 2022MNRAS.510.4767L
#STARFORGE_FEEDBACK_TRACERS=3       # adds tracer fields to recycled gas in the SINGLE_STAR modules to explicitly track gas coming from jets versus main-sequence winds versus supernovae. tracer fields added to metal fields; 0 for jets, 1 for winds, 2 for SNe, 3 for all.
# --------------------
# ----- Dust grain/particulate/aerosol module special options
#GRAIN_RDI_TESTPROBLEM             # top-level flag to enable a variety of test problem behaviors, customized for the idealized studies of dust dynamics in Moseley et al 2019MNRAS.489..325M, Seligman et al 2019MNRAS.485.3991S, Steinwandel et al arXiv:2111.09335, Ji et al arXiv:2112.00752, Hopkins et al 2020MNRAS.496.2123H and arXiv:2107.04608, Squire et al 2022MNRAS.510..110S. Cite these if used.
#GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE    # Make the idealized external grain acceleration grain size-dependent; equivalent to assuming an absorption efficiency Q=1. Cite GRAIN_RDI_TESTPROBLEM papers.
#GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION # Enables idealized radiation injection by a source population designed to set up outflow test problems with live radiation-hydrodynamics as in Hopkins et al., arXiv:2107.04608. Cite that paper if this module is used.
#IO_DUST_NOT_IN_ICFILE              # special flag needed if restarting from a snapshot (flag=2) with no dust. Will set dust species via params arguments and assume an MRN size distribution, if GALSF_ISMDUSTCHEM_MODEL is active 
# --------------------
# ----- MPI & Parallel-FFTW De-Bugging
#DOUBLEPRECISION_FFTW               # FFTW in double precision to match libraries
#DISABLE_ALIGNED_ALLOC              # disable calls to 'aligned_alloc', needed for older C99-only versions of GCC compilers [everything C11+ -should- be compatible and not need this]
# --------------------
# ----- Load-Balancing
#ALLOW_IMBALANCED_GASPARTICLELOAD   # increases All.MaxPartGas to All.MaxPart: can allow better load-balancing in some cases, but uses more memory. But use me if you run into errors where it can't fit the domain (where you would increase PartAllocFac, but can't for some reason)
####################################################################################################




####################################################################################################-
####################################################################################################-
##-
##- LEGACY CODE & PARTIALLY-IMPLEMENTED FEATURES: BEWARE EVERYTHING BELOW THIS LINE !
##-
####################################################################################################-
####################################################################################################-

####################################################################################################-
#-------------------------------- misc dev flags needing to be incorporated here (in testing)
#GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF=(2.0) #- star particle formation of O-stars is spread over this multiple of the free-fall time; requires GALSF_SFR_IMF_SAMPLING. developed by PFH
#GALSF_MERGER_STARCLUSTER_PARTICLES         #- module which merges star particles together meeting certain core-collapse conditions, so they can be cosmologically evolved. developed by PFH
#SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM           #- module for special nuclear zoom-in simulations. currently entirely custom behavior, not designed for wide use.
#CBE_INTEGRATOR_WITHGRADIENTS               #- module being developed. early stage
#CRFLUID_INJECTION_AT_SHOCKS=(0.1)          #- inject CRs (using standard spectra) at resolved shocks throughout simulation. value=maximum fraction of shock energy flux to convert to accelerated CRs. currently using hard threshold of >Mach 5, >100 km/s shock velocity, along with spurious shock detection. coded by PFH, in testing
#SINK_CR_INJECTION_AT_TERMINATION=(0.25)    #- inject CRs (requires SINK_WIND_SPAWN and SINK_COSMIC_RAYS) approximately at spawned-cell termination shocks. developed by Kung-Yi Su, this version testing by PFH, currently uses very simple deceleration to fraction of launch velocity (=value set here) to determine when to inject
#SINK_TEST_WIND_MIXED_FASTSLOW=(1.e5)       #- BH has a slow outflow and fast jet, where fast jet has 1/50th the mass-loading and this speed in kms by default
#USE_TIMESTEP_DILATION_FOR_ZOOMS            #- enable time dilation modules, need to customize for applications, cannot be simply generically turned on without coding how they will work
#DILATION_FOR_STELLAR_KINEMATICS_ONLY       #- special version of time dilation designed for stellar kinematics in e.g. dense star clusters or galaxy centers
#SINK_RIAF_SUBEDDINGTON_MODEL=(0.01)        #- enable an arbitrary modular variation in the radiative efficiency of BHs as a function of eddington ratio or other particle properties, with the critical transition to the jet mode at this eddington ratio (defined in terms of mdot/mdot_crit)
#OUTPUT_SHOCK_MACH_NUMBER                   #- compute and output the shock Mach number for each gas cell, using the information in the Riemann problem and reconstruction, plus additional converging flow and spurious compression checks.
####################################################################################################-

############################################################################################################################-
#------------------ ISM Dust Chemical Evolution Models (follow growth, destruction, and size evolution of different grain types)
#----------------- Users of any of these modules should cite Choban et al., 2022/25 for the methods/implementation in GIZMO and FIRE
############################################################################################################################-
#GALSF_ISMDUSTCHEM_MODEL=(2+4+8)    #- enable live dust evolution model (must select either elemental or species or other model codes as well)
                                    #- model = 1: "dust by element" dust evolution model based off Bekki(2013)/McKinnon+(2016). Track generalized silicates and carbonaceous dust.
                                    #- model = 2: "dust by species" dust evolution model based off Zhukovska+(2008/2016/2018). Tracks silicates (set composition), carbonaceous, SiC, and metallic iron dust along with optional iron nanoparticles and/or O reservoir dust species.
                                    #- model = 4: additional metallic iron dust nano-particles with set fraction assumed to be locked in silicate dust as inclusions based on Zhukovska+(2018)
                                    #- model = 8: additional oxygen bearing dust species which is a simple match to observations of MW oxygen depletion since they cannot be explained with purely silicate dust
                                    #- model = 16: modified "dust by species" model to be used with GALSF_ISMDUSTCHEM_GRAINSIZEEVO. Tracks size evolution of silicates (set composition), carbonaceous, and metallic iron
                                    #- model = 32: modified "dust by species" model to be used with GALSF_ISMDUSTCHEM_GRAINSIZEEVO. Tracks size evolution of silicates w/ extra O and Fe (set composition) and carbonaceous
#GALSF_ISMDUSTCHEM_GRAINSIZEEVO=16  #- enable grain size evolution model w/ N number of logarithmically spaced bins (must also turn on GALSF_ISMDUSTCHEM_MODEL= 2 + (16 or 32) only)
############################################################################################################################-

#DISABLE_HERMITE_INTEGRATION      # turn off the 4th-order Hermite integration of sinks, leaving them on plain KDK leapfrog. Diagnostic: leapfrog's order is known (2nd), so it is a control for convergence measurements. See precompiler_logic.h, where it negates the HERMITE_INTEGRATION default
