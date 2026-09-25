/* This structure holds all the information that is stored for each particle of the simulation. */

extern ALIGN(32) struct particle_data
{
    short int Type;                 /*!< flags particle type.  0=gas, 1=halo/high-res dm, 2=alt dm/disk/collisionless, 3=pic/dust/bulge/alt dm, 4=new stars, 5=sink */
    short int TimeBin;
    MyIDType ID;                    /*! < unique ID of particle (assigned at beginning of the simulation) */
    MyIDType ID_child_number;       /*! < child number for particles 'split' from main (retain ID, get new child number) */
#ifndef SINK_WIND_SPAWN
    int ID_generation;              /*! < generation (need to track for particle-splitting to ensure each 'child' gets a unique child number */
#else
    MyIDType ID_generation;
#endif
    
    integertime Ti_begstep;         /*!< marks start of current timestep of particle on integer timeline */
    integertime Ti_current;         /*!< current time of the particle */
    
    ALIGN(32) Vec3<MyDouble> Pos;   /*!< particle position at its current time */
    MyDouble Mass;                  /*!< particle mass */

    Vec3<MyDouble> Vel;             /*!< particle velocity at its current time */
    Vec3<MyDouble> dp;
    MyFloat Particle_DivVel;        /*!< velocity divergence of neighbors (for predict step) */
    
    Vec3<MyDouble> GravAccel;       /*!< particle acceleration due to gravity */
#ifdef PMGRID
    Vec3<MyFloat> GravPM;           /*!< particle acceleration due to long-range PM gravity force */
#endif
    MyFloat OldAcc;                    /*!< magnitude of old gravitational force. Used in relative opening criterion */
#ifdef SPECIAL_POINT_MOTION
    Vec3<MyFloat> Acc_Total_PrevStep;  /*!< old total acceleration on a given cell/particle */
#endif
#ifdef HERMITE_INTEGRATION
    Vec3<MyFloat> Hermite_OldAcc;
    Vec3<MyFloat> OldPos;
    Vec3<MyFloat> OldVel;
    Vec3<MyFloat> OldJerk;
#ifdef HERMITE_SYMMETRIC_TIMESTEP_ACTIVE
    MyDouble dt_crit_last;  /*!< step criterion (code time, all factors and caps, not quantised) from the latest get_timestep call; 0 if not Hermite-eligible then */
    MyDouble dt_crit_prev;  /*!< dt_crit_last as it stood when the current step started */
#endif
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    MyFloat TreeMass;  /*!< Mass seen by the particle as it sums up the gravitational force from the tree - should be equal to total mass, a useful debug diagnostic  */
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL)
    MyFloat Potential;        /*!< gravitational potential */
#if defined(EVALPOTENTIAL) && defined(PMGRID)
    MyFloat PM_Potential;
#endif
#endif
#if defined(GALSF_SFR_TIDAL_HILL_CRITERION) || defined(TIDAL_TIMESTEP_CRITERION) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(OUTPUT_TIDAL_TENSOR) || (defined(SINGLE_STAR_TIMESTEPPING) && (SINGLE_STAR_TIMESTEPPING > 0)) || defined(ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)
#define COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    SymmetricTensor2<MyFloat> tidal_tensorps;            /*!< tidal tensor (=second derivatives of grav. potential) */
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    MyFloat tidal_tensor_mag_prev;                       /*!< saved frobenius norm of the tidal tensor, from the previous timestep >*/
    SymmetricTensor2<MyFloat> tidal_tensorps_prevstep;   /*!< save the entire tensor if this is active >*/
    MyFloat tidal_zeta;                                  /*!< also need to calculate an analog of the ags zeta variable here >*/
#endif
#ifdef PMGRID
    SymmetricTensor2<MyFloat> tidal_tensorpsPM;          /*!< for TreePM simulations, long range tidal field */
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
    Vec3<double> GravJerk;
#endif
    
#ifdef GALSF
    MyFloat StellarAge;        /*!< formation time of star particle */
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
#ifdef GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF
    MyFloat TimeDistribOfStarFormation; /*!< free-fall time at the moment of star formation, which defines for this particle the delay distribution for forming the relevant O-stars */
    MyFloat IMF_WeightedMeanStellarFormationTime; /*!< weighted mean stellar formation time, to use instead of the normal stellarage parameter on-the-fly */
#endif
#endif
    
    MyFloat KernelRadius;           /*!< search radius around particle for neighbors/interactions */
    MyFloat NumNgb;                 /*!< neighbor number around particle */
    MyFloat DrkernNgbFactor;        /*!< correction factor needed for varying kernel lengths */
#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    MyFloat DensityAroundParticle;         /*!< gas density in the neighborhood of the collisionless particle (evaluated from neighbors) */
#endif
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES) || defined(COOLING)
    Vec3<MyFloat> GradRho;          /*!< gas density gradient evaluated simply from the neighboring particles, for collisionless centers */
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
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
    MyFloat MassReturn_ThisTimeStep; /* gas return from stellar winds */
#ifdef GALSF_FB_FIRE_RPROCESS
    MyFloat RProcessEvent_ThisTimeStep; /* R-process event tracker */
#endif
#ifdef GALSF_FB_FIRE_AGE_TRACERS
    MyFloat AgeDeposition_ThisTimeStep; /* age-tracer deposition */
#endif
#endif
#endif
#ifdef GALSF_FB_MECHANICAL
#define AREA_WEIGHTED_SUM_ELEMENTS 12 /* number of weights needed for full momentum-and-energy conserving system */
    MyFloat Area_weighted_sum[AREA_WEIGHTED_SUM_ELEMENTS]; /* normalized weights for particles in kernel weighted by area, not mass */
#endif
#ifdef GALSF_FB_FIRE_RT_LOCALRP
    MyFloat NewStar_Momentum_For_JetFeedback; /* amount of momentum to return from protostellar jet sub-grid model */
#endif
    
#if defined(GRAIN_FLUID)
    MyFloat Grain_Size;
    MyFloat Gas_Density;
    MyFloat Gas_InternalEnergy;
    Vec3<MyFloat> Gas_Velocity;
    MyFloat Grain_AccelTimeMin;
#if defined(GRAIN_BACKREACTION)
    Vec3<MyFloat> Grain_DeltaMomentum;
#endif
#if defined(GRAIN_LORENTZFORCE)
    Vec3<MyFloat> Gas_B;
#endif
#endif
#if defined(PIC_MHD)
    short int MHD_PIC_SubType;
#endif
    
#if defined(SINK_PARTICLES)
    MyIDType SwallowID;
    int IndexMapToTempStruc;   /*!< allows for mapping to SinkTempInfo struc */
#ifdef SINK_WIND_SPAWN
    MyFloat unspawned_wind_mass;    /*!< tabulates the wind mass which has not yet been spawned */
#ifdef SINGLE_STAR_FB_JETS
    MyFloat unspawned_jet_mass;    /*!< separate reservoir for jet (accretion) mass not yet spawned; jets and MS winds bank into their own reservoirs and only whichever currently dominates (P.wind_mode) actually drains via discrete spawning, see sink.cc */
#endif
#endif
#ifdef SINK_COUNTPROGS
    int Sink_CountProgs;
#endif
    MyFloat Sink_Mass;
    MyFloat Sink_Formation_Mass; /* initial mass of sink (total particle) when it formed */
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
    MyFloat Sink_Mdot_ROI;
    MyFloat Sink_ROI;
#endif
#if defined(SINK_GRAVCAPTURE_FIXEDSINKRADIUS)
    MyFloat SinkRadius;
#endif
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
    MyFloat dt_since_last_gas_search; /* keep track of time since the sink's last neighbor search and gas interaction (for feedback/accretion) */
    short int do_gas_search_this_timestep; /* flag for deciding whether to do gas stuff for a given timestep */
#endif
#ifdef GRAIN_FLUID
    MyFloat Sink_Dust_Mass;
#endif
#ifdef RT_REINJECT_ACCRETED_PHOTONS
    MyFloat Sink_accreted_photon_energy;
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
    MyFloat SwallowTime; /* freefall time of a particle onto a sink particle  */
#endif
#if defined(SINGLE_STAR_TIMESTEPPING)
    MyFloat Sink_SurroundingGasVel; /* Relative speed of sink to surrounding gas  */
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
    int Sink_Ngb_Flag; /* whether or not the gas lives in a sink's hydro stencil */
#endif
#ifdef SINK_ALPHADISK_ACCRETION
    MyFloat Sink_Mass_Reservoir;
#endif
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
    MyFloat Sink_AccretionDeficit; /* difference between continuously-accreted and discretely-accreted masses, needs to be evolved to ensure exact conservation with some modules */
#endif
#ifdef SINK_FOLLOW_ACCRETED_ANGMOM
    Vec3<MyFloat> Sink_Specific_AngMom;
#endif
#ifdef SINK_RETURN_BFLUX
    Vec3<MyDouble> B;
#endif
#ifdef JET_DIRECTION_FROM_KERNEL_AND_SINK
    MyFloat Mgas_in_Kernel;
    Vec3<MyFloat> Jgas_in_Kernel;
#endif
    MyFloat Sink_Mdot;
    int Sink_TimeBinGasNeighbor;
#if defined(SINGLE_STAR_TIMESTEPPING)
    MyFloat Sink_dr_to_NearestGasNeighbor;
#endif
#ifdef SINK_REPOSITION_ON_POTMIN
    Vec3<MyFloat> Sink_PotentialMinimumOfNeighborsPos;
    MyFloat Sink_PotentialMinimumOfNeighbors;
#endif
#endif  /* if defined(SINK_PARTICLES) */
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    MyFloat MencInRcrit;
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
#if defined(SINGLE_STAR_FIND_BINARIES) || (SINGLE_STAR_TIMESTEPPING > 0)
    MyDouble Min_Sink_OrbitalTime; //orbital time for binary
    Vec3<MyDouble> comp_dx; //position offset of binary companion - this will be evolved in the Kepler solution while we use the Pos attribute to track the binary COM
    Vec3<MyDouble> comp_dv; //velocity offset of binary companion - this will be evolved in the Kepler solution while we use the Vel attribute to track the binary COM velocity
    MyDouble comp_Mass; //mass of binary companion
    int is_in_a_binary; // flag whether star is in a binary or not
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    MyFloat Min_Sink_Freefall_time;
    MyFloat Min_Sink_Approach_Time;
#if (SINGLE_STAR_TIMESTEPPING > 0)
    int SuperTimestepFlag; // >=2 if allowed to super-timestep (increases with each drift/kick), 1 if a candidate for super-timestepping, 0 otherwise
    MyDouble COM_dt_tidal; //timescale from tidal tensor evaluated at the center of mass without contribution from the companion
    Vec3<MyDouble> COM_GravAccel; //gravitational acceleration evaluated at the center of mass without contribution from the companion
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    MyFloat MaxFeedbackVel; // maximum signal velocity of any feedback mechanism emanating from the star
    MyFloat Min_Sink_FeedbackTime;  // minimum time for feedback to arrive from a star
#endif
#endif
#endif
    
    
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
    MyFloat ProtoStellarAge; /*!< record the proto-stellar age instead of age */
    MyFloat ProtoStellarRadius_inSolar; /*!< protostellar radius (also tracks evolution from protostar to ZAMS star) */
    int ProtoStellarStage; /* Track the stage of protostellar evolution, 0: pre collapse, 1: no burning, 2: fixed Tc burning, 3: variable Tc burning, 4: shell burning, 5: main sequence, 6: supernova, see Offner 2009 Appendix B*/ //IO flag IO_STAGE_PROTOSTAR
    MyFloat Mass_D; /* Mass of gas in the protostar that still contains D to burn */ // IO flag IO_MASS_D_PROTOSTAR
    MyFloat StarLuminosity_Solar; /* the total luminosity of the star in L_solar units*/ //IO flag IO_LUM_SINGLESTAR
    MyFloat ZAMS_Mass; /* The mass the star has when reaching the main sequence */ //IO flag IO_ZAMS_MASS
#ifdef SINGLE_STAR_FB_WINDS
    MyFloat Wind_direction[6]; // direction of wind launches, to reduce anisotropy launches go along a random axis then a random perpendicular one, then one perpendicular to both.
    int wind_mode; // tells what kind of wind model to use, 1 for particle spawning and 2 for using the FIRE wind module
    double wind_mode_time; /* All.Time of the last wind mode change, read by the spawning hysteresis. Double because it stores
                              All.Time: at large code times a float's quantum exceeds the interval between evaluations, the stored
                              value stops moving, and the hysteresis clock reads zero elapsed forever -- pinning the mode. */
#endif
#ifdef  SINGLE_STAR_FB_SNE
    MyFloat Mass_final; //final mass of the star before going SN (Since this is not saved to snapshots, hard restarts in the middle of spawning an SN will do weird things)
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
    MyFloat DM_KernelRadius;
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
    
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    MyFloat Time_Of_Last_MergeSplit;
#endif
    
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    MyFloat Time_Of_Last_SmoothedVelUpdate;
#endif
    
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
    MyFloat AGS_zeta;               /*!< correction term for adaptive gravitational softening lengths */
#endif
    
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    MyDouble AGS_KernelRadius;          /*!< smoothing length (for gravitational forces) */
    MyDouble AGS_vsig;          /*!< signal velocity of particle approach, to properly time-step */
#endif
    
#if defined(WAKEUP)
    short int wakeup;                     /*!< flag to wake up particle */
#endif
    
#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
    MyFloat StarParticleEffectiveSize;   /*!< effective 'size' of a star particle at formation */
#endif
    
#ifdef DM_FUZZY
    MyFloat AGS_Density;                /*!< density calculated corresponding to AGS routine (over interacting DM neighbors) */
    Vec3<MyFloat> AGS_Gradients_Density;   /*!< density gradient calculated corresponding to AGS routine (over interacting DM neighbors) */
    Mat3<MyFloat> AGS_Gradients2_Density;   /*!< density gradient calculated corresponding to AGS routine (over interacting DM neighbors) */
    MyFloat AGS_Numerical_QuantumPotential; /*!< additional potential terms 'generated' by un-resolved compression [numerical diffusivity] */
    MyFloat AGS_Dt_Numerical_QuantumPotential; /*!< time derivative of the above */
#if (DM_FUZZY > 0)
    MyFloat AGS_Psi_Re;
    MyFloat AGS_Psi_Re_Pred;
    MyFloat AGS_Dt_Psi_Re;
    Vec3<MyFloat> AGS_Gradients_Psi_Re;
    Mat3<MyFloat> AGS_Gradients2_Psi_Re;
    MyFloat AGS_Psi_Im;
    MyFloat AGS_Psi_Im_Pred;
    MyFloat AGS_Dt_Psi_Im;
    Vec3<MyFloat> AGS_Gradients_Psi_Im;
    Mat3<MyFloat> AGS_Gradients2_Psi_Im;
    MyFloat AGS_Dt_Psi_Mass;
#endif
#endif
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    Mat3<MyDouble> NV_T;                                           /*!< holds the tensor used for gradient estimation */
#endif
#ifdef CBE_INTEGRATOR
    double CBE_basis_moments[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];         /* moments per basis function */
    double CBE_basis_moments_dt[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];      /* time-derivative of moments per basis function */
#ifdef CBE_INTEGRATOR_WITHGRADIENTS
    CBE_basis_moments_Gradients[CBE_INTEGRATOR_NBASIS][3]; /* gradients of the scalar weight of each basis function */
#endif
#endif

    /* member functions */
    inline integertime integertime_step() const { /*!< integer timestep for this particle */
#ifndef WAKEUP
        return GET_INTEGERTIME_FROM_TIMEBIN(TimeBin);
#else
        return dt_step;
#endif
    }

    inline double Get_Particle_Size() const { /*!< effective particle/cell size from kernel radius and neighbor number */
#if (NUMDIMS == 1)
        return 2.00000 * KernelRadius / NumNgb; // (2)^(1/1)
#elif (NUMDIMS == 2)
        return 1.77245 * KernelRadius / NumNgb; // (pi)^(1/2)
#else
        return 1.61199 * KernelRadius / NumNgb; // (4pi/3)^(1/3)
#endif
    }
}
*P,                /*!< holds particle data on local processor */
*DomainPartBuf;        /*!< buffer for particle data used in domain decomposition */

