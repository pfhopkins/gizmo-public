#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"


/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */
/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified heavily
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO: initializing
 * new/modified variables, rewriting order-of-operations, standardizing
 * where some calls occur, and adding new subroutines that need to be called
 * on startup now.
 */

/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If gas cells are present,
 *  the initial gas kernel lengths are determined.
 */
void init(void)
{
    int i, j; double a3, atime, a2_fac;

#ifdef MAGNETIC
    double gauss2gizmo = All.UnitMagneticField_in_gauss / UNIT_B_IN_GAUSS;
    /* NOTE: we will always work -internally- in code units where MU_0 = 1; hence the 4pi here; [much simpler, but be sure of your conversions!] */
#endif

#ifdef SINK_PARTICLES
    int count_holes = 0;
#endif

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();

    if(RestartFlag != 1) {All.MinMassForParticleMerger = 0; All.MaxMassForParticleSplit = 0;}

    if(RestartFlag == 3 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");}
        endrun(0);
    }

    if(RestartFlag == 4 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if snapshot should be converted\n");}
        endrun(0);
    }

    if(RestartFlag == 5 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if power spectrum and two-point correlation function should be calculated\n");}
        endrun(0);
    }

    if(RestartFlag == 6 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if velocity power spectrum for the gas cells should be calculated\n");}
        endrun(0);
    }


    switch (All.ICFormat)
    {
        case 1:
        case 2:
        case 3:
        case 4:
            if(RestartFlag >= 2 && RestartSnapNum >= 0)
            {
                char fname[MAX_PATH_BUFFERSIZE_TOUSE];
                if(All.NumFilesPerSnapshot > 1) {snprintf(fname, MAX_PATH_BUFFERSIZE_TOUSE, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase, RestartSnapNum);}
                    else {snprintf(fname, MAX_PATH_BUFFERSIZE_TOUSE, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);}
                read_ic(fname);

            }
            else {read_ic(All.InitCondFile);}
            break;

        default:
            if(ThisTask == 0) {printf("ICFormat=%d not supported.\n", All.ICFormat);}
            endrun(0);
    }

#ifdef CHIMES_INITIALISE_IN_EQM
    for (i = 0; i < N_gas; i++) {allocate_gas_abundances_memory(&(ChimesGasVars[i]), &ChimesGlobalVars);}
#endif

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();


#if defined(COOLING) && !defined(CHIMES)
    IonizeParams();
#endif

    All.Ti_Current = 0;
    if(All.ComovingIntegrationOn)
    {
        All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
        a3 = All.Time * All.Time * All.Time; atime = All.Time; a2_fac = (All.Time * All.Time);
    }
    else
    {
        All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
        a3 = atime = a2_fac = 1;
    }

    set_softenings();

    All.NumCurrentTiStep = 0;	/* setup some counters */
    All.DomainFacWork = 0.5; All.DomainFacLoad = 0.5; All.DomainFacWorkGas = 0.0; /* seed adaptive load-balance weights (restarts restore them from All) */
    All.SnapshotFileCount = 0;
    if(RestartFlag == 2)
    {
        if(RestartSnapNum < 0)
        {
            char *underscore = strrchr(All.InitCondFile, '_');
            if(!underscore)
            {
                char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
                snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n", All.InitCondFile);
                terminate(buf);
            }
            else {All.SnapshotFileCount = atoi(underscore + 1) + 1;}
        }
        else {All.SnapshotFileCount = RestartSnapNum + 1;}
    }

#ifdef OUTPUT_LINEOFSIGHT
    All.Ti_nextlineofsight = (int) (log(All.TimeFirstLineOfSight / All.TimeBegin) / All.Timebase_interval);
    if(RestartFlag == 2) {endrun(78787);}
#endif

    All.TotNumOfForces = 0;
    All.TopNodeAllocFactor = 0.008; /* this will start from a low value and be iteratively increased until it is well-behaved */
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    All.TopNodeAllocFactor = 0.1; /* for optimization on startup this needs to be increased for these extreme dynamic range runs */
#endif
    All.TreeAllocFactor = 0.45; /* this will also iteratively increase to fit the particle distribution */
    /* To construct the BH-tree for N particles, somewhat less than N
     internal tree-nodes are necessary for ‘normal’ particle distributions.
     TreeAllocFactor sets the number of internal tree-nodes allocated in units of the particle number.
     By experience, space for ~0.65N internal nodes is usually fully sufficient for typical clustered
     particle distributions, so a value of 0.7 should put you on the safe side. If the employed particle
     number per processor is very small (less than a thousand or so), or if there are many particle pairs
     with identical or nearly identical coordinates, a higher value may be required. Since the number of
     particles on a given processor may be higher by a factor PartAllocFactor than the average particle
     number, the total amount of memory requested for the Barnes-Hut tree on a single processor scales proportional
     to PartAllocFactor*TreeAllocFactor. */

#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    for(i = 0; i < SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM; i++)
    {
        All.SpecialParticle_Position_ForRefinement[i] = {};
        All.Mass_Accreted_By_SpecialParticle[i]=0; All.Mass_of_SpecialParticle[i]=0;
    }
#endif
    
#ifdef GALSF_LIMIT_FBTIMESTEPS_FROM_BELOW
    if(RestartFlag != 1) {All.Dt_Since_LastFBCalc_Gyr=0; All.Dt_Min_Between_FBCalc_Gyr=((double)(GALSF_LIMIT_FBTIMESTEPS_FROM_BELOW))/1.e9;}
#endif

#ifdef BOX_PERIODIC
    if(All.ComovingIntegrationOn) {check_omega();}
#endif
    All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#if (defined(SINK_PARTICLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
    All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif

    for(i = 0; i < GRAVCOSTLEVELS; i++) {All.LevelToTimeBin[i] = 0;}

    for(i = 0; i < NumPart; i++) {for(j = 0; j < GRAVCOSTLEVELS; j++) {P[i].GravCost[j] = 0;}}

    if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
        {for(i=0;i<NumPart;i++) {P[i].Vel *= sqrt(All.Time)*All.Time;}}

#ifdef DM_SIDM
    init_self_interactions();
#endif

#ifdef METALS
    for(j=0;j<NUM_METAL_SPECIES;j++) {All.SolarAbundances[j]=0;} // initialize all to zero
    All.SolarAbundances[0]=0.02;        // all metals (by mass); present photospheric abundances from Asplund et al. 2009 (Z=0.0134, proto-solar=0.0142) in notes;
    //   also Anders+Grevesse 1989 (older, but hugely-cited compilation; their Z=0.0201, proto-solar=0.0213)
#ifdef COOL_METAL_LINES_BY_SPECIES
    All.SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314), with proto-solar Y=0.27
    All.SolarAbundances[2]=3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3); proto-solar from Asplund=8.47 -> 2.53e-3
    All.SolarAbundances[3]=1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3); PS=7.87->7.41e-4
    All.SolarAbundances[4]=8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3); PS=8.73->6.13e-3
    All.SolarAbundances[5]=2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3); PS=7.97->1.34e-3
    All.SolarAbundances[6]=9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4); PS=7.64->7.57e-4
    All.SolarAbundances[7]=1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4); PS=7.55->7.12e-4
    All.SolarAbundances[8]=6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4); PS=7.16->3.31e-4
    All.SolarAbundances[9]=1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4); PS=6.38->6.87e-5
    All.SolarAbundances[10]=1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3); PS=7.54->1.38e-3
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // new default abundances; using Asplund et al. 2009 proto-solar abundances
    All.SolarAbundances[0]=0.0142; if(NUM_METAL_SPECIES>=10) {
        All.SolarAbundances[1]=0.27030; All.SolarAbundances[2]=2.53e-3; All.SolarAbundances[3]=7.41e-4; All.SolarAbundances[4]=6.13e-3; All.SolarAbundances[5]=1.34e-3;
        All.SolarAbundances[6]=7.57e-4; All.SolarAbundances[7]=7.12e-4; All.SolarAbundances[8]=3.31e-4; All.SolarAbundances[9]=6.87e-5; All.SolarAbundances[10]=1.38e-3;}
#endif
#endif
#if defined(GALSF_ISMDUSTCHEM_MODEL)
    Initialize_ISMDustChem_Global_Variables();
#endif
#endif


    for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
        P[i].GravAccel = {};

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* init tidal tensor for first output (not used for calculation) */
        for(j=0;j<3;j++) {int kt; for(kt=0;kt<3;kt++) {P[i].tidal_tensorps[j][kt]=0;}}
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        P[i].tidal_tensor_mag_prev = 0; P[i].tidal_zeta=0; for(j=0;j<3;j++) {int kt; for(kt=0;kt<3;kt++) {P[i].tidal_tensorps_prevstep[j][kt]=0;}}
#endif
#ifdef PMGRID
        for(j=0;j<3;j++) {int kt; for(kt=0;kt<3;kt++) {P[i].tidal_tensorpsPM[j][kt]=0;}}
#endif
#endif

#ifdef ADAPTIVE_TREEFORCE_UPDATE
        P[i].time_since_last_treeforce = 0;
        P[i].tdyn_step_for_treeforce = 0;
#endif        
        

#ifdef PMGRID
        P[i].GravPM = {};
#endif
        P[i].Ti_begstep = 0;
        P[i].Ti_current = (integertime)0;
        P[i].TimeBin = 0;
        if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS) {P[i].OldAcc = 0;}	/* Do not zero in 2lpt case as masses are stored here */

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
        P[i].Potential = 0;
#endif
#ifdef GALSF
        if(RestartFlag == 0) {
#if defined(INPUT_READ_SINKPROPS) && defined(SINGLE_STAR_SINK_DYNAMICS)
            if(P[i].Type != 5) /* sink particles have their age read from the IC */
#endif
            {P[i].StellarAge = 0;}
        }
#ifdef GALSF_SFR_IMF_VARIATION
        if(RestartFlag == 0) {P[i].IMF_Mturnover = 2.0;} /* gives a solar-type IMF for our calculations in current code */
#endif
#ifdef GALSF_SFR_IMF_SAMPLING
        if(RestartFlag == 0) {P[i].IMF_NumMassiveStars = 0;}
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL) && defined(SINGLE_STAR_RESTART_FROM_FIRESIM)
        if(RestartFlag == 2) {P[i].IMF_NumMassiveStars = 0;}
#endif
#endif
#endif
        
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
        if(RestartFlag == 0 && P[i].ProtoStellarStage <= 0) {P[i].ProtoStellarStage = 0;} /* only reset if not already set from IC (via INPUT_READ_SINKPROPS) */
#endif

        if(RestartFlag != 1)
        {
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES)
            P[i].DensityAroundParticle = 0;
            P[i].GradRho[0]=0;
            P[i].GradRho[1]=0;
            P[i].GradRho[2]=1;
#endif
#ifdef SINK_WIND_SPAWN
            /* spawn reservoirs for sinks read from a fresh IC rather than formed in-run, which would
               otherwise hold whatever was in the malloc slot: fresh-IC reads skip every block above
               IO_U, so nothing has assigned these. Any snapshot-reading mode (RestartFlag >= 2) reads
               both reservoirs back -- or zero-fills when the datasets are absent -- so zeroing there
               would clobber restored mass already debited from the sink. */
            if(RestartFlag == 0)
            {
                P[i].unspawned_wind_mass = 0;
#ifdef SINGLE_STAR_FB_JETS
                P[i].unspawned_jet_mass = 0;
#endif
            }
#endif
#ifdef HERMITE_INTEGRATION
            /* the retained start-of-step state, normally written by the first half-step kick; seed it
               here so anything that reads it before that kick -- the IO_HERMITE_SYNC output blocks of a
               snapshot written straight after the ICs are read -- predicts from the actual state rather
               than from whatever the allocator left behind */
            P[i].OldPos = P[i].Pos; P[i].OldVel = P[i].Vel;
            P[i].Hermite_OldAcc = {}; P[i].OldJerk = {};
#ifdef HERMITE_SYMMETRIC_TIMESTEP_ACTIVE
            P[i].dt_crit_last = P[i].dt_crit_prev = 0; /* no criterion history: the default step rule until there is one */
#endif
#endif
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
#if defined(SINGLE_STAR_FB_SNE)
            P[i].Mass_final = P[i].Mass; // best guess, only matters if we restart in the middle of spawning an SN
#endif
#if defined(SINGLE_STAR_FB_WINDS)
            P[i].wind_mode = 0; P[i].wind_mode_time = All.Time; // wind_mode=0 triggers reset on first call; time initialised for hysteresis
            Vec3<double> nx,ny,nz; int kw; get_random_orthonormal_basis(P[i].ID,nx,ny,nz); for(kw=0;kw<3;kw++) {P[i].Wind_direction[kw] = nx[kw]; P[i].Wind_direction[kw+3] = ny[kw];}
#endif
#endif
#if defined(GALSF_FB_FIRE_RT_LOCALRP)
            P[i].NewStar_Momentum_For_JetFeedback = 0;
#endif
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
            P[i].SNe_ThisTimeStep = 0;
#endif
#ifdef GALSF_FB_MECHANICAL
            int k; for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {P[i].Area_weighted_sum[k] = 0;}
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
            P[i].MassReturn_ThisTimeStep = 0;
#ifdef GALSF_FB_FIRE_RPROCESS
            P[i].RProcessEvent_ThisTimeStep = 0;
#endif
#ifdef GALSF_FB_FIRE_AGE_TRACERS
            P[i].AgeDeposition_ThisTimeStep = 0;
#endif
#endif
#endif
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
            P[i].Time_Of_Last_MergeSplit = All.TimeBegin;
#endif
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
            P[i].Time_Of_Last_SmoothedVelUpdate = All.TimeBegin;
#endif
        }

#if defined(INIT_STELLAR_METALS_AGES_DEFINED) && defined(GALSF)
        if(RestartFlag == 0) {
#if defined(INPUT_READ_SINKPROPS) && defined(SINGLE_STAR_SINK_DYNAMICS)
            if(P[i].Type != 5) /* sink particles have their age read from the IC */
#endif
            {P[i].StellarAge = -2.0 * All.InitStellarAgeinGyr / (UNIT_TIME_IN_GYR) * get_random_number(P[i].ID + 3);}
        }
#endif
        
#ifdef GRAIN_FLUID
        if((RestartFlag == 0) && ((1 << P[i].Type) & (GRAIN_PTYPES)))
        {
            int grain_subtype = 1; P[i].Grain_Size = 0; /* default assumption about particulate sub-type for operations below */
#if defined(PIC_MHD)
            grain_subtype = P[i].MHD_PIC_SubType; /* check if the 'grains' are really PIC elements */
#endif
            /* Change grain mass to change the distribution of sizes.  Grain_Size_Spectrum_Powerlaw parameter sets d\mu/dln(R_d) ~ R_d^Grain_Size_Spectrum_Powerlaw */
            if(grain_subtype <= 2)
            {
                P[i].Grain_Size = All.Grain_Size_Min * exp( gsl_rng_uniform(random_generator) * log(All.Grain_Size_Max/All.Grain_Size_Min) );
                if(All.Grain_Size_Max > All.Grain_Size_Min*1.0001 && fabs(All.Grain_Size_Spectrum_Powerlaw) != 0) {
                    P[i].Mass *= (All.Grain_Size_Spectrum_Powerlaw/(pow(All.Grain_Size_Max/All.Grain_Size_Min,All.Grain_Size_Spectrum_Powerlaw)-1.)) *
                    pow(P[i].Grain_Size/All.Grain_Size_Min,All.Grain_Size_Spectrum_Powerlaw) * log(All.Grain_Size_Max/All.Grain_Size_Min);}
#ifdef GRAIN_RDI_TESTPROBLEM /* initialize various quantities for test problems from parameters set in the ICs */
                P[i].Mass *= All.Dust_to_Gas_Mass_Ratio;
                int k, non_gdir=1; Vec3<double> A={}, B={}, A_cross_B={}; double amag, rho_gas_expected, cs_gas_expected, acc_ang=All.Vertical_Grain_Accel_Angle * M_PI/180., tS0, a0, ct=1, tau2=0, ct2=0, w0, agamma=9.*M_PI/128.; B[2]=1; if(GRAV_DIRECTION_RDI==1) {non_gdir=2;}
                rho_gas_expected = 1.*UNIT_DENSITY_IN_CGS, cs_gas_expected = 1.*UNIT_VEL_IN_CGS; /* guess for the gas density here [set custom for e.g. stratified problems */
                tS0 = ((0.626657 * sqrt(GAMMA_DEFAULT) * P[i].Grain_Size * All.Grain_Internal_Density) / (rho_gas_expected * cs_gas_expected)) / UNIT_TIME_IN_CGS; /* stopping time [Epstein] for driftvel->0 [cgs->code units] */
                A[GRAV_DIRECTION_RDI]=cos(acc_ang)*All.Vertical_Grain_Accel - All.Vertical_Gravity_Strength; A[0]=sin(acc_ang)*All.Vertical_Grain_Accel; /* define angles/direction of external acceleration */
                amag=sqrt(A.norm_sq()+MIN_REAL_NUMBER); A /= amag;
                a0 = (tS0 * amag / (1.+All.Dust_to_Gas_Mass_Ratio)) / (cs_gas_expected/UNIT_VEL_IN_CGS) ; /* acc * tS0 / (1+mu) -- we're assuming that the code unit velocity equals the sound speed, for simplicity here */
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
                a0 *= All.Grain_Size_Max / P[i].Grain_Size;
#endif
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                double q_a = (0.75*All.Grain_Q_at_MaxGrainSize) / (All.Grain_Internal_Density*All.Grain_Size_Max), kappa_0 = All.Grain_Absorbed_Fraction_vs_Total_Extinction * q_a * All.Dust_to_Gas_Mass_Ratio; // this will be in cgs here
                double rho_base_setup = 1.*UNIT_DENSITY_IN_CGS, H_scale_setup = 1.*boxSize_X*UNIT_LENGTH_IN_CGS; // define in code units the -assumed- initial scaling of the base gas density and vertical scale-length (PROBLEM SPECIFIC HERE!)
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
                kappa_0 *= sqrt(All.Grain_Size_Max / All.Grain_Size_Min); // opacity must be corrected for dependence of Q on grainsize or lack thereof
#endif
                a0 *= exp(-kappa_0*rho_base_setup*H_scale_setup*(1.-exp(-P[i].Pos[2]/H_scale_setup))); // attenuate incident flux (and reduce acceleration) according to equilibrium expectation, if we're using single-scattering radiation pressure [otherwise comment this line out] //
#endif
                w0=sqrt((sqrt(1.+4.*agamma*a0*a0)-1.)/(2.*agamma)); // exact solution if no Lorentz forces and Epstein drag //
#ifdef GRAIN_LORENTZFORCE
                double Bmag, tL_i=0, tau2_0=0, f_tau_guess2=0; B = {All.BiniX, All.BiniY, All.BiniZ}; Bmag=B.norm(); B /= Bmag;
                tL_i = (All.Grain_Charge_Parameter*sqrt(rho_gas_expected/UNIT_DENSITY_IN_CGS)/((All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS)*(All.Grain_Size_Max/UNIT_LENGTH_IN_CGS))) * pow(All.Grain_Size_Max/P[i].Grain_Size,2) * Bmag; // 1/t_Lorentz in code units
                ct=dot(A,B); ct2=ct*ct; tau2_0=pow(tS0*tL_i,2); // variables for below //
                for(k=0;k<20;k++)
                {
                   tau2 = tau2_0 / (1. + agamma*w0*w0); // guess tau [including velocity dependence] //
                   f_tau_guess2 = (1.+tau2*ct2) / (1.+tau2); // what the projection factor (reduction in w from projection) would be //
                   w0 = sqrt((sqrt(1.+4.*agamma*a0*a0*f_tau_guess2)-1.)/(2.*agamma)); // re-calculate w0 with this //
                }
#endif
                w0 /= sqrt((1.+tau2)*(1.+tau2*ct2)); // ensures normalization to unity with convention below //
                A_cross_B = cross(A, B);
                P[i].Vel = (A + A_cross_B*sqrt(tau2) + B*(tau2*ct)) * (w0 * (cs_gas_expected/UNIT_VEL_IN_CGS));
#ifdef BOX_SHEARING
                // now add linearly the NHS drift solution for our shearing box setup
                double v00 = -All.Pressure_Gradient_Accel / (2. * BOX_SHEARING_OMEGA_BOX_CENTER);
                double v_K = -(P[i].Pos[0]-boxHalf_X) * BOX_SHEARING_Q*BOX_SHEARING_OMEGA_BOX_CENTER;
                double tau_s = tS0 * BOX_SHEARING_OMEGA_BOX_CENTER;
                v00 /= (1. + tau_s*tau_s); // appears in both terms here //
                P[i].Vel[0] += v00 * 2.*tau_s; // radial drift
                P[i].Vel[BOX_SHEARING_PHI_COORDINATE] = v_K + v00; // azimuthal drift relative to keplerian frame
#endif
#endif // closes rdi_testproblem
            }
            P[i].Gas_Density = P[i].Gas_InternalEnergy = 0; P[i].Gas_Velocity = {}; P[i].Grain_AccelTimeMin = MAX_REAL_NUMBER;
#if defined(GRAIN_BACKREACTION)
            P[i].Grain_DeltaMomentum = {};
#endif
#if defined(GRAIN_LORENTZFORCE)
            P[i].Gas_B = {};
#endif
        } // closes check on restartflag and particle type
#endif // closes grain_fluid


#ifdef METALS
        if(RestartFlag == 0) {
#if defined(INIT_STELLAR_METALS_AGES_DEFINED)
            P[i].Metallicity[0] = All.InitMetallicityinSolar*All.SolarAbundances[0];
#else
            P[i].Metallicity[0] = 0;
#endif
            /* initialize abundance ratios. for now, assume solar */
            for(j=0;j<NUM_METAL_SPECIES;j++) {P[i].Metallicity[j]=All.SolarAbundances[j]*(P[i].Metallicity[0]/All.SolarAbundances[0]);}
            /* need to allow for a primordial He abundance */
            if(NUM_LIVE_SPECIES_FOR_COOLTABLES>=10) P[i].Metallicity[1]=(1.-HYDROGEN_MASSFRAC)+(All.SolarAbundances[1]-(1.-HYDROGEN_MASSFRAC))*P[i].Metallicity[0]/All.SolarAbundances[0];
        } // if(RestartFlag == 0)

#if defined(GALSF_ISMDUSTCHEM_MODEL)
    if (P[i].Type == 0) {Initialize_ISMDustChem_Particle_Variables(i);}
#endif

#ifdef CHIMES
#ifdef COOL_METAL_LINES_BY_SPECIES
	if (P[i].Type == 0)
	  {
	    double H_mass_fraction = 1.0 - (P[i].Metallicity[0] + P[i].Metallicity[1]);
	    ChimesGasVars[i].element_abundances[0] = (ChimesFloat) (P[i].Metallicity[1] / (4.0 * H_mass_fraction));   // He
	    ChimesGasVars[i].element_abundances[1] = (ChimesFloat) (P[i].Metallicity[2] / (12.0 * H_mass_fraction));  // C
	    ChimesGasVars[i].element_abundances[2] = (ChimesFloat) (P[i].Metallicity[3] / (14.0 * H_mass_fraction));  // N
	    ChimesGasVars[i].element_abundances[3] = (ChimesFloat) (P[i].Metallicity[4] / (16.0 * H_mass_fraction));  // O
	    ChimesGasVars[i].element_abundances[4] = (ChimesFloat) (P[i].Metallicity[5] / (20.0 * H_mass_fraction));  // Ne
	    ChimesGasVars[i].element_abundances[5] = (ChimesFloat) (P[i].Metallicity[6] / (24.0 * H_mass_fraction));  // Mg
	    ChimesGasVars[i].element_abundances[6] = (ChimesFloat) (P[i].Metallicity[7] / (28.0 * H_mass_fraction));  // Si
	    ChimesGasVars[i].element_abundances[7] = (ChimesFloat) (P[i].Metallicity[8] / (32.0 * H_mass_fraction));  // S
	    ChimesGasVars[i].element_abundances[8] = (ChimesFloat) (P[i].Metallicity[9] / (40.0 * H_mass_fraction));  // Ca
	    ChimesGasVars[i].element_abundances[9] = (ChimesFloat) (P[i].Metallicity[10] / (56.0 * H_mass_fraction)); // Fe

	    ChimesGasVars[i].metallicity = (ChimesFloat) (P[i].Metallicity[0] / 0.0129);  // In Zsol. CHIMES uses Zsol = 0.0129.
	    ChimesGasVars[i].dust_ratio = ChimesGasVars[i].metallicity;
	  }
#else
	if (ThisTask == 0)
	  {
	    printf("ERROR: Config flags CHIMES and METALS are switched on, but COOL_METAL_LINES_BY_SPECIES is switched off. \n");
	    printf("If you want to include metals with CHIMES, you will also need to switch on COOL_METAL_LINES_BY_SPECIES. Aborting. \n");
	    endrun(202);
	  }
#endif // COOL_METAL_LINES_BY_SPECIES
#endif // CHIMES
#else
#ifdef CHIMES
	if (P[i].Type == 0)
	  {
	    double H_mass_fraction = HYDROGEN_MASSFRAC;
	    ChimesGasVars[i].element_abundances[0] = (ChimesFloat) ((1.0 - H_mass_fraction) / (4.0 * H_mass_fraction));  // He
	    for (j = 1; j < 10; j++) {ChimesGasVars[i].element_abundances[j] = 0.0;}
	    ChimesGasVars[i].metallicity = 0.0;
	    ChimesGasVars[i].dust_ratio = 0.0;
	  }
#endif // CHIMES
#endif // METALS



#ifdef SINK_PARTICLES
#if (SINGLE_STAR_SINK_FORMATION & 8)
        P[i].Sink_Ngb_Flag = 0;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
 // start with a large value (> plausible values v_ejecta or v_wind) as a conservative choice when starting up a simulation with an active feedback-emmiting star - this will get updated to a more reasonable value once the particle walks the gravity tree, but need this to ensure the first timestep is stable.
        P[i].MaxFeedbackVel = 1e4 / UNIT_VEL_IN_KMS;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
	    P[i].Min_Sink_Approach_Time = P[i].Min_Sink_Freefall_time = MAX_REAL_NUMBER;
#if (SINGLE_STAR_TIMESTEPPING > 0)
	    P[i].SuperTimestepFlag = 0;
#endif
#endif
        if(P[i].Type == 5)
        {
            count_holes++;
            if(RestartFlag == 0)
            {
                P[i].Sink_Mass = All.SeedSinkMass;
                P[i].Sink_Formation_Mass = P[i].Mass;
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
                P[i].Sink_Mdot_ROI = 0;
                P[i].Sink_ROI = 0;
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
                P[i].Sink_Mass = P[i].Mass;
#endif
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION // properly initialize luminosity
                if(P[i].ProtoStellarRadius_inSolar <= 0) {P[i].ProtoStellarRadius_inSolar = ps_radius_MS_in_solar(P[i].Sink_Mass);} /* initialize protostellar radius if not set */
                singlestar_subgrid_protostellar_evolution_update_track(i,0,0);
#if (SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION == 2)
                P[i].ZAMS_Mass = P[i].Sink_Mass;
                calculate_individual_stellar_luminosity(P[i].Sink_Mdot, P[i].Sink_Mass, i);
#endif                
#endif
#ifdef GRAIN_FLUID
                P[i].Sink_Dust_Mass = 0;
#endif
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
                P[i].SinkRadius = KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER * SinkParticle_GravityKernelRadius;
#endif
#ifdef SINK_ALPHADISK_ACCRETION
                P[i].Sink_Mass_Reservoir = All.SeedReservoirMass;
#endif
#ifdef SINK_FOLLOW_ACCRETED_ANGMOM
                double sink_mu=2*get_random_number(P[i].ID+3)-1, sink_phi=2*M_PI*get_random_number(P[i].ID+4), sink_sin=sqrt(1-sink_mu*sink_mu);
                double spin_prefac = All.G * P[i].Sink_Mass / C_LIGHT_CODE; // assume initially maximally-spinning black hole with random orientation
                P[i].Sink_Specific_AngMom[0]=spin_prefac*sink_sin*cos(sink_phi); P[i].Sink_Specific_AngMom[1]=spin_prefac*sink_sin*sin(sink_phi); P[i].Sink_Specific_AngMom[2]=spin_prefac*sink_mu;
#endif
#ifdef SINK_COUNTPROGS
                P[i].Sink_CountProgs = 1;
#endif
            }
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
            P[i].dt_since_last_gas_search = 0;
            P[i].do_gas_search_this_timestep = 1;
#endif 
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
            if(RestartFlag != 1) {P[i].Sink_AccretionDeficit = 0;}
#endif
        }
#endif
    }

#ifdef SINK_PARTICLES
    MPI_Allreduce(&count_holes, &All.TotSinks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    for(i = 0; i < TIMEBINS; i++) {TimeBinActive[i] = 1;}

    reconstruct_timebins();

#ifdef PMGRID
    All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

    for(i = 0; i < N_gas; i++)	/* initialize gas/fluid cell properties */
    {
        CellP[i].InternalEnergyPred = CellP[i].InternalEnergy;

        CellP[i].VelPred = P[i].Vel;
        CellP[i].HydroAccel = {};

        P[i].Particle_DivVel = 0;
        CellP[i].ConditionNumber = 1;
        CellP[i].DtInternalEnergy = 0;
        if(RestartFlag == 0)
        {
            /* first-guess temperature: gamma_eos_value() -- and through it the H2 partition series --
               reads this before the first EOS solve has written it, and a fresh IC read leaves it as
               malloc garbage (a binary IC cannot even probe for the block). A huge garbage value is
               exactly what lets the rotational series run away. Neutral-atomic ideal-gas estimate;
               the first set_eos_pressure() replaces it with the solved value. Snapshot restarts read
               the written Temperature back, so they are left alone. */
            CellP[i].Temperature = 1.22 * (GAMMA_DEFAULT-1.) * U_TO_TEMP_UNITS * CellP[i].InternalEnergy;
        }
        CellP[i].FaceClosureError = 0;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        CellP[i].MaxKineticEnergyNgb = 0;
#endif
        CellP[i].Mass = P[i].Mass;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        CellP[i].dMass = 0;
        CellP[i].DtMass = 0;
        CellP[i].MassTrue = P[i].Mass;
        CellP[i].GravWorkTerm = {};
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
        P[i].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        if(1 & ADAPTIVE_GRAVSOFT_FORALL) {P[i].AGS_KernelRadius = P[i].KernelRadius;} else {P[i].AGS_KernelRadius = All.ForceSoftening[P[i].Type];}
#endif
#endif

#ifdef CONDUCTION
        CellP[i].Kappa_Conduction = 0;
#endif
#ifdef MHD_NON_IDEAL
        CellP[i].Eta_MHD_OhmicResistivity_Coeff = 0;
        CellP[i].Eta_MHD_HallEffect_Coeff = 0;
        CellP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
#ifdef VISCOSITY
        CellP[i].Eta_ShearViscosity = 0;
        CellP[i].Zeta_BulkViscosity = 0;
#endif


#ifdef TURB_DIFFUSION
        CellP[i].TD_DiffCoeff = 0;

#ifdef TURB_DIFF_DYNAMIC
        int u, v; /* start with the standard Smagorinsky-Lilly constant from Kolmogorov theory */
        CellP[i].TD_DynDiffCoeff = 0.01;
        CellP[i].h_turb = 0;
        CellP[i].FilterWidth_bar = 0;
        CellP[i].MagShear_bar = 0;
        CellP[i].MagShear = 0;
        CellP[i].Norm_hat = 0;
        CellP[i].Dynamic_numerator = 0;
        CellP[i].Dynamic_denominator = 0;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        CellP[i].TD_DynDiffCoeff_error = 0;
#endif
        if (RestartFlag != 7) {
            CellP[i].Velocity_bar = {};
            CellP[i].Velocity_hat = {};
        }
        for (u = 0; u < 3; u++) {
            for (v = 0; v < 3; v++) {
                CellP[i].VelShear_bar[u][v] = 0;
            }
        }
#endif
#endif

        if(RestartFlag == 0)
        {
#ifndef INPUT_READ_KERNELRADIUS
            P[i].KernelRadius = 0;
#endif
            CellP[i].Density = -1;
#ifdef COOLING
#ifndef CHIMES
            CellP[i].Ne = 1.0;
#endif
#if defined(COOL_MOLECFRAC_NONEQM)
            CellP[i].MolecularMassFraction = 0.0; CellP[i].MolecularMassFraction_perNeutralH = 0.0; // start atomic
#endif
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
            CellP[i].Rad_Flux_UV = 0;
            CellP[i].Rad_Flux_EUV = 0;
#endif
#ifdef CHIMES_STELLAR_FLUXES
	    int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {CellP[i].Chimes_fluxPhotIon[kc] = 0; CellP[i].Chimes_G0[kc] = 0;}
#endif
#ifdef SINK_COMPTON_HEATING
            CellP[i].Rad_Flux_AGN = 0;
#endif
        }
#ifdef GALSF_SUBGRID_WINDS
        if(RestartFlag == 0) {CellP[i].DelayTime = 0;}
#if (GALSF_SUBGRID_WIND_SCALING==1)
        CellP[i].HostHaloMass = 0;
#endif
#endif
#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
        CellP[i].DelayTimeHII = 0;
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        CellP[i].DelayTimeCoolingSNe = 0;
#endif
#ifdef GALSF
        CellP[i].Sfr = 0;
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED)
        CellP[i].AlphaVirial_SF_TimeSmoothed = 0;
#endif
#endif
#ifdef COSMIC_RAY_FLUID
        if(RestartFlag == 0) {for(j=0;j<N_CR_PARTICLE_BINS;j++) {CellP[i].CosmicRayEnergy[j] = 0;}}
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
        if(RestartFlag != 1) {CellP[i].DtCREgyNewInjectionFromShocks = 0;}
#endif
#if defined(SINK_CR_INJECTION_AT_TERMINATION)
        if(RestartFlag != 1) {CellP[i].Sink_CR_Energy_Available_For_Injection = 0;}
#endif
#if defined(CRFLUID_EVOLVE_SPECTRUM)
        //if(RestartFlag == 0) {for(j=0;j<N_CR_PARTICLE_BINS;j++) {CellP[i].CosmicRay_PwrLaw_Slopes_in_Bin[j] = -2.5; CellP[i].CosmicRay_Number_in_Bin[j] = 0; CellP[i].DtCosmicRay_Number_in_Bin[j] = 0;}} // initialize a flat spectrum in each bin
        if(RestartFlag == 0) {for(j=0;j<N_CR_PARTICLE_BINS;j++) {CellP[i].CosmicRay_Number_in_Bin[j] = 0; CellP[i].DtCosmicRay_Number_in_Bin[j] = 0;}} // initialize the number in each bin
        if(RestartFlag == 2) { // if we don't directly evolve the slopes, we do write them out and read them in, so need to re-construct the correct number in bin from what we actually read, which was the -slope- information
#if defined(CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART)
            double e0 = 0.5 * P[i].Mass * CellP[i].InternalEnergy; // snapshot from which we read does not have the full CR info, so we need to initialize it from something //
#if (CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART==1)
            e0 = CellP[i].CosmicRayEnergy[0]; // we had one value of energy (total) available to read in //
#endif
            // now define the desired spectrum //
            for(j=0;j<N_CR_PARTICLE_BINS;j++) {
                int species = CR_species_ID_in_bin[j];
                double f_norm = 1.e-20;
                double f_elec = 0.05; // fraction of the energy to put into e- as opposed to p+ at injection [early experiments with 'observed'  fraction ~ 1% give lower e-/p+ actually observed in the end, so tentative favoring closer to equal at injection? but not run to z=0, so U_rad high from CMB; still experimenting here]
                if(species == -1) {f_norm = f_elec;} // e-
                if(species == +1) {f_norm = 1.-f_elec;} // p
                if(species == -2) {f_norm = 1.e-10 * f_elec;} // e+ (initialize to negligible since want to start with primary)
                if(species > 1)
                {
                    double Zfac = P[i].Metallicity[0]/All.SolarAbundances[0]; // scale heavier elements to the metallicity of the gas into which CRs are being accelerated
                    Zfac *= pow(return_CRbin_CRmass_in_mp(-1,j) / fabs(return_CRbin_CR_charge_in_e(-1,j)) , 1.5); // approximate injection factor for a constant-beta distribution at a given R_GV needed below
                    if(species == 2) {f_norm = 3.7e-9 * Zfac;} // B (for standard elements initialize to solar ratios assuming similar energy/nucleon)
                    if(species == 3) {f_norm = 2.4e-3 * Zfac;} // C
                    if(species == 4) {f_norm = 1.4e-10 * Zfac;} // Be7+9 (stable)
                    if(species == 5) {f_norm = 1.4e-20 * Zfac;} // Be10 (radioactive)
                    if(species == 6) {f_norm = 0.0094 * Zfac;} // CNO (combined bin)
                }
                double e_tmp = 1.0 * f_norm * e0, x_RGV = CR_global_rigidity_at_bin_center[j], a_0=0.715197, x_0=1.7, gamma_0=-2.3, gamma_1=0.8, fac, slope; // hadrons
                if(species < 0) {a_0=0.506309; x_0=1.0; gamma_0=-0.6; gamma_1=1.3;} // leptons
                double dlnR = log(CR_global_max_rigidity_in_bin[j]/CR_global_min_rigidity_in_bin[j]); // bin-width needed below for assignment
                double qx0=pow(x_RGV/x_0,gamma_0), qx1=pow(x_RGV/x_0,gamma_1); fac=a_0*dlnR/(qx0+qx1); slope=-(gamma_0*qx0 + gamma_1*qx1)/(qx0+qx1)-2.; // adopt an extremely simple two-power law spectrum, identical in E space for everything, except normalization, to initialize
                //double fac = 2.3 / (pow(x_RGV,-0.6) + pow(x_RGV,0.8)), slope = (3. - 4.*pow(x_RGV,1.4)) / (5. + 5.*pow(x_RGV,1.4)) - 2.;
                if(CR_check_if_bin_is_nonrelativistic(j)) {slope -= 1.;} // correct for NR terms
                CellP[i].CosmicRayEnergy[j] = e_tmp * fac; CellP[i].CosmicRay_Number_in_Bin[j] = slope; // actually assign the energy and power-law slope
            }
#endif
            /* now initialize the number in each bin from the slopes that we have either read in or assumed */
            for(j=0;j<N_CR_PARTICLE_BINS;j++) {
                double slope_from_snapshot = CellP[i].CosmicRay_Number_in_Bin[j];
                CellP[i].CosmicRay_Number_in_Bin[j] = CR_get_number_in_bin_from_slope(i,j,CellP[i].CosmicRayEnergy[j],slope_from_snapshot);
            }
        }
#endif
#endif
#ifdef MAGNETIC
#if defined MHD_B_SET_IN_PARAMS
        if(RestartFlag == 0)
        {			/* Set only when starting from ICs */
            CellP[i].B[0]=CellP[i].BPred[0] = All.BiniX;
            CellP[i].B[1]=CellP[i].BPred[1] = All.BiniY;
            CellP[i].B[2]=CellP[i].BPred[2] = All.BiniZ;
        }
#endif /*MHD_B_SET_IN_PARAMS*/
        CellP[i].BPred *= a2_fac * gauss2gizmo;
        CellP[i].B = CellP[i].BPred;
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
        CellP[i].Balpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
        CellP[i].Phi = CellP[i].PhiPred = CellP[i].DtPhi = 0;
#endif
#ifdef SINK_RETURN_BFLUX
        P[i].B[0] = P[i].B[1] = P[i].B[2] = 0;
#endif
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        CellP[i].alpha = 0.0;
#endif
#if defined(SINK_THERMALFEEDBACK)
        CellP[i].Injected_Sink_Energy = 0;
#endif
    }

#ifndef BOX_SHEARING
#if (NUMDIMS==2)
    for(i = 0; i < NumPart; i++)
    {
        P[i].Pos[2] = 0;
        //P[i].Vel[2] = 0; // this should be set in the ICs, not here //

        P[i].GravAccel[2] = 0;

        if(P[i].Type == 0)
        {
            CellP[i].VelPred[2] = 0;
            CellP[i].HydroAccel[2] = 0;
        }
    }
#endif
#endif

#if (NUMDIMS==1)
    for(i = 0; i < NumPart; i++)
    {
        P[i].Pos[1] = P[i].Pos[2] = 0;
        //P[i].Vel[1] = P[i].Vel[2] = 0; // this should be set in the ICs, not here //

        P[i].GravAccel[1] = P[i].GravAccel[2] = 0;

        if(P[i].Type == 0)
        {
            CellP[i].VelPred[1] = CellP[i].VelPred[2] = 0;
            CellP[i].HydroAccel[1] = CellP[i].HydroAccel[2] = 0;
        }
    }
#endif

#ifdef ASSIGN_NEW_IDS
    assign_unique_ids();
#endif
    /* assign other ID parameters needed */

    if(RestartFlag==0) {for(i = 0; i < NumPart; i++) {P[i].ID_child_number = 0; P[i].ID_generation = 0;}}
#ifdef NO_CHILD_IDS_IN_ICS
    if(RestartFlag != 1) {for(i = 0; i < NumPart; i++) {P[i].ID_child_number = 0; P[i].ID_generation = 0;}}
#endif

#ifdef TEST_FOR_IDUNIQUENESS
    test_id_uniqueness();
#endif

    Flag_FullStep = 1;		/* to ensure that Peano-Hilbert order is done */
    TreeReconstructFlag = 1;
    DomainReconstructFlag = 1;	/* fresh start: both the decomposition and the build are owed */

#ifdef SINK_WIND_SPAWN
    Max_Unspawned_MassUnits_fromSink = 0;
#endif

#ifdef SHIFT_BY_HALF_BOX
    for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
        {
            double boxtmp = 0;
            if(j==0) {boxtmp = boxSize_X;}
            if(j==1) {boxtmp = boxSize_Y;}
            if(j==2) {boxtmp = boxSize_Z;}
            P[i].Pos[j] += 0.5 * boxtmp;
        }
#endif


    Gas_split = 0;
#ifdef GALSF
    Stars_converted = 0;
#endif
    domain_Decomposition(0, 0, 0);	/* do initial domain decomposition (gives equal numbers of particles) */

    set_softenings();

    /* will build tree */
    ngb_treebuild();

    All.Ti_Current = 0;

    if(RestartFlag != 3 && RestartFlag != 5) {setup_smoothinglengths();}

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(RestartFlag != 3 && RestartFlag != 5) {ags_setup_smoothinglengths();}
#endif
#ifdef CBE_INTEGRATOR
    do_cbe_initialization();
#endif

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)
    if(RestartFlag != 3 && RestartFlag != 5) {disp_setup_smoothinglengths();}
#endif
#endif

#if defined GALSF_SFR_IMF_VARIATION
    for(i = 0; i < NumPart; i++) {P[i].IMF_Mturnover = 2.0;} // reset to normal IMF
#endif

#if defined(WAKEUP) && defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
    for(i=0;i<NumPart;i++) {P[i].wakeup=0;}
#endif


    /* HELLO! This here is where you should insert custom code for hard-wiring the ICs of various test problems */



    density();
    for(i = 0; i < N_gas; i++)	/* initialize gas/fluid cell properties */
    {
        int k __attribute__((unused)) = 0;
        CellP[i].InternalEnergyPred = CellP[i].InternalEnergy;
        CellP[i].recent_refinement_flag = 0; // always initialize to zero, no recent refinements

        // re-match the predicted and initial velocities and B-field values, just to be sure //
        CellP[i].VelPred = P[i].Vel;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==0)
        CellP[i].ParticleVel = {}; // set these to zero and forget them, for the rest of the run //
#endif

#ifdef MAGNETIC
        CellP[i].B = CellP[i].BPred * (P[i].Mass / CellP[i].Density); // convert to the conserved unit V*B //
        CellP[i].BPred = CellP[i].B; CellP[i].DtB = {}; CellP[i].BField_prerefinement = {};
#endif
#ifdef COSMIC_RAY_FLUID
        for(k=0;k<N_CR_PARTICLE_BINS;k++)
        {
            CellP[i].CosmicRayEnergyPred[k]=CellP[i].CosmicRayEnergy[k]; CellP[i].CosmicRayDiffusionCoeff[k]=0; CellP[i].DtCosmicRayEnergy[k]=0;
            CellP[i].CosmicRayFlux[k] = {}; CellP[i].CosmicRayFluxPred[k] = {};
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(j=0;j<2;j++) {CellP[i].CosmicRayAlfvenEnergy[k][j]=0; CellP[i].CosmicRayAlfvenEnergyPred[k][j]=0; CellP[i].DtCosmicRayAlfvenEnergy[k][j]=0;}
#endif
        }
#endif
#if defined(EOS_ELASTIC)
        if(RestartFlag != 1)
        {
            CellP[i].Dt_Elastic_Stress_Tensor = CellP[i].Elastic_Stress_Tensor_Pred = CellP[i].Elastic_Stress_Tensor = {};
        } else {
            CellP[i].Elastic_Stress_Tensor_Pred = CellP[i].Elastic_Stress_Tensor; CellP[i].Dt_Elastic_Stress_Tensor = {};
        }
#endif
        CellP[i].DtInternalEnergy = 0;
#if defined(COOLING) && !defined(COOLING_OPERATOR_SPLIT)
        CellP[i].CoolingIsOperatorSplitThisTimestep = 1; /* default to more conservative split */
#endif
        CellP[i].Mass = P[i].Mass;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        CellP[i].dMass = 0;
        CellP[i].DtMass = 0;
        CellP[i].MassTrue = P[i].Mass;
        CellP[i].GravWorkTerm = {};
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
        P[i].AGS_zeta = 0;
#endif
#ifdef WAKEUP
        if(RestartFlag!=0) {P[i].wakeup=0;}
        NeedToWakeupParticles = 0;
        NeedToWakeupParticles_local = 0;
#endif
#ifdef SUPER_TIMESTEP_DIFFUSION
        CellP[i].Super_Timestep_Dt_Explicit = 0;
        CellP[i].Super_Timestep_j = 0;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        CellP[i].Rad_Flux_UV = 0;
        CellP[i].Rad_Flux_EUV = 0;
#endif
#ifdef SINK_COMPTON_HEATING
        CellP[i].Rad_Flux_AGN = 0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[i].Rad_E_gamma[kf]=0;}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[i].Rad_Flux[kf] = {};}}
#endif
#if defined(COSMIC_RAY_SUBGRID_LEBRON)
        CellP[i].SubGrid_CosmicRayEnergyDensity = 0;
#endif
        
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES>=2)
        if(RestartFlag != 1) {if(P[i].ID == All.SpawnedWindCellID) {P[i].ID += 1;}} // reset any of these so can obey desired merge-split rules
#endif

#ifdef COOL_GRACKLE
        if(RestartFlag == 0)
        {
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            CellP[i].grHI    = HYDROGEN_MASSFRAC;
            CellP[i].grHII   = 1.0e-20;
            CellP[i].grHM    = 1.0e-20;
            CellP[i].grHeI   = 1.0 - HYDROGEN_MASSFRAC;
            CellP[i].grHeII  = 1.0e-20;
            CellP[i].grHeIII = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            CellP[i].grH2I   = 1.0e-20;
            CellP[i].grH2II  = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            CellP[i].grDI    = 2.0 * 3.4e-5;
            CellP[i].grDII   = 1.0e-20;
            CellP[i].grHDI   = 1.0e-20;
#endif
        }
#endif

    }


    /* we should define the maximum and minimum particle masses
        below/above which particles are merged/split */
    if(RestartFlag != 1)
    {
        double mass_min = MAX_REAL_NUMBER;
        double mass_max = -MAX_REAL_NUMBER;
        double mass_tot = 0;
        for(i = 0; i < N_gas; i++)	/* initialize gas/fluid cell properties */
        {
            mass_tot += P[i].Mass;
            if(P[i].Mass > mass_max) mass_max = P[i].Mass;
            if(P[i].Mass < mass_min) mass_min = P[i].Mass;
        }
        /* broadcast this and get the min and max values over all processors */
        double mpi_mass_min, mpi_mass_max;
        MPI_Allreduce(&mass_min, &mpi_mass_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&mass_max, &mpi_mass_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double mpi_splitmerge_readmin, mpi_splitmerge_readmax; /* check if this has been initialized by broadcasting to all processors */
        MPI_Allreduce(&All.MinMassForParticleMerger, &mpi_splitmerge_readmin, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&All.MaxMassForParticleSplit, &mpi_splitmerge_readmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if(mpi_splitmerge_readmin <= 0) { /* initialize if this isn't saved in the ICs */
            All.MinMassForParticleMerger = 0.49 * mpi_mass_min;
#ifdef GALSF_GENERATIONS
            All.MinMassForParticleMerger /= (float)GALSF_GENERATIONS;
#endif
        } else {All.MinMassForParticleMerger = mpi_splitmerge_readmin;} /* use the version from the ICs */
        if(mpi_splitmerge_readmax <= 0) {All.MaxMassForParticleSplit  = 3.01 * mpi_mass_max;} else {All.MaxMassForParticleSplit = mpi_splitmerge_readmax;}
#ifdef MERGESPLIT_HARDCODE_MAX_MASS
        All.MaxMassForParticleSplit = MERGESPLIT_HARDCODE_MAX_MASS;
#endif
#ifdef MERGESPLIT_HARDCODE_MIN_MASS
        All.MinMassForParticleMerger = MERGESPLIT_HARDCODE_MIN_MASS;
#endif

#ifdef SINGLE_STAR_SINK_DYNAMICS /* Get mean gas mass, used in various subroutines */
        double mpi_mass_tot; long mpi_Ngas; long Ngas_l = (long) N_gas;
        MPI_Allreduce(&mass_tot, &mpi_mass_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Ngas_l, &mpi_Ngas, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        All.MeanGasParticleMass = mpi_mass_tot/( (double)mpi_Ngas );
        if(RestartFlag==0){
            for(i=0; i<NumPart; i++){
                if(P[i].Type==5){P[i].Sink_Formation_Mass = All.MeanGasParticleMass;} // will behave as if this sink formed from a gas cell with the average mass
            }
        }
#endif
    }

    if(RestartFlag == 3)
    {
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
        if(ThisTask == 0) {printf("*AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE* Computation of softening lengths... \n");}
        ags_setup_smoothinglengths();
        if(ThisTask == 0) {printf("*AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE* Computation of softening lengths done. \n");}
#endif

#ifdef FOF
        fof_fof(RestartSnapNum);
#endif
        endrun(0);
    }

#ifdef OUTPUT_TWOPOINT_ENABLED
    if(RestartFlag == 5)
    {
        /* calculating powerspec and twopoint function */
#ifdef PMGRID
        long_range_init_regionsize();
#ifdef BOX_PERIODIC
        /* determine global and local particle numbers */
        int n, n_type[6]; long long ntot_type_all[6];
        for(n = 0; n < 6; n++) {n_type[n] = 0;}
        for(n = 0; n < NumPart; n++) {n_type[P[n].Type]++;}
        sumup_large_ints(6, n_type, ntot_type_all);
        calculate_power_spectra(RestartSnapNum, ntot_type_all);
#endif
#endif
        force_treebuild(NumPart, NULL);
        twopoint();
        endrun(0);
    }
#endif


    if(RestartFlag == 4)
    {
        All.Time = All.TimeBegin = header.time;
        snprintf(All.SnapshotFileBase, 100, "%s_converted", All.SnapshotFileBase);
        if(ThisTask == 0) {printf("Start writing file %s\n", All.SnapshotFileBase);}
        printf("RestartSnapNum %d\n", RestartSnapNum);

        All.TopNodeAllocFactor = 0.008;

        savepositions(RestartSnapNum);
        endrun(0);
    }
    

#if defined(COOL_MOLECFRAC_NONEQM)
    if(RestartFlag == 2) // should have read in CellP[i].MolecularMassFraction_perNeutralH
    {
        CellP[i].MolecularMassFraction_perNeutralH = DMIN(1,DMAX(0,CellP[i].MolecularMassFraction_perNeutralH));
        CellP[i].MolecularMassFraction = DMIN(1,DMAX(0, 1.-CellP[i].Ne/1.25)) * CellP[i].MolecularMassFraction_perNeutralH;
    }
#endif
    

#ifdef CHIMES_INITIALISE_IN_EQM
    if (RestartFlag != 1)
      {
	/* Note that stellar fluxes computed through the
	 * gravity tree are all zero at this stage,
	 * because the gravitational forces have not yet
	 * been computed. So the equilibrium abundances
	 * computed here include only the extragalactic UVB. */
	if (ThisTask == 0)
	  printf("Computing equilibrium CHIMES abundances. \n");

	int iter_number;

#ifdef _OPENMP
	int ThisThread;

#pragma omp parallel private(i, iter_number, ThisThread)
	{
	  ThisThread = omp_get_thread_num();

#pragma omp for schedule(dynamic)
#endif
	  for(i = 0; i < N_gas; i++)
	    {
	      initialise_gas_abundances(&(ChimesGasVars[i]), &ChimesGlobalVars);

#ifdef CHIMES_TURB_DIFF_IONS
	      chimes_update_turbulent_abundances(i, 1, P, CellP);
#endif

	      chimes_update_gas_vars(i, P, CellP);

	      // Evolve the chemistry for (1 / nH) Myr (limited to 1 Gyr) ten times at fixed temperature.
	      ChimesGasVars[i].hydro_timestep = (ChimesFloat) DMIN(3.16e13 / ChimesGasVars[i].nH_tot, 3.16e16);
	      ChimesGasVars[i].ThermEvolOn = 0;

	      for (iter_number = 0; iter_number < 10; iter_number++) chimes_network(&(ChimesGasVars[i]), &ChimesGlobalVars);


#ifdef CHIMES_TURB_DIFF_IONS
	      chimes_update_turbulent_abundances(i, 1, P, CellP);
#endif
	    }
#ifdef _OPENMP
	} // End of parallel block
#endif
      } // RestartFlag != 1
#endif // CHIMES_INITIALISE_IN_EQM
}



/*! This routine computes the mass content of the box and compares it to the specified value of Omega-matter.  If discrepant, the run is terminated. */
#ifdef BOX_PERIODIC
void check_omega(void)
{
    double mass = 0, masstot, omega; int i;
    for(i = 0; i < NumPart; i++) {mass += P[i].Mass;}
    MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    omega = masstot / (boxSize_X*boxSize_Y*boxSize_Z) / (3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G));
#ifdef GR_TABULATED_COSMOLOGY_G
    omega *= All.Gini / All.G;
#endif
    if(fabs(omega - All.OmegaMatter) > 1.0e-2) // look for a 1% tolerance of omega-matter
        {PRINT_WARNING("\n\nMass content in the ICs accounts only for Omega_M=%g,\nbut you specified Omega_M=%g in the parameterfile.\nRun will stop.\n",omega, All.OmegaMatter); endrun(1);}
}
#endif


/*! This function is used to find an initial kernel length (what used to be called the
 *  'smoothing length' for SPH, but is just the kernel size for the mesh-free methods) for each gas
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the kernel length is provided to the function density(), which will
 *  then iterate if needed to find the right kernel length.
 */
void setup_smoothinglengths(void)
{
    int i, no, p;
    if((RestartFlag == 0)||(RestartFlag==2)) // best for stability if we re-calc KernelRadius for snapshot restarts //
    {
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES) || defined(GRAIN_FLUID)
        for(i = 0; i < NumPart; i++)
#else
        for(i = 0; i < N_gas; i++)
#endif
        {
                no = Father[i];
                /* Negative means the treebuild never reached this particle (see
                 * force_treebuild_single), so there is no parent node to size the kernel from.
                 * Start from the root rather than indexing Nodes with -1; the guess is crude but
                 * finite, and gets refined iteratively. */
                if(no < 0) {no = All.MaxPart;}
                while(2 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0) {break;}
                    no = p;
                }

                if((RestartFlag == 0)||(P[i].Type != 0)) // if Restartflag==2, use the saved KernelRadius of the gas as initial guess //
                {
#ifndef INPUT_READ_KERNELRADIUS
#if NUMDIMS == 3
                    P[i].KernelRadius = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.333333) * Nodes[no].len;
#endif
#if NUMDIMS == 2
                    P[i].KernelRadius = pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.5) * Nodes[no].len;
#endif
#if NUMDIMS == 1
                    P[i].KernelRadius = All.DesNumNgb * (P[i].Mass / Nodes[no].u.d.mass) * Nodes[no].len;
#endif
#ifndef SELFGRAVITY_OFF
                    double soft = All.ForceSoftening[P[i].Type];
                    int ags_active_flag = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                    if(P[i].Type == 0) {ags_active_flag = 1;}
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                    if((ADAPTIVE_GRAVSOFT_FORALL) & (1 << P[i].Type)) {ags_active_flag = 1;}
#endif
                    if(ags_active_flag) {
                        /* adaptive softening: fixed softening is just a floor and can be arbitrarily small,
                           so only reset KernelRadius if the tree data itself is degenerate */
                        if((Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {P[i].KernelRadius = 0.05 * All.BoxSize;}
                    } else {
                        if(soft != 0) {if((P[i].KernelRadius>100.*soft)||(P[i].KernelRadius<=0.01*soft)||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {P[i].KernelRadius = soft;}}
                    }
#else
                    if((Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {
#if (defined(BOX_PERIODIC) || defined(BOX_SHEARING) || defined(BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE) || defined(BOX_LONG_X) || defined(BOX_LONG_Y) || defined(BOX_LONG_Z))
                        P[i].KernelRadius = 0.05 * All.BoxSize;
#else
                        P[i].KernelRadius = 1;
#endif
                    }
#endif
#endif // INPUT_READ_KERNELRADIUS
                } // closes if((RestartFlag == 0)||(P[i].Type != 0))
            }
    }
    if((RestartFlag==0 || RestartFlag==2) && All.ComovingIntegrationOn) {for(i=0;i<N_gas;i++) {P[i].KernelRadius *= pow(All.OmegaMatter/All.OmegaBaryon,1./NUMDIMS);}} /* correct (crudely) for baryon fraction, used in the estimate above for KernelRadius */

#ifdef SINK_PARTICLES
    if(RestartFlag==0 || RestartFlag==2) {for(i=0;i<NumPart;i++) {if(P[i].Type == 5) {P[i].KernelRadius = All.ForceSoftening[P[i].Type];}}}
#endif

#ifdef GRAIN_FLUID
    if(RestartFlag==0 || RestartFlag==2) {for(i=0;i<NumPart;i++) {P[i].KernelRadius *= pow(2.,1./NUMDIMS);}} /* very rough correction assuming comparable numbers of dust and gas elements */
#endif

    density();
}


void assign_unique_ids(void)
{
    int i, *numpartlist;
    MyIDType idfirst;

    numpartlist = (int *) mymalloc("numpartlist", NTask * sizeof(int));

    MPI_Allgather(&NumPart, 1, MPI_INT, numpartlist, 1, MPI_INT, MPI_COMM_WORLD);

    idfirst = 1;

    for(i = 0; i < ThisTask; i++)
        idfirst += numpartlist[i];

    for(i = 0; i < NumPart; i++)
    {
        P[i].ID = idfirst;
        idfirst++;
    }

    myfree(numpartlist);
}


#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
void ags_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            P[i].Particle_DivVel = 0;
            P[i].AGS_zeta = 0;
            if(ags_density_isactive(i) || P[i].Type==0) // type is AGS-active //
            {
                if(P[i].Type > 0)
                {
                    no = Father[i];
                    if(no < 0) {no = All.MaxPart;} /* orphaned by the treebuild; see init note above */
                    while(10 * All.AGS_DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                    {
                        p = Nodes[no].u.d.father;
                        if(p < 0) break;
                        no = p;
                    }
                    P[i].AGS_KernelRadius = 2. * pow(1.0/VOLUME_NORM_COEFF_FOR_NDIMS * All.AGS_DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                    double soft = All.ForceSoftening[P[i].Type];
                    if(soft != 0)
                    {
                        if((P[i].AGS_KernelRadius>1e6*soft)||(P[i].AGS_KernelRadius<=1e-3*soft)||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {P[i].AGS_KernelRadius = 1.e2*soft;} /* random guess to get things started here, thats all */
                    }
                } else {
                    P[i].AGS_KernelRadius = P[i].KernelRadius;
                }
            } else {
                P[i].AGS_KernelRadius = All.ForceSoftening[P[i].Type]; /* not AGS-active, use fixed softening */
            }
        }
    }
    ags_density();
#ifdef DM_FUZZY
    do_dm_fuzzy_initialization();
#endif
}
#endif // AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE


#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==2)
void disp_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            if(P[i].Type == 0)
            {
                no = Father[i];
                if(no < 0) {no = All.MaxPart;} /* orphaned by the treebuild; see init note above */
                while(10 * 2.0 * 64 * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0) {break;}
                    no = p;
                }
                CellP[i].KernelRadiusDM = pow(1.0/VOLUME_NORM_COEFF_FOR_NDIMS * 2.0 * 64 * P[i].Mass / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                double soft = All.ForceSoftening[P[i].Type];
                if(soft != 0) {if((CellP[i].KernelRadiusDM >1000.*soft)||(CellP[i].KernelRadiusDM<=0.01*soft)||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {CellP[i].KernelRadiusDM = soft;}} /* lower clause tests the DM kernel it resets -- it read the gas kernel, so a degenerate small KernelRadiusDM was never clamped */
            }
        }
    }
    if(ThisTask == 0) {printf("computing DM Vel_disp around gas particles.\n");}
    disp_density();
}
#endif
#endif


void test_id_uniqueness(void)
{
    double t0, t1;
#ifndef BOX_BND_PARTICLES
    int i;
    MyIDType *ids, *ids_first;
#endif

    if(ThisTask == 0)
    {
        printf("Testing ID uniqueness...\n");
    }

    if(NumPart == 0)
    {
        printf("need at least one particle per cpu\n");
        endrun(8);
    }

    t0 = my_second();

#ifndef BOX_BND_PARTICLES
    ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
    ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));

    for(i = 0; i < NumPart; i++)
        ids[i] = P[i].ID;

    parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

    for(i = 1; i < NumPart; i++)
        if(ids[i] == ids[i - 1])
        {
            printf("non-unique ID=%llu found on task=%d   (i=%d NumPart=%d)\n", (unsigned long long) ids[i], ThisTask, i, NumPart);
            endrun(12);
        }

    MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

    if(ThisTask < NTask - 1)
        if(ids[NumPart - 1] == ids_first[ThisTask + 1])
        {
            printf("non-unique ID=%llu found on task=%d\n", (unsigned long long) ids[NumPart - 1], ThisTask);
            endrun(13);
        }

    myfree(ids_first);
    myfree(ids);
#endif

    t1 = my_second();

    if(ThisTask == 0)
    {
        printf("success.  took=%g sec\n", timediff(t0, t1));
    }
}

int compare_IDs(const void *a, const void *b)
{
    if(*((MyIDType *) a) < *((MyIDType *) b)) {return -1;}
    if(*((MyIDType *) a) > *((MyIDType *) b)) {return +1;}
    return 0;
}
