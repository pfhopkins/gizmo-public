#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly initializing
 * new/modified variables, as needed)
 */

/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the initial gas kernel lengths are determined.
 */
void init(void)
{
    int i, j; double a3, atime, a2_fac;

#ifdef MAGNETIC
    double gauss2gizmo = All.UnitMagneticField_in_gauss / UNIT_B_IN_GAUSS;
    /* NOTE: we will always work -internally- in code units where MU_0 = 1; hence the 4pi here; [much simpler, but be sure of your conversions!] */
#endif

#ifdef BLACK_HOLES
    int count_holes = 0;
#endif

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();

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
                char fname[1000];
                if(All.NumFilesPerSnapshot > 1) {sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase, RestartSnapNum);}
                    else {sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);}
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
    All.SnapshotFileCount = 0;
    if(RestartFlag == 2)
    {
        if(RestartSnapNum < 0)
        {
            char *underscore = strrchr(All.InitCondFile, '_');
            if(!underscore)
            {
                char buf[1000];
                sprintf(buf, "Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n", All.InitCondFile);
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
    All.TreeAllocFactor = 0.45; /* this will also iteratively increase to fit the particle distribution */
    /* To construct the BH-tree for N particles, somewhat less than N
     internal tree-nodes are necessary for ‘normal’ particle distributions.
     TreeAllocFactor sets the number of internal tree-nodes allocated in units of the particle number.
     By experience, space for ~0.65N internal nodes is usually fully sufficient for typical clustered
     particle distributions, so a value of 0.7 should put you on the safe side. If the employed particle
     number per processor is very small (less than a thousand or so), or if there are many particle pairs
     with identical or nearly identical coordinates, a higher value may be required. Since the number of
     particles on a given processor may be higher by a factor PartAllocFactor than the average particle
     number, the total amount of memory requested for the BH tree on a single processor scales proportional
     to PartAllocFactor*TreeAllocFactor. */



#ifdef BOX_PERIODIC
    if(All.ComovingIntegrationOn) {check_omega();}
#endif
    All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#if (defined(BLACK_HOLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
    All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif

    for(i = 0; i < GRAVCOSTLEVELS; i++) {All.LevelToTimeBin[i] = 0;}

    for(i = 0; i < NumPart; i++) {for(j = 0; j < GRAVCOSTLEVELS; j++) {P[i].GravCost[j] = 0;}}

    if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
        {for(i=0;i<NumPart;i++) {for(j=0;j<3;j++) {P[i].Vel[j] *= sqrt(All.Time)*All.Time;}}}

#ifdef DM_SIDM
    init_self_interactions();
#endif


    for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
        for(j = 0; j < 3; j++) {P[i].GravAccel[j] = 0;}

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* init tidal tensor for first output (not used for calculation) */
        P[i].tidal_tensorps[0][0]=P[i].tidal_tensorps[0][1]=P[i].tidal_tensorps[0][2]=0;
        P[i].tidal_tensorps[1][0]=P[i].tidal_tensorps[1][1]=P[i].tidal_tensorps[1][2]=0;
        P[i].tidal_tensorps[2][0]=P[i].tidal_tensorps[2][1]=P[i].tidal_tensorps[2][2]=0;
#ifdef PMGRID
        P[i].tidal_tensorpsPM[0][0]=P[i].tidal_tensorpsPM[0][1]=P[i].tidal_tensorpsPM[0][2]=0;
        P[i].tidal_tensorpsPM[1][0]=P[i].tidal_tensorpsPM[1][1]=P[i].tidal_tensorpsPM[1][2]=0;
        P[i].tidal_tensorpsPM[2][0]=P[i].tidal_tensorpsPM[2][1]=P[i].tidal_tensorpsPM[2][2]=0;
#endif
#endif
#ifdef GDE_DISTORTIONTENSOR
        /* find caustics by sign analysis of configuration space distortion */
        P[i].last_determinant = 1.0;
#ifdef OUTPUT_GDE_LASTCAUSTIC
        P[i].lc_Time = 0.0; /* all entries zero -> no caustic yet */
        P[i].lc_Pos[0] = 0.0; P[i].lc_Pos[1] = 0.0; P[i].lc_Pos[2] = 0.0;
        P[i].lc_Vel[0] = 0.0; P[i].lc_Vel[1] = 0.0; P[i].lc_Vel[2] = 0.0;
        P[i].lc_rho_normed_cutoff = 0.0;
        P[i].lc_Dir_x[0] = 0.0; P[i].lc_Dir_x[1] = 0.0; P[i].lc_Dir_x[2] = 0.0;
        P[i].lc_Dir_y[0] = 0.0; P[i].lc_Dir_y[1] = 0.0; P[i].lc_Dir_y[2] = 0.0;
        P[i].lc_Dir_z[0] = 0.0; P[i].lc_Dir_z[1] = 0.0; P[i].lc_Dir_z[2] = 0.0;
        P[i].lc_smear_x = 0.0; P[i].lc_smear_y = 0.0; P[i].lc_smear_z = 0.0;
#endif
        for(i1 = 0; i1 < 6; i1++) {for(i2 = 0; i2 < 6; i2++) {if(i1 == i2) {P[i].distortion_tensorps[i1][i2] = 1.0;} else {P[i].distortion_tensorps[i1][i2] = 0.0;}}}
        if(All.ComovingIntegrationOn) /* for cosmological simulations we do init here, not read from ICs */
        {
#ifndef GDE_READIC
            P[i].caustic_counter = 0.0; /* no caustic passages in the beginning */
#ifndef GDE_LEAN
            P[i].a0 = All.TimeBegin; /* Lagrange time of particle */
            /* approximation: perfect Hubble Flow -> peculiar sheet orientation is exactly zero */
            for(i1 = 0; i1 < 3; i1++) {for(i2 = 0; i2 < 3; i2++) {GDE_VMATRIX(i,i1,i2) = 0.0;}}
            /* approximation: initial stream density equals background density */
            P[i].init_density = All.OmegaMatter * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
#else
            All.GDEInitStreamDensity = All.OmegaMatter * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
#endif
#endif
        }
#ifndef GDE_LEAN
        /* annihilation stuff */
        P[i].s_1_last = 1.0; P[i].s_2_last = 1.0; P[i].s_3_last = 1.0; P[i].second_deriv_last = 0.0; P[i].rho_normed_cutoff_last = 1.0;
        P[i].s_1_current = 1.0; P[i].s_2_current = 1.0; P[i].s_3_current = 1.0; P[i].second_deriv_current = 0.0; P[i].rho_normed_cutoff_current = 1.0;
        P[i].annihilation = 0.0; P[i].analytic_caustics = 0.0; P[i].analytic_annihilation = 0.0;
#endif
        if(All.ComovingIntegrationOn) {P[i].stream_density = GDE_INITDENSITY(i) / (All.TimeBegin * All.TimeBegin * All.TimeBegin);} else {P[i].stream_density = GDE_INITDENSITY(i);}
#endif /* GDE_DISTORTIONTENSOR */

#ifdef ADAPTIVE_TREEFORCE_UPDATE
        P[i].time_since_last_treeforce = 0;
        P[i].tdyn_step_for_treeforce = 0;
#endif        
        

#ifdef KEEP_DM_HSML_AS_GUESS
        if(RestartFlag != 1)
            P[i].DM_Hsml = -1;
#endif

#ifdef PMGRID
        for(j = 0; j < 3; j++) {P[i].GravPM[j] = 0;}
#endif
        P[i].Ti_begstep = 0;
        P[i].Ti_current = (integertime)0;
        P[i].TimeBin = 0;

        if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS) {P[i].OldAcc = 0;}	/* Do not zero in 2lpt case as masses are stored here */

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
        P[i].Potential = 0;
#endif
#ifdef GALSF
        if(RestartFlag == 0)
        {
            P[i].StellarAge = 0;
#ifdef GALSF_SFR_IMF_VARIATION
            P[i].IMF_Mturnover = 2.0; /* gives a solar-type IMF for our calculations in current code */
#endif
#ifdef GALSF_SFR_IMF_SAMPLING
            P[i].IMF_NumMassiveStars = 0;
#endif
        }
#endif

        if(RestartFlag != 1)
        {
#if defined(DO_DENSITY_AROUND_STAR_PARTICLES)
            P[i].DensAroundStar = 0;
            P[i].GradRho[0]=0;
            P[i].GradRho[1]=0;
            P[i].GradRho[2]=1;
#endif
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
            P[i].SNe_ThisTimeStep = 0;
#endif
#ifdef GALSF_FB_MECHANICAL
            int k; for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {P[i].Area_weighted_sum[k] = 0;}
#endif
        }

#if defined(INIT_STELLAR_METALS_AGES_DEFINED) && defined(GALSF)
        if(RestartFlag == 0) {P[i].StellarAge = -2.0 * All.InitStellarAgeinGyr / (UNIT_TIME_IN_GYR) * get_random_number(P[i].ID + 3);}
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
                int k, non_gdir=1; double A[3]={0}, B[3]={0}, A_cross_B[3]={0}, amag, rho_gas_expected, acc_ang=All.Vertical_Grain_Accel_Angle * M_PI/180., tS0, a0, ct=1, tau2=0, ct2=0, w0, agamma=9.*M_PI/128.; B[2]=1; if(GRAV_DIRECTION_RDI==1) {non_gdir=2;}
                rho_gas_expected = 1; /* guess for the gas density here [set custom for e.g. stratified problems */
                tS0 = 0.626657 * P[i].Grain_Size * sqrt(GAMMA_DEFAULT) / rho_gas_expected; /* stopping time [Epstein] for driftvel->0 */
                A[GRAV_DIRECTION_RDI]=cos(acc_ang)*All.Vertical_Grain_Accel - All.Vertical_Gravity_Strength; A[0]=sin(acc_ang)*All.Vertical_Grain_Accel; /* define angles/direction of external acceleration */
                amag=sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+MIN_REAL_NUMBER); A[0]/=amag; A[1]/=amag; A[2]/=amag;
                a0 = tS0 * amag / (1.+All.Dust_to_Gas_Mass_Ratio); /* acc * tS0 / (1+mu) */
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
                a0 *= All.Grain_Size_Max / P[i].Grain_Size;
#endif
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                double q_a = (0.75*All.Grain_Q_at_MaxGrainSize) / (All.Grain_Internal_Density*All.Grain_Size_Max), kappa_0 = All.Grain_Absorbed_Fraction_vs_Total_Extinction * q_a * All.Dust_to_Gas_Mass_Ratio;
                double rho_base_setup = 1., H_scale_setup = 1.; // define in code units the -assumed- initial scaling of the base gas density and vertical scale-length (PROBLEM SPECIFIC HERE!)
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
                kappa_0 *= sqrt(All.Grain_Size_Max / All.Grain_Size_Min); // opacity must be corrected for dependence of Q on grainsize or lack thereof
#endif
                a0 *= exp(-kappa_0*rho_base_setup*H_scale_setup*(1.-exp(-P[i].Pos[2]/H_scale_setup))); // attenuate incident flux (and reduce acceleration) according to equilibrium expectation, if we're using single-scattering radiation pressure [otherwise comment this line out] //
#endif
                w0=sqrt((sqrt(1.+4.*agamma*a0*a0)-1.)/(2.*agamma)); // exact solution if no Lorentz forces and Epstein drag //
#ifdef GRAIN_LORENTZFORCE
                double Bmag, tL_i=0, tau2_0=0, f_tau_guess2=0; B[0]=All.BiniX; B[1]=All.BiniY; B[2]=All.BiniZ; Bmag=sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]); B[0]/=Bmag; B[1]/=Bmag; B[2]/=Bmag;
                tL_i = (All.Grain_Charge_Parameter/All.Grain_Size_Max) * pow(All.Grain_Size_Max/P[i].Grain_Size,2) * Bmag; // 1/Lorentz in code units
                ct=A[0]*B[0]+A[1]*B[1]+A[2]*B[2]; ct2=ct*ct; tau2_0=pow(tS0*tL_i,2); // variables for below //
                for(k=0;k<20;k++)
                {
                   tau2 = tau2_0 / (1. + agamma*w0*w0); // guess tau [including velocity dependence] //
                   f_tau_guess2 = (1.+tau2*ct2) / (1.+tau2); // what the projection factor (reduction in w from projection) would be //
                   w0 = sqrt((sqrt(1.+4.*agamma*a0*a0*f_tau_guess2)-1.)/(2.*agamma)); // re-calculate w0 with this //
                }
#endif
                w0 /= sqrt((1.+tau2)*(1.+tau2*ct2)); // ensures normalization to unity with convention below //
                A_cross_B[0]=A[1]*B[2]-A[2]*B[1]; A_cross_B[1]=A[2]*B[0]-A[0]*B[2]; A_cross_B[2]=A[0]*B[1]-A[1]*B[0];
                for(k=0;k<3;k++) {P[i].Vel[k]=w0*(A[k] + sqrt(tau2)*A_cross_B[k] + tau2*ct*B[k]);}
#ifdef BOX_SHEARING
                // now add linearly the NHS drift solution for our shearing box setup
                double v00 = -All.Pressure_Gradient_Accel / (2. * BOX_SHEARING_OMEGA_BOX_CENTER);
                double v_K = -(P[i].Pos[0]-boxHalf_X) * BOX_SHEARING_Q*BOX_SHEARING_OMEGA_BOX_CENTER;
                double tau_s = tS0 * P[i].Grain_Size * BOX_SHEARING_OMEGA_BOX_CENTER;
                v00 /= (1. + tau_s*tau_s); // appears in both terms here //
                P[i].Vel[0] += v00 * 2.*tau_s; // radial drift
                P[i].Vel[BOX_SHEARING_PHI_COORDINATE] = v_K + v00; // azimuthal drift relative to keplerian frame
#endif
#endif // closes rdi_testproblem
            }
            P[i].Gas_Density = P[i].Gas_InternalEnergy = P[i].Gas_Velocity[0]=P[i].Gas_Velocity[1]=P[i].Gas_Velocity[2]=0; P[i].Grain_AccelTimeMin = MAX_REAL_NUMBER;
#if defined(GRAIN_BACKREACTION)
            P[i].Grain_DeltaMomentum[0]=P[i].Grain_DeltaMomentum[1]=P[i].Grain_DeltaMomentum[2]=0;
#endif
#if defined(GRAIN_LORENTZFORCE)
            P[i].Gas_B[0]=P[i].Gas_B[1]=P[i].Gas_B[2]=0;
#endif
        } // closes check on restartflag and particle type
#endif // closes grain_fluid



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
#endif

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



#ifdef BLACK_HOLES
#ifdef BH_WAKEUP_GAS
	    if(P[i].Type == 0) {P[i].LowestBHTimeBin = TIMEBINS;}
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
        P[i].BH_Ngb_Flag = 0;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
	    P[i].min_bh_approach_time = P[i].min_bh_freefall_time = MAX_REAL_NUMBER;
#if (SINGLE_STAR_TIMESTEPPING > 0)
	    P[i].SuperTimestepFlag = 0;
#endif
#endif
        if(P[i].Type == 5)
        {
            count_holes++;
            if(RestartFlag == 0)
            {
                BPP(i).BH_Mass = All.SeedBlackHoleMass;
#ifdef SINGLE_STAR_SINK_DYNAMICS
                BPP(i).BH_Mass = P[i].Mass;
#endif
#ifdef GRAIN_FLUID
                BPP(i).BH_Dust_Mass = 0;
#endif
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                BPP(i).SinkRadius = All.ForceSoftening[5];
#endif
#ifdef BH_ALPHADISK_ACCRETION
                BPP(i).BH_Mass_AlphaDisk = All.SeedAlphaDiskMass;
#endif
#ifdef BH_FOLLOW_ACCRETED_ANGMOM
                double bh_mu=2*get_random_number(P[i].ID+3)-1, bh_phi=2*M_PI*get_random_number(P[i].ID+4), bh_sin=sqrt(1-bh_mu*bh_mu);
                double spin_prefac = All.G * P[i].BH_Mass / C_LIGHT_CODE; // assume initially maximally-spinning BH with random orientation
                P[i].BH_Specific_AngMom[0]=spin_prefac*bh_sin*cos(bh_phi); P[i].BH_Specific_AngMom[1]=spin_prefac*bh_sin*sin(bh_phi); P[i].BH_Specific_AngMom[2]=spin_prefac*bh_mu;
#endif
#ifdef BH_COUNTPROGS
                BPP(i).BH_CountProgs = 1;
#endif
            }
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
	    P[i].dt_since_last_gas_search = 0;
	    P[i].do_gas_search_this_timestep = 1;
#endif 
        }
#endif
    }

#ifdef BLACK_HOLES
    MPI_Allreduce(&count_holes, &All.TotBHs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    for(i = 0; i < TIMEBINS; i++) {TimeBinActive[i] = 1;}

    reconstruct_timebins();

#ifdef PMGRID
    All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

    for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;

        for(j = 0; j < 3; j++)
        {
            SphP[i].VelPred[j] = P[i].Vel[j];
            SphP[i].HydroAccel[j] = 0;
            //SphP[i].dMomentum[j] = 0;//manifest-indiv-timestep-debug//
        }

        //SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
        P[i].Particle_DivVel = 0;
        SphP[i].ConditionNumber = 1;
        SphP[i].DtInternalEnergy = 0;
        SphP[i].FaceClosureError = 0;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        SphP[i].MaxKineticEnergyNgb = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[i].dMass = 0;
        SphP[i].DtMass = 0;
        SphP[i].MassTrue = P[i].Mass;
        for(j=0;j<3;j++) SphP[i].GravWorkTerm[j] = 0;
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_HSML_CALCULATION_IS_ACTIVE)
        PPPZ[i].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        if(1 & ADAPTIVE_GRAVSOFT_FORALL) {PPP[i].AGS_Hsml = PPP[i].Hsml;} else {PPP[i].AGS_Hsml = All.ForceSoftening[0];}
#endif
#endif

#ifdef CONDUCTION
        SphP[i].Kappa_Conduction = 0;
#endif
#ifdef MHD_NON_IDEAL
        SphP[i].Eta_MHD_OhmicResistivity_Coeff = 0;
        SphP[i].Eta_MHD_HallEffect_Coeff = 0;
        SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
#ifdef VISCOSITY
        SphP[i].Eta_ShearViscosity = 0;
        SphP[i].Zeta_BulkViscosity = 0;
#endif


#ifdef TURB_DIFFUSION
        SphP[i].TD_DiffCoeff = 0;

#ifdef TURB_DIFF_DYNAMIC
        int u, v; /* start with the standard Smagorinsky-Lilly constant from Kolmogorov theory */
        SphP[i].TD_DynDiffCoeff = 0.01;
        SphP[i].h_turb = 0;
        SphP[i].FilterWidth_bar = 0;
        SphP[i].MagShear_bar = 0;
        SphP[i].MagShear = 0;
        SphP[i].Norm_hat = 0;
        SphP[i].Dynamic_numerator = 0;
        SphP[i].Dynamic_denominator = 0;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        SphP[i].TD_DynDiffCoeff_error = 0;
#endif
        for (u = 0; u < 3; u++) {
            if (RestartFlag != 7) {
                SphP[i].Velocity_bar[u] = 0;
                SphP[i].Velocity_hat[u] = 0;
            }
            for (v = 0; v < 3; v++) {
                SphP[i].VelShear_bar[u][v] = 0;
            }
        }
#endif
#endif

        if(RestartFlag == 0)
        {
#ifndef INPUT_READ_HSML
            PPP[i].Hsml = 0;
#endif
            SphP[i].Density = -1;
#ifdef COOLING
#ifndef CHIMES
            SphP[i].Ne = 1.0;
#endif
#if defined(COOL_MOLECFRAC_NONEQM)
            SphP[i].MolecularMassFraction = 0.0; SphP[i].MolecularMassFraction_perNeutralH = 0.0; // start atomic
#endif
#endif
#ifdef CHIMES_STELLAR_FLUXES
	    int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {SphP[i].Chimes_fluxPhotIon[kc] = 0; SphP[i].Chimes_G0[kc] = 0;}
#endif
#ifdef BH_COMPTON_HEATING
            SphP[i].Rad_Flux_AGN = 0;
#endif
        }
#ifdef GALSF_SUBGRID_WINDS
        if(RestartFlag == 0) {SphP[i].DelayTime = 0;}
#if (GALSF_SUBGRID_WIND_SCALING==1)
        SphP[i].HostHaloMass = 0;
#endif
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        SphP[i].DelayTimeCoolingSNe = 0;
#endif
#ifdef GALSF
        SphP[i].Sfr = 0;
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED)
        SphP[i].AlphaVirial_SF_TimeSmoothed = 0;
#endif
#endif
#ifdef COSMIC_RAY_FLUID
        if(RestartFlag == 0) {for(j=0;j<N_CR_PARTICLE_BINS;j++) {SphP[i].CosmicRayEnergy[j] = 0;}}
#endif
#ifdef MAGNETIC
#if defined MHD_B_SET_IN_PARAMS
        if(RestartFlag == 0)
        {			/* Set only when starting from ICs */
            SphP[i].B[0]=SphP[i].BPred[0] = All.BiniX;
            SphP[i].B[1]=SphP[i].BPred[1] = All.BiniY;
            SphP[i].B[2]=SphP[i].BPred[2] = All.BiniZ;
        }
#endif /*MHD_B_SET_IN_PARAMS*/
        for(j = 0; j < 3; j++)
        {
            SphP[i].BPred[j] *= a2_fac * gauss2gizmo;
            SphP[i].B[j] = SphP[i].BPred[j];
        }
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
        SphP[i].Balpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
        SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#ifdef BH_RETURN_BFLUX
        P[i].B[0] = P[i].B[1] = P[i].B[2] = 0;
#endif
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        SphP[i].alpha = 0.0;
#endif
#if defined(BH_THERMALFEEDBACK)
        SphP[i].Injected_BH_Energy = 0;
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
            SphP[i].VelPred[2] = 0;
            SphP[i].HydroAccel[2] = 0;
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
            SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0;
            SphP[i].HydroAccel[1] = SphP[i].HydroAccel[2] = 0;
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

#ifdef BH_WIND_SPAWN
    MaxUnSpanMassBH     = 0;
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

#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    if(RestartFlag != 3 && RestartFlag != 5) {ags_setup_smoothinglengths();}
#endif

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)
    if(RestartFlag != 3 && RestartFlag != 5) {disp_setup_smoothinglengths();}
#endif
#endif

#if defined GALSF_SFR_IMF_VARIATION
    for(i = 0; i < NumPart; i++) {P[i].IMF_Mturnover = 2.0;} // reset to normal IMF
#endif

#if defined(WAKEUP) && defined(AGS_HSML_CALCULATION_IS_ACTIVE)
    for(i=0;i<NumPart;i++) {P[i].wakeup=0;}
#endif


    /* HELLO! This here is where you should insert custom code for hard-wiring the ICs of various test problems */



    density();
    for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
        int k; k=0;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;

        // re-match the predicted and initial velocities and B-field values, just to be sure //
        for(j=0;j<3;j++) SphP[i].VelPred[j]=P[i].Vel[j];
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==0)
        for(j=0;j<3;j++) {SphP[i].ParticleVel[j] = 0;} // set these to zero and forget them, for the rest of the run //
#endif

#ifdef MAGNETIC
        for(j=0;j<3;j++) {SphP[i].B[j] = SphP[i].BPred[j] * P[i].Mass / SphP[i].Density;} // convert to the conserved unit V*B //
        for(j=0;j<3;j++) {SphP[i].BPred[j]=SphP[i].B[j]; SphP[i].DtB[j]=0;}
#endif
#ifdef COSMIC_RAY_FLUID
        for(k=0;k<N_CR_PARTICLE_BINS;k++)
        {
            SphP[i].CosmicRayEnergyPred[k]=SphP[i].CosmicRayEnergy[k]; SphP[i].CosmicRayDiffusionCoeff[k]=0; SphP[i].DtCosmicRayEnergy[k]=0;
#ifdef CRFLUID_M1
            for(j=0;j<3;j++) {SphP[i].CosmicRayFlux[k][j]=0; SphP[i].CosmicRayFluxPred[k][j]=0;}
#endif
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(j=0;j<2;j++) {SphP[i].CosmicRayAlfvenEnergy[k][j]=0; SphP[i].CosmicRayAlfvenEnergyPred[k][j]=0; SphP[i].DtCosmicRayAlfvenEnergy[k][j]=0;}
#endif
        }
#endif
#if defined(EOS_ELASTIC)
        if(RestartFlag != 1)
        {
            for(k=0;k<3;k++) {for(j=0;j<3;j++) {SphP[i].Dt_Elastic_Stress_Tensor[j][k] = SphP[i].Elastic_Stress_Tensor_Pred[j][k] = SphP[i].Elastic_Stress_Tensor[j][k] = 0;}}
        } else {
            for(k=0;k<3;k++) {for(j=0;j<3;j++) {SphP[i].Elastic_Stress_Tensor_Pred[j][k] = SphP[i].Elastic_Stress_Tensor[j][k]; SphP[i].Dt_Elastic_Stress_Tensor[j][k] = 0;}}
        }
#endif
        //SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
        SphP[i].DtInternalEnergy = 0;
#if defined(COOLING) && !defined(COOLING_OPERATOR_SPLIT)
        SphP[i].CoolingIsOperatorSplitThisTimestep = 1; /* default to more conservative split */
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[i].dMass = 0;
        SphP[i].DtMass = 0;
        SphP[i].MassTrue = P[i].Mass;
        for(j=0;j<3;j++) SphP[i].GravWorkTerm[j] = 0;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_HSML_CALCULATION_IS_ACTIVE)
        PPPZ[i].AGS_zeta = 0;
#endif
#ifdef WAKEUP
        if(RestartFlag!=0) {PPPZ[i].wakeup=0;}
        NeedToWakeupParticles = 0;
        NeedToWakeupParticles_local = 0;
#endif
#ifdef SUPER_TIMESTEP_DIFFUSION
        SphP[i].Super_Timestep_Dt_Explicit = 0;
        SphP[i].Super_Timestep_j = 0;
#endif
#ifdef BH_COMPTON_HEATING
        SphP[i].Rad_Flux_AGN = 0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {SphP[i].Rad_E_gamma[kf]=0;}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(j=0;j<3;j++) {SphP[i].Rad_Flux[kf][j]=0;}}}
#endif

#ifdef COOL_GRACKLE
        if(RestartFlag == 0)
        {
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            SphP[i].grHI    = HYDROGEN_MASSFRAC;
            SphP[i].grHII   = 1.0e-20;
            SphP[i].grHM    = 1.0e-20;
            SphP[i].grHeI   = 1.0 - HYDROGEN_MASSFRAC;
            SphP[i].grHeII  = 1.0e-20;
            SphP[i].grHeIII = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            SphP[i].grH2I   = 1.0e-20;
            SphP[i].grH2II  = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            SphP[i].grDI    = 2.0 * 3.4e-5;
            SphP[i].grDII   = 1.0e-20;
            SphP[i].grHDI   = 1.0e-20;
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
        for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
        {
            mass_tot += P[i].Mass;
            if(P[i].Mass > mass_max) mass_max = P[i].Mass;
            if(P[i].Mass < mass_min) mass_min = P[i].Mass;
        }
        /* broadcast this and get the min and max values over all processors */
        double mpi_mass_min,mpi_mass_max;
        MPI_Allreduce(&mass_min, &mpi_mass_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&mass_max, &mpi_mass_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        All.MinMassForParticleMerger = 0.49 * mpi_mass_min;
#ifdef SINGLE_STAR_SINK_DYNAMICS /* Get mean gas mass, used in various subroutiens */
        double mpi_mass_tot; long mpi_Ngas; long Ngas_l = (long) N_gas;
        MPI_Allreduce(&mass_tot, &mpi_mass_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Ngas_l, &mpi_Ngas, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        All.MeanGasParticleMass = mpi_mass_tot/( (double)mpi_Ngas );
#endif
#ifdef GALSF_GENERATIONS
        All.MinMassForParticleMerger /= (float)GALSF_GENERATIONS;
#endif
        All.MaxMassForParticleSplit  = 3.01 * mpi_mass_max;
#ifdef MERGESPLIT_HARDCODE_MAX_MASS
        All.MaxMassForParticleSplit = MERGESPLIT_HARDCODE_MAX_MASS;
#endif
#ifdef MERGESPLIT_HARDCODE_MIN_MASS
        All.MinMassForParticleMerger = MERGESPLIT_HARDCODE_MIN_MASS;
#endif
    }


#ifdef PM_HIRES_REGION_CLIPDM
    if(RestartFlag != 1)
    {
        double mpi_m_hires_max, m_hires_max=0.0;
        for(i=0; i<NumPart; i++) {if(P[i].Type==1) {if(P[i].Mass > m_hires_max) {m_hires_max=P[i].Mass;}}}
        MPI_Allreduce(&m_hires_max, &mpi_m_hires_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        All.MassOfClippedDMParticles = mpi_m_hires_max;
    }
#endif


    if(RestartFlag == 3)
    {
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
        if(ThisTask == 0) {printf("*AGS_HSML_CALCULATION_IS_ACTIVE* Computation of softening lengths... \n");}
        ags_setup_smoothinglengths();
        if(ThisTask == 0) {printf("*AGS_HSML_CALCULATION_IS_ACTIVE* Computation of softening lengths done. \n");}
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
        sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
        if(ThisTask == 0) {printf("Start writing file %s\n", All.SnapshotFileBase);}
        printf("RestartSnapNum %d\n", RestartSnapNum);

        All.TopNodeAllocFactor = 0.008;

        savepositions(RestartSnapNum);
        endrun(0);
    }
    

#if defined(COOL_MOLECFRAC_NONEQM)
    if(RestartFlag == 2) // should have read in SphP[i].MolecularMassFraction_perNeutralH
    {
        SphP[i].MolecularMassFraction_perNeutralH = DMIN(1,DMAX(0,SphP[i].MolecularMassFraction_perNeutralH));
        SphP[i].MolecularMassFraction = DMIN(1,DMAX(0, 1.-SphP[i].Ne/1.25)) * SphP[i].MolecularMassFraction_perNeutralH;
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
	      chimes_update_turbulent_abundances(i, 1);
#endif

	      chimes_update_gas_vars(i);

	      // Evolve the chemistry for (1 / nH) Myr (limited to 1 Gyr) ten times at fixed temperature.
	      ChimesGasVars[i].hydro_timestep = (ChimesFloat) DMIN(3.16e13 / ChimesGasVars[i].nH_tot, 3.16e16);
	      ChimesGasVars[i].ThermEvolOn = 0;

	      for (iter_number = 0; iter_number < 10; iter_number++) chimes_network(&(ChimesGasVars[i]), &ChimesGlobalVars);


#ifdef CHIMES_TURB_DIFF_IONS
	      chimes_update_turbulent_abundances(i, 1);
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
    if((RestartFlag == 0)||(RestartFlag==2)) // best for stability if we re-calc Hsml for snapshot restarts //
    {
#if defined(DO_DENSITY_AROUND_STAR_PARTICLES) || defined(GRAIN_FLUID)
        for(i = 0; i < NumPart; i++)
#else
        for(i = 0; i < N_gas; i++)
#endif
        {
                no = Father[i];
                while(2 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0) {break;}
                    no = p;
                }

                if((RestartFlag == 0)||(P[i].Type != 0)) // if Restartflag==2, use the saved Hsml of the gas as initial guess //
                {
#ifndef INPUT_READ_HSML
#if NUMDIMS == 3
                    PPP[i].Hsml = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.333333) * Nodes[no].len;
#endif
#if NUMDIMS == 2
                    PPP[i].Hsml = pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.5) * Nodes[no].len;
#endif
#if NUMDIMS == 1
                    PPP[i].Hsml = All.DesNumNgb * (P[i].Mass / Nodes[no].u.d.mass) * Nodes[no].len;
#endif
#ifndef SELFGRAVITY_OFF
                    if(All.SofteningTable[P[i].Type] != 0)
                    {
                        if((PPP[i].Hsml>100.*All.SofteningTable[P[i].Type])||(PPP[i].Hsml<=0.01*All.SofteningTable[P[i].Type])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0))
                            {PPP[i].Hsml = All.SofteningTable[P[i].Type];}
                    }
#else
                    if((Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {PPP[i].Hsml = All.SofteningTable[P[i].Type];}
#endif
#endif // INPUT_READ_HSML
                } // closes if((RestartFlag == 0)||(P[i].Type != 0))
            }
    }
    if((RestartFlag==0 || RestartFlag==2) && All.ComovingIntegrationOn) {for(i=0;i<N_gas;i++) {PPP[i].Hsml *= pow(All.OmegaMatter/All.OmegaBaryon,1./NUMDIMS);}} /* correct (crudely) for baryon fraction, used in the estimate above for Hsml */

#ifdef BLACK_HOLES
    if(RestartFlag==0 || RestartFlag==2) {for(i=0;i<NumPart;i++) {if(P[i].Type == 5) {PPP[i].Hsml = All.SofteningTable[P[i].Type];}}}
#endif

#ifdef GRAIN_FLUID
    if(RestartFlag==0 || RestartFlag==2) {for(i=0;i<NumPart;i++) {PPP[i].Hsml *= pow(2.,1./NUMDIMS);}} /* very rough correction assuming comparable numbers of dust and gas elements */
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


#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
void ags_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            P[i].Particle_DivVel = 0;
            PPPZ[i].AGS_zeta = 0;
            if(ags_density_isactive(i) || P[i].Type==0) // type is AGS-active //
            {
                if(P[i].Type > 0)
                {
                    no = Father[i];
                    while(10 * All.AGS_DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                    {
                        p = Nodes[no].u.d.father;
                        if(p < 0) break;
                        no = p;
                    }
                    PPP[i].AGS_Hsml = 2. * pow(1.0/NORM_COEFF * All.AGS_DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                    if(All.SofteningTable[P[i].Type] != 0)
                    {
                        if((PPP[i].AGS_Hsml>1e6*All.ForceSoftening[P[i].Type])||(PPP[i].AGS_Hsml<=1e-3*All.ForceSoftening[P[i].Type])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0))
                            PPP[i].AGS_Hsml = 1e2 * All.ForceSoftening[P[i].Type]; /* random guess to get things started here, thats all */
                    }
                } else {
                    PPP[i].AGS_Hsml = PPP[i].Hsml;
                }
            } else {
                PPP[i].AGS_Hsml = All.ForceSoftening[P[i].Type]; /* not AGS-active, use fixed softening */
            }
        }
    }
    ags_density();
#ifdef DM_FUZZY
    do_dm_fuzzy_initialization();
#endif
}
#endif // AGS_HSML_CALCULATION_IS_ACTIVE


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
                while(10 * 2.0 * 64 * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0) {break;}
                    no = p;
                }
                SphP[i].HsmlDM = pow(1.0/NORM_COEFF * 2.0 * 64 * P[i].Mass / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                if(All.SofteningTable[P[i].Type] != 0)
                {
                    if((SphP[i].HsmlDM >1000.*All.SofteningTable[P[i].Type])||(PPP[i].Hsml<=0.01*All.SofteningTable[P[i].Type])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {SphP[i].HsmlDM = All.SofteningTable[P[i].Type];}
                }
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
