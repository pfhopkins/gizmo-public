/*! \file analytic_gravity.h
 *  \brief externally-specified (analytic) gravity goes here
 *
 *  This file contains supplemental code if you want to add an 
 *   -analytic- potential or gravitational force in the code, 
 *   rather than solely relying on the calculated self-gravity. 
 *   Note that the terms here are added at the end of the self-gravity
 *   loop, so if you want to keep self-gravity, but add these, you need
 *   to make sure that your routine -adds to- the GravAccel values, rather 
 *   than re-setting them entirely.
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void add_analytic_gravitational_forces(void);
void GravAccel_set_zeros_if_needed(void);
void GravAccel_StaticPlummerSphere(void);
void GravAccel_StaticHernquist(void);
void GravAccel_StaticIsothermalSphere(void);
void GravAccel_KeplerianOrbit(void);
void GravAccel_KeplerianTestProblem(void);
void GravAccel_GrowingDiskPotential(void);
void GravAccel_StaticNFW(void);
void GravAccel_RayleighTaylorTest(void);
void GravAccel_ShearingSheet(void);
void GravAccel_PaczynskyWiita(void);
void GravAccel_RDITestProblem(void);
void GravAccel_GMCTurbInit(void);

/* parent routine which decides which (if any) analytic gravitational forces are applied */
void add_analytic_gravitational_forces()
{
    GravAccel_set_zeros_if_needed();    // initial book-keeping: make sure relevant terms are initialized/reset if needed

    /* now add the appropriate [if any] analytic gravitational forces */
#ifdef GRAVITY_ANALYTIC
    //GravAccel_RayleighTaylorTest();     // uniform vertical force for Rayleigh-Taylor-type tests
    //GravAccel_StaticPlummerSphere();    // plummer sphere
    //GravAccel_StaticHernquist();        // hernquist-profile sphere
    //GravAccel_StaticIsothermalSphere(); // singular (but finite) isothermal sphere
    //GravAccel_KeplerianOrbit();         // keplerian disk (2D or 3D)
    //GravAccel_KeplerianTestProblem();   // keplerian disk with special boundaries for test problem
    //GravAccel_GrowingDiskPotential();   // time-dependent (adiabatically growing) disk
    //GravAccel_StaticNFW();              // NFW profile sphere
    //GravAccel_PaczynskyWiita();         // Paczynsky-Wiita pseudo-Newtonian potential
#endif
#ifdef STARFORGE_GMC_TURBINIT
    GravAccel_GMCTurbInit();              // uniform sphere harmonic potential + r^-3 halo outside to confine stirred turbulent gas
#endif

#ifdef BOX_SHEARING
    GravAccel_ShearingSheet();            // adds coriolis and centrifugal terms for shearing-sheet approximation
#endif
#ifdef GRAIN_RDI_TESTPROBLEM
    GravAccel_RDITestProblem();           // vertical gravity+external acceleration for grain-RDI-wind tests
#endif
}



/* first if the 'self gravity off' options are enabled, we need to ensure the appropriate terms are disabled here */
void GravAccel_set_zeros_if_needed()
{
#if defined(SELFGRAVITY_OFF) || defined(RT_SELFGRAVITY_OFF) /* zero gravaccel [difference is that RT_SELFGRAVITY_OFF... option still computes everything above ]*/
    int i; for(i=FirstActiveParticle; i>=0; i=NextActiveParticle[i]) {P[i].GravAccel[0]=P[i].GravAccel[1]=P[i].GravAccel[2]=0;}
#if defined(COMPUTE_TIDAL_TENSOR_IN_GRAVTREE)
    for(i=FirstActiveParticle; i>=0; i=NextActiveParticle[i]) {P[i].tidal_tensorps[0][0]=P[i].tidal_tensorps[0][1]=P[i].tidal_tensorps[0][2]=P[i].tidal_tensorps[1][0]=P[i].tidal_tensorps[1][1]=P[i].tidal_tensorps[1][2]=P[i].tidal_tensorps[2][0]=P[i].tidal_tensorps[2][1]=P[i].tidal_tensorps[2][2]=0;}
#endif
#endif
}



/* external forces for dusty-box problem */
void GravAccel_RDITestProblem()
{
#ifdef GRAIN_RDI_TESTPROBLEM
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {   /* add the relevant vertical field for non-anchored particles */
        if(P[i].ID > 0 && (P[i].Type==0 || ((1 << P[i].Type) & (GRAIN_PTYPES))))
        {
#if defined(BOX_SHEARING) && (BOX_SHEARING != 4)
            double mu_g = All.Vertical_Gravity_Strength/(1.+All.Dust_to_Gas_Mass_Ratio); // unstratified box, work in compensated/free-falling frame here, with respect to vertical gravity //
            if(P[i].Type==0) {P[i].GravAccel[GRAV_DIRECTION_RDI]+=All.Dust_to_Gas_Mass_Ratio*mu_g;} else {P[i].GravAccel[GRAV_DIRECTION_RDI]-=mu_g;}
#else
            P[i].GravAccel[GRAV_DIRECTION_RDI] -= All.Vertical_Gravity_Strength; /* everything feels same vertical gravity */
#if defined(GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION) /* this is a hack for this problem to prevent the bottom boundary layer of gas from detachinng in a spurious way that prevents numerical flux from propagating in the z-direction */
            if(P[i].Type==0) {
                double h_exp = 0.05*(0.5+get_random_number(2*i+10+P[i].ID)) + 0.5 * (0.8*All.Vertical_Grain_Accel*All.Dust_to_Gas_Mass_Ratio) * All.Time*All.Time;
                if(P[i].Pos[GRAV_DIRECTION_RDI] < h_exp) {P[i].GravAccel[GRAV_DIRECTION_RDI] = All.Vertical_Gravity_Strength * pow(1. - P[i].Pos[GRAV_DIRECTION_RDI]/h_exp, 8) + 1.;}
                h_exp = 0.05*(0.5+get_random_number(2*i+10+P[i].ID)) + 0.5 * (1.0*All.Vertical_Grain_Accel*All.Dust_to_Gas_Mass_Ratio) * All.Time*All.Time;
                if(P[i].Pos[GRAV_DIRECTION_RDI] < h_exp) {P[i].GravAccel[GRAV_DIRECTION_RDI] = All.Vertical_Gravity_Strength * (1.+2.*get_random_number(i+2+P[i].ID) -  P[i].Pos[2]/(0.5*h_exp));}
            }
#endif
#endif
#ifdef BOX_SHEARING
            if(P[i].Type==0) {P[i].GravAccel[0] += All.Pressure_Gradient_Accel;} /* gas feels pressure gradient force in radial direction as well */
#endif
            double acc = All.Vertical_Grain_Accel;
#ifdef RT_OPACITY_FROM_EXPLICIT_GRAINS
            acc = 0; /* this is calculated separately, if this flag is on, from the explicitly-evolved radiation field */
#endif
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
            acc *= All.Grain_Size_Max / P[i].Grain_Size;
#endif
            if((1 << P[i].Type) & (GRAIN_PTYPES))
            {
                P[i].GravAccel[GRAV_DIRECTION_RDI] += acc * cos(All.Vertical_Grain_Accel_Angle * M_PI/180.);
                P[i].GravAccel[0] += acc * sin(All.Vertical_Grain_Accel_Angle * M_PI/180.);
            }
        }
    }
#endif
}



/* adds coriolis and centrifugal terms for shearing-sheet approximation */
void GravAccel_ShearingSheet()
{
#ifdef BOX_SHEARING
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        /* centrifugal force term (depends on distance from box center) */
        P[i].GravAccel[0] += 2.*(P[i].Pos[0]-boxHalf_X) * BOX_SHEARING_Q*BOX_SHEARING_OMEGA_BOX_CENTER*BOX_SHEARING_OMEGA_BOX_CENTER;
        /* coriolis force terms */
        double vp=0; if(P[i].Type==0) {vp=SphP[i].VelPred[BOX_SHEARING_PHI_COORDINATE];} else {vp=P[i].Vel[BOX_SHEARING_PHI_COORDINATE];}
        P[i].GravAccel[0] += 2.*vp * BOX_SHEARING_OMEGA_BOX_CENTER;
        if(P[i].Type==0) {vp=SphP[i].VelPred[0];} else {vp=P[i].Vel[0];}
        P[i].GravAccel[BOX_SHEARING_PHI_COORDINATE] -= 2.*vp * BOX_SHEARING_OMEGA_BOX_CENTER;
#if (BOX_SHEARING==4) /* add vertical gravity to the force law */
        P[i].GravAccel[2] -= BOX_SHEARING_OMEGA_BOX_CENTER * BOX_SHEARING_OMEGA_BOX_CENTER * (P[i].Pos[2]-boxHalf_Z);
#endif
    }
#endif
}



/* constant vertical acceleration for Rayleigh-Taylor test problem */
void GravAccel_RayleighTaylorTest()
{
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        {if(P[i].ID != 0) {P[i].GravAccel[1]=-0.5;}} /* now add the constant vertical field */
}



/* static unit Plummer Sphere (assumes G=M=a=1) */
void GravAccel_StaticPlummerSphere()
{
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#ifdef GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE
        for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        double r2, r; r2 = dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2]; r = sqrt(r2);
        for(k=0;k<3;k++) {P[i].GravAccel[k] += -dp[k] / pow(r2 + 1, 1.5);}
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        double f=pow(1+r2, 1.5), f2=pow(1+r2, 2.5);
        for(k=0;k<3;k++) {P[i].tidal_tensorps[k][k]-=1/f; int j; for(j=0;j<3;j++) {P[i].tidal_tensorps[k][j]+=3*dp[k]*dp[j]/f2;}}
#endif
    }
}



/* static Hernquist Profile (parameters specified in the routine below) */
void GravAccel_StaticHernquist()
{
    double HQ_Mtot=100, HQ_a=20; /* total mass and scale-length "a" [both in code units] */
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#ifdef GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE
        for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        double r2 = dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2], r = sqrt(r2), f = r+HQ_a, m = HQ_Mtot*(r/f)*(r/f);
        for(k=0;k<3;k++) {P[i].GravAccel[k] += -All.G * m * dp[k]/(r2*r);}
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        double f0=All.G*HQ_Mtot, fa=f0*(2/f+1/r)/(r2*f*f), fxx=-f0/(r*f*f);
        for(k=0;k<3;k++) {P[i].tidal_tensorps[k][k]+=fxx; int j; for(j=0;j<3;j++) {P[i].tidal_tensorps[k][j]+=fa*dp[k]*dp[j];}}
#endif
    }
}



/* static singular Isothermal Sphere Profile (parameters specified in the routine below) */
void GravAccel_StaticIsothermalSphere()
{
    double ISO_Mmax=100, ISO_Rmax=200; /* total mass inside rmax, the maximum radius with mass (outside of which density=0, just set Rmax very large if you want an infinite SIS) */
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#ifdef GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE
        for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        double r2 = dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2], r = sqrt(r2);
        double m = ISO_Mmax; if(r < ISO_Rmax) {m *= r/ISO_Rmax;} /* mass enclosed ~r, until Rmax, where it cuts off and remains constant */
        for(k=0;k<3;k++) {P[i].GravAccel[k] += -All.G * m * dp[k]/(r2*r);}
    }
}

/* potential of a uniform sphere of mass M and radius R, plus a r^-3 density profile outside for a gentle infinite confining potential - used for initializing turbulence in isolated spheres */
void GravAccel_GMCTurbInit()
{
#ifdef STARFORGE_GMC_TURBINIT
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k] - 0.5*All.BoxSize;}
        double r2 = dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2], r = sqrt(r2);
	double M = 0.808 * All.TotN_gas * All.MeanGasParticleMass, R=All.BoxSize/10; // these are for the default settings of MakeCloud's uniform sphere IC, adjust for your problem!
	double menc = DMIN(M,M*pow(r/R,3)) + DMAX(0,3*M*log(r/R)); // uniform sphere plus a r^-3 surrounding halo with density matched at the sphere radius
        for(k=0;k<3;k++) {P[i].GravAccel[k] += -All.G * STARFORGE_GMC_ALPHA * menc * dp[k]/(r2*r);}
    }
#endif
}



/* time-dependent potential of an adiabatically-growing disk */
void GravAccel_GrowingDiskPotential()
{
    int n_table = 14; // number of table entries below (must match!)
    // scale factor for cosmological runs (must be in monotonic increasing order!)
    double t_disk_table[14] = {0.2, 0.250, 0.266, 0.285, 0.308, 0.333, 0.363, 0.400, 0.444, 0.500, 0.572, 0.667, 0.800, 1.000};
    // m12i parameters: from Shea's fits:
    double m_disk_table[14] = {0.0, 0.061, 0.088, 0.117, 0.153, 0.223, 0.348, 0.429, 0.581, 1.118, 2.004, 3.008, 4.403, 6.001}; // disk mass in code units
    double r_disk_table[14] = {1.0, 5.071, 7.513, 6.787, 6.162, 3.277, 4.772, 3.964, 3.418, 2.511, 2.463, 1.503, 1.005, 1.150}; // disk scale length in code units
    double z_disk_table[14] = {1.0, 4.185, 8.971, 5.089, 3.532, 3.057, 4.557, 2.117, 1.828, 0.809, 0.217, 0.148, 0.335, 0.404}; // disk scale height in code units
    /* before the particle loop, interpolate the relevant quantities to the simulation time */
    double t=All.Time, dt=0, r2, Zterm, Rterm, Rterm2, myfacR; int i, i0=0, i1=0, k;
    if(t<=t_disk_table[0])
    {
        i0=i1=0;
    } else if(t>=t_disk_table[n_table-1]) {
        i0=i1=n_table-1;
    } else {
        for(k=1;k<n_table;k++) {if(t_disk_table[k] > t) {i1=k; break;}}
        i0=i1-1; dt=(t - t_disk_table[i0])/(t_disk_table[i1]-t_disk_table[i0]);
    }
    double m_disk = m_disk_table[i0] + dt * (m_disk_table[i1]-m_disk_table[i0]);
    double r_disk = r_disk_table[i0] + dt * (r_disk_table[i1]-r_disk_table[i0]);
    double z_disk = z_disk_table[i0] + dt * (z_disk_table[i1]-z_disk_table[i0]);
    /* ok now we can assign actual accelerations */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#ifdef GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE
        for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1];
        Zterm = sqrt(z_disk*z_disk + dp[2]*dp[2]); /* sqrt((Zdisk^2 + dZ^2); appears several times  */
        Rterm = r_disk + Zterm; Rterm2 = sqrt(r2 + Rterm*Rterm); Rterm2 = Rterm2*Rterm2*Rterm2;
        myfacR = -All.G * m_disk / Rterm2; /* has units s^-2, so  multiply by length to get accel.  no sign; handle that in min_xyz_to_bh */
        /* remember, min_xyz_to_bh = x_BH - myx => positive if x_BH > myx => acceleration is in positive x if x_BH > myx, which is correct (attractive) */
        P[i].GravAccel[0] += myfacR * dp[0]; P[i].GravAccel[1] += myfacR * dp[1];
        P[i].GravAccel[2] += myfacR * dp[2] * Rterm/Zterm; // this has units of:  M*L^3*M^-1*T^-2*L^2*L^-1*L^-3 = L/T^2
    }
}



/* Keplerian forces (G=M=1): useful for orbit, MRI, planetary disk problems */
void GravAccel_KeplerianOrbit()
{
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE)
        for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#elif defined(BOX_PERIODIC)
        dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y; dp[2] -= boxHalf_Z;
#endif
        double r2 = dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2], r = sqrt(r2);
        for(k=0;k<3;k++) {P[i].GravAccel[k] -= dp[k] / (r2*r);}
    }
}



/* Keplerian forces (G=M=1): this is a specific (bounded and softened) version 
 used just for the Keplerian disk test problem */
void GravAccel_KeplerianTestProblem()
{
    double x00=4.0, y00=4.0; /* 2D center of orbit: the is hard-coded for the relevant test problem */
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double r = pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),0.5);
        if((r > 0.35)&(r < 2.1))
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
        if(r <= 0.35)
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00)*pow(r/0.35,2) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[0] += +(P[i].Pos[0]-x00)*(0.35-r)/0.35 / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00)*pow(r/0.35,2) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] += +(P[i].Pos[1]-y00)*(0.35-r)/0.35 / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
        if(r >= 2.1)
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00)*(1+(r-2.1)/0.1) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00)*(1+(r-2.1)/0.1) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
    }
}



/* static NFW potential (parameters set below) */
void GravAccel_StaticNFW()
{
    double NFW_M200=100, NFW_C=10; /* NFW mass inside R200 (in code units), and concentration =R200/Rs */
    double R200 = pow(NFW_M200*All.G/(100.*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits), 1./3.), Rs=R200/NFW_C; /* using R200 = R where mean density = 200x critical density, and Rs=R200/c200 */
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#ifdef GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE
        for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#elif defined(BOX_PERIODIC)
        dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y; dp[2] -= boxHalf_Z;
#endif
        double r2=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2], r=sqrt(r2), x=r/Rs, cfac=log(1+NFW_C)-NFW_C/(1+NFW_C);
        if(r>0) {
            double mfac = (log(1+x)-x/(1+x)) / (x*x); if(x<=0.04) {mfac=0.5-2.*x/3.+0.75*x*x;} /* expression works well for larger x, small-x leads to potential numerical errors */
            for(k=0;k<3;k++) {P[i].GravAccel[k] += -All.G * mfac * NFW_M200/(cfac*Rs*Rs) * (dp[k]/r);}}
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //
}



/* Paczysnky Wiita pseudo-Newtonian potential, G = M_sol = c = 1 */
void GravAccel_PaczynskyWiita()
{
    double PACZYNSKY_WIITA_MASS = 1.0; // Mass to use for the Paczynksy-Wiita analytic gravity pseudo-Newtonian potential (in solar masses)
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#ifdef GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE
        for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        double r2 = dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2], r = sqrt(r2), r_g = 2*PACZYNSKY_WIITA_MASS;
        if(r > r_g)
        {
            double q = PACZYNSKY_WIITA_MASS/((r - r_g)*(r - r_g));
            for(k=0;k<3;k++) {P[i].GravAccel[k] -= q * P[i].Pos[k]/r;}
        }
    }
}



#ifdef PARTICLE_EXCISION
void apply_excision(void)
{
    double EXCISION_MASS = 0; // mass of the excised object. Used to move the excision boundary so as to capture bound objects. If zero the excision boundary will not move
    double EXCISION_INIT_RADIUS = 0; // initial excision radius
    double EXCISION_ETA = 1; // remove particles with radius < EXCISION_ETA R_excision
    double excision_radius = EXCISION_ETA * pow(EXCISION_INIT_RADIUS*EXCISION_INIT_RADIUS*EXCISION_INIT_RADIUS +
                                                3.*sqrt(2. * All.G * EXCISION_MASS) * pow(EXCISION_INIT_RADIUS, 3./2.) * All.Time +
                                                9./2. * All.G * EXCISION_MASS * All.Time*All.Time, 1./3.);
    int i,k; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0)
        {
            double dp[3]; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];}
#ifdef GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE
            for(k=0;k<3;k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
            double r2 = dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2], r = sqrt(r2);
            if(r < excision_radius) {P[i].Mass = 0;}
        }
    }
}
#endif

