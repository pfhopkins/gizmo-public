/*! \file mesh_motion.h
 *  \brief externally-specified mesh motion goes here
 *
 *  This file contains supplemental code if you want to add an
 *   arbitrary mesh motion to the code, where the mesh-generating
 *   points move with velocities set according to any user-specified function
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#ifdef HYDRO_MESHLESS_FINITE_VOLUME

// define the routines which are given below //
void set_mesh_motion(int i);
void MeshMotion_FixedGrid(int i);
void MeshMotion_KeplerianOrbit(int i);
void MeshMotion_CircularOrbitExternalGravity(int i);
void MeshMotion_FreeFallExternalGravity(int i);
void MeshMotion_UniformExpansion(int i);
void MeshMotion_UniformCollapse(int i);
void MeshMotion_ShearingSheet(int i);



void set_mesh_motion(int i)
{
#if (HYDRO_FIX_MESH_MOTION==4)
    MeshMotion_FixedGrid(i);      // non-moving
    //MeshMotion_KeplerianOrbit(i); // circular (cylindrical) Keplerian orbit, e.g. disk-type geometries (as DISCO code)
    //MeshMotion_CircularOrbitExternalGravity(i); // mesh moves on circular orbits according to the gravitational acceleration (can be specified analytically) //
    //MeshMotion_FreeFallExternalGravity(i); // alternatively, assume the mesh is free-accelerating in the potential //
    //MeshMotion_UniformExpansion(i); // uniformly-expanding mesh (for e.g. hubble-flow, outflow, or explosion-type geometries) //
    //MeshMotion_UniformCollapse(i); // uniformly-contracting mesh (for e.g. gravitational collapse or accretion-type geometries) //
#ifdef BOX_SHEARING
    MeshMotion_ShearingSheet(i); // mesh moves with equilibrium shear-flow in the shearing-sheet approximation //
#endif
#endif
    return;
}



/* no mesh motion (simple but here anyways) */
void MeshMotion_FixedGrid(int i)
{
    int k; for(k=0;k<3;k++) {SphP[i].ParticleVel[k]=0;}
}



/* Mesh moves with Keplerian speed (G=M=1) in a cylindrical mid-plane configuration:
    This is particularly useful for orbit, MRI, and disk problems.
    This is similar to the functionality in the DISCO code by Paul Duffel.
    The shearing motion of the background is folded into the smooth mesh motion. */
void MeshMotion_KeplerianOrbit(int i)
{
    double dp[3]; dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE)
    int k; for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#elif defined(BOX_PERIODIC)
    dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y; dp[2] -= boxHalf_Z;
#endif
    double r2 = dp[0]*dp[0] + dp[1]*dp[1], r = sqrt(r2), Omega=sqrt(1./(r2*r)); // orbital frequency
    SphP[i].ParticleVel[0] = Omega * dp[1];
    SphP[i].ParticleVel[1] = -Omega * dp[0];
    SphP[i].ParticleVel[2] = 0;
}



/* Mesh moves with the circular velocity according to whatever the background gravitational acceleration is. This is
    useful for setting the mesh to correspond to circular orbits in arbitrary, complicated background potentials.
    The variable 'GravAccel' is used for this, which is calculated self-consistently by the code under normal self-gravity,
    or can be set by the analytic gravity module. So setting this, with the analytic Keplerian potential, for example, will
    accomplish the same thing as the KeplerianOrbit routine above. Note that for a given acceleration, any the circular velocity
    can point anywhere in a plane perpendicular to the acceleration vector. To break this degeneracy, we assume the circular
    orbits point in the r x z direction (the same convention as used for the Keplerian disk above). */
void MeshMotion_CircularOrbitExternalGravity(int i)
{
    double a, dp[3], r, r2, Omega; int k;
    dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE)
    for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#elif defined(BOX_PERIODIC)
    dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y;
#endif
    r2 = dp[0]*dp[0] + dp[1]*dp[1]; r = sqrt(r2);
    a=0; for(k=0;k<3;k++) {a+=P[i].GravAccel[k]*P[i].GravAccel[k];}
    if(a > 0)
    {
        Omega = sqrt(a/r); // orbital frequency
        SphP[i].ParticleVel[0] = Omega * dp[1];
        SphP[i].ParticleVel[1] = -Omega * dp[0];
        SphP[i].ParticleVel[2] = 0;
    }
}



/* Mesh moves assuming pure gravitational free-fall, according to the 'GravAccel' quantity. This is useful for
    simulations of gravitational collapse in e.g. collapsing cores, at least up to a specific point (the
    routine below can trivially be modified to halt, at a certain point, if things become too dense) */
void MeshMotion_FreeFallExternalGravity(int i)
{
    double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); int k;
    for(k=0;k<3;k++) {SphP[i].ParticleVel[k]+=P[i].GravAccel[k]*dt;}
}


/* Mesh moves assuming a uniform expansion (e.g. a uniform free-expanding wind or Hubble-flow type
    expansion, homologously expanding away from a center). This is particularly useful for
    simpulations of winds, explosions, and rapidly-expanding media. Obviously not intended for
    cosmological expansion (just the mesh itself expanding, here) */
void MeshMotion_UniformExpansion(int i)
{
    double dvdr = 1; // velocity divergence (or effective "Hubble constant") of the flow, in code units //
    int k; double dp[3]; dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE)
    for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#elif defined(BOX_PERIODIC)
    dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y;
#endif
    for(k=0;k<3;k++) {SphP[i].ParticleVel[k] = dvdr * dp[k];}
}



/* Mesh moves assuming a uniform contraction (e.g. a uniform compression or gravitational free-fall).
 This is useful for some types of implosion simulations or self-gravitating collapse (although the free-fall
 example above may be more appropriate for gravitational collapse) */
void MeshMotion_UniformCollapse(int i)
{
    double dvdr = -1; // velocity divergence (or effective "Hubble constant") of the flow, in code units //
    int k; double dp[3]; dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE)
    for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#elif defined(BOX_PERIODIC)
    dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y;
#endif
    for(k=0;k<3;k++) {SphP[i].ParticleVel[k] = dvdr * dp[k];}
}



/* Mesh moves according to the equilibrium shear in the shearing-box coordinate system. Useful for
    shearing-sheet and shearing-box simulations, especially when highly sub-sonic local phenomena
    (e.g. sub-sonic turbulence or grain drift) need to be followed */
void MeshMotion_ShearingSheet(int i)
{
#ifdef BOX_SHEARING
    SphP[i].ParticleVel[0] = SphP[i].ParticleVel[1] = SphP[i].ParticleVel[2] = 0;
    SphP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] = -BOX_SHEARING_Q * (P[i].Pos[0]-boxHalf_X) * BOX_SHEARING_OMEGA_BOX_CENTER; // equilibrium motion is purely in phi
#ifdef GRAIN_RDI_TESTPROBLEM
     SphP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] -= All.Pressure_Gradient_Accel / (2. * BOX_SHEARING_OMEGA_BOX_CENTER); // equilibrium motion is purely in phi
#endif
#endif
}







#endif




