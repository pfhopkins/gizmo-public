/* --------------------------------------------------------------------------------- */
/* this is the sub-routine where we actually extrapolate quantities to the faces 
    and set up, then solve the pair-wise Riemann problem for the method */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
//#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && !(defined(GALSF) || defined(COOLING))
//#define DO_HALFSTEP_FOR_MESHLESS_METHODS 1
//#endif
#if (SLOPE_LIMITER_TOLERANCE==0)
#define HYDRO_FACE_AREA_LIMITER // use more restrictive face-area limiter in the simulations [some applications this is useful, but unclear if we can generally apply it] //
#endif
//#if !defined(PROTECT_FROZEN_FIRE) && !defined(HYDRO_FACE_AREA_LIMITER)
//#define HYDRO_FACE_VOLUME_RECONSTRUCTION_CORRECTION
//#endif
    
    double s_star_ij,s_i,s_j,v_frame[3],dummy_pressure,distance_from_i[3],distance_from_j[3],leak_vs_tol=0;
#if !(defined(HYDRO_KERNEL_SURFACE_VOLCORR) || defined(EOS_ELASTIC))
    leak_vs_tol = 0.5 * (local.FaceClosureError+SphP[j].FaceClosureError);
#endif
    dummy_pressure=face_area_dot_vel=face_vel_i=face_vel_j=Face_Area_Norm=0;
    double Pressure_i = local.Pressure, Pressure_j = SphP[j].Pressure;
#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)
    if((Pressure_i<0)||(Pressure_j<0)) /* negative pressures are allowed, but dealt with below by a constant shift and re-shift, which should be invariant for HLLC with the MFM method */
    {
        dummy_pressure = -DMIN(Pressure_i,Pressure_j);
        Pressure_i += dummy_pressure; Pressure_j += dummy_pressure;
        dummy_pressure *= 1. - tensile_correction_factor; /* we still need to include an effective stress for large negative pressures when elements are too close, to prevent tensile instability */
    }
#endif
#ifdef COSMIC_RAY_FLUID
    for(k=0;k<N_CR_PARTICLE_BINS;k++)
    {
        Fluxes.CosmicRayPressure[k] = 0;
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
        Fluxes.CosmicRayAlfvenEnergy[k][0] = Fluxes.CosmicRayAlfvenEnergy[k][1] = 0;
#endif
    }
#endif
    
    /* --------------------------------------------------------------------------------- */
    /* define volume elements and interface position */
    /* --------------------------------------------------------------------------------- */
    V_j = P[j].Mass / SphP[j].Density;
    s_star_ij = 0;
    //
#if !defined(MHD_CONSTRAINED_GRADIENT)
     //s_star_ij = 0.5 * kernel.r * (PPP[j].Hsml - local.Hsml) / (local.Hsml + PPP[j].Hsml); // old test, doesn't account for Hsml changing for condition number reasons
     //s_star_ij = 0.5 * kernel.r * (local.Density - SphP[j].Density) / (local.Density + SphP[j].Density); // frame with zero mass flux in a first-order reconstruction //
#endif
    //
    /* ------------------------------------------------------------------------------------------------------------------- */
    /* now we're ready to compute the volume integral of the fluxes (or equivalently an 'effective area'/face orientation) */
    /* ------------------------------------------------------------------------------------------------------------------- */
#include "compute_finitevol_faces.h" /* insert code block for computing Face_Area_Vec, Face_Area_Norm, etc. */

    if(Face_Area_Norm == 0)
    {
        memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifdef DIVBCLEANING_DEDNER
        Riemann_out.phi_normal_mean=Riemann_out.phi_normal_db=0;
#endif
    } else {
        if((Face_Area_Norm<=0)||(isnan(Face_Area_Norm))) {PRINT_WARNING("PANIC! Face_Area_Norm=%g Mij=%g/%g wk_ij=%g/%g Vij=%g/%g dx/dy/dz=%g/%g/%g NVT=%g/%g/%g NVT_j=%g/%g/%g \n",Face_Area_Norm,local.Mass,P[j].Mass,kernel.wk_i,kernel.wk_j,V_i,V_j,kernel.dp[0],kernel.dp[1],kernel.dp[2],local.NV_T[0][0],local.NV_T[0][1],local.NV_T[0][2],SphP[j].NV_T[0][0],SphP[j].NV_T[0][1],SphP[j].NV_T[0][2]); fflush(stdout);}
        double n_unit[3]; for(k=0;k<3;k++) {n_unit[k] = Face_Area_Vec[k] / Face_Area_Norm;} /* define useful unit vector for below */

        /* --------------------------------------------------------------------------------- */
        /* extrapolate the conserved quantities to the interaction face between the particles */
        /* first we define some useful variables for the extrapolation */
        /* --------------------------------------------------------------------------------- */
        s_i =  0.5 * kernel.r;
        s_j = -0.5 * kernel.r;
#ifdef DO_HALFSTEP_FOR_MESHLESS_METHODS
        /* advance the faces a half-step forward in time (given our leapfrog scheme, this actually has
            very, very weak effects on the errors. nonetheless it does help a small amount in reducing
            certain types of noise and oscillations (but not always!) */
        s_i += 0.5 * DMIN(dt_hydrostep_i,dt_hydrostep_j) * (local.Vel[0]*kernel.dp[0] + local.Vel[1]*kernel.dp[1] + local.Vel[2]*kernel.dp[2]) * rinv;
        s_j += 0.5 * DMIN(dt_hydrostep_i,dt_hydrostep_j) * (VelPred_j[0]*kernel.dp[0] + VelPred_j[1]*kernel.dp[1] + VelPred_j[2]*kernel.dp[2]) * rinv;
#endif
#ifdef DO_UPWIND_TIME_CENTERING
        //(simple up-winding formulation: use if desired instead of time-centered fluxes)//
        double delta_halfstep_i = kernel.sound_i*0.5*dt_hydrostep_i*(All.cf_afac3/All.cf_atime); if(delta_halfstep_i>s_i) {delta_halfstep_i=s_i;}
        double delta_halfstep_j = kernel.sound_j*0.5*dt_hydrostep_j*(All.cf_afac3/All.cf_atime); if(delta_halfstep_j>-s_j) {delta_halfstep_j=-s_j;}
        s_i = s_star_ij - s_i + delta_halfstep_i; /* projection element for gradients */
        s_j = s_star_ij - s_j - delta_halfstep_j;
#else
        s_i = s_star_ij - s_i; /* projection element for gradients */
        s_j = s_star_ij - s_j;
#endif
        distance_from_i[0]=kernel.dp[0]*rinv; distance_from_i[1]=kernel.dp[1]*rinv; distance_from_i[2]=kernel.dp[2]*rinv;
        for(k=0;k<3;k++) {distance_from_j[k] = distance_from_i[k] * s_j; distance_from_i[k] *= s_i;}
        //for(k=0;k<3;k++) {v_frame[k] = 0.5 * (VelPred_j[k] + local.Vel[k]);}
        for(k=0;k<3;k++) {v_frame[k] = rinv * (-s_i*VelPred_j[k] + s_j*local.Vel[k]);} // allows for face to be off-center (to second-order)
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
        for(k=0;k<3;k++) {v_frame[k] = rinv * (-s_i*ParticleVel_j[k] + s_j*local.ParticleVel[k]);}
#endif
        // (note that in the above, the s_i/s_j terms are crossed with the opposing velocity terms: this is because the face is closer to the
        //   particle with the smaller smoothing length; so it's values are slightly up-weighted //
        
        /* we need the face velocities, dotted into the face vector, for correction back to the lab frame */
        for(k=0;k<3;k++) {face_vel_i+=local.Vel[k]*n_unit[k]; face_vel_j+=VelPred_j[k]*n_unit[k];}
        face_vel_i /= All.cf_atime; face_vel_j /= All.cf_atime;
        face_area_dot_vel = rinv*(-s_i*face_vel_j + s_j*face_vel_i);
        
        /* also will need approach velocities to determine maximum upwind pressure */
        double v2_approach = 0;
        double vdotr2_phys = kernel.vdotr2;
        if(All.ComovingIntegrationOn) {vdotr2_phys -= All.cf_hubble_a2 * r2;}
        vdotr2_phys *= 1/(kernel.r * All.cf_atime);
        if(vdotr2_phys < 0) {v2_approach = vdotr2_phys*vdotr2_phys;}
        double vdotf2_phys = face_vel_i - face_vel_j; // need to be careful of sign here //
        if(vdotf2_phys < 0) {v2_approach = DMAX( v2_approach , vdotf2_phys*vdotf2_phys );}
        
        
        /* now we do the reconstruction (second-order reconstruction at the face) */
        int recon_mode = 1; // default to 'normal' reconstruction: some special physics will set this to zero for low-order reconstructions
#ifdef BH_WIND_SPAWN
        if((P[j].ID==All.AGNWindID)||(local.ConditionNumber<0)) {recon_mode = 0;} // one of the particles is a wind particle: use a low-order reconstruction for safety
#endif
#if defined(GALSF) || defined(COOLING)
        if(fabs(vdotr2_phys)*UNIT_VEL_IN_KMS > 1000.) {recon_mode = 0;} // particle approach/recession velocity > 1000 km/s: be extra careful here!
#endif
        //if(kernel.r > local.Hsml || kernel.r > PPP[j].Hsml) {recon_mode = 0;} // some extrapolation: this is more conservative but does help preserve contact discontinuities (perhaps too well?)
        if(leak_vs_tol > 1) {recon_mode = 0;}
        
        double rho_i=local.Density, rho_j=SphP[j].Density, P_i=Pressure_i, P_j=Pressure_j; // initialize for below
#if defined(HYDRO_FACE_VOLUME_RECONSTRUCTION_CORRECTION)
        if(Vi_inv_corr*Vj_inv_corr < 0.999) {recon_mode = 0;} // if this correction is needed, reconstruction is unsafe, must revert to lower-order
        rho_i*=Vi_inv_corr; P_i*=Vi_inv_corr; rho_j*=Vj_inv_corr; P_j*=Vj_inv_corr; // do scalar correction
#endif

        reconstruct_face_states(rho_i, local.Gradients.Density, rho_j, SphP[j].Gradients.Density,
                                distance_from_i, distance_from_j, &Riemann_vec.L.rho, &Riemann_vec.R.rho, recon_mode);
        reconstruct_face_states(P_i, local.Gradients.Pressure, P_j, SphP[j].Gradients.Pressure,
                                distance_from_i, distance_from_j, &Riemann_vec.L.p, &Riemann_vec.R.p, recon_mode);
#ifdef EOS_GENERAL
        reconstruct_face_states(local.InternalEnergyPred, local.Gradients.InternalEnergy, SphP[j].InternalEnergyPred, SphP[j].Gradients.InternalEnergy,
                                distance_from_i, distance_from_j, &Riemann_vec.L.u, &Riemann_vec.R.u, recon_mode);
        reconstruct_face_states(kernel.sound_i, local.Gradients.SoundSpeed, kernel.sound_j, SphP[j].Gradients.SoundSpeed,
                                distance_from_i, distance_from_j, &Riemann_vec.L.cs, &Riemann_vec.R.cs, recon_mode);
#endif
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.Vel[k], local.Gradients.Velocity[k], VelPred_j[k], SphP[j].Gradients.Velocity[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.v[k], &Riemann_vec.R.v[k], recon_mode);
            Riemann_vec.L.v[k] -= v_frame[k]; Riemann_vec.R.v[k] -= v_frame[k];
        }
#ifdef MAGNETIC
        int slim_mode = 1;
#ifdef MHD_CONSTRAINED_GRADIENT
        if((local.ConditionNumber < 0) || (SphP[j].FlagForConstrainedGradients == 0)) {slim_mode = 1;} else {slim_mode = -1;}
#endif
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.BPred[k], local.Gradients.B[k], BPred_j[k], SphP[j].Gradients.B[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.B[k], &Riemann_vec.R.B[k], slim_mode);
        }
#ifdef DIVBCLEANING_DEDNER
        reconstruct_face_states(local.PhiPred, local.Gradients.Phi, PhiPred_j, SphP[j].Gradients.Phi,
                                distance_from_i, distance_from_j, &Riemann_vec.L.phi, &Riemann_vec.R.phi, 2);
#endif
#endif
        

#ifdef DO_HALFSTEP_FOR_MESHLESS_METHODS
        /* advance the faces a half-step forward in time (given our leapfrog scheme, this actually has
            very, very weak effects on the errors. nonetheless it does help a small amount in reducing 
            certain types of noise and oscillations */
        double dt_half = 0.5*DMIN(dt_hydrostep_i,dt_hydrostep_j);
        for(k=0;k<3;k++)
        {
            Riemann_vec.R.rho -= dt_half * local.Density * local.Gradients.Velocity[k][k];
            Riemann_vec.L.rho -= dt_half * SphP[j].Density * SphP[j].Gradients.Velocity[k][k];
            Riemann_vec.R.p -= dt_half * GAMMA(j) * Pressure_i * local.Gradients.Velocity[k][k];
            Riemann_vec.L.p -= dt_half * GAMMA(j) * Pressure_j * SphP[j].Gradients.Velocity[k][k];
            double dv_l_half = -dt_half * local.Gradients.Pressure[k] / local.Density;
            double dv_r_half = -dt_half * SphP[j].Gradients.Pressure[k] / SphP[j].Density;
            Riemann_vec.R.v[k] += 0.5 * (dv_l_half - dv_r_half);
            Riemann_vec.L.v[k] += 0.5 * (dv_r_half - dv_l_half);
            v_frame[k] += 0.5*(dv_l_half + dv_r_half);
        }
#endif
        
       
        /* estimate maximum upwind pressure */
        double press_i_tot = Pressure_i + local.Density * v2_approach;
        double press_j_tot = Pressure_j + SphP[j].Density * v2_approach;
#ifdef MAGNETIC
        press_i_tot += 0.5 * kernel.b2_i * fac_magnetic_pressure;
        press_j_tot += 0.5 * kernel.b2_j * fac_magnetic_pressure;
#endif
        double press_tot_limiter;
#ifdef MAGNETIC
        press_tot_limiter = 2.0 * 1.1 * All.cf_a3inv * (press_i_tot + press_j_tot);
#else 
        press_tot_limiter = 1.1 * All.cf_a3inv * DMAX( press_i_tot , press_j_tot );
#endif
#if defined(EOS_GENERAL) || defined(HYDRO_MESHLESS_FINITE_VOLUME)
        press_tot_limiter *= 2.0;
#endif
#if (SLOPE_LIMITER_TOLERANCE==2)
        press_tot_limiter *= 100.0; // large number
#endif
        if(recon_mode==0) {press_tot_limiter = DMAX(press_tot_limiter , DMAX(DMAX(Pressure_i,Pressure_j),2.*DMAX(local.Density,SphP[j].Density)*v2_approach));}
#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)
        press_tot_limiter = 1.e10*(press_tot_limiter+1.); // it is unclear how this particular limiter behaves for solid-body EOS's, so for now, disable it in these cases
#endif

        
        /* --------------------------------------------------------------------------------- */
        /* Alright! Now we're actually ready to solve the Riemann problem at the particle interface */
        /* --------------------------------------------------------------------------------- */
        Riemann_solver(Riemann_vec, &Riemann_out, n_unit, press_tot_limiter);
        /* before going on, check to make sure we have a valid Riemann solution */
        if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M))||(Riemann_out.P_M>1.4*press_tot_limiter))
        {
            /* go to a linear reconstruction of P, rho, and v, and re-try */
            Riemann_vec.R.p = Pressure_i; Riemann_vec.L.p = Pressure_j;
            Riemann_vec.R.rho = local.Density; Riemann_vec.L.rho = SphP[j].Density;
            for(k=0;k<3;k++) {Riemann_vec.R.v[k]=local.Vel[k]-v_frame[k]; Riemann_vec.L.v[k]=VelPred_j[k]-v_frame[k];}
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.BPred[k]; Riemann_vec.L.B[k]=BPred_j[k];}
#ifdef DIVBCLEANING_DEDNER
            Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = PhiPred_j;
#endif
#endif
#ifdef EOS_GENERAL
            Riemann_vec.R.u = local.InternalEnergyPred; Riemann_vec.L.u = SphP[j].InternalEnergyPred;
            Riemann_vec.R.cs = kernel.sound_i; Riemann_vec.L.cs = kernel.sound_j;
#endif
            Riemann_solver(Riemann_vec, &Riemann_out, n_unit, 1.4*press_tot_limiter);
            if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
            {
                /* ignore any velocity difference between the particles: this should gaurantee we have a positive pressure! */
                Riemann_vec.R.p = Pressure_i; Riemann_vec.L.p = Pressure_j;
                Riemann_vec.R.rho = local.Density; Riemann_vec.L.rho = SphP[j].Density;
                for(k=0;k<3;k++) {Riemann_vec.R.v[k]=0; Riemann_vec.L.v[k]=0;}
#ifdef MAGNETIC
                for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.BPred[k]; Riemann_vec.L.B[k]=BPred_j[k];}
#ifdef DIVBCLEANING_DEDNER
                Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = PhiPred_j;
#endif
#endif
#ifdef EOS_GENERAL
                Riemann_vec.R.u = local.InternalEnergyPred; Riemann_vec.L.u = SphP[j].InternalEnergyPred;
                Riemann_vec.R.cs = kernel.sound_i; Riemann_vec.L.cs = kernel.sound_j;
#endif
                Riemann_solver(Riemann_vec, &Riemann_out, n_unit, 2.0*press_tot_limiter);
                if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
                {
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
                    printf("Riemann Solver Failed to Find Positive Pressure!: Pmax=%g PL/M/R=%g/%g/%g Mi/j=%g/%g rhoL/R=%g/%g H_ij=%g/%g vL=%g/%g/%g vR=%g/%g/%g n_unit=%g/%g/%g BL=%g/%g/%g BR=%g/%g/%g phiL/R=%g/%g \n",
                           press_tot_limiter,Riemann_vec.L.p,Riemann_out.P_M,Riemann_vec.R.p,local.Mass,P[j].Mass,Riemann_vec.L.rho,Riemann_vec.R.rho,local.Hsml,PPP[j].Hsml,
                           local.Vel[0]-v_frame[0],local.Vel[1]-v_frame[1],local.Vel[2]-v_frame[2],
                           VelPred_j[0]-v_frame[0],VelPred_j[1]-v_frame[1],VelPred_j[2]-v_frame[2],
                           n_unit[0],n_unit[1],n_unit[2],
                           Riemann_vec.L.B[0],Riemann_vec.L.B[1],Riemann_vec.L.B[2],
                           Riemann_vec.R.B[0],Riemann_vec.R.B[1],Riemann_vec.R.B[2],
                           Riemann_vec.L.phi,Riemann_vec.R.phi);
#else
                    printf("Riemann Solver Failed to Find Positive Pressure!: Pmax=%g PL/M/R=%g/%g/%g Mi/j=%g/%g rhoL/R=%g/%g vL=%g/%g/%g vR=%g/%g/%g n_unit=%g/%g/%g \n",
                           press_tot_limiter,Riemann_vec.L.p,Riemann_out.P_M,Riemann_vec.R.p,local.Mass,P[j].Mass,Riemann_vec.L.rho,Riemann_vec.R.rho,
                           Riemann_vec.L.v[0],Riemann_vec.L.v[1],Riemann_vec.L.v[2],
                           Riemann_vec.R.v[0],Riemann_vec.R.v[1],Riemann_vec.R.v[2],n_unit[0],n_unit[1],n_unit[2]);
#endif
                    endrun(1234);
                }
            }
        } // closes loop of alternative reconstructions if invalid pressures are found //
        
        /* --------------------------------------------------------------------------------- */
        /* Calculate the fluxes (EQUATION OF MOTION) -- all in physical units -- */
        /* --------------------------------------------------------------------------------- */
        if((Riemann_out.P_M>0)&&(!isnan(Riemann_out.P_M)))
        {
            if(All.ComovingIntegrationOn) {for(k=0;k<3;k++) v_frame[k] /= All.cf_atime;} // correct units
#ifdef TURB_DIFF_METALS
            mdot_estimated = Riemann_out.Mdot_estimated * Face_Area_Norm; // needed below
#endif
            if(dummy_pressure != 0) // if we had to shift the pressure zero-point, use this point to correct back to the (allowed) negative pressures //
            {
                for(k=0;k<3;k++) {Riemann_out.Fluxes.v[k] -= dummy_pressure * n_unit[k];} /* total momentum flux */
                Riemann_out.Fluxes.p -= dummy_pressure * Riemann_out.S_M; // default: total energy flux = v_frame.dot.mom_flux. note this is in the frame here, will correct for frame motion below. //
            }
            
            /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame' which we do following Pakmor et al. 2011 */
            for(k=0;k<3;k++)
            {
                Riemann_out.Fluxes.p += v_frame[k] * Riemann_out.Fluxes.v[k];
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
                /* Riemann_out->Fluxes.rho is un-modified */
                Riemann_out.Fluxes.p += (0.5*v_frame[k]*v_frame[k])*Riemann_out.Fluxes.rho;
                Riemann_out.Fluxes.v[k] += v_frame[k] * Riemann_out.Fluxes.rho; /* just boost by frame vel (as we would in non-moving frame) */
#endif
            }
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Riemann_out.Fluxes.B[k] += -v_frame[k] * Riemann_out.B_normal_corrected;} /* v dotted into B along the normal to the face (careful of sign here) */
#endif
            
            /* ok now we can actually apply this to the EOM */
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
            Fluxes.rho = Face_Area_Norm * Riemann_out.Fluxes.rho;
#endif
            Fluxes.p = Face_Area_Norm * Riemann_out.Fluxes.p; // this is really Dt of --total-- energy, need to subtract KE component for e */
            for(k=0;k<3;k++) {Fluxes.v[k] = Face_Area_Norm * Riemann_out.Fluxes.v[k];} // momentum flux (need to divide by mass) //
#if defined(COSMIC_RAY_FLUID) && defined(HYDRO_MESHLESS_FINITE_VOLUME)
            /* here we simply assume that if there is mass flux, the cosmic ray fluid is advected -with the mass flux-, taking an
             implicit constant (zeroth-order) reconstruction of the CR energy density at the face (we could reconstruct the CR
             properties at the face, and calculate a more accurate advection term; however at that stage we should actually be
             including them self-consistently in the Riemann problem */
            for(k=0;k<N_CR_PARTICLE_BINS;k++)
            {
                if(Fluxes.rho < 0)
                {
                    Fluxes.CosmicRayPressure[k] = Fluxes.rho * (local.CosmicRayPressure[k]*V_i/((GAMMA_COSMICRAY(k)-1.)*local.Mass)); /* note: CosmicRayPressure and V_i have comoving units, their product has physical units */
                } else {
                    Fluxes.CosmicRayPressure[k] = Fluxes.rho * (CosmicRayPressure_j[k]*V_j/((GAMMA_COSMICRAY(k)-1.)*P[j].Mass));
                }
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
                int kAlf=0; for(kAlf=0;kAlf<2;kAlf++) {if(Fluxes.rho<0) {Fluxes.CosmicRayAlfvenEnergy[k][kAlf]+=local.CosmicRayAlfvenEnergy[k][kAlf]*Fluxes.rho/local.Mass;} else {Fluxes.CosmicRayAlfvenEnergy[k][kAlf]+=SphP[j].CosmicRayAlfvenEnergy[k][kAlf]*Fluxes.rho/local.Mass;}}
#endif
            }
#endif
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Fluxes.B[k] = Face_Area_Norm * Riemann_out.Fluxes.B[k];} // magnetic flux (B*V) //
            Fluxes.B_normal_corrected = -Riemann_out.B_normal_corrected * Face_Area_Norm;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME)
            //Fluxes.phi = Riemann_out.Fluxes.phi * Face_Area_Norm; // after testing, we now prefer the mass-based fluxes, now //
            if(Fluxes.rho < 0)
            {
                Fluxes.phi = Fluxes.rho * Riemann_vec.R.phi;
            } else {
                Fluxes.phi = Fluxes.rho * Riemann_vec.L.phi; // phi_ij, phi_L/R, or midpoint phi?
            }
#endif
#endif // MAGNETIC

#if defined(HYDRO_MESHLESS_FINITE_MASS) && (SLOPE_LIMITER_TOLERANCE < 2) && !(defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)) // below is defined for adiabatic ideal fluids, don't use for materials
            /* for MFM, do the face correction for adiabatic flows here */
            double SM_over_ceff = fabs(Riemann_out.S_M) / DMIN(kernel.sound_i,kernel.sound_j); // for now use sound speed here (more conservative) vs magnetosonic speed //
            /* if SM is sufficiently large, we do nothing to the equations */
            if((SM_over_ceff < epsilon_entropic_eos_big && All.ComovingIntegrationOn >= 0) || (leak_vs_tol>1))
            {
                /* ok SM is small, we should use adiabatic equations instead */
#ifdef MAGNETIC
                // convert the face pressure P_M to a gas pressure alone (subtract the magnetic term) //
                for(k=0;k<3;k++) {Riemann_out.P_M -= 0.5*Riemann_out.Face_B[k]*Riemann_out.Face_B[k];}
#endif
                int use_entropic_energy_equation = 1;
                double facenorm_pm = Riemann_out.P_M * Face_Area_Norm;
                double PdV_fac = Riemann_out.P_M * vdotr2_phys / All.cf_a2inv;
                double PdV_i = kernel.dwk_i * V_i*V_i * local.DhsmlNgbFactor * PdV_fac;
                double PdV_j = kernel.dwk_j * V_j*V_j * PPP[j].DhsmlNgbFactor * PdV_fac;
                double du_old = facenorm_pm * (Riemann_out.S_M + face_area_dot_vel);
                double du_new = 0.5 * (PdV_i - PdV_j + facenorm_pm * (face_vel_i+face_vel_j));
                // more detailed check for intermediate cases //
                double cnum2 = SphP[j].ConditionNumber*SphP[j].ConditionNumber;
                if(SM_over_ceff > epsilon_entropic_eos_small && cnum2 < cnumcrit2)
                {
                    if(Pressure_i/local.Density != Pressure_j/SphP[j].Density)
                    {
                        if(Pressure_i/local.Density > Pressure_j/SphP[j].Density)
                        {
                            double dtoj = -du_old + facenorm_pm * face_vel_j;
                            if(dtoj > 0) {use_entropic_energy_equation=0;} else
                                {if(dtoj < 0) {if(dtoj > -du_new+facenorm_pm*face_vel_j) {use_entropic_energy_equation=0;}}}
                        } else {
                            double dtoi = +du_old - facenorm_pm * face_vel_i;
                            if(dtoi > 0) {use_entropic_energy_equation=0;} else
                                {if(dtoi < 0) {if(dtoi > +du_new-facenorm_pm*face_vel_i) {use_entropic_energy_equation=0;}}}
                        }
                    }
                }
                if(cnum2 >= cnumcrit2) {use_entropic_energy_equation=1;}
                // alright, if we've come this far, we need to subtract -off- the thermal energy part of the flux, and replace it //
                if(use_entropic_energy_equation) {Fluxes.p += du_new - du_old;}
            }
#endif // closes adiabatic flow face correction check //

        } else {
            /* nothing but bad riemann solutions found! */
            memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifdef DIVBCLEANING_DEDNER
            Riemann_out.phi_normal_mean=Riemann_out.phi_normal_db=0;
#endif
        }
    } // Face_Area_Norm != 0
}
