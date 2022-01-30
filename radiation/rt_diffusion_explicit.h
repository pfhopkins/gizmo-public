/* --------------------------------------------------------------------------------- */
/* ... explicit radiation diffusion/flux transport evaluation ...
 *
 * For SPH, this relys on the (noisy and zeroth-order inconsistent) SPH second-derivative
 *  operator. So a large kernel is especially useful to minimize the systematic errors.
 *  For MFM/MFV methods, the consistent finite-volume formulation is used.
 *  In either case, since we solve the conduction equations explicitly, a stronger timestep
 *  restriction is necessary (since the equations are parabolic); this is in timestep.c
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
#define RT_ENHANCED_NUMERICAL_DIFFUSION /* option which increases numerical diffusion, to get smoother solutions, if desired; akin to slopelimiters~0 model */
{
    // first define some variables needed regardless //
    double c_light_eff = C_LIGHT_CODE_REDUCED, rsol_corr = RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    double v_frame[3]={0}; for(k=0;k<3;k++) {v_frame[k]=0.5*(ParticleVel_j[k] + local.ParticleVel[k])/All.cf_atime;} // frame velocity, not fluid velocity, is what appears here
#else
    double v_frame[3]={0}; for(k=0;k<3;k++) {v_frame[k]=0.5*(VelPred_j[k] + local.Vel[k])/All.cf_atime;} // variable to use below //
#endif
#if defined(RT_INFRARED)
    double Fluxes_Rad_E_gamma_T_weighted_IR = 0;
#endif
    double Fluxes_Rad_E_gamma[N_RT_FREQ_BINS];
    
    
#if !defined(RT_EVOLVE_FLUX) /* this means we just solve the diffusion equation for the eddington tensor, done in the loop below */
    int k_freq;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        Fluxes_Rad_E_gamma[k_freq] = 0;
        double kappa_ij = 0.5 * (local.RT_DiffusionCoeff[k_freq] + rt_diffusion_coefficient(j,k_freq)); // physical
        if((kappa_ij>0)&&(local.Mass>0)&&(P[j].Mass>0))
        {
            double scalar_i = local.Rad_E_gamma[k_freq] / V_i; // volumetric photon number density in this frequency bin (1/code volume) //
            double scalar_j = SphP[j].Rad_E_gamma_Pred[k_freq] / V_j;
            
            double d_scalar = (scalar_i - scalar_j); // units (1/code volume)
            double conduction_wt = kappa_ij * All.cf_a3inv/All.cf_atime;  // weight factor and conversion to physical units
            double cmag=0., grad_norm=0, grad_dot_x_ij=0.0;
            for(k=0;k<3;k++)
            {
                /* the flux is determined by the energy density gradient */
                double grad = 0.5 * (local.Gradients.Rad_E_gamma_ET[k_freq][k] + SphP[j].Gradients.Rad_E_gamma_ET[k_freq][k]); // (1/(code volume*code length))
                double grad_direct = d_scalar * kernel.dp[k] * rinv*rinv; // (1/(code volume*code length))
                grad_dot_x_ij += grad * kernel.dp[k]; // dp = local - j
                grad = MINMOD_G( grad , grad_direct );
#if defined(GALSF) || defined(COOLING) || defined(BLACK_HOLES)
                double grad_direct_vs_abs_fac = 2.0;
#else
                double grad_direct_vs_abs_fac = 5.0;
#endif
                if(grad*grad_direct < 0) {if(fabs(grad_direct) > grad_direct_vs_abs_fac*fabs(grad)) {grad = 0.0;}}
                cmag += Face_Area_Vec[k] * grad;
                grad_norm += grad*grad;
            }
            double A_dot_grad_alignment = cmag*cmag / (Face_Area_Norm*Face_Area_Norm * grad_norm); // dimensionless
            cmag *= -conduction_wt; // multiplies through the coefficient to get actual flux (physical) //

            /* here we add the HLL-like correction term. this greatly reduces noise and improves the stability of the diffusion.
            	however it comes at the cost of (significant) additional numerical diffusion */
            double v_eff_touse = DMIN(c_light_eff , kappa_ij / Particle_Size_j); // physical
            double c_hll = DMIN(0.5*fabs(face_vel_i-face_vel_j) + v_eff_touse , c_light_eff);
            double q = 0.5 * c_hll * kernel.r * All.cf_atime / fabs(1.e-37 + kappa_ij); q = (0.2 + q) / (0.2 + q + q*q); // physical
            double d_scalar_tmp = d_scalar - grad_dot_x_ij; // (1/code volume)
            double d_scalar_hll = MINMOD(d_scalar , d_scalar_tmp) * All.cf_a3inv; // physical
            double hll_tmp = -A_dot_grad_alignment * q * Face_Area_Norm * c_hll * d_scalar_hll; // physical
            
            /* add asymptotic-preserving correction so that numerical flux doesn't dominate in optically thick limit */
            double tau_c_j = Particle_Size_j * SphP[j].Rad_Kappa[k_freq]*SphP[j].Density*All.cf_a3inv; // = L_particle / (lambda_mean_free_path) = L*kappa*rho (physical) //
            double hll_corr = 1./(1. + 1.5*DMAX(tau_c_i[k_freq],tau_c_j));
            hll_tmp *= hll_corr;
            
            double thold_hll = 2.0*fabs(cmag);
            if(fabs(hll_tmp)>thold_hll) {hll_tmp*=thold_hll/fabs(hll_tmp);}
            double cmag_corr = cmag + hll_tmp;
            cmag = MINMOD(HLL_DIFFUSION_COMPROMISE_FACTOR*cmag, cmag_corr);
            /* flux-limiter to ensure flow is always down the local gradient */
            double f_direct = -conduction_wt * (1./9.) * Face_Area_Norm*d_scalar*rinv; // physical
            double check_for_stability_sign = f_direct*cmag;
            if((check_for_stability_sign < 0) && (fabs(f_direct) > HLL_DIFFUSION_OVERSHOOT_FACTOR*fabs(cmag))) {cmag = 0;}

            // prevent super-luminal local fluxes //
            double R_flux = fabs(cmag) / (3. * Face_Area_Norm * c_light_eff * (fabs(d_scalar)*All.cf_a3inv) + 1.e-37); // physical
            R_flux = (1. + 12.*R_flux) / (1. + 12.*R_flux*(1.+R_flux)); // 12 arbitrary but >>1 gives good behavior here //
#ifndef FREEZE_HYDRO
            cmag *= R_flux;
#endif
            
            // now we need to add the advective flux. note we do this after the limiters above, since those are designed for the diffusive terms, and this is simpler and more stable. we do this zeroth order (super-diffusive, but that's fine for our purposes)
            double cmag_adv=0, fluxlim_ij=1, v_Area_dot_rt=0, v_Fluid_dot_rt=0; for(k=0;k<3;k++) {v_Area_dot_rt += v_frame[k] * Face_Area_Vec[k]; v_Fluid_dot_rt += (0.5*(local.Vel[k]+VelPred_j[k])/All.cf_atime) * Face_Area_Vec[k];}
            double scalar_ij_phys = 2.*scalar_i*scalar_j/(scalar_i+scalar_j) * All.cf_a3inv; // use harmonic mean here, to weight lower value
#ifdef RT_FLUXLIMITER
            double fluxlim_j = return_flux_limiter(j,k_freq);
            double fluxlim_i = local.RT_DiffusionCoeff[k_freq] * local.Density * local.Rad_Kappa[k_freq] / c_light_eff; /* figure this out by what we've passed to save an extra variable being sent, here */
            fluxlim_ij = 0.5 * (fluxlim_i+fluxlim_j);
#endif
            double fac_fluxlim = rsol_corr*(4./3.)*fluxlim_ij*v_Fluid_dot_rt - v_Area_dot_rt; /* also need to account for where the RSOL terms appear in the asymptotic flux with an RSOL */
#ifdef RT_RADPRESSURE_IN_HYDRO
            fac_fluxlim = rsol_corr*fluxlim_ij*v_Fluid_dot_rt - v_Area_dot_rt + (2.*rt_absorb_frac_albedo(j,k_freq))*(fluxlim_ij/3.)*v_Fluid_dot_rt; // when P is included in hydro solver, some terms are automatically included, some not, so the pre-factor here needs to be revised
#endif
            cmag_adv += scalar_ij_phys * fac_fluxlim; // need to be careful with the sign here. since this is an oriented area and A points from j to i, need to watch the sign here (recently flipped in push - need to test). the 1/3 sometimes above owes to the fact that this is really the --pressure-- term for FLD-like methods, the energy term is implicitly part of the flux already if we're actually doing this correctly //
            cmag_adv += fabs(v_Area_dot_rt) * (scalar_j-scalar_i); // hll (rusanov) flux to stabilize and smooth flow
            if(fabs(cmag_adv)>0) // now limit the advective flux like we limit other fluxes
            {
                double cmag_max = 0.5 * DMAX(1,0.5*(tau_c_i[k_freq]+tau_c_j)) * (fabs(v_Area_dot_rt/Face_Area_Norm)/c_light_eff) * fabs(cmag); /* the 0.5 factor gives better result for the outflow test */
                if(fabs(cmag_adv) > cmag_max) {cmag_adv *= cmag_max/fabs(cmag_adv);}
                cmag += cmag_adv;
            }
            
            cmag *= dt_hydrostep; // all in physical units //
            if(fabs(cmag) > 0)
            {
                // enforce a flux limiter for stability (to prevent overshoot) //
                double thold_hll = 0.25 * DMIN(fabs(scalar_i*V_i-scalar_j*V_j),DMAX(fabs(scalar_i*V_i),fabs(scalar_j*V_j))); // physical
#ifdef RT_ENHANCED_NUMERICAL_DIFFUSION
                thold_hll *= 2.0;
#endif
                if(check_for_stability_sign<0) {thold_hll *= 1.e-2;}
                if(fabs(cmag)>thold_hll) {cmag *= thold_hll/fabs(cmag);}
                cmag /= dt_hydrostep;
                Fluxes_Rad_E_gamma[k_freq] += cmag;
#ifdef RT_INFRARED // define advected radiation temperature based on direction of net radiation flow //
                if(k_freq==RT_FREQ_BIN_INFRARED) {if(cmag > 0) {Fluxes_Rad_E_gamma_T_weighted_IR = cmag/(MIN_REAL_NUMBER+SphP[j].Radiation_Temperature);} else {Fluxes_Rad_E_gamma_T_weighted_IR = cmag/(MIN_REAL_NUMBER+local.Radiation_Temperature);}}
#endif
            } // if(conduction_wt > 0)
            
        } // close check that kappa and particle masses are positive
    }


#else /* RT_EVOLVE_FLUX is ON, so we don't solve a diffusion equation, but a system of two standard advection-like equations */


    int k_freq; double Fluxes_Rad_Flux[N_RT_FREQ_BINS][3], V_i_phys = V_i / All.cf_a3inv, V_j_phys = V_j / All.cf_a3inv;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        Fluxes_Rad_E_gamma[k_freq] = 0;
        Fluxes_Rad_Flux[k_freq][0]=Fluxes_Rad_Flux[k_freq][1]=Fluxes_Rad_Flux[k_freq][2]=0;
        double scalar_i = local.Rad_E_gamma[k_freq] / V_i_phys; // volumetric photon number density in this frequency bin (E_phys/L_phys^3)//
        double scalar_j = SphP[j].Rad_E_gamma_Pred[k_freq] / V_j_phys;
        if((scalar_i+scalar_j>0)&&(local.Mass>0)&&(P[j].Mass>0)&&(dt_hydrostep>0)&&(Face_Area_Norm>0))
        {
            double d_scalar = scalar_i - scalar_j;
            double face_dot_flux=0., cmag=0., cmag_flux[3]={0}, grad_norm=0, flux_i[3]={0}, flux_j[3]={0}, thold_hll;
            double kappa_i = local.RT_DiffusionCoeff[k_freq], kappa_j = rt_diffusion_coefficient(j,k_freq), kappa_ij = 0.5*(kappa_i+kappa_j); // physical units

            /* calculate the eigenvalues for the HLLE flux-weighting */
            for(k=0;k<3;k++)
            {
                flux_i[k] = local.Rad_Flux[k_freq][k]/V_i_phys - rsol_corr*v_frame[k]*scalar_i;; flux_j[k] = SphP[j].Rad_Flux_Pred[k_freq][k]/V_j_phys - rsol_corr*v_frame[k]*scalar_j; // units (E_phys/[t_phys*L_phys^2]) [physical]. include advective flux terms here
                double grad = 0.5*(flux_i[k] + flux_j[k]);
                grad_norm += grad*grad;
                face_dot_flux += Face_Area_Vec[k] * grad; /* remember, our 'flux' variable is a volume-integral */
            }
            grad_norm = sqrt(grad_norm) + MIN_REAL_NUMBER;
            double cos_theta_face_flux = face_dot_flux / (Face_Area_Norm * grad_norm); // angle between flux and face vector normal
            if(cos_theta_face_flux < -1) {cos_theta_face_flux=-1;} else {if(cos_theta_face_flux > 1) {cos_theta_face_flux=1;}}
#if 0
            double sthreeinv = 1./sqrt(3.);
            double reduced_flux = grad_norm / (c_light_eff * 0.5*(scalar_i+scalar_j)); // |f|/(c*E): ratio of flux to optically thin limit
            if(reduced_flux > 1) {reduced_flux=1;} else {if(reduced_flux < 0) {reduced_flux=0;}}
            double lam_m, lam_p, wt, y_f=1.-reduced_flux, y_f_h=sqrt(y_f), y_f_h2=sqrt(y_f_h), cth=cos_theta_face_flux/2.;
            wt = (1. + cos_theta_face_flux)*(1. + cos_theta_face_flux) / 4.;
            lam_m = sthreeinv - reduced_flux*(+cth + (1.-(y_f_h2+wt*y_f_h)/(1.+wt))*(+cth+sthreeinv));
            wt = (1. - cos_theta_face_flux)*(1. - cos_theta_face_flux) / 4.;
            lam_p = sthreeinv - reduced_flux*(-cth + (1.-(y_f_h2+wt*y_f_h)/(1.+wt))*(-cth+sthreeinv));
            /* alternative expression not accounting for rotation from Audit et al. astro-ph 0206281 */
            lam_p=lam_m=1; double f=DMIN(1.,DMAX(0.,reduced_flux)), eps_f=DMAX(f,DMIN(1.,1./(MIN_REAL_NUMBER+tau_c_i[k_freq]))), xd=sqrt(4.-3.*f*f), x=(3.+4.*f*f)/(5.+2.*xd), xp=2.*f/xd, xpe=1.*xp, xq=sqrt(DMAX(xpe*xpe+4.*(x-xp*f),0));
            lam_m=(xq-xpe)/2.; lam_p=(xq+xpe)/2.; if(lam_p>1 || lam_m <0) {lam_p=1.; lam_m=0.;}
            if(lam_p < 0) {lam_p=0;} else {if(lam_p>1) {lam_p=1;}}
            if(lam_m < 0) {lam_m=0;} else {if(lam_m>1) {lam_m=1;}}
            double hlle_wtfac_f, hlle_wtfac_u, eps_wtfac_f = 1.0e-10; // minimum weight
            if((lam_m==0)&&(lam_p==0)) {hlle_wtfac_f=0.5;} else {hlle_wtfac_f=lam_p/(lam_p+lam_m);}
            if(hlle_wtfac_f < eps_wtfac_f) {hlle_wtfac_f=eps_wtfac_f;} else {if(hlle_wtfac_f > 1.-eps_wtfac_f) {hlle_wtfac_f=1.-eps_wtfac_f;}}
            hlle_wtfac_u = hlle_wtfac_f * (1.-hlle_wtfac_f) * (lam_p + lam_m); // weight for addition of diffusion term
#else
            double hlle_wtfac_u=0.5, hlle_wtfac_f=0.5; // this corresponds to the Global-Lax-Friedrichs (GLF) flux function
#endif
            /* the flux is already known (its explicitly evolved, rather than determined by the gradient of the energy density */
            for(k=0;k<3;k++) {cmag += Face_Area_Vec[k] * (hlle_wtfac_f*flux_i[k] + (1.-hlle_wtfac_f)*flux_j[k]);} /* remember, our 'flux' variable is a volume-integral [all physical units here] */

            /* now compute the 'flux source term' - divergence of the radiation pressure tensor */
            double ET_dot_Face_i[3]={0}, ET_dot_Face_j[3]={0};
            eddington_tensor_dot_vector(local.ET[k_freq],Face_Area_Vec,ET_dot_Face_i); /* compute face dotted into eddington tensors for both sides */
            eddington_tensor_dot_vector(SphP[j].ET[k_freq],Face_Area_Vec,ET_dot_Face_j); /* compute face dotted into eddington tensors for both sides */
#ifdef RT_ENHANCED_NUMERICAL_DIFFUSION
            for(k=0;k<3;k++) {cmag_flux[k] += c_light_eff*c_light_eff * Face_Area_Vec[k] * 0.5*(scalar_i + scalar_j)/3.;}
#else
            for(k=0;k<3;k++) {cmag_flux[k] += c_light_eff*c_light_eff * (hlle_wtfac_f*scalar_i*ET_dot_Face_i[k] + (1.-hlle_wtfac_f)*scalar_j*ET_dot_Face_j[k]);}
#endif
            
            /* add asymptotic-preserving correction so that numerical flux doesn't unphysically dominate in optically thick limit */
            double v_eff_touse = DMIN(c_light_eff , kappa_ij / Particle_Size_j); // physical
            double c_hll = DMIN( 0.5*fabs(face_vel_i-face_vel_j) + DMAX(1.,hlle_wtfac_u) * v_eff_touse , c_light_eff ); // physical
            double tau_c_j = Particle_Size_j * SphP[j].Rad_Kappa[k_freq]*(SphP[j].Density*All.cf_a3inv); // = L_particle / (lambda_mean_free_path) = L*kappa*rho [physical units] //
            double hll_corr = 1./(1. + 1.5*DMAX(tau_c_i[k_freq],tau_c_j));
            /* q below is a limiter to try and make sure the diffusion speed given by the hll flux doesn't exceed the diffusion speed in the diffusion limit */
            double q = 0.5 * c_hll * (kernel.r*All.cf_atime) / fabs(MIN_REAL_NUMBER + kappa_ij); q = (0.2 + q) / (0.2 + q + q*q); // physical
            double renormerFAC = DMIN(1.,fabs(cos_theta_face_flux*cos_theta_face_flux * q * hll_corr));            
#ifdef RT_ENHANCED_NUMERICAL_DIFFUSION
            renormerFAC = 1;
#endif
#if !defined(RT_ENHANCED_NUMERICAL_DIFFUSION)
            double RTopticaldepth = DMIN(Particle_Size_i,Particle_Size_j) * c_light_eff / kappa_ij;
            double reductionfactor = sqrt((1.0-exp(-RTopticaldepth*RTopticaldepth))) / RTopticaldepth;
            double reducedcM1 = reductionfactor * c_light_eff / sqrt(3.0);
            v_eff_touse = DMAX(reducedcM1 , 2.*All.cf_afac3*kernel.vsig);
            c_hll = DMIN( 0.5*fabs(face_vel_i-face_vel_j) + v_eff_touse , c_light_eff );
            hll_corr = 1. / (1. + 1.5 * v_eff_touse * DMAX(Particle_Size_j/kappa_j , Particle_Size_i/kappa_i));
            q = 0.5 * c_hll * (kernel.r * All.cf_atime) / fabs(MIN_REAL_NUMBER + kappa_ij); q = (0.2 + q) / (0.2 + q + q*q);
            renormerFAC = DMIN(1.,fabs(cos_theta_face_flux*cos_theta_face_flux * q * hll_corr));

            double scalar_jr=scalar_j, scalar_ir=scalar_i, d_scalar_hll=d_scalar, d_scalar_ij=0;
            for(k=0;k<3;k++) {scalar_jr+=0.5*kernel.dp[k]*local.Gradients.Rad_E_gamma_ET[k_freq][k]; scalar_ir-=0.5*kernel.dp[k]*SphP[j].Gradients.Rad_E_gamma_ET[k_freq][k];}
            d_scalar_ij=scalar_ir-scalar_jr; if((d_scalar_ij*d_scalar>0)&&(fabs(d_scalar_ij)<fabs(d_scalar))) {d_scalar_hll=d_scalar_ij;}
            d_scalar = d_scalar_hll;
#endif

            /* flux-limiter to ensure flow is always down the local gradient [no 'uphill' flow] */
            double f_direct = -Face_Area_Norm * c_hll * d_scalar * renormerFAC; // simple HLL term for frame moving at 1/2 inter-particle velocity: here not limited [physical units] //
            if(dt_hydrostep>0) {
                double f_direct_max = 0.25 * fabs(d_scalar) * DMIN(V_i_phys , V_j_phys) / dt_hydrostep;
                if(fabs(f_direct) > f_direct_max) {f_direct *= f_direct_max / fabs(f_direct);}} else {f_direct=0;}
            double sign_c0 = f_direct*cmag;
            if((sign_c0 < 0) && (fabs(f_direct) > fabs(cmag))) {cmag = 0;}
#ifdef RT_ENHANCED_NUMERICAL_DIFFUSION
            if(cmag != 0)
#endif
            {
                if(f_direct != 0)
                {
                    thold_hll = (0.5*hlle_wtfac_u) * fabs(cmag); // add hll term but flux-limited //
//#ifndef RT_ENHANCED_NUMERICAL_DIFFUSION
                    if(fabs(f_direct) > thold_hll) {f_direct *= thold_hll/fabs(f_direct);}
//#endif
                    cmag += f_direct;
                }
                // enforce a flux limiter for stability (to prevent overshoot) //
                cmag *= dt_hydrostep; // all in physical units //
                double sVi = scalar_i*V_i_phys, sVj = scalar_j*V_j_phys; // physical units //
                thold_hll = 0.25 * DMAX(fabs(sVi-sVj), DMAX(fabs(sVi), fabs(sVj)));
#ifdef RT_ENHANCED_NUMERICAL_DIFFUSION
                thold_hll *= 2.0; // allow this term to be more generous //
#ifdef BH_WIND_SPAWN // 
		if(local.ConditionNumber < 0 || P[j].ID == All.AGNWindID) {thold_hll *= 0.25;}  // be extra conservative if dealing with fluxes involving jet cells - won't be particularly accurate anyway
#endif
#endif

//#ifndef RT_ENHANCED_NUMERICAL_DIFFUSION
                if(sign_c0 < 0) {thold_hll *= 1.e-2;} // if opposing signs, restrict this term //
//#endif
                if(fabs(cmag)>thold_hll) {cmag *= thold_hll/fabs(cmag);}
                cmag /= dt_hydrostep;
                Fluxes_Rad_E_gamma[k_freq] += cmag; // returned in physical units //
#ifdef RT_INFRARED // define advected radiation temperature based on direction of net radiation flow //
                if(k_freq==RT_FREQ_BIN_INFRARED) {if(cmag > 0) {Fluxes_Rad_E_gamma_T_weighted_IR = cmag/(MIN_REAL_NUMBER+SphP[j].Radiation_Temperature);} else {Fluxes_Rad_E_gamma_T_weighted_IR = cmag/(MIN_REAL_NUMBER+local.Radiation_Temperature);}}
#endif
            } // cmag != 0
            
            /* alright, now we calculate and limit the HLL diffusive fluxes for the radiative fluxes */
            double hll_mult_dmin = 1;
            for(k=0;k<3;k++)
            {
                // going to 0.5 instead of 1 here strongly increases diffusion: still get squared HII region, but very bad shadowing
                double f_direct = -Face_Area_Norm * c_hll * (flux_i[k]-flux_j[k]); // * renormerFAC; [all physical units]
                if(f_direct != 0)
                {
                    thold_hll = 1.0 * fabs(cmag_flux[k]) / fabs(f_direct); // coefficient of 0.5-1: 0.5=longer-lasting shadows, more 'memory' effect/shape distortion of HII regions; 1=fill in shadows faster, less HII distortion
                    if(thold_hll < hll_mult_dmin) {hll_mult_dmin = thold_hll;}
                }
            }
                        
            // now we need to add the advective flux. note we do this after the limiters above, since those are designed for the diffusive terms, and this is simpler and more stable. we do this zeroth order (super-diffusive, but that's fine for our purposes)
            double flux_dot_face=0,v_Area_dot_rt=0; for(k=0;k<3;k++) {flux_dot_face += 0.5*(flux_i[k]+flux_j[k])*Face_Area_Vec[k]; v_Area_dot_rt += rsol_corr*v_frame[k]*Face_Area_Vec[k];} // order of operations needs care here. following through on the divergence theorem carefully, the transport follows v, with A dotted into flux here, rather than A.v or some other law
            //for(k=0;k<3;k++) {cmag_flux[k] += -v_frame[k]*flux_dot_face + 0.5*fabs(v_Area_dot_rt)*(flux_j[k]-flux_i[k]);} // need to be careful with the sign here. since this is an oriented area and A points from j to i, need to flip the sign here //
            for(k=0;k<3;k++) {cmag_flux[k] += -0.5*(flux_i[k]+flux_j[k])*v_Area_dot_rt + 0.5*fabs(v_Area_dot_rt)*(flux_j[k]-flux_i[k]);} // need to be careful with the sign here. since this is an oriented area and A points from j to i, need to flip the sign here //

            for(k=0;k<3;k++)
            {
                double f_direct = -Face_Area_Norm * c_hll * (flux_i[k] - flux_j[k]); // [physical units]
                double sign_agreement = f_direct * cmag_flux[k];
                if((sign_agreement < 0) && (fabs(f_direct) > fabs(cmag_flux[k]))) {cmag_flux[k] = 0;}
                if(cmag_flux[k] != 0)
                {
                    cmag_flux[k] += hll_mult_dmin * f_direct; // add diffusive flux //
                    /* flux-limiter to prevent overshoot */
                    cmag_flux[k] *= dt_hydrostep;
                    thold_hll = DMIN( DMAX( DMIN(fabs(local.Rad_Flux[k_freq][k]), fabs(SphP[j].Rad_Flux_Pred[k_freq][k])) , fabs(local.Rad_Flux[k_freq][k]-SphP[j].Rad_Flux_Pred[k_freq][k]) ) , DMAX(fabs(local.Rad_Flux[k_freq][k]), fabs(SphP[j].Rad_Flux_Pred[k_freq][k])) );
                    double fii=V_i_phys*scalar_i*c_light_eff, fjj=V_j_phys*scalar_j*c_light_eff; // physical units //
                    double tij = DMIN( DMAX( DMIN(fabs(fii),fabs(fjj)) , fabs(fii-fjj) ) , DMAX(fabs(fii),fabs(fjj)) );
                    thold_hll = 0.25 * DMAX( thold_hll , 0.5*tij );
                    if(fabs(cmag_flux[k])>thold_hll) {cmag_flux[k] *= thold_hll/fabs(cmag_flux[k]);}
                    Fluxes_Rad_Flux[k_freq][k] += cmag_flux[k] / dt_hydrostep; // all returned in physical units
                }
            }
        } // close check that energy and masses are positive
    }

#endif
    
    // assign the actual fluxes //
    for(k=0;k<N_RT_FREQ_BINS;k++) {out.Dt_Rad_E_gamma[k] += Fluxes_Rad_E_gamma[k];}
    if(j_is_active_for_fluxes) {for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[j].Dt_Rad_E_gamma[k] -= Fluxes_Rad_E_gamma[k];}}
#if defined(RT_INFRARED)
    out.Dt_Rad_E_gamma_T_weighted_IR += Fluxes_Rad_E_gamma_T_weighted_IR;
    if(j_is_active_for_fluxes) {SphP[j].Dt_Rad_E_gamma_T_weighted_IR -= Fluxes_Rad_E_gamma_T_weighted_IR;}
#endif
#ifdef RT_EVOLVE_FLUX
    for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {out.Dt_Rad_Flux[k][k_dir] += Fluxes_Rad_Flux[k][k_dir];}}
    if(j_is_active_for_fluxes) {for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {SphP[j].Dt_Rad_Flux[k][k_dir] -= Fluxes_Rad_Flux[k][k_dir];}}}
#endif
    
}
