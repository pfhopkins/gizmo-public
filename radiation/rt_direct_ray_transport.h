/* --------------------------------------------------------------------------------- */
/* ... explicit radiation transport evaluation along rays (for Jiang et al. module) ...
 *
 * Solve the explicit transport problem for intensity along rays on the ray-grid
 *   carried by elements, following Jiang et al. 2014's method.
 *   Currently a piecewise-constant reconstruction (adding higher-order shortly,
 *   this is temporary but being used for testing before doing so)
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */

#ifdef RT_EVOLVE_INTENSITIES
if(local.Mass>0 && P[j].Mass>0 && dt_hydrostep>0 && Face_Area_Norm>0)
{
    double c_light = RT_SPEEDOFLIGHT_REDUCTION * (C/All.UnitVelocity_in_cm_per_s); // RSOL
    double c_hll = 0.5*fabs(face_vel_i-face_vel_j) + c_light, V_i_phys = V_i / All.cf_a3inv, V_j_phys = V_j / All.cf_a3inv, sthreeinv = 1./sqrt(3.); // physical units
    int k_freq, k_angle;
    for(k_freq=0; k_freq<N_RT_FREQ_BINS; k_freq++)
    {
        double kappa_ij = 0.5 * (local.RT_DiffusionCoeff[k_freq] + rt_diffusion_coefficient(j,k_freq)); // physical
        for(k_angle=0; k_angle<N_RT_INTENSITY_BINS; k_angle++)
        {
            double scalar_i=local.Intensity_Pred[k_freq][k_angle], scalar_j=SphP[j].Intensity_Pred[k_freq][k_angle], d_scalar=scalar_i-scalar_j, scalar_ij=0.5*(scalar_i+scalar_j);
            if(scalar_i+scalar_j > 0)
            {
                double Face_Area_Proj = All.RT_Intensity_Direction[k_angle][0]*Face_Area_Vec[0]+All.RT_Intensity_Direction[k_angle][1]*Face_Area_Vec[1]+All.RT_Intensity_Direction[k_angle][2]*Face_Area_Vec[2];
                cmag = c_light * Face_Area_Proj * scalar_ij;
                
#if 0
                double face_dot_flux=0., cmag=0., cmag_flux[3]={0}, grad_norm=0, flux_i[3]={0}, flux_j[3]={0}, thold_hll;
                /* calculate the eigenvalues for the HLLE flux-weighting */
                for(k=0;k<3;k++)
                {
                    flux_i[k] = local.Flux[k_freq][k]/V_i_phys; flux_j[k] = SphP[j].Flux_Pred[k_freq][k]/V_j_phys; // units (E_phys/[t_phys*L_phys^2])
                    double grad = 0.5*(flux_i[k] + flux_j[k]);
                    grad_norm += grad*grad;
                    face_dot_flux += Face_Area_Vec[k] * grad; /* remember, our 'flux' variable is a volume-integral */
                }
                grad_norm = sqrt(grad_norm) + MIN_REAL_NUMBER;
                double reduced_flux = grad_norm / ((C/All.UnitVelocity_in_cm_per_s) * 0.5*(scalar_i+scalar_j)); // |F|/(c*E): ratio of flux to optically thin limit
                if(reduced_flux > 1) {reduced_flux=1;} else {if(reduced_flux < 0) {reduced_flux=0;}}
                double cos_theta_face_flux = face_dot_flux / (Face_Area_Norm * grad_norm); // angle between flux and face vector normal
                if(cos_theta_face_flux < -1) {cos_theta_face_flux=-1;} else {if(cos_theta_face_flux > 1) {cos_theta_face_flux=1;}}
                double lam_m, lam_p, wt, y_f=1.-reduced_flux, y_f_h=sqrt(y_f), y_f_h2=sqrt(y_f_h), cth=cos_theta_face_flux/2.;
                wt = (1. + cos_theta_face_flux)*(1. + cos_theta_face_flux) / 4.;
                lam_m = sthreeinv - reduced_flux*(+cth + (1.-(y_f_h2+wt*y_f_h)/(1.+wt))*(+cth+sthreeinv));
                wt = (1. - cos_theta_face_flux)*(1. - cos_theta_face_flux) / 4.;
                lam_p = sthreeinv - reduced_flux*(-cth + (1.-(y_f_h2+wt*y_f_h)/(1.+wt))*(-cth+sthreeinv));
                
                if(lam_p < 0) {lam_p=0;} else {if(lam_p>1) {lam_p=1;}}
                if(lam_m < 0) {lam_m=0;} else {if(lam_m>1) {lam_m=1;}}
                double hlle_wtfac_f, hlle_wtfac_u, eps_wtfac_f = 1.0e-10; // minimum weight
                if((lam_m==0)&&(lam_p==0)) {hlle_wtfac_f=0.5;} else {hlle_wtfac_f=lam_p/(lam_p+lam_m);}
                if(hlle_wtfac_f < eps_wtfac_f) {hlle_wtfac_f=eps_wtfac_f;} else {if(hlle_wtfac_f > 1.-eps_wtfac_f) {hlle_wtfac_f=1.-eps_wtfac_f;}}
                hlle_wtfac_u = hlle_wtfac_f * (1.-hlle_wtfac_f) * (lam_p + lam_m); // weight for addition of diffusion term
                
                /* add asymptotic-preserving correction so that numerical flux doesn't unphysically dominate in optically thick limit */
                double v_eff_light = DMIN(c_light , kappa_ij / (Get_Particle_Size(j)*All.cf_atime)); // physical
                c_hll = 0.5*fabs(face_vel_i-face_vel_j) + DMAX(1.,hlle_wtfac_u) * v_eff_light; // physical
                double tau_c_j = Get_Particle_Size(j)*All.cf_atime * SphP[j].Kappa_RT[k_freq]*(SphP[j].Density*All.cf_a3inv); // = L_particle / (lambda_mean_free_path) = L*kappa*rho [physical units] //
                double hll_corr = 1./(1. + 1.5*DMAX(tau_c_i[k_freq],tau_c_j));
                /* q below is a limiter to try and make sure the diffusion speed given by the hll flux doesn't exceed the diffusion speed in the diffusion limit */
                double q = 0.5 * c_hll * (kernel.r*All.cf_atime) / fabs(MIN_REAL_NUMBER + kappa_ij); q = (0.2 + q) / (0.2 + q + q*q); // physical
                double renormerFAC = DMIN(1.,fabs(cos_theta_face_flux*cos_theta_face_flux * q * hll_corr));
#endif

                
                /* flux-limiter to ensure flow is always down the local gradient [no 'uphill' flow] */
                double f_direct = -fabs(Face_Area_Proj) * c_hll*renormerFAC * d_scalar, sign_c0=f_direct*cmag; // simple HLL term for frame moving at 1/2 inter-particle velocity: here not limited [physical units] //
                if((sign_c0 < 0) && (fabs(f_direct) > fabs(cmag))) {cmag = 0;}
                if(cmag != 0)
                {
                    if(f_direct != 0)
                    {
                        double thold_hll = (0.5*hlle_wtfac_u) * fabs(cmag); // add hll term but flux-limited //
                        if(fabs(f_direct) > thold_hll) {f_direct *= thold_hll/fabs(f_direct);}
                        cmag += f_direct;
                    }
                    // enforce a flux limiter for stability (to prevent overshoot) //
                    cmag *= dt_hydrostep; // all in physical units //
                    double sVi = scalar_i*V_i_phys, sVj = scalar_j*V_j_phys; // physical units //
                    thold_hll = 0.25 * DMIN(fabs(sVi-sVj), DMAX(fabs(sVi), fabs(sVj)));
                    if(sign_c0 < 0) {thold_hll *= 1.e-2;} // if opposing signs, restrict this term //
                    if(fabs(cmag)>thold_hll) {cmag *= thold_hll/fabs(cmag);}
                    cmag /= dt_hydrostep;
                    out.Dt_Intensity[k_freq][k_angle] += cmag;
                    if(j_is_active_for_fluxes) {SphP[j].Dt_Intensity[k_freq][k_angle] -= cmag;}
#ifdef RT_INFRARED
                    if(k_freq==RT_FREQ_BIN_INFRARED) // define advected radiation temperature based on direction of net radiation flow //
                    {
                        double Fluxes_E_gamma_T_weighted_IR=0;
                        if(cmag > 0) {Fluxes_E_gamma_T_weighted_IR = cmag/(MIN_REAL_NUMBER+SphP[j].Radiation_Temperature);} else {Fluxes_E_gamma_T_weighted_IR = cmag/(MIN_REAL_NUMBER+local.Radiation_Temperature);}
                        out.Dt_E_gamma_T_weighted_IR += Fluxes_E_gamma_T_weighted_IR;
                        if(j_is_active_for_fluxes) {SphP[j].Dt_E_gamma_T_weighted_IR -= Fluxes_E_gamma_T_weighted_IR;}
                    }
#endif
                } // cmag != 0
            } // scalar_i+scalar_j > 0
        } // k_angle loop
    } // k_freq loop
}
#endif
