/* --------------------------------------------------------------------------------- */
/* ... real conduction evaluation ...
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
{
    double scalar_i = local.InternalEnergyPred; // physical units
    double scalar_j = SphP[j].InternalEnergyPred; // physical units
    double kappa_i = local.Kappa_Conduction; // physical units
    double kappa_j = SphP[j].Kappa_Conduction; // physical units
    
    if((kappa_i>MIN_REAL_NUMBER)&&(kappa_j>MIN_REAL_NUMBER)&&(local.Mass>0)&&(P[j].Mass>0))
    {
        double d_scalar = scalar_i - scalar_j;
        double rho_i, rho_j, rho_ij;
        rho_i = local.Density*All.cf_a3inv; rho_j = SphP[j].Density*All.cf_a3inv; rho_ij = 0.5*(rho_i+rho_j); // physical units
        
        // NOT SPH: Now we use the more accurate finite-volume formulation, with the effective faces we have already calculated //
        double *grad_i = local.Gradients.InternalEnergy;  // physical u / code length
        double *grad_j = SphP[j].Gradients.InternalEnergy;
        
        double flux_wt = rho_ij;
        double diffusion_wt = 0.5*(kappa_i+kappa_j);
#ifdef COOLING
        diffusion_wt = kappa_i*kappa_j/diffusion_wt;
#endif
        int do_isotropic = 1;
        double b_hll=1, cmag=0, wt_i=0.5, wt_j=0.5, grad_dot_x_ij = 0;
        double grad_ij[3];
        for(k=0;k<3;k++)
        {
            double q_grad = wt_i*grad_i[k] + wt_j*grad_j[k];
            double q_direct = d_scalar * kernel.dp[k] * rinv*rinv; // physical u / code length
            grad_dot_x_ij += q_grad * kernel.dp[k]; // physical u
#ifdef MAGNETIC
            grad_ij[k] = MINMOD_G(q_grad , q_direct);
            if(q_grad*q_direct < 0) {if(fabs(q_direct) > 5.*fabs(q_grad)) {grad_ij[k] = 0.0;}}
#else
            grad_ij[k] = MINMOD(q_grad , q_direct);
#endif
        }
#ifdef MAGNETIC
        if(bhat_mag > 0)
        {
            do_isotropic = 0;
            double B_interface_dot_grad_T = 0.0, grad_mag = 0.0;
            for(k=0;k<3;k++)
            {
                B_interface_dot_grad_T += bhat[k] * grad_ij[k];
                grad_mag += grad_ij[k]*grad_ij[k];
            }
            for(k=0;k<3;k++) {cmag += bhat[k] * Face_Area_Vec[k];} // physical
            cmag *= B_interface_dot_grad_T; // physical / code length
            if(grad_mag > 0) {grad_mag = sqrt(grad_mag);} else {grad_mag=1;}
            b_hll = B_interface_dot_grad_T / grad_mag; // physical
            b_hll *= b_hll; // physical
        }
#endif
        if(do_isotropic) {for(k=0;k<3;k++) {cmag += Face_Area_Vec[k] * grad_ij[k];}}
        cmag /= All.cf_atime; // cmag has units of u/r -- convert to physical
        
        /* obtain HLL correction terms for Reimann problem solution */
        double d_scalar_tmp = d_scalar - grad_dot_x_ij; // physical
        double d_scalar_hll = MINMOD(d_scalar , d_scalar_tmp);
        double hll_corr = b_hll * flux_wt * HLL_correction(d_scalar_hll, 0, flux_wt, diffusion_wt) / (-diffusion_wt); // physical
        double cmag_corr = cmag + hll_corr;
        cmag = MINMOD(HLL_DIFFUSION_COMPROMISE_FACTOR*cmag, cmag_corr);
        /* slope-limiter to ensure heat always flows from hot to cold */
        double d_scalar_b = b_hll * d_scalar; // physical
        double f_direct = Face_Area_Norm*d_scalar_b*rinv/All.cf_atime; // physical
        double check_for_stability_sign = f_direct*cmag; // physical
        if((check_for_stability_sign < 0) && (fabs(f_direct) > HLL_DIFFUSION_OVERSHOOT_FACTOR*fabs(cmag))) {cmag = 0;}
        cmag *= -diffusion_wt; /* multiply through coefficient to get flux (physical units) */
        
        
        /* follow that with a flux limiter as well */
        diffusion_wt = dt_hydrostep * cmag; // all in physical units //
        if(fabs(diffusion_wt) > 0)
        {
            // enforce a flux limiter for stability (to prevent overshoot) //
            //double du_ij_cond = 1.0*DMIN(local.Mass*scalar_i, P[j].Mass*scalar_j);
            double du_ij_cond = DMIN( 0.25*fabs(local.Mass*scalar_i-P[j].Mass*scalar_j) , DMAX(local.Mass*scalar_i , P[j].Mass*scalar_j));
            if(check_for_stability_sign<0) {du_ij_cond *= 1.e-2;}
            if(fabs(diffusion_wt)>du_ij_cond) {diffusion_wt *= du_ij_cond/fabs(diffusion_wt);}
            Fluxes.p += diffusion_wt / dt_hydrostep;
        } // if(diffusion_wt > 0)
        
    } // close check that kappa and particle masses are positive
}
