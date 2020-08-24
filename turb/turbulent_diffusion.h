/* --------------------------------------------------------------------------------- */
/* ... turbulent diffusion (sub-grid) models ...
 *
 *  The basic equations here follow the Smagorinky eddy diffusion model. For SPH, the
 *    discretization comes from Wadsley 2008 & Shen 2010. However, some caution is needed,
 *    for SPH, this relys on the (noisy and zeroth-order inconsistent) SPH second-derivative
 *    operator. So a large kernel is especially useful to minimize the systematic errors.
 *  For MFM/MFV methods, the consistent finite-volume formulation is used, which
 *    greatly minimizes artificial (numerical) diffusion.
 *  In either case, since we solve the diffusion equations explicitly, a stronger timestep
 *    restriction is necessary (since the equations are parabolic); this is in timestep.c.
 *    This is very important (many implementations of these equations in the literature
 *    do not include the appropriate timestep and flux-limiters; that makes the equations
 *    numerically unstable (you can get an answer, but it might be wrong, independent of resolution)
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
#ifdef TURB_DIFF_METALS // turbulent diffusion of metals (passive scalar mixing) //
        
    if((local.Mass>0)&&(P[j].Mass>0)&&((local.TD_DiffCoeff>MIN_REAL_NUMBER)||(SphP[j].TD_DiffCoeff>MIN_REAL_NUMBER)))
    {
        double wt_i=0.5, wt_j=0.5, cmag, d_scalar;
        double diffusion_wt = wt_i*local.TD_DiffCoeff + wt_j*SphP[j].TD_DiffCoeff; // physical
#ifdef HYDRO_SPH
        diffusion_wt *= 0.5*(local.Density + SphP[j].Density)*All.cf_a3inv; // physical
#else
        diffusion_wt *= Riemann_out.Face_Density; // physical
#endif
        /* calculate implied mass flux 'across the boundary' to prevent excessively large coefficients */
        double massflux = fabs( Face_Area_Norm * diffusion_wt / (DMIN(kernel.h_i,kernel.h_j)*All.cf_atime) * dt_hydrostep / (DMIN(local.Mass,P[j].Mass)) );
        if(massflux > 0.25) {diffusion_wt *= 0.25/massflux;}
        
        int k_species;
        double rho_i = local.Density*All.cf_a3inv, rho_j = SphP[j].Density*All.cf_a3inv, rho_ij = 0.5*(rho_i+rho_j); // physical
        for(k_species=0;k_species<NUM_METAL_SPECIES;k_species++)
        {
            cmag = 0.0;
            double grad_dot_x_ij = 0.0;
            d_scalar = local.Metallicity[k_species]-P[j].Metallicity[k_species]; // physical
            for(k=0;k<3;k++)
            {
                double grad_direct = d_scalar * kernel.dp[k] * rinv*rinv; // 1/code length
#ifdef TURB_DIFF_METALS_LOWORDER
                double grad_ij = grad_direct;
#else
                double grad_ij = wt_i*local.Gradients.Metallicity[k_species][k] + wt_j*SphP[j].Gradients.Metallicity[k_species][k]; // 1/code length
#endif
                grad_dot_x_ij += grad_ij * kernel.dp[k]; // physical
                grad_ij = MINMOD(grad_ij , grad_direct);
                cmag += Face_Area_Vec[k] * grad_ij; // 1/code length
            }
            cmag /= All.cf_atime; // cmag has units of 1/r -- convert to physical

            double d_scalar_tmp = d_scalar - grad_dot_x_ij; // physical
            double d_scalar_hll = MINMOD(d_scalar , d_scalar_tmp);
            double hll_corr = rho_ij * HLL_correction(d_scalar_hll, 0, rho_ij, diffusion_wt) / (-diffusion_wt); // physical
            double cmag_corr = cmag + hll_corr;
            cmag = MINMOD(1.5*cmag, cmag_corr);
            double f_direct = Face_Area_Norm*d_scalar*rinv/All.cf_atime; // physical
            if((f_direct*cmag < 0) && (fabs(f_direct) > HLL_DIFFUSION_OVERSHOOT_FACTOR*fabs(cmag))) {cmag = 0;}
            
            cmag *= -diffusion_wt * dt_hydrostep; // physical
            if(fabs(cmag) > 0)
            {
                double zlim = 0.25 * DMIN( DMIN(local.Mass,P[j].Mass)*fabs(d_scalar) , DMAX(local.Mass*local.Metallicity[k_species] , P[j].Mass*P[j].Metallicity[k_species]) );
                if(fabs(cmag)>zlim) {cmag*=zlim/fabs(cmag);}
#ifndef HYDRO_SPH
                double dmet = (P[j].Metallicity[k_species]-local.Metallicity[k_species]) * fabs(mdot_estimated) * dt_hydrostep;
                cmag = MINMOD(dmet,cmag); // limiter based on mass exchange from MFV HLLC solver //
#endif
                out.Dyield[k_species] += cmag;
                if(j_is_active_for_fluxes) {SphP[j].Dyield[k_species] -= cmag;}
            }
        }
    }
    
#endif
}
