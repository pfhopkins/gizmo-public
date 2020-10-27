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
 * This file was taken almost entirely from the turbulent_diffusion.h file, which was written by 
 * Phil Hopkins (phopkins@caltech.edu) for GIZMO. It has been modified by Alex Richings to apply 
 * those same methods to the diffusion of ions and molecules in the CHIMES module. 
 */
/* --------------------------------------------------------------------------------- */
{
#ifdef CHIMES_TURB_DIFF_IONS  // turbulent diffusion of CHIMES ions and molecules (passive scalar mixing) //
        
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

        // Loop through all CHIMES species to compute the change in
        // the total number of ions/molecules of that species.
        for(k_species = 0; k_species < ChimesGlobalVars.totalNumberOfSpecies; k_species++)
        {
            cmag = 0.0;
            double grad_dot_x_ij = 0.0;

            // For d_scalar, we want the difference in the number of ions per unit mass.
            // This is analogous to turbulent_diffusion.h, which uses the difference in
            // Metallicity (i.e. metal mass over total mass) here.
            
            double local_abundance_times_mass = local.ChimesNIons[k_species], external_abundance_times_mass; // values we will need
            #pragma omp atomic read
            external_abundance_times_mass = SphP[j].ChimesNIons[k_species]; // this variable can be reset owing to how we code this, needs to be carefully read for thread safety here.
            double external_abundance_times_mass_0 = external_abundance_times_mass; // save initial value for below
            
            d_scalar = (local_abundance_times_mass / local.Mass) - (external_abundance_times_mass / P[j].Mass); // physical
            for(k = 0; k < 3; k++)
            {
                double grad_direct = d_scalar * kernel.dp[k] * rinv * rinv; // 1/code length

                // NOTE: For CHIMES, we ONLY do this with the LOWORDER method 
                double grad_ij = grad_direct;
                grad_dot_x_ij += grad_ij * kernel.dp[k]; // physical
                cmag += Face_Area_Vec[k] * grad_ij; // 1/code length
            }
            cmag /= All.cf_atime; // cmag has units of 1/r -- convert to physical

            double d_scalar_tmp = d_scalar - grad_dot_x_ij; // physical
            double d_scalar_hll = MINMOD(d_scalar , d_scalar_tmp);
            double hll_corr = rho_ij * HLL_correction(d_scalar_hll, 0, rho_ij, diffusion_wt) / (-diffusion_wt); // physical
            double cmag_corr = cmag + hll_corr;
            cmag = MINMOD(1.5 * cmag, cmag_corr);
            double f_direct = Face_Area_Norm * d_scalar * rinv / All.cf_atime; // physical
            if((f_direct * cmag < 0) && (fabs(f_direct) > HLL_DIFFUSION_OVERSHOOT_FACTOR * fabs(cmag))) {cmag = 0;}
            
            cmag *= -diffusion_wt * dt_hydrostep; // physical
            if(fabs(cmag) > 0)
            {
                double zlim = 0.25 * DMIN( DMIN(local.Mass, P[j].Mass) * fabs(d_scalar) , DMAX(local_abundance_times_mass , external_abundance_times_mass) );
                if(fabs(cmag) > zlim) {cmag *= zlim / fabs(cmag);}
#ifndef HYDRO_SPH
                double dmet = ((external_abundance_times_mass / P[j].Mass) - (local_abundance_times_mass / local.Mass)) * fabs(mdot_estimated) * dt_hydrostep;
                cmag = MINMOD(dmet,cmag); // limiter based on mass exchange from MFV HLLC solver //
#endif
                out.ChimesIonsYield[k_species] += cmag;
                external_abundance_times_mass -= cmag; // that metal mass must come out of the neighbor element
                #pragma omp atomic
                SphP[j].ChimesNIons[k_species] += external_abundance_times_mass - external_abundance_times_mass_0; // here we enforce machine-accurate conservation by swapping right here. that means we need to be very careful to do this in a thread-safe manner
            }
        }
    } 
#endif // CHIMES_TURB_DIFF_IONS 
}
