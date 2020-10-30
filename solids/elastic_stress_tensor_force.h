/* --------------------------------------------------------------------------------- */
/* ... elastic stress evaluation. currently a first-order solver. HLL artificial diffusivity
 *      terms not included here because with current formulation what is in the Reimann solver
 *      should already be sufficient, up to tensile check included. But this needs to be
 *      directly tested.
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    if( ((local.Mass>0)&&(P[j].Mass>0)) )
    {
        int k_v, j_v;
        double FNormT=local.Mass * P[j].Mass * fabs(kernel.dwk_i+kernel.dwk_j) / (local.Density * SphP[j].Density) * All.cf_atime*All.cf_atime;
        double FVec[3]={0}; for(k=0;k<3;k++) {FVec[k] = FNormT * kernel.dp[k]/kernel.r;} // use this particular face formulation for fluxes here to avoid tensile instability //

        double cmag[3]={0}, v_interface[3], wt_i=-0.5, wt_j=-0.5; // we need a minus sign at some point; its handy to just include it now in the weights //
        for(k=0;k<3;k++) {v_interface[k] = (wt_i*local.Vel[k] + wt_j*VelPred_j[k]) / All.cf_atime;} // physical units //
        double rho_i = local.Density * All.cf_a3inv, rho_j = SphP[j].Density * All.cf_a3inv;

        // terms for HLL fluxes
        double wt_r = rho_i*local.SoundSpeed * rho_j*SphP[j].SoundSpeed / (rho_i*local.SoundSpeed + rho_j*SphP[j].SoundSpeed) * FNormT; // physical
        double cT_i=sqrt(All.Tillotson_EOS_params[local.CompositionType][10]/rho_i), cT_j=sqrt(All.Tillotson_EOS_params[SphP[j].CompositionType][10]/rho_j); // physical
        double wt_t = rho_i*cT_i * rho_j*cT_j / (rho_i*cT_i + rho_j*cT_j) * FNormT; // physical
        double wt_rt = (wt_r - wt_t) * kernel.vdotr2 * rinv*rinv*All.cf_atime*All.cf_atime;
        // evaluate force from deviatoric stress tensor //
#if 0   // simple version for tensile correction based on which pressures are negative or not
        for(j_v=0;j_v<3;j_v++) {for(k_v=0;k_v<3;k_v++) {
            cmag[j_v] += (wt_i*local.Elastic_Stress_Tensor[j_v][k_v] + wt_j*SphP[j].Elastic_Stress_Tensor[j_v][k_v]) * FVec[k_v]*All.cf_a2inv; // direct flux
            if(local.Pressure < 0) {cmag[j_v] -= tensile_correction_factor *   local.Elastic_Stress_Tensor[j_v][k_v] * wt_i * FVec[k_v]*All.cf_a2inv;} // tensile correction
            if(SphP[j].Pressure < 0) {cmag[j_v] -= tensile_correction_factor * SphP[j].Elastic_Stress_Tensor[j_v][k_v] * wt_j * FVec[k_v]*All.cf_a2inv;} // tensile correction
        }}
#else   // fancier model where tensile correction is always applied to compressive principle component of the stress tensor [have to solve for eigenvalues of the stress tensor]
        int ij_switch; for(ij_switch=0;ij_switch<2;ij_switch++) {
            double nvt[NUMDIMS*NUMDIMS]={0}, norm_m=0, wtfac=0; if(ij_switch==0) {wtfac=wt_i;} else {wtfac=wt_j;}
            for(k_v=0;k_v<NUMDIMS;k_v++) {for(j_v=0;j_v<NUMDIMS;j_v++) {
            if(ij_switch==0) {nvt[NUMDIMS*k_v + j_v] = local.Elastic_Stress_Tensor[k_v][j_v];} else {nvt[NUMDIMS*k_v + j_v] = SphP[j].Elastic_Stress_Tensor_Pred[k_v][j_v];}
            norm_m += nvt[NUMDIMS*k_v + j_v]*nvt[NUMDIMS*k_v + j_v];}} // initialize auxiliary array to store for feeding to GSL eigen routine
            if(norm_m > MIN_REAL_NUMBER) {
                gsl_matrix_view M = gsl_matrix_view_array(nvt,NUMDIMS,NUMDIMS); gsl_vector *eigvals = gsl_vector_alloc(NUMDIMS); gsl_matrix *eigvecs = gsl_matrix_alloc(NUMDIMS,NUMDIMS);
                gsl_eigen_symmv_workspace *v = gsl_eigen_symmv_alloc(NUMDIMS); gsl_eigen_symmv(&M.matrix, eigvals, eigvecs, v);
                double eigenval_k=0, eigenvec_k[NUMDIMS]={0}, A_dot_v=0, prefac=0;
                for(k_v=0;k_v<NUMDIMS;k_v++) {eigenval_k = gsl_vector_get(eigvals, k_v); A_dot_v = 0;
                    for(j_v=0;j_v<NUMDIMS;j_v++) {eigenvec_k[j_v] = gsl_matrix_get(eigvecs, j_v, k_v); A_dot_v += FVec[j_v]*eigenvec_k[j_v];}
                    prefac = wtfac * eigenval_k * (A_dot_v*All.cf_a2inv); if(eigenval_k > 0) {prefac *= 1. - tensile_correction_factor;}
                    if(!isnan(eigenval_k)) {for(j_v=0;j_v<NUMDIMS;j_v++) {cmag[j_v] += prefac * eigenvec_k[j_v];}} // evaluate S.Face = S.[sum of eigenvalues times projection on Face onto each eigenvector]
                }
                gsl_eigen_symmv_free(v); gsl_vector_free(eigvals); gsl_matrix_free(eigvecs);} // free memory
        }
#endif
        for(j_v=0;j_v<3;j_v++) {cmag[j_v] -= wt_rt * kernel.dp[j_v]*All.cf_atime + wt_t * kernel.dv[j_v]/All.cf_atime;} // HLL-type fluxes
        /* // below limiter doesn't actually appear necessary in our tests, thus far; worth more detailed testing the future //
        if(kernel.vdotr2 < 0) // check if particles are approaching
        {
            double astress_dot_r=0, a_dot_r=0; for(k_v=0;k_v<3;k_v++) {astress_dot_r += cmag[k_v]*kernel.dp[k_v]; a_dot_r += Fluxes.v[k_v]*kernel.dp[k_v];}
            if((astress_dot_r < 0) && (a_dot_r+astress_dot_r < 0)) // check if this makes them approach faster (and if -net- term is also making them approach faster)
            {
                double cmag_dir[3]={0}, a_dot_r_alt=0; // check if the face area orientation is the difference
                for(j_v=0;j_v<3;j_v++)
                {
                    for(k_v=0;k_v<3;k_v++) {cmag_dir[j_v] += (wt_i*local.Elastic_Stress_Tensor[j_v][k_v] + wt_j*SphP[j].Elastic_Stress_Tensor[j_v][k_v]) * FNormT*kernel.dp[k_v]/kernel.r;}
                    a_dot_r_alt += kernel.dp[j_v]*cmag[j_v];
                }
                if((a_dot_r_alt >= 0)||(a_dot_r_alt+astress_dot_r>=0)) {for(k_v=0;k_v<3;k_v++) {cmag[k_v]=cmag_dir[k_v];}} else // if this solves it, use this, done!
                {
                    double h_eff = 0.5*(Particle_Size_i + Particle_Size_j); // effective inter-particle spacing around these elements
                    double corrfac = kernel.r/h_eff; corrfac=corrfac*corrfac; corrfac=corrfac/(1+corrfac); // factor which interpolates between 1 if r>>mean spacing, 0 if r<<mean spacing
                    for(k_v=0;k_v<3;k_v++) {cmag[k_v] *= corrfac;} // applies 'artificial stress' term which behaves similarly to Monaghan 2000
                    a_dot_r *= corrfac;
                }
            }
        }
        */
        
        double cmag_E = cmag[0]*v_interface[0] + cmag[1]*v_interface[1] + cmag[2]*v_interface[2];
        /* now add a flux-limiter to prevent overshoot (even when the directions are correct) */
        /*
        if(dt_hydrostep > 0)
        {
            double cmag_lim = 0.25 * (local.Mass + P[j].Mass) * (v_interface[0]*v_interface[0]+v_interface[1]*v_interface[1]+v_interface[2]*v_interface[2]) / dt_hydrostep;
            if(fabs(cmag_E) > cmag_lim)
            {
                double corr_visc = cmag_lim / fabs(cmag_E);
                cmag[0]*=corr_visc; cmag[1]*=corr_visc; cmag[2]*=corr_visc; cmag_E*=corr_visc;
            }
        }
        */
        /* ok now we can finally add this to the numerical fluxes */
        for(k=0;k<3;k++) {Fluxes.v[k] += cmag[k];}
        Fluxes.p += cmag_E;
    } // close check that kappa and particle masses are positive
}
