/* --------------------------------------------------------------------------------- */
/* ... non-ideal MHD term evaluation ...
 *
 * For SPH, this relys on the anisoptropic SPH second-derivative
 *  operator. So a large kernel is especially useful to minimize the systematic errors.
 *  For MFM/MFV methods, the consistent finite-volume formulation is used.
 *  In either case, since we solve the diffusion equations explicitly, a stronger timestep
 *  restriction is necessary (since the equations are not strictly hyperbolic); this is in timestep.c
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
double bflux_from_nonideal_effects[3]={0};
if((local.Mass > 0) && (P[j].Mass > 0))
{
    // set the effective scalar coefficients //
    double eta_i, eta_j, eta_ohmic, eta_hall, eta_ad;
    eta_i = local.Eta_MHD_OhmicResistivity_Coeff; eta_j = SphP[j].Eta_MHD_OhmicResistivity_Coeff;
    if(eta_i>0 && eta_j>0) {eta_ohmic = 2*eta_i*eta_j / (eta_i+eta_j);} else {eta_ohmic = 0;}
    eta_i = local.Eta_MHD_HallEffect_Coeff; eta_j = SphP[j].Eta_MHD_HallEffect_Coeff;
    if(eta_i>0 && eta_j>0) {eta_hall = 2*eta_i*eta_j / (eta_i+eta_j);} else {eta_hall = 0;}
    eta_i = local.Eta_MHD_AmbiPolarDiffusion_Coeff; eta_j = SphP[j].Eta_MHD_AmbiPolarDiffusion_Coeff;
    if(eta_i>0 && eta_j>0) {eta_ad = 2*eta_i*eta_j / (eta_i+eta_j);} else {eta_ad = 0;}
    
    // only go further if these are non-zero //
    double eta_max = DMAX(eta_ohmic , DMAX(eta_hall, eta_ad));
    int k_xyz_A=0, k_xyz_B=0;
    if(eta_max > 0)
    {
        // define the current J //
        double J_current[3], d_scalar[3], rinv2 = rinv*rinv;
        double J_direct[3]={0}, Jmag=0, grad_dot_x_ij[3]={0};
        for(k=0;k<3;k++) {d_scalar[k] = local.BPred[k] - BPred_j[k];}
        for(k=0;k<3;k++)
        {
            int k2;
            if(k==0) {k_xyz_A=2; k_xyz_B=1;}
            if(k==1) {k_xyz_A=0; k_xyz_B=2;}
            if(k==2) {k_xyz_A=1; k_xyz_B=0;}
            double tmp_grad_A = 0.5*(local.Gradients.B[k_xyz_A][k_xyz_B] + SphP[j].Gradients.B[k_xyz_A][k_xyz_B]); // construct averaged slopes //
            double tmp_grad_B = 0.5*(local.Gradients.B[k_xyz_B][k_xyz_A] + SphP[j].Gradients.B[k_xyz_B][k_xyz_A]);
            J_current[k] = tmp_grad_B - tmp_grad_A; // determine contribution to J //
            // calculate the 'direct' J needed for stabilizing numerical diffusion terms //
            J_direct[k] = rinv2*(kernel.dp[k_xyz_A]*d_scalar[k_xyz_B]-kernel.dp[k_xyz_B]*d_scalar[k_xyz_A]);
            if(J_current[k]*J_direct[k] < 0) {if(fabs(J_direct[k]) > 5.*fabs(J_current[k])) {J_current[k] = 0.0;}}
            Jmag += J_current[k]*J_current[k];
            for(k2=0;k2<3;k2++) {grad_dot_x_ij[k] += 0.5*(local.Gradients.B[k][k2]+SphP[j].Gradients.B[k][k2]) * kernel.dp[k2];}
        }
        
        // calculate the actual fluxes //
        double b_flux[3]={0}, JcrossB[3], JcrossBcrossB[3];
        JcrossB[0] = (bhat[2]*J_current[1]-bhat[1]*J_current[2]);
        JcrossB[1] = (bhat[0]*J_current[2]-bhat[2]*J_current[0]);
        JcrossB[2] = (bhat[1]*J_current[0]-bhat[0]*J_current[1]);
        JcrossBcrossB[0] = (bhat[2]*JcrossB[1]-bhat[1]*JcrossB[2]);
        JcrossBcrossB[1] = (bhat[0]*JcrossB[2]-bhat[2]*JcrossB[0]);
        JcrossBcrossB[2] = (bhat[1]*JcrossB[0]-bhat[0]*JcrossB[1]);
        if(eta_ohmic>0) {for(k=0;k<3;k++) {b_flux[k] += -eta_ohmic * J_current[k];}} // ohmic ~ J
        if(eta_ad>0) {for(k=0;k<3;k++) {b_flux[k] += eta_ad * JcrossBcrossB[k];}} // a.d. ~ (JxB)xB
        if(eta_hall>0) {for(k=0;k<3;k++) {b_flux[k] += -eta_hall * JcrossB[k];}} // hall ~ (JxB)
        
        // calculate dB/dt = Area.cross.Flux //
        bflux_from_nonideal_effects[0] = Face_Area_Vec[1]*b_flux[2] - Face_Area_Vec[2]*b_flux[1];
        bflux_from_nonideal_effects[1] = Face_Area_Vec[2]*b_flux[0] - Face_Area_Vec[0]*b_flux[2];
        bflux_from_nonideal_effects[2] = Face_Area_Vec[0]*b_flux[1] - Face_Area_Vec[1]*b_flux[0];
        
        // ok now construct the numerical fluxes based on the numerical diffusion coefficients and the direct-difference values between elements
        double eta_0 = v_hll * kernel.r * All.cf_atime; // standard numerical diffusivity needed for stabilizing fluxes
        double eta_ohmic_0=0,eta_ad_0=0,eta_hall_0=0,q=0; // now limit the numerical diffusivity to avoid unphysically fast diffusion
        if(eta_ohmic>0) {q=eta_0/eta_ohmic; eta_ohmic_0=eta_0*(0.2 + q)/(0.2 + q + q*q);}
        if(eta_ad>0) {q=eta_0/eta_ad; eta_ad_0=eta_0*(0.2 + q)/(0.2 + q + q*q);}
        if(eta_hall>0) {q=eta_0/eta_hall; eta_hall_0=eta_0*(0.2 + q)/(0.2 + q + q*q);}
        double b_flux_direct[3]={0}, JcrossB_direct[3], JcrossBcrossB_direct[3], db_direct[3];
        JcrossB_direct[0] = (bhat[2]*J_direct[1]-bhat[1]*J_direct[2]);
        JcrossB_direct[1] = (bhat[0]*J_direct[2]-bhat[2]*J_direct[0]);
        JcrossB_direct[2] = (bhat[1]*J_direct[0]-bhat[0]*J_direct[1]);
        JcrossBcrossB_direct[0] = (bhat[2]*JcrossB_direct[1]-bhat[1]*JcrossB_direct[2]);
        JcrossBcrossB_direct[1] = (bhat[0]*JcrossB_direct[2]-bhat[2]*JcrossB_direct[0]);
        JcrossBcrossB_direct[2] = (bhat[1]*JcrossB_direct[0]-bhat[0]*JcrossB_direct[1]);
        if(eta_ohmic_0>0) {for(k=0;k<3;k++) {b_flux_direct[k] += -eta_ohmic_0 * J_direct[k];}}
        if(eta_ad_0>0) {for(k=0;k<3;k++) {b_flux_direct[k] += eta_ad_0 * JcrossBcrossB_direct[k];}}
        if(eta_hall_0>0) {for(k=0;k<3;k++) {b_flux_direct[k] += -eta_hall_0 * JcrossB_direct[k];}}
        db_direct[0] = Face_Area_Vec[1]*b_flux_direct[2] - Face_Area_Vec[2]*b_flux_direct[1];
        db_direct[1] = Face_Area_Vec[2]*b_flux_direct[0] - Face_Area_Vec[0]*b_flux_direct[2];
        db_direct[2] = Face_Area_Vec[0]*b_flux_direct[1] - Face_Area_Vec[1]*b_flux_direct[0];
        
        double bfluxmag = 0;
        for(k=0;k<3;k++) {bfluxmag += bflux_from_nonideal_effects[k]*bflux_from_nonideal_effects[k];}
        bfluxmag /= (1.e-37 + Jmag * eta_max*eta_max * Face_Area_Norm*Face_Area_Norm);
        
        for(k=0;k<3;k++)
        {
            double d_scalar_tmp = d_scalar[k] - grad_dot_x_ij[k];
            double d_scalar_hll = MINMOD(d_scalar[k] , d_scalar_tmp);
            double hll_corr = bfluxmag * HLL_correction(d_scalar_hll,0.,1.,eta_max);
            double db_corr = bflux_from_nonideal_effects[k] + hll_corr;
            bflux_from_nonideal_effects[k] = MINMOD(1.1*bflux_from_nonideal_effects[k], db_corr);
            if(hll_corr!=0)
            {
                double db_direct_tmp = fabs(db_direct[k]) * hll_corr / fabs(hll_corr);
                if((bflux_from_nonideal_effects[k]*db_direct_tmp < 0) && (fabs(db_direct_tmp) > 1.0*fabs(bflux_from_nonideal_effects[k]))) {bflux_from_nonideal_effects[k]=0;}
            }
            if((bflux_from_nonideal_effects[k]*db_direct[k] < 0) && (fabs(db_direct[k]) > 10.0*fabs(bflux_from_nonideal_effects[k]))) {bflux_from_nonideal_effects[k]=0;}
        }
        
        // finally add one more flux-limiter to prevent the change in B from exceeding a large threshold in a single timestep //
        double cmag_B2 = bhat_mag * (bhat[0]*bflux_from_nonideal_effects[0]+bhat[1]*bflux_from_nonideal_effects[1]+bhat[2]*bflux_from_nonideal_effects[2]);
        if(dt_hydrostep > 0)
        {
            double cmag_lim = 0.01 * bhat_mag*bhat_mag / dt_hydrostep;
            if(fabs(cmag_B2) > cmag_lim)
            {
                double corr_visc = cmag_lim / fabs(cmag_B2);
                bflux_from_nonideal_effects[0]*=corr_visc; bflux_from_nonideal_effects[1]*=corr_visc; bflux_from_nonideal_effects[2]*=corr_visc;
            }
        }
        // -now- we can finally add this to the numerical fluxes //
        for(k=0;k<3;k++) {Fluxes.B[k] += bflux_from_nonideal_effects[k];}
    }
}
