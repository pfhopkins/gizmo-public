/* --------------------------------------------------------------------------------- */
/* ... real anisotropic viscosity evaluation ...
 *
 * For SPH, this relys on the (noisy and zeroth-order inconsistent) SPH second-derivative
 *  operator. So a large kernel is especially useful to minimize the systematic errors.
 *  For MFM/MFV methods, the consistent finite-volume formulation is used.
 *  In either case, since we solve the viscous equations explicitly, a stronger timestep
 *  restriction is necessary (since the equations are not strictly hyperbolic); this is in timestep.c
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    if( (((local.Eta_ShearViscosity>MIN_REAL_NUMBER)&&(SphP[j].Eta_ShearViscosity>MIN_REAL_NUMBER)) ||
         ((local.Zeta_BulkViscosity>MIN_REAL_NUMBER)&&(SphP[j].Zeta_BulkViscosity>MIN_REAL_NUMBER))) &&
       ((local.Mass>0)&&(P[j].Mass>0)) )
    {
        int k_v;
        double a_inv = 1 / All.cf_atime;
        double v_interface[3];
        double cmag[3];
        double wt_i,wt_j;
        wt_i = wt_j = 0.5;
        for(k=0;k<3;k++) {v_interface[k] = (wt_i*local.Vel[k] + wt_j*VelPred_j[k]) * a_inv;} // physical units //
        
        // use a geometric average, since we want to weight the smaller of the two coefficients //
        double eta = 0.5 * (local.Eta_ShearViscosity + SphP[j].Eta_ShearViscosity);
        if(eta > 0) {eta = local.Eta_ShearViscosity * SphP[j].Eta_ShearViscosity / eta;} else {eta = 0;} // also converts to physical units
        double zeta = 0.5 * (local.Zeta_BulkViscosity + SphP[j].Zeta_BulkViscosity);
        if(zeta > 0) {zeta = local.Zeta_BulkViscosity * SphP[j].Zeta_BulkViscosity / zeta;} else {zeta = 0;} // also converts to physical units
        double viscous_wt_physical = DMAX(eta,zeta);
        // we need a minus sign at some point; its handy to just include it now in the weights //
        wt_i *= -1.; wt_j *= -1.;
        int do_non_mhd = 1;
        
        
#ifdef MAGNETIC
        /* should use the solution in the appropriate face of the Riemann problem for interface values */
        double B_interface[3],B_interface_mag2=0;
        for(k=0;k<3;k++)
        {
            B_interface[k] = bhat[k] * bhat_mag;
            B_interface_mag2 += B_interface[k]*B_interface[k];
        }
        double bhat_dot_gradvhat = 1;
        double bhat_dot_gradvhat_direct = 1;
        double grad_v_mag = 0;
        double Bi_proj = 1;
        if(B_interface_mag2 >= 0)
        {
            /* Braginskii formulation of leading-order anisotropic viscosity */
            do_non_mhd = 0;
            double rhs = 0, rhs_direct = 0;
            double one_third = 1./3. * B_interface_mag2;
            Bi_proj = (kernel.dp[0]*B_interface[0]+kernel.dp[1]*B_interface[1]+kernel.dp[2]*B_interface[2]) / B_interface_mag2;
            for(k_v=0;k_v<3;k_v++)
            {
                double tmp;
                double tmp_kv = B_interface[k_v]*B_interface[k_v] - one_third;
                double b_tensor_dot_face = 0;
                // note that because of the symmetry of the tensors here, the order of k,k_v doesn't matter //
                for(k=0;k<3;k++)
                {
                    double grad_v = wt_i*local.Gradients.Velocity[k_v][k]+wt_j*SphP[j].Gradients.Velocity[k_v][k];
                    double grad_direct = -kernel.dv[k_v] * kernel.dp[k] * rinv*rinv;
                    grad_v = MINMOD_G( grad_v, grad_direct);
                    if(grad_v*grad_direct < 0) {if(fabs(grad_direct) > 2.*fabs(grad_v)) {grad_v = 0.0;}}
                    
                    if(k==k_v) {tmp=tmp_kv;} else {tmp=B_interface[k_v]*B_interface[k];}
                    rhs += tmp * grad_v; // units vcode/rcode
                    rhs_direct += tmp * kernel.dv[k] * kernel.dp[k_v]; // units vcode*rcode
                    grad_v_mag += grad_v*grad_v; // units (vcode/rcode)^2
                    b_tensor_dot_face += tmp * Face_Area_Vec[k]; // physical
                }
                cmag[k_v] = b_tensor_dot_face / B_interface_mag2; // physical
            }
            rhs /= B_interface_mag2;
            grad_v_mag = sqrt(grad_v_mag);
            bhat_dot_gradvhat = rhs / grad_v_mag; // physical (dimensionless)
            bhat_dot_gradvhat_direct = 3.*rinv*rinv * (rhs_direct / B_interface_mag2); // units vcode/rcode
            /* ok now just multipy the scalar contraction of the B tensor and shear tensor to get the fluxes */
            for(k_v=0;k_v<3;k_v++) {cmag[k_v] *= 3.*eta*rhs;} // units vcode/rcode
        }
#endif
        
        /* standard Navier-Stokes equations: viscosity is decomposed into the shear and bulk viscosity terms */
        double cmag_dir[3];
        if(do_non_mhd==1)
        {
            double divv_i=0,divv_j=0;
            for(k=0;k<3;k++)
            {
                divv_i += local.Gradients.Velocity[k][k];
                divv_j += SphP[j].Gradients.Velocity[k][k];
            }
            double divv = (wt_i*divv_i + wt_j*divv_j) * (zeta - eta*(2./3.));
            wt_i*=eta; wt_j*=eta;
            
            double Pxx = 2*(wt_i*local.Gradients.Velocity[0][0]+wt_j*SphP[j].Gradients.Velocity[0][0]) + divv;
            double Pyy = 2*(wt_i*local.Gradients.Velocity[1][1]+wt_j*SphP[j].Gradients.Velocity[1][1]) + divv;
            double Pzz = 2*(wt_i*local.Gradients.Velocity[2][2]+wt_j*SphP[j].Gradients.Velocity[2][2]) + divv;
            double Pxy = wt_i*(local.Gradients.Velocity[0][1]+local.Gradients.Velocity[1][0]) +
                         wt_j*(SphP[j].Gradients.Velocity[0][1]+SphP[j].Gradients.Velocity[1][0]);
            double Pxz = wt_i*(local.Gradients.Velocity[0][2]+local.Gradients.Velocity[2][0]) +
                         wt_j*(SphP[j].Gradients.Velocity[0][2]+SphP[j].Gradients.Velocity[2][0]);
            double Pyz = wt_i*(local.Gradients.Velocity[1][2]+local.Gradients.Velocity[2][1]) +
                         wt_j*(SphP[j].Gradients.Velocity[1][2]+SphP[j].Gradients.Velocity[2][1]);
            
            cmag[0] = Pxx*Face_Area_Vec[0] + Pxy*Face_Area_Vec[1] + Pxz*Face_Area_Vec[2]; // units vcode/rcode
            cmag[1] = Pxy*Face_Area_Vec[0] + Pyy*Face_Area_Vec[1] + Pyz*Face_Area_Vec[2];
            cmag[2] = Pxz*Face_Area_Vec[0] + Pyz*Face_Area_Vec[1] + Pzz*Face_Area_Vec[2];
            
            double dv_dir = (zeta-eta*2./3.)*(kernel.dv[0]*kernel.dp[0]+kernel.dv[1]*kernel.dp[1]+kernel.dv[2]*kernel.dp[2]);
            double Pxx_direct = eta*2.*kernel.dv[0]*kernel.dp[0] + dv_dir;
            double Pyy_direct = eta*2.*kernel.dv[1]*kernel.dp[1] + dv_dir;
            double Pzz_direct = eta*2.*kernel.dv[2]*kernel.dp[2] + dv_dir;
            double Pxy_direct = eta*(kernel.dv[0]*kernel.dp[1]+kernel.dv[1]*kernel.dp[0]);
            double Pxz_direct = eta*(kernel.dv[0]*kernel.dp[2]+kernel.dv[2]*kernel.dp[0]);
            double Pyz_direct = eta*(kernel.dv[2]*kernel.dp[1]+kernel.dv[1]*kernel.dp[2]);
            cmag_dir[0] = (Pxx_direct*Face_Area_Vec[0] + Pxy_direct*Face_Area_Vec[1] + Pxz_direct*Face_Area_Vec[2])/Face_Area_Norm*kernel.r;
            cmag_dir[1] = (Pxy_direct*Face_Area_Vec[0] + Pyy_direct*Face_Area_Vec[1] + Pyz_direct*Face_Area_Vec[2])/Face_Area_Norm*kernel.r;
            cmag_dir[2] = (Pxz_direct*Face_Area_Vec[0] + Pyz_direct*Face_Area_Vec[1] + Pzz_direct*Face_Area_Vec[2])/Face_Area_Norm*kernel.r;
        }
        
        /* slope-limit this to be sure that viscosity always acts in the proper direction when there is local noise */
        double rho_i = local.Density*All.cf_a3inv, rho_j = SphP[j].Density*All.cf_a3inv, rho_ij=0.5*(rho_i+rho_j);

        /* convert to physical units */
        for(k_v=0;k_v<3;k_v++) {cmag[k_v] *= All.cf_a2inv;}
        
        for(k_v=0;k_v<3;k_v++)
        {
#ifdef MAGNETIC
            double b_hll_eff = DMAX(DMIN(1. , 3.*bhat_dot_gradvhat*bhat_dot_gradvhat) , 0.01);
            double dv_visc = 3. * (B_interface[k_v]*Bi_proj - kernel.dp[k_v]/3.) * bhat_dot_gradvhat_direct * a_inv; // physical
            double thold_ptot_hll = 0.1 * exp(-2. * grad_v_mag * kernel.r / (1.e-30 + fabs(kernel.dv[k_v]))); // units ok
#else
            double b_hll_eff=1;
            double dv_visc = cmag_dir[k_v] * rinv*rinv / (All.cf_atime * DMAX(eta,zeta)); // physical
            double thold_ptot_hll = 0.03;
#endif
            /* obtain HLL correction terms for Reimann problem solution */
            double hll_tmp = rho_ij * HLL_correction(dv_visc,-dv_visc,rho_ij,viscous_wt_physical); // physical
            double fluxlimiter_absnorm = -DMAX(eta,zeta) * sqrt(b_hll_eff) * Face_Area_Norm * dv_visc*rinv*a_inv; // physical
            double ptot = DMIN(local.Mass,P[j].Mass)*sqrt(kernel.dv[0]*kernel.dv[0]+
                                                          kernel.dv[1]*kernel.dv[1]+
                                                          kernel.dv[2]*kernel.dv[2]) / (1.e-37 + dt_hydrostep) * a_inv; // physical
            double thold_hll = 0.1*ptot;
            if(fabs(hll_tmp) > thold_hll) {hll_tmp *= thold_hll/fabs(hll_tmp);}
            hll_tmp *= b_hll_eff;
            if(cmag[k_v] * hll_tmp <= 0)
            {
                thold_hll = b_hll_eff * thold_ptot_hll * ptot;
                if(fabs(hll_tmp) > thold_hll) {hll_tmp *= thold_hll / fabs(hll_tmp);}
            } else {
                thold_hll = b_hll_eff * DMAX(fabs(0.5*cmag[k_v]) , 0.3*thold_ptot_hll*ptot);
                if(fabs(hll_tmp) > thold_hll) {hll_tmp *= thold_hll / fabs(hll_tmp);}
            }
            double cmag_corr = cmag[k_v] + hll_tmp;
            cmag[k_v] = MINMOD(cmag[k_v], cmag_corr);
            double check_for_stability_sign = fluxlimiter_absnorm*cmag[k_v];
            if((check_for_stability_sign < 0) && (fabs(fluxlimiter_absnorm) > HLL_DIFFUSION_OVERSHOOT_FACTOR*fabs(cmag[k_v]))) {cmag[k_v] = 0;}
            
#if defined(GALSF)
            if(check_for_stability_sign < 0) {cmag[k_v]=0;}
#endif
        }
        
        double v_dot_dv=0; for(k=0;k<3;k++) {v_dot_dv += kernel.dv[k] * cmag[k] * a_inv;} // physical
        if(v_dot_dv>0)
        {
            for(k=0;k<3;k++) {cmag[k] = 0;}
        } else {
            double KE_com=0; for(k=0;k<3;k++) {KE_com += kernel.dv[k]*kernel.dv[k] * All.cf_a2inv;} // physical
            KE_com *= 0.25 * (local.Mass + P[j].Mass);
            double dKE_q = fabs(v_dot_dv) * dt_hydrostep / (1.e-40 + KE_com);
            double threshold_tmp = 1.0;
            double lim_corr=1; if(dKE_q > threshold_tmp) {lim_corr = threshold_tmp/dKE_q;}
            for(k=0;k<3;k++) {cmag[k] *= lim_corr;}
        }
        
        
        /* now add a flux-limiter to prevent overshoot (even when the directions are correct) */
        double cmag_E = cmag[0]*v_interface[0] + cmag[1]*v_interface[1] + cmag[2]*v_interface[2];
        if(dt_hydrostep > 0)
        {
            double cmag_lim = 0.25 * (local.Mass + P[j].Mass) * (v_interface[0]*v_interface[0]+v_interface[1]*v_interface[1]+v_interface[2]*v_interface[2]) / dt_hydrostep;
            if(fabs(cmag_E) > cmag_lim)
            {
                double corr_visc = cmag_lim / fabs(cmag_E);
                cmag[0]*=corr_visc; cmag[1]*=corr_visc; cmag[2]*=corr_visc; cmag_E*=corr_visc;
            }
        }
        /* ok now we can finally add this to the numerical fluxes */
        for(k=0;k<3;k++) {Fluxes.v[k] += cmag[k];}
        Fluxes.p += cmag_E;
        
    } // close check that kappa and particle masses are positive
}
