/* quantum pressure-tensor computation to calculate the relevant fluxes between faces of particles, within the gravity routine */
#ifdef DM_FUZZY

#define DM_FUZZY_USE_SIMPLER_HLL_SOLVER 0    /* determines which solver will be used for DM_FUZZY=0; =1 is the newer, simpler, but more diffusive solver */


if((local.Type==1) && (P[j].Type==1)) // only acts between DM particles of type 1 (can easily change if desired)
{
    /* since this isn't a super-expensive calculation, and we need to walk the gravity tree for both 'sides' anyways, we will effectively do this twice each timestep */
    double V_i=local.V_i, V_j=get_particle_volume_ags(j), wt_i=V_i, wt_j=V_j, Face_Area_Vec[3], Face_Area_Norm=0, vface_i_minus_j=0, dv[3]; // calculate densities (in physical units)
    // calculate effective faces and face velocity between elements //
    if((fabs(V_i-V_j)/DMIN(V_i,V_j))/NUMDIMS > 1.25) {wt_i=wt_j=2.*V_i*V_j/(V_i+V_j);} else {wt_i=V_i; wt_j=V_j;} // limiter to prevent faces from going geometrically crazy in disparate limits
    for(k=0;k<3;k++)
    {
        Face_Area_Vec[k] = (kernel.wk_i*wt_i * (local.NV_T[k][0]*kernel.dp[0] + local.NV_T[k][1]*kernel.dp[1] + local.NV_T[k][2]*kernel.dp[2]) +
                            kernel.wk_j*wt_j * (P[j].NV_T[k][0]*kernel.dp[0] + P[j].NV_T[k][1]*kernel.dp[1] + P[j].NV_T[k][2]*kernel.dp[2])) * All.cf_atime*All.cf_atime; // physical units
        Face_Area_Norm += Face_Area_Vec[k]*Face_Area_Vec[k]; // physical units
        dv[k] = kernel.dv[k] / All.cf_atime; // physical units: dp and dv = local - j = R - L, always //
        vface_i_minus_j += dv[k] * Face_Area_Vec[k]; // physical units
    }
    Face_Area_Norm = sqrt(Face_Area_Norm); vface_i_minus_j /= Face_Area_Norm;
    // convert everything needed below into physical units //
    double igrad[3], jgrad[3], i2grad[3][3], j2grad[3][3], fac_g = All.cf_a3inv/All.cf_atime, fac_g2 = All.cf_a3inv*All.cf_a2inv, dp[3]; int m;
    for(k=0;k<3;k++)
    {
        dp[k] = kernel.dp[k] * All.cf_atime;
        igrad[k] = fac_g * local.AGS_Gradients_Density[k];
        jgrad[k] = fac_g * P[j].AGS_Gradients_Density[k];
        for(m=0;m<3;m++)
        {
            i2grad[k][m] = fac_g2 * local.AGS_Gradients2_Density[k][m];
            j2grad[k][m] = fac_g2 * P[j].AGS_Gradients2_Density[k][m];
        }
    }
    
#if (DM_FUZZY == 0)
    double prev_acc = All.G*All.cf_a2inv * P[j].Mass * P[j].OldAcc, flux[3]={0}, dt_egy_Numerical_QuantumPotential=0, m_mean = 0.5*(local.Mass+P[j].Mass), rho_i=local.Mass/V_i*All.cf_a3inv, rho_j=P[j].Mass/V_j*All.cf_a3inv, dt = local.dtime, AGS_Numerical_QuantumPotential = 0.5*(local.AGS_Numerical_QuantumPotential/V_i + P[j].AGS_Numerical_QuantumPotential/V_j)*All.cf_a3inv;
    double HLLwt = (0.5*(kernel.wk_i/kernel.hinv3_i + kernel.wk_j/kernel.hinv3_j)) * (0.5*(kernel.h_i+kernel.h_j)/kernel.r); HLLwt = 10.*HLLwt*HLLwt; // strong dissipation terms allowed for very-close particles, where second-derivative diverges, otherwise weak (no diffusion) //
    // actually compute the fluxes now, this is the key routine, below //

#if (DM_FUZZY_USE_SIMPLER_HLL_SOLVER == 0)
    do_dm_fuzzy_flux_computation_old(HLLwt, dt, m_mean, prev_acc, dp, dv, jgrad, igrad, j2grad, i2grad, rho_j, rho_i, vface_i_minus_j, Face_Area_Vec, flux, AGS_Numerical_QuantumPotential, &dt_egy_Numerical_QuantumPotential);
#else
    do_dm_fuzzy_flux_computation(HLLwt, dt, prev_acc, dv, jgrad, igrad, j2grad, i2grad, rho_j, rho_i, vface_i_minus_j, Face_Area_Vec, flux, P[j].AGS_Numerical_QuantumPotential/V_j*All.cf_a3inv, local.AGS_Numerical_QuantumPotential/V_i*All.cf_a3inv, &dt_egy_Numerical_QuantumPotential);
#endif
    out.AGS_Dt_Numerical_QuantumPotential += dt_egy_Numerical_QuantumPotential; for(k=0;k<3;k++) {out.acc[k] += flux[k] / (local.Mass * All.cf_a2inv);} // assign back to particles

#else

    double h_2m = 0.5*All.ScalarField_hbar_over_mass;
    double Psi_Re_R, Psi_Re_L, d_Psi_Re_R[3], d_Psi_Re_L[3], v_face[3],
           Psi_Im_R, Psi_Im_L, d_Psi_Im_R[3], d_Psi_Im_L[3], Flux_Re=0, Flux_Im=0, Flux_M=0;
    
    for(k=0;k<3;k++) {v_face[k] = 0.5*(local.Vel[k]+P[j].Vel[k]) / All.cf_atime;}
    
    dm_fuzzy_reconstruct_and_slopelimit(&Psi_Re_R, d_Psi_Re_R, &Psi_Re_L, d_Psi_Re_L,
                                        local.AGS_Psi_Re, local.AGS_Gradients_Psi_Re, local.AGS_Gradients2_Psi_Re,
                                        P[j].AGS_Psi_Re_Pred * P[j].AGS_Density / P[j].Mass,
                                        P[j].AGS_Gradients_Psi_Re, P[j].AGS_Gradients2_Psi_Re, dp);

    dm_fuzzy_reconstruct_and_slopelimit(&Psi_Im_R, d_Psi_Im_R, &Psi_Im_L, d_Psi_Im_L,
                                        local.AGS_Psi_Im, local.AGS_Gradients_Psi_Im, local.AGS_Gradients2_Psi_Im,
                                        P[j].AGS_Psi_Im_Pred * P[j].AGS_Density / P[j].Mass,
                                        P[j].AGS_Gradients_Psi_Im, P[j].AGS_Gradients2_Psi_Im, dp);

    double psi2_L = Psi_Re_L*Psi_Re_L + Psi_Im_L*Psi_Im_L, psi2_R = Psi_Re_R*Psi_Re_R + Psi_Im_R*Psi_Im_R;
    /*
    double Adotv = 0; for(k=0;k<3;k++) {Adotv += Face_Area_Vec[k]*v_face[k];}
    double AdotDp_Re_R=0; for(k=0;k<3;k++) {AdotDp_Re_R += Face_Area_Vec[k]*d_Psi_Re_R[k];}
    double AdotDp_Re_L=0; for(k=0;k<3;k++) {AdotDp_Re_L += Face_Area_Vec[k]*d_Psi_Re_L[k];}
    double AdotDp_Im_R=0; for(k=0;k<3;k++) {AdotDp_Im_R += Face_Area_Vec[k]*d_Psi_Im_R[k];}
    double AdotDp_Im_L=0; for(k=0;k<3;k++) {AdotDp_Im_L += Face_Area_Vec[k]*d_Psi_Im_L[k];}
     */
    
    for(k=0;k<3;k++)
    {
        Flux_Re += 0.5 * Face_Area_Vec[k] * ( -h_2m*(d_Psi_Im_R[k]+d_Psi_Im_L[k]) + v_face[k]*(Psi_Re_R+Psi_Re_L) );
        Flux_Im += 0.5 * Face_Area_Vec[k] * (  h_2m*(d_Psi_Re_R[k]+d_Psi_Re_L[k]) + v_face[k]*(Psi_Im_R+Psi_Im_L) );
        Flux_M  += 0.5 * Face_Area_Vec[k] * ( 2.*h_2m * (Psi_Im_R * d_Psi_Re_R[k] - Psi_Re_R * d_Psi_Im_R[k] +
                                                         Psi_Im_L * d_Psi_Re_L[k] - Psi_Re_L * d_Psi_Im_L[k])
                                             + v_face[k]*(psi2_L + psi2_R) );
    }
    double k_eff = 0.3 / kernel.r, cs_eff = 0.5 * h_2m * k_eff, prefac = 0.5 * Face_Area_Norm * cs_eff;
    Flux_Re += prefac * (Psi_Re_R-Psi_Re_L);
    Flux_Im += prefac * (Psi_Im_R-Psi_Im_L);
    Flux_M  += prefac * (psi2_R - psi2_L);

    
    out.AGS_Dt_Psi_Re += Flux_Re;
    out.AGS_Dt_Psi_Im += Flux_Im;
    out.AGS_Dt_Psi_Mass += Flux_M;
#endif

} // total bracket (for variable protection)
#endif


