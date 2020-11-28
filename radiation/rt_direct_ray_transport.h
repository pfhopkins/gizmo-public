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

#if defined(RT_EVOLVE_INTENSITIES)
if(local.Mass>0 && P[j].Mass>0 && dt_hydrostep>0 && Face_Area_Norm>0)
{
    double c_light_eff=C_LIGHT_CODE_REDUCED, rsol_fac=c_light_eff/C_LIGHT_CODE, V_i_invphys=All.cf_a3inv/V_i, V_j_invphys=All.cf_a3inv/V_j, vfluid_minus_vface_dotA=0, cminusv_n_dotA[N_RT_INTENSITY_BINS]={0}, sigma_j=Particle_Size_j*(SphP[j].Density*All.cf_a3inv); int k_freq, k_angle;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION<5)
    double v_frame[3]={0}; for(k=0;k<3;k++) {vfluid_minus_vface_dotA+=0.5*((ParticleVel_j[k]+local.ParticleVel[k])-(local.Vel[k]+VelPred_j[k]))/All.cf_atime * Face_Area_Vec[k];} // frame velocity, not fluid velocity, is what appears here. physical units
#endif
    for(k_angle=0; k_angle<N_RT_INTENSITY_BINS; k_angle++) {for(k=0;k<3;k++) {cminusv_n_dotA[k_angle]+=(c_light_eff*All.Rad_Intensity_Direction[k_angle][k]-rsol_fac*0.5*(local.Vel[k]+VelPred_j[k])/All.cf_atime) * Face_Area_Vec[k];}} // face area dotted into 'residual' speed here [not actual velocity addition, just a mathematical trick here, hence no relativistic terms]. physical units
    
    for(k_freq=0; k_freq<N_RT_FREQ_BINS; k_freq++)
    {
        // following Jiang et al. we reduce the advection speed here when the cell optical depth is large
        double tau_c_j = SphP[j].Rad_Kappa[k_freq] * sigma_j; // = L_particle / (lambda_mean_free_path) = L*kappa*rho (physical) //
        double q_tau = 10. * 0.5*(tau_c_i[k_freq]+tau_c_j), a_tau=1; if(q_tau>3.5) {a_tau=1./q_tau;} else {if(q_tau<0.1) {a_tau=1.-0.25*q_tau*q_tau;} else {a_tau=sqrt(1.-exp(-q_tau*q_tau))/q_tau;}}

        for(k_angle=0; k_angle<N_RT_INTENSITY_BINS; k_angle++)
        {
            double scalar_ij = 0.5*(local.Rad_Intensity_Pred[k_freq][k_angle]*V_i_invphys + SphP[j].Rad_Intensity_Pred[k_freq][k_angle]*V_j_invphys); // physical
            double cmag = scalar_ij * (vfluid_minus_vface_dotA + a_tau*cminusv_n_dotA[k_angle]); // 0th-order flux
            out.Dt_Rad_Intensity[k_freq][k_angle] += cmag;
            if(j_is_active_for_fluxes) {SphP[j].Dt_Rad_Intensity[k_freq][k_angle] -= cmag;}
        }
    }
}
#endif
