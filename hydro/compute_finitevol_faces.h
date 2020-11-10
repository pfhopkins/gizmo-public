    Face_Area_Norm=0; // need this to be set here for any module below


/*  if we are using fixed-grid mode for the code, we manually set the areas to the correct geometric areas, which is trivial */
#ifdef HYDRO_REGULAR_GRID
    Face_Area_Norm = calculate_face_area_for_cartesian_mesh(kernel.dp, rinv, Particle_Size_i, Face_Area_Vec);
#endif


/* this is where we have the voronoi-reconstruction faces in that experimental branch of the code */


/* mesh-free reconstructions (MFM/MFV). if the effective gradient matrix is well-conditioned: we can safely use the consistent EOM here.  */
#if (defined(HYDRO_MESHLESS_FINITE_MASS) || defined(HYDRO_MESHLESS_FINITE_VOLUME)) && !defined(HYDRO_REGULAR_GRID)
    double wt_i,wt_j; wt_i=V_i; wt_j=V_j;
    /* note the 'default' formulation from Lanson and Vila takes wt_i=V_i, wt_j=V_j; but this assumes negligible variation in h between particles;
        it is more accurate to use a centered wt (centered face area), which we get by linear interpolation, in extreme discontinuities of particle separation */
#if (SLOPE_LIMITER_TOLERANCE != 2) && !((defined(HYDRO_FACE_AREA_LIMITER) || !defined(PROTECT_FROZEN_FIRE)) && (HYDRO_FIX_MESH_MOTION >= 5)) // unless using the most aggressive reconstruction, we will limit face-area disparity here //
#if defined(COOLING) || (SLOPE_LIMITER_TOLERANCE==0)
    if((fabs(V_i-V_j)/DMIN(V_i,V_j))/NUMDIMS > 1.25) {wt_i=wt_j=2.*V_i*V_j/(V_i+V_j);} else {wt_i=V_i; wt_j=V_j;} //wt_i=wt_j = 2.*V_i*V_j / (V_i + V_j); // more conservatively, could use DMIN(V_i,V_j), but that is less accurate
#else
    if((fabs(V_i-V_j)/DMIN(V_i,V_j))/NUMDIMS > 1.50) {wt_i=wt_j=(V_i*Particle_Size_j+V_j*Particle_Size_i)/(Particle_Size_i+Particle_Size_j);} else {wt_i=V_i; wt_j=V_j;} //wt_i=wt_j = (V_i*PPP[j].Hsml + V_j*local.Hsml) / (local.Hsml+PPP[j].Hsml); // should these be H, or be -effective sizes- //
#endif
#elif defined(GALSF)
    if( (fabs(log(V_i/V_j)/NUMDIMS) > 1.25) && (kernel.r > local.Hsml || kernel.r > PPP[j].Hsml) ) {wt_i=wt_j=(V_i*Particle_Size_j+V_j*Particle_Size_i)/(Particle_Size_i+Particle_Size_j);}
#endif

    double facenormal_dot_dp = 0;
    for(k=0;k<3;k++)
    {
        Face_Area_Vec[k] = kernel.wk_i * wt_i * (local.NV_T[k][0]*kernel.dp[0] + local.NV_T[k][1]*kernel.dp[1] + local.NV_T[k][2]*kernel.dp[2])
                         + kernel.wk_j * wt_j * (SphP[j].NV_T[k][0]*kernel.dp[0] + SphP[j].NV_T[k][1]*kernel.dp[1] + SphP[j].NV_T[k][2]*kernel.dp[2]);
        Face_Area_Vec[k] *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
        Face_Area_Norm += Face_Area_Vec[k]*Face_Area_Vec[k];
        facenormal_dot_dp += Face_Area_Vec[k] * kernel.dp[k]; /* check that face points same direction as vector normal: should be true for positive-definite (well-conditioned) NV_T */
    }


#if defined(HYDRO_TENSOR_FACE_CORRECTIONS) /* experimental! not really ready for full use yet! */
    {
        double wkwt_i=V_i*kernel.wk_i/All.cf_a2inv, wkwt_j=V_j*kernel.wk_j/All.cf_a2inv, Bi_dot_dx=0, Bj_dot_dx=0;
        for(k=0;k<3;k++) {Bi_dot_dx += local.Tensor_MFM_Face_Corrections[k] * kernel.dp[k]; Bj_dot_dx -= SphP[j].Tensor_MFM_Face_Corrections[k] * kernel.dp[k];}
        wkwt_i *= 1. + Bi_dot_dx; wkwt_j *= 1. + Bj_dot_dx;
        
        Face_Area_Vec[0] = wkwt_i * (local.Tensor_MFM_Face_Corrections[3]*kernel.dp[0] + local.Tensor_MFM_Face_Corrections[4]*kernel.dp[1] + local.Tensor_MFM_Face_Corrections[5]*kernel.dp[2]);
        Face_Area_Vec[1] = wkwt_i * (local.Tensor_MFM_Face_Corrections[4]*kernel.dp[0] + local.Tensor_MFM_Face_Corrections[6]*kernel.dp[1] + local.Tensor_MFM_Face_Corrections[7]*kernel.dp[2]);
        Face_Area_Vec[2] = wkwt_i * (local.Tensor_MFM_Face_Corrections[5]*kernel.dp[0] + local.Tensor_MFM_Face_Corrections[7]*kernel.dp[1] + local.Tensor_MFM_Face_Corrections[8]*kernel.dp[2]);

        Face_Area_Vec[0]+= wkwt_j * (SphP[j].Tensor_MFM_Face_Corrections[3]*kernel.dp[0] + SphP[j].Tensor_MFM_Face_Corrections[4]*kernel.dp[1] + SphP[j].Tensor_MFM_Face_Corrections[5]*kernel.dp[2]);
        Face_Area_Vec[1]+= wkwt_j * (SphP[j].Tensor_MFM_Face_Corrections[4]*kernel.dp[0] + SphP[j].Tensor_MFM_Face_Corrections[6]*kernel.dp[1] + SphP[j].Tensor_MFM_Face_Corrections[7]*kernel.dp[2]);
        Face_Area_Vec[2]+= wkwt_j * (SphP[j].Tensor_MFM_Face_Corrections[5]*kernel.dp[0] + SphP[j].Tensor_MFM_Face_Corrections[7]*kernel.dp[1] + SphP[j].Tensor_MFM_Face_Corrections[8]*kernel.dp[2]);

        Face_Area_Norm=0; facenormal_dot_dp=0; for(k=0;k<3;k++) {Face_Area_Norm+=Face_Area_Vec[k]*Face_Area_Vec[k]; facenormal_dot_dp+=Face_Area_Vec[k]*kernel.dp[k];} /* check that face points same direction as vector normal: should be true for positive-definite (well-conditioned) NV_T */
    }
#endif

    
#if defined(KERNEL_CRK_FACES) /* well-tested at this point */
    {
        // order of Tensor_CRK_Face_Corrections: A, B[3], (dA+A*B)[3], (dA.B+A.dB)[3][3] //
        double wk_ij = 0.5*(kernel.wk_i+kernel.wk_j), dwk_ij = 0.5*(kernel.dwk_i+kernel.dwk_j) / (MIN_REAL_NUMBER + kernel.r);
        double Bi_dot_dx = 0, Bj_dot_dx = 0, dAi_etc_dot_dx[3]={0}, dAj_etc_dot_dx[3]={0};
        for(k=0;k<3;k++)
        {
            Bi_dot_dx +=   local.Tensor_CRK_Face_Corrections[k+1] * kernel.dp[k];
            Bj_dot_dx -= SphP[j].Tensor_CRK_Face_Corrections[k+1] * kernel.dp[k];
            int k_x;
            for(k_x=0;k_x<3;k_x++)
            {
                dAi_etc_dot_dx[k_x] +=   local.Tensor_CRK_Face_Corrections[7+3*k+k_x] * kernel.dp[k];
                dAj_etc_dot_dx[k_x] -= SphP[j].Tensor_CRK_Face_Corrections[7+3*k+k_x] * kernel.dp[k];
            }
        }
        Face_Area_Norm = 0; facenormal_dot_dp = 0;
        for(k=0;k<3;k++)
        {
            double Ai = -V_i*V_j*(    local.Tensor_CRK_Face_Corrections[0] * (1. + Bi_dot_dx) * dwk_ij * (+kernel.dp[k])
                                 + (  local.Tensor_CRK_Face_Corrections[4+k] + dAi_etc_dot_dx[k]) * wk_ij);
            double Aj = -V_i*V_j*(  SphP[j].Tensor_CRK_Face_Corrections[0] * (1. + Bj_dot_dx) * dwk_ij * (-kernel.dp[k])
                                 + (SphP[j].Tensor_CRK_Face_Corrections[4+k] + dAj_etc_dot_dx[k]) * wk_ij);
            Face_Area_Vec[k] = (Ai - Aj) * All.cf_atime*All.cf_atime;
            Face_Area_Norm += Face_Area_Vec[k]*Face_Area_Vec[k];
            facenormal_dot_dp += Face_Area_Vec[k] * kernel.dp[k]; /* check that face points same direction as vector normal: should be true for positive-definite (well-conditioned) NV_T */
        }
    }
#endif


    /* if the effective gradient matrix is ill-conditioned (or not positive-definite!): for stability, we revert to the "RSPH" EOM */
    if((SphP[j].ConditionNumber*SphP[j].ConditionNumber > 1.0e12 + cnumcrit2) || (facenormal_dot_dp < 0))
    {
        Face_Area_Norm = -(wt_i*V_i*kernel.dwk_i + wt_j*V_j*kernel.dwk_j) / kernel.r;
        Face_Area_Norm *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
        Face_Area_Vec[0] = Face_Area_Norm * kernel.dp[0];
        Face_Area_Vec[1] = Face_Area_Norm * kernel.dp[1];
        Face_Area_Vec[2] = Face_Area_Norm * kernel.dp[2];
        Face_Area_Norm = Face_Area_Norm * Face_Area_Norm * r2;
    }

    Face_Area_Norm = sqrt(Face_Area_Norm); /* for use in what follows */


/* this corrects certain reconstructed quantities in the reimann solver, if active, and the reconstruction implies an artificially large volume due to large discontinuity in the h-distribution. still experimental. */
#if defined(HYDRO_FACE_VOLUME_RECONSTRUCTION_CORRECTION)
    double Vi_phys=V_i/All.cf_a3inv, Vj_phys=V_j/All.cf_a3inv, Vi_inv_corr=1., Vj_inv_corr=1., Vol_min = fabs(0.5*facenormal_dot_dp) * All.cf_atime / 3.; // area of an oblique pyramid, origin at one origin point, face at midpoint
    if(Vol_min > DMIN(Vi_phys,Vj_phys)) // minimum volume extrapolated to face is larger than one of the particle volumes, so we will apply a correction to avoid too-large a flux
    {
        if(Vol_min > Vi_phys) {Vi_inv_corr = Vi_phys/Vol_min;}
        if(Vol_min > Vj_phys) {Vj_inv_corr = Vj_phys/Vol_min;}
    }
#endif


/* check if face area exceeds maximum geometric allowed limit (can occur when particles with -very- different Hsml interact at the edge of the kernel, limited to geometric max to prevent numerical instability */
#if (defined(HYDRO_FACE_AREA_LIMITER) || !defined(PROTECT_FROZEN_FIRE)) && (HYDRO_FIX_MESH_MOTION >= 5)
    double Amax = DMIN(Get_Particle_Expected_Area(Particle_Size_i) , Get_Particle_Expected_Area(Particle_Size_j)); // minimum of area "i" or area "j": this subroutine takes care of dimensionality, etc. note inputs are all in -physical- units here
    if(Face_Area_Norm > Amax) {for(k=0;k<3;k++) {Face_Area_Vec[k] *= (Amax/Face_Area_Norm);} Face_Area_Norm = Amax;} /* set the face area to the maximum limit, and reset the face vector as well [ direction is preserved, just area changes] */
#endif



#endif
