#define GAMMA_G0 (GAMMA_DEFAULT)
#define GAMMA_G1 ((GAMMA_G0-1.0)/(2.0*GAMMA_G0))
#define GAMMA_G2 ((GAMMA_G0+1.0)/(2.0*GAMMA_G0))
#define GAMMA_G3 ((2.0*GAMMA_G0/(GAMMA_G0-1.0)))
#define GAMMA_G4 (2.0/(GAMMA_G0-1.0))
#define GAMMA_G5 (2.0/(GAMMA_G0+1.0))
#define GAMMA_G6 ((GAMMA_G0-1.0)/(GAMMA_G0+1.0))
#define GAMMA_G7 (0.5*(GAMMA_G0-1.0))
#define GAMMA_G8 (1.0/GAMMA_G0)
#define GAMMA_G9 (GAMMA_G0-1.0)

#define TOL_ITER 1.e-6
#define NMAX_ITER 1000

#ifdef MHD_CONSTRAINED_GRADIENT
#define DEDNER_IMPLICIT_LIMITER 0.25
#else 
#define DEDNER_IMPLICIT_LIMITER 0.75
#endif

#if defined(TURB_DIFFUSION)
#define SAVE_FACE_DENSITY 1
#endif
#if defined(MAGNETIC)
#define SAVE_FACE_BFIELD 1
#endif
#if defined(VISCOSITY) && defined(MAGNETIC)
//#define SAVE_FACE_VFIELD 1 // for now we use a simpler midpoint velocity; this should be updated for hydro solvers (not just mhd) //
#endif


/*!
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. 
 *   However some of the sub-routines here are adopted from other codes (with permission), in particular
 *   AREPO by Volker Springel and ATHENA by Jim Stone. These sections are identified explicitly in the
 *   code below (anything not identified as such was written by PFH). 
 */


/* --------------------------------------------------------------------------------- */
/* some structures with the conserved variables to pass to/from the Riemann solver */
/* --------------------------------------------------------------------------------- */
struct Input_vec_Riemann
{
    struct Conserved_var_Riemann L;
    struct Conserved_var_Riemann R;
};
struct Riemann_outputs
{
    MyDouble P_M;
    MyDouble S_M;
#ifdef MAGNETIC
    MyDouble B_normal_corrected;
    MyDouble cfast_L;
    MyDouble cfast_R;
#ifdef DIVBCLEANING_DEDNER
    MyDouble phi_normal_mean;
    MyDouble phi_normal_db;
#endif
#endif
#ifdef SAVE_FACE_DENSITY
    MyDouble Face_Density;
#endif
#ifdef SAVE_FACE_BFIELD
    MyDouble Face_B[3];
#endif
#ifdef SAVE_FACE_VFIELD
    MyDouble Face_Vel[3];
#endif
#ifdef TURB_DIFF_METALS
    MyDouble Mdot_estimated;
#endif
    struct Conserved_var_Riemann Fluxes;
};
struct rotation_matrix
{
    MyDouble n[3];
    MyDouble m[3];
    MyDouble p[3];
};


/* --------------------------------------------------------------------------------- */
/* function definitions */
/* --------------------------------------------------------------------------------- */
static inline double actual_slopelimiter(double dQ_1, double dQ_2);
static inline double get_dQ_from_slopelimiter(double dQ_1, MyFloat grad[3], struct kernel_hydra kernel, double rinv);
void Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3], double press_tot_limiter);
double guess_for_pressure(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                          double v_line_L, double v_line_R, double cs_L, double cs_R);
void sample_reimann_standard(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                             double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);
void sample_reimann_vaccum_internal(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                                    double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);
void sample_reimann_vaccum_right(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                                 double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);
void sample_reimann_vaccum_left(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                                double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);
void Riemann_solver_exact(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                            double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R);
void HLLC_fluxes(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                 double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R, double S_L, double S_R);
void get_wavespeeds_and_pressure_star(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                                      double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R,
                                      double *S_L_out, double *S_R_out, double press_tot_limiter);
void Riemann_solver_Rusanov(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                            double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R);
void Riemann_solver_KurganovTadmor_PWK(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                            double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R);
void HLLC_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                         double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R, double press_tot_limiter);
void convert_face_to_flux(struct Riemann_outputs *Riemann_out, double n_unit[3]);
int iterative_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                              double v_line_L, double v_line_R, double cs_L, double cs_R);
void reconstruct_face_states(double Q_i, MyFloat Grad_Q_i[3], double Q_j, MyFloat Grad_Q_j[3],
                             double distance_from_i[3], double distance_from_j[3], double *Q_L, double *Q_R, int mode);
#ifdef MAGNETIC
void HLLD_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double press_tot_limiter);
void rotate_states_to_face(struct Input_vec_Riemann *Riemann_vec, double n_unit[3], struct rotation_matrix *rot_matrix);
void rotate_fluxes_back_to_lab(struct Riemann_outputs *Riemann_out, struct rotation_matrix rot_matrix);
#endif



/* --------------------------------------------------------------------------------- */
/* reconstruction procedure (use to extrapolate from cell/particle centered quantities to faces) */
/*  (reconstruction and slope-limiter from P. Hopkins) */
/* --------------------------------------------------------------------------------- */
void reconstruct_face_states(double Q_i, MyFloat Grad_Q_i[3], double Q_j, MyFloat Grad_Q_j[3],
                             double distance_from_i[3], double distance_from_j[3], double *Q_L, double *Q_R, int mode)
{
    if(mode == 0)
    {
        /* zeroth order reconstruction: this case is trivial */
        *Q_R = Q_i;
        *Q_L = Q_j;
        return;
    }
    /* check for the (also) trivial case of equal values on both sides */
    if(Q_i==Q_j) {*Q_L=*Q_R=Q_i; return;}
    
    /* first order reconstruction */
    *Q_R = Q_i + Grad_Q_i[0]*distance_from_i[0] + Grad_Q_i[1]*distance_from_i[1] + Grad_Q_i[2]*distance_from_i[2];
    *Q_L = Q_j + Grad_Q_j[0]*distance_from_j[0] + Grad_Q_j[1]*distance_from_j[1] + Grad_Q_j[2]*distance_from_j[2];

    /* here we do our slightly-fancy slope-limiting */
    double Qmin,Qmax,Qmed,Qmax_eff,Qmin_eff,fac,Qmed_max,Qmed_min;
#ifdef MAGNETIC
    double fac_minmax;
    double fac_meddev;
    if(mode == 1)
    {
        fac_minmax = 0.5; //0.375; /* 0.5, 0.1 works; 1.0 unstable; 0.75 is stable but begins to 'creep' */
        fac_meddev = 0.375; //0.25; /* 0.25,0.375 work well; 0.5 unstable; 0.44 is on-edge */
    }
    if(mode == 2)
    {
#ifdef MHD_CONSTRAINED_GRADIENT
        fac_minmax = 0.0;
        fac_meddev = 0.0;
#else
        fac_minmax = 0.0;
        fac_meddev = 0.25;
#endif
    }
    if(mode == -1)
    {
        fac_minmax = MHD_CONSTRAINED_GRADIENT_FAC_MINMAX;
        fac_meddev = MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV;
    }
#else
    double fac_minmax = 0.5; /* 0.5, 0.1 works; 1.0 unstable; 0.75 is stable but begins to 'creep' */
    double fac_meddev = 0.375; /* 0.25,0.375 work well; 0.5 unstable; 0.44 is on-edge */
#if (SLOPE_LIMITER_TOLERANCE == 2)
    fac_minmax=0.75;
    fac_meddev=0.40;
#endif
#endif
#if (SLOPE_LIMITER_TOLERANCE == 0)
    fac_minmax=0.0;
    fac_meddev=0.0;
#endif
#if defined(KERNEL_CRK_FACES) && (SLOPE_LIMITER_TOLERANCE > 0)
    fac_minmax=0.75; fac_meddev=0.5; // default to aggressive limiters, but with additional limiter below //
#endif

    
    /* get the max/min vals, difference, and midpoint value */
    Qmed = 0.5*(Q_i+Q_j);
    if(Q_i<Q_j) {Qmax=Q_j; Qmin=Q_i;} else {Qmax=Q_i; Qmin=Q_j;}
    fac = fac_minmax * (Qmax-Qmin);
    if(mode == -1) {fac += MHD_CONSTRAINED_GRADIENT_FAC_MAX_PM * fabs(Qmed);}
    Qmax_eff = Qmax + fac; /* 'overshoot tolerance' */
    Qmin_eff = Qmin - fac; /* 'undershoot tolerance' */
    /* check if this implies a sign from the min/max values: if so, we re-interpret the derivative as a
     logarithmic derivative to prevent sign changes from occurring */
    if(mode > 0)
    {
        if(Qmax<0) {if(Qmax_eff>0) Qmax_eff=Qmax*Qmax/(Qmax-(Qmax_eff-Qmax));} // works with 0.5,0.1 //
        if(Qmin>0) {if(Qmin_eff<0) Qmin_eff=Qmin*Qmin/(Qmin+(Qmin-Qmin_eff));}
    }
    /* also allow tolerance to over/undershoot the exact midpoint value in the reconstruction */
    fac = fac_meddev * (Qmax-Qmin);
    if(mode == -1) {fac += MHD_CONSTRAINED_GRADIENT_FAC_MED_PM * fabs(Qmed);}
    Qmed_max = Qmed + fac;
    Qmed_min = Qmed - fac;
    if(Qmed_max>Qmax_eff) Qmed_max=Qmax_eff;
    if(Qmed_min<Qmin_eff) Qmed_min=Qmin_eff;
    /* now check which side is which and apply these limiters */
    if(Q_i<Q_j)
    {
        if(*Q_R<Qmin_eff) *Q_R=Qmin_eff;
        if(*Q_R>Qmed_max) *Q_R=Qmed_max;
        if(*Q_L>Qmax_eff) *Q_L=Qmax_eff;
        if(*Q_L<Qmed_min) *Q_L=Qmed_min;
#if defined(KERNEL_CRK_FACES)
        if(*Q_R > *Q_L)
        {
            double Q0L = *Q_L, Q0R = *Q_R, Qh = 0.5*(Q0L+Q0R);
            if(Q0R > Q_j)  {if(Qh > Q_j) {*Q_R=Qh; *Q_L=Qh;} else {*Q_R=Q_j; if(Q0L < Q_i) {*Q_L=Q_i;} else {*Q_L=Q0L;}}}
            if(Q0L < Q_i)  {if(Qh < Q_i) {*Q_L=Qh; *Q_R=Qh;} else {*Q_L=Q_i; if(Q0R > Q_j) {*Q_R=Q_j;} else {*Q_R=Q0R;}}}
        }
#endif
    } else {
        if(*Q_R>Qmax_eff) *Q_R=Qmax_eff;
        if(*Q_R<Qmed_min) *Q_R=Qmed_min;
        if(*Q_L<Qmin_eff) *Q_L=Qmin_eff;
        if(*Q_L>Qmed_max) *Q_L=Qmed_max;
#if defined(KERNEL_CRK_FACES)
        if(*Q_R < *Q_L)
        {
            double Q0L = *Q_L, Q0R = *Q_R, Qh = 0.5*(Q0L+Q0R);
            if(Q0L > Q_i)  {if(Qh > Q_i) {*Q_L=Qh; *Q_R=Qh;} else {*Q_L=Q_i; if(Q0R < Q_j) {*Q_R=Q_j;} else {*Q_R=Q0R;}}}
            if(Q0R < Q_j)  {if(Qh < Q_j) {*Q_R=Qh; *Q_L=Qh;} else {*Q_R=Q_j; if(Q0L > Q_i) {*Q_L=Q_i;} else {*Q_L=Q0L;}}}
        }
#endif
    }
#if defined(KERNEL_CRK_FACES) && defined(KERNEL_CRK_FACES_EXPERIMENTAL_SLOPELIMITERS)
    double Q0L=*Q_L, Q0R=*Q_R, dQ_ij=fabs(Q_j-Q_i), dQ_LR=fabs(Q0L-Q0R);
    if(dQ_LR > dQ_ij)
    {
        double Q0=0.5*(Q0L+Q0R), dQ0=0.5*(Q0R-Q0L), alpha = dQ_ij/dQ_LR;
        *Q_R = Q0 + alpha*dQ0; *Q_L = Q0 - alpha*dQ0;
    }
#endif

    /* done! */
}


/* --------------------------------------------------------------------------------- */
/* slope limiter: put other limiters here, will replace all calculations with this */
/*  (not currently used, but optional if we want to use other limiters, cited as noted below) */
/* --------------------------------------------------------------------------------- */
static inline double actual_slopelimiter(double dQ_1, double dQ_2)
{
    //return dQ_2; /* no limiter (unstable, use for tests) */
    if(((dQ_1<0)&&(dQ_2<0))||((dQ_1>0)&&(dQ_2>0)))
    {
        if(dQ_1<0) {return DMAX(dQ_1,dQ_2);} else {return DMIN(dQ_1,dQ_2);} /* Sweby: minmod */
        //return dQ_2 * DMIN(DMIN(2,0.5*(1+dQ_1/dQ_2)),2*dQ_1/dQ_2); /* monotonized central slope limiter */
        //return 2*dQ_1*dQ_2/(dQ_1+dQ_2); /* Van Leer (1974, 2006) limiter */
        //return dQ_2 * DMAX(DMIN(1,2*dQ_1/dQ_2),DMIN(2,dQ_1/dQ_2)); /* Sweby: superbee */
    } else {
        return 0;
    }
}



/* --------------------------------------------------------------------------------- */
/* simple function to get the slope limiter from the difference and gradient vectors */
/*  (not currently used, but optional if we want to use other limiters, cited as noted below) */
/* --------------------------------------------------------------------------------- */
static inline double get_dQ_from_slopelimiter(double dQ_1, MyFloat grad[3], struct kernel_hydra kernel, double rinv)
{
    double dQ_2 = (grad[0]*kernel.dp[0] + grad[1]*kernel.dp[1] + grad[2]*kernel.dp[2]) * rinv;
    return actual_slopelimiter(dQ_1,dQ_2);
}


/* --------------------------------------------------------------------------------- */
/* Top-level Riemann solver routine: call this, it will call sub-routines */
/*  (written by P. Hopkins, this is just a wrapper though for the various sub-routines) */
/* --------------------------------------------------------------------------------- */
void Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3], double press_tot_limiter)
{
    if((Riemann_vec.L.p < 0 && Riemann_vec.R.p < 0)||(Riemann_vec.L.rho < 0)||(Riemann_vec.R.rho < 0))
    {
        printf("FAILURE: Unphysical Inputs to Reimann Solver: Left P/rho=%g/%g, Right P/rho=%g/%g \n",
               Riemann_vec.L.p,Riemann_vec.L.rho,Riemann_vec.R.p,Riemann_vec.R.rho); fflush(stdout);
        Riemann_out->P_M = 0;
        return;
    }
    
    if(All.ComovingIntegrationOn)
    {
        /* first convert the input variables to -PHYSICAL- units so the answer makes sense:
         note that we don't require a Hubble-flow correction, because we're solving at the face */
        int k;
        for(k=0;k<3;k++)
        {
            Riemann_vec.L.v[k] /= All.cf_atime;
            Riemann_vec.R.v[k] /= All.cf_atime;
#ifdef MAGNETIC
            Riemann_vec.L.B[k] *= All.cf_a2inv;
            Riemann_vec.R.B[k] *= All.cf_a2inv;
#endif
        }
        Riemann_vec.L.rho *= All.cf_a3inv;
        Riemann_vec.R.rho *= All.cf_a3inv;
        Riemann_vec.L.p *= All.cf_a3inv / All.cf_afac1;
        Riemann_vec.R.p *= All.cf_a3inv / All.cf_afac1;
#ifdef DIVBCLEANING_DEDNER
        Riemann_vec.L.phi *= All.cf_a3inv;
        Riemann_vec.R.phi *= All.cf_a3inv;
#endif
#ifdef EOS_GENERAL
        Riemann_vec.L.cs *= All.cf_afac3;
        Riemann_vec.R.cs *= All.cf_afac3;
        Riemann_vec.L.u /= All.cf_afac1;
        Riemann_vec.R.u /= All.cf_afac1;
#endif
    }
#ifndef EOS_GENERAL
    /* here we haven't reconstructed the sound speeds and internal energies explicitly, so need to do it from pressure, density */
    Riemann_vec.L.cs = sqrt(GAMMA_G0 * Riemann_vec.L.p / Riemann_vec.L.rho);
    Riemann_vec.R.cs = sqrt(GAMMA_G0 * Riemann_vec.R.p / Riemann_vec.R.rho);
    Riemann_vec.L.u  = Riemann_vec.L.p / (GAMMA_G9 * Riemann_vec.L.rho);
    Riemann_vec.R.u  = Riemann_vec.R.p / (GAMMA_G9 * Riemann_vec.R.rho);
#endif
    
#ifdef MAGNETIC
    struct rotation_matrix rot_matrix;
    /* rotate the state to the face, so we can use the 1D riemann solver */
    rotate_states_to_face(&Riemann_vec, n_unit, &rot_matrix);
    /* send it to our 1D MHD solver, which will try several solutions until it gets a positive pressure */
    HLLD_Riemann_solver(Riemann_vec, Riemann_out, press_tot_limiter);
    /* rotate the returned fluxes back to the lab frame */
    rotate_fluxes_back_to_lab(Riemann_out, rot_matrix);

#else
    
    double cs_L = Riemann_vec.L.cs;
    double cs_R = Riemann_vec.R.cs;
    double h_L = Riemann_vec.L.p/Riemann_vec.L.rho + Riemann_vec.L.u + 0.5*(Riemann_vec.L.v[0]*Riemann_vec.L.v[0]+Riemann_vec.L.v[1]*Riemann_vec.L.v[1]+Riemann_vec.L.v[2]*Riemann_vec.L.v[2]);
    double h_R = Riemann_vec.R.p/Riemann_vec.R.rho + Riemann_vec.R.u + 0.5*(Riemann_vec.R.v[0]*Riemann_vec.R.v[0]+Riemann_vec.R.v[1]*Riemann_vec.R.v[1]+Riemann_vec.R.v[2]*Riemann_vec.R.v[2]);
    double v_line_L = Riemann_vec.L.v[0]*n_unit[0] + Riemann_vec.L.v[1]*n_unit[1] + Riemann_vec.L.v[2]*n_unit[2];
    double v_line_R = Riemann_vec.R.v[0]*n_unit[0] + Riemann_vec.R.v[1]*n_unit[1] + Riemann_vec.R.v[2]*n_unit[2];

#ifdef HYDRO_REPLACE_RIEMANN_KT
    Riemann_solver_KurganovTadmor_PWK(Riemann_vec, Riemann_out, n_unit, v_line_L, v_line_R, cs_L, cs_R, h_L, h_R);
    if((Riemann_out->P_M<=0)||(isnan(Riemann_out->P_M))) /* check if it failed, if so compute HLLC instead */
#endif
    HLLC_Riemann_solver(Riemann_vec, Riemann_out, n_unit, v_line_L, v_line_R, cs_L, cs_R, h_L, h_R, press_tot_limiter);
    /* check if HLLC failed: if so, compute the KT flux instead */
    if((Riemann_out->P_M<0)||(isnan(Riemann_out->P_M)))
       Riemann_solver_KurganovTadmor_PWK(Riemann_vec, Riemann_out, n_unit, v_line_L, v_line_R, cs_L, cs_R, h_L, h_R);
#ifdef EOS_GENERAL
    /* check if HLLC+KT failed: if so, compute the Rusanov flux instead */
    if((Riemann_out->P_M<0)||(isnan(Riemann_out->P_M)))
        Riemann_solver_Rusanov(Riemann_vec, Riemann_out, n_unit, v_line_L, v_line_R, cs_L, cs_R, h_L, h_R);
#else
#if !defined(COOLING) && !defined(GALSF)
    /* go straight to the expensive but exact solver (only for hydro with polytropic eos!) */
    if((Riemann_out->P_M<0)||(isnan(Riemann_out->P_M))||(Riemann_out->P_M>press_tot_limiter))
    {
        Riemann_solver_exact(Riemann_vec, Riemann_out, n_unit, v_line_L, v_line_R, cs_L, cs_R, h_L, h_R);
#ifdef SAVE_FACE_DENSITY
        Riemann_out->Face_Density = 0.5*(Riemann_vec.L.rho + Riemann_vec.R.rho);
#endif
    }
#endif // cooling/galsf
#endif // eos_general
#endif // magnetic

#ifdef TURB_DIFF_METALS
    {
        double rho_k,s_k,v_k,cs_m = DMAX(Riemann_vec.L.cs,Riemann_vec.R.cs);
#ifdef MAGNETIC
        double v_line_L = Riemann_vec.L.v[0]*n_unit[0] + Riemann_vec.L.v[1]*n_unit[1] + Riemann_vec.L.v[2]*n_unit[2];
        double v_line_R = Riemann_vec.R.v[0]*n_unit[0] + Riemann_vec.R.v[1]*n_unit[1] + Riemann_vec.R.v[2]*n_unit[2];
#endif
        if(Riemann_out->S_M>0) {rho_k=Riemann_vec.L.rho; v_k=v_line_L; s_k=DMIN(v_line_L,v_line_R)-cs_m;} else {rho_k=Riemann_vec.R.rho; v_k=v_line_R; s_k=DMAX(v_line_L,v_line_R)+cs_m;}
        if(s_k != Riemann_out->S_M) {Riemann_out->Mdot_estimated=rho_k*Riemann_out->S_M*(s_k-v_k)/(s_k-Riemann_out->S_M);} else {Riemann_out->Mdot_estimated=0;}
    }
#endif // turb_diff_metals
}






/* -------------------------------------------------------------------------------------------------------------- */
/* the HLLC Riemann solver: try this first - it's approximate, but fast, and accurate for our purposes */
/*  (wrapper for sub-routines to evaluate hydro reimann problem) */
/* -------------------------------------------------------------------------------------------------------------- */
/* HLLC: hydro (no MHD) */
void HLLC_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                        double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R, double press_tot_limiter)
{
    double S_L,S_R;
    get_wavespeeds_and_pressure_star(Riemann_vec, Riemann_out, n_unit,  v_line_L, v_line_R, cs_L, cs_R, h_L, h_R, &S_L, &S_R, press_tot_limiter);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    /* check if we have a valid solution, and if so, compute the fluxes */
    if((Riemann_out->P_M>=0)&&(!isnan(Riemann_out->P_M)))
        HLLC_fluxes(Riemann_vec, Riemann_out, n_unit,  v_line_L, v_line_R, cs_L, cs_R, h_L, h_R, S_L, S_R);
#else
    Riemann_out->Fluxes.rho = 0; /* vanishes by definition in this frame */
    int k; for(k=0;k<3;k++) {Riemann_out->Fluxes.v[k] = Riemann_out->P_M * n_unit[k];} /* becomes extremely simple for MFM in this frame */
    Riemann_out->Fluxes.p = Riemann_out->P_M * Riemann_out->S_M; /* becomes extremely simple for MFM in this frame */
#ifdef SAVE_FACE_DENSITY
    if((Riemann_out->S_M==0) || ((Riemann_out->P_M<=0)&&(!isnan(Riemann_out->P_M))))
    {
        Riemann_out->Face_Density = 0.5*(Riemann_vec.L.rho+Riemann_vec.R.rho);
    } else {
        Riemann_out->Face_Density = 0.5 * (Riemann_vec.L.rho*(S_L-v_line_L)/(S_L-Riemann_out->S_M) +
                                           Riemann_vec.R.rho*(S_R-v_line_R)/(S_R-Riemann_out->S_M));
        if(Riemann_out->S_M==S_L) {Riemann_out->Face_Density=Riemann_vec.L.rho;}
        if(Riemann_out->S_M==S_R) {Riemann_out->Face_Density=Riemann_vec.R.rho;}
    }
#endif
#endif
}



/*  Rusanov flux: generally not used, because it's too diffusive. But here if we need it. 
        (this implementation written by P. Hopkins) */
void Riemann_solver_Rusanov(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                            double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R)
{
    /* estimate wave speed and simplest-average intermediate state (Primitive Variable Riemann Solvers approximate Riemann solver from Toro) */
    double S_L, S_R, S_plus, P_M, S_M;
    if((v_line_R - v_line_L) > DMAX(cs_L,cs_R)) // first check for vacuum solution, which is not accounted for in this
    {
        Riemann_out->P_M = P_M = MIN_REAL_NUMBER; Riemann_out->S_M = S_L = S_R = S_plus = S_M = 0;
    } else {
        S_plus = DMAX(DMAX(fabs(v_line_L - cs_L), fabs(v_line_R - cs_R)), DMAX(fabs(v_line_L + cs_L), fabs(v_line_R + cs_R)));
        S_L=-S_plus; S_R=S_plus;
        //double rho_csnd_hat=0.5*(Riemann_vec.L.rho+Riemann_vec.R.rho) * 0.5*(cs_L+cs_R);
        //P_M = 0.5 * ((Riemann_vec.L.p + Riemann_vec.R.p) + (v_line_L-v_line_R) * rho_csnd_hat);
        //S_M = 0.5 * ((v_line_R+v_line_L) + (Riemann_vec.L.p-Riemann_vec.R.p) / (rho_csnd_hat));
        S_M = ((v_line_L*Riemann_vec.L.rho + v_line_R*Riemann_vec.R.rho) + S_plus*(Riemann_vec.L.rho - Riemann_vec.R.rho)) / (Riemann_vec.L.rho + Riemann_vec.R.rho);
        P_M = 0.5 * (Riemann_vec.L.p + Riemann_vec.R.p + (v_line_L - v_line_R) * ((2.*S_plus + v_line_L - v_line_R) * Riemann_vec.L.rho * Riemann_vec.R.rho / (Riemann_vec.L.rho + Riemann_vec.R.rho)));
        Riemann_out->P_M = P_M;
        Riemann_out->S_M = S_M;
    }
    int k;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    if((P_M > 0)&&(!isnan(P_M)))
    {
        /* flux = (1/2) * ( F_L + F_R ) - (S_plus/2) * (Q_R - Q_L) */
        double f_rho_left = Riemann_vec.L.rho * (v_line_L + S_plus);
        double f_rho_right = Riemann_vec.R.rho * (v_line_R - S_plus);
        Riemann_out->Fluxes.rho = 0.5 * (f_rho_left + f_rho_right);
        for(k=0;k<3;k++) {Riemann_out->Fluxes.v[k] = 0.5 * (f_rho_left * Riemann_vec.L.v[k] + f_rho_right * Riemann_vec.R.v[k] +
                                                            (Riemann_vec.L.p + Riemann_vec.R.p) * n_unit[k]);}
        Riemann_out->Fluxes.p = 0.5 * (f_rho_left * h_L + f_rho_right * h_R + S_plus * (Riemann_vec.R.p - Riemann_vec.L.p));
    }
#else
    double f_rho_left = (2.*S_plus + v_line_L - v_line_R) * Riemann_vec.L.rho * Riemann_vec.R.rho / (Riemann_vec.L.rho + Riemann_vec.R.rho); // using solution for frame co-moving w contact discontinuity
    double f_rho_right = -f_rho_left; Riemann_out->Fluxes.rho = 0; // by definition in this frame
    for(k=0;k<3;k++) {Riemann_out->Fluxes.v[k] = 0.5 * (f_rho_left * Riemann_vec.L.v[k] + f_rho_right * Riemann_vec.R.v[k] +
                                                        (Riemann_vec.L.p + Riemann_vec.R.p) * n_unit[k]);}
    Riemann_out->Fluxes.p = 0.5 * (f_rho_left * h_L + f_rho_right * h_R + S_plus * (Riemann_vec.R.p - Riemann_vec.L.p));
#endif
#ifdef SAVE_FACE_DENSITY
    Riemann_out->Face_Density = 0.5 * (Riemann_vec.L.rho*(S_L-v_line_L)/(S_L-Riemann_out->S_M) +
                                       Riemann_vec.R.rho*(S_R-v_line_R)/(S_R-Riemann_out->S_M));
    if(Riemann_out->S_M==S_L) {Riemann_out->Face_Density=Riemann_vec.L.rho;}
    if(Riemann_out->S_M==S_R) {Riemann_out->Face_Density=Riemann_vec.R.rho;}
    if((Riemann_out->P_M<=0)&&(!isnan(Riemann_out->P_M))) {Riemann_out->Face_Density = 0.5*(Riemann_vec.L.rho+Riemann_vec.R.rho);}
#endif
    return;
}




/* generalized Kurganov-Tadmor flux as derived in Panuelos, Wadsley, and Kevlahan [PWK] 'Low Shear Diffusion Central Schemes for Particle Methods' (this implementation written by P. Hopkins).
    This is qualitatively similar to a Rusanov flux, as they are essentially Lax-Friedrichs schemes with some improvements */
void Riemann_solver_KurganovTadmor_PWK(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                            double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R)
{
    /* estimate wave speed using the PWK 'switch' alpha */
    double S_L, S_R, S_M, nu, alpha; int k; nu=0;
#if (SLOPE_LIMITER_TOLERANCE==0) || defined(HYDRO_FACE_AREA_LIMITER) || defined(HYDRO_RIEMANN_KT_UNLIMITED)
    alpha = 1; /* default to the more dissipative but smoother limiter function */
#else
    double delta_threshold = 0.001 * 0.5*(Riemann_vec.L.rho+Riemann_vec.R.rho) * 0.5*(cs_L+cs_R); /* alpha is non-zero only if relative momentum is appreciable fraction of this value */
    double dv2=0,dv[3]; for(k=0;k<3;k++) {dv[k]=Riemann_vec.R.rho*Riemann_vec.R.v[k] - Riemann_vec.L.rho*Riemann_vec.L.v[k]; dv2+=dv[k]*dv[k];} /* calculate relative momentum */
    double dvL=fabs(Riemann_vec.R.rho*v_line_R-Riemann_vec.L.rho*v_line_L); /* difference along normal [switch here is designed to reduce diffusion in shear flows] */
    alpha = dvL / sqrt(delta_threshold*delta_threshold + dv2); /* for either very weak velocity differences, or for almost entirely shear-flows, this cuts off the diffusivity */
#endif
    S_L=alpha*cs_L; if(v_line_L>0) {S_L+=v_line_L;} else {S_L-=v_line_L;} /* fastest left-side wavespeed */
    S_R=alpha*cs_R; if(v_line_R>0) {S_R+=v_line_R;} else {S_R-=v_line_R;} /* fastest right-side wavespeed */
    S_M=DMAX(S_L,S_R);
    Riemann_out->S_M = DMAX(cs_L+fabs(v_line_L) , cs_R+fabs(v_line_R)); /* note this does not have the limiter, otherwise corrections in hydro_core_meshless will be incorrect and unstable */
    double f_rho_left = Riemann_vec.L.rho * (v_line_L + S_M);
    double f_rho_right = Riemann_vec.R.rho * (v_line_R - S_M);
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) /* MFV face vel is null, as we are in frame co-moving with face already */
    Riemann_out->P_M = 0.5 * (Riemann_vec.R.p + Riemann_vec.L.p);
    Riemann_out->Fluxes.rho = 0.5 * (f_rho_left + f_rho_right);
    for(k=0;k<3;k++) {Riemann_out->Fluxes.v[k] = 0.5*(f_rho_left*Riemann_vec.L.v[k] + f_rho_right*Riemann_vec.R.v[k]) + (Riemann_out->P_M)*n_unit[k];}
    Riemann_out->Fluxes.p = 0.5 * (f_rho_left*h_L + f_rho_right*h_R + S_M*(Riemann_vec.R.p - Riemann_vec.L.p));
#else
    double denom = 1/(Riemann_vec.L.rho*v_line_L - Riemann_vec.R.rho*v_line_R + S_M*(Riemann_vec.L.rho+Riemann_vec.R.rho));
    Riemann_out->P_M = (f_rho_left*Riemann_vec.R.p - f_rho_right*Riemann_vec.L.p) * denom;
    Riemann_out->Fluxes.rho = 0; /* vanishes by definition in this frame */
    for(k=0;k<3;k++) {Riemann_out->Fluxes.v[k] = (Riemann_vec.R.v[k]-Riemann_vec.L.v[k])*f_rho_left*f_rho_right*denom + (Riemann_out->P_M)*n_unit[k];}
    Riemann_out->Fluxes.p = (S_M*(f_rho_left*Riemann_vec.R.p + f_rho_right*Riemann_vec.L.p) + (h_R-h_L)*f_rho_left*f_rho_right) * denom;
#endif
#ifdef SAVE_FACE_DENSITY
    Riemann_out->Face_Density = 0.5*(Riemann_vec.L.rho+Riemann_vec.R.rho);
#endif
    return;
}



/* here we obtain wave-speed and pressure estimates for the 'star' region for the HLLC solver or the 
    Lagrangian (contact-wave) method; note we keep trying several methods here in the hopes of eventually getting a
    valid (positive-pressure) solution */
/*  (written by P. Hopkins) */
void get_wavespeeds_and_pressure_star(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                                      double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R,
                                      double *S_L_out, double *S_R_out, double press_tot_limiter)
{
    double S_L, S_R;
    
    /* first, check for vacuum conditions, not accounted for in the standard HLLC scheme */
    if((v_line_R - v_line_L) > DMAX(cs_L,cs_R))
    {
        Riemann_out->P_M = MIN_REAL_NUMBER; Riemann_out->S_M = S_L = S_R = 0;
    } else {
        /* Gaburov: 'simplest' HLLC discretization with weighting scheme */
        double PT_L = Riemann_vec.L.p;
        double PT_R = Riemann_vec.R.p;
        S_L = DMIN(v_line_L,v_line_R) - DMAX(cs_L,cs_R);
        S_R = DMAX(v_line_L,v_line_R) + DMAX(cs_L,cs_R);
        double rho_wt_L = Riemann_vec.L.rho*(S_L-v_line_L);
        double rho_wt_R = Riemann_vec.R.rho*(S_R-v_line_R);
        Riemann_out->S_M = ((PT_R-PT_L) + rho_wt_L*v_line_L - rho_wt_R*v_line_R) / (rho_wt_L - rho_wt_R);
        Riemann_out->P_M = (PT_L*rho_wt_R - PT_R*rho_wt_L + rho_wt_L*rho_wt_R*(v_line_R - v_line_L)) / (rho_wt_R - rho_wt_L);
        if(Riemann_out->P_M <= MIN_REAL_NUMBER) {Riemann_out->P_M = MIN_REAL_NUMBER; Riemann_out->S_M = S_L = S_R = 0;}
        
        if((Riemann_out->P_M <= 0)||(isnan(Riemann_out->P_M))||(Riemann_out->P_M>press_tot_limiter))
        {
            /* failed: compute Roe-averaged values (Roe 1981) [roe-averaging not strictly necessary for HLLC, though it improves accuracy */
            /* note that enthalpy H=(Etotal+P)/d is averaged, with Etotal=Ekinetic+Einternal (Einternal=P/(gamma-1)) */
            double sqrt_rho_L = sqrt(Riemann_vec.L.rho);
            double sqrt_rho_R = sqrt(Riemann_vec.R.rho);
            double sqrt_rho_inv = 1 / (sqrt_rho_L + sqrt_rho_R);
            double vx_roe = (sqrt_rho_L*Riemann_vec.L.v[0] + sqrt_rho_R*Riemann_vec.R.v[0]) * sqrt_rho_inv;
            double vy_roe = (sqrt_rho_L*Riemann_vec.L.v[1] + sqrt_rho_R*Riemann_vec.R.v[1]) * sqrt_rho_inv;
            double vz_roe = (sqrt_rho_L*Riemann_vec.L.v[2] + sqrt_rho_R*Riemann_vec.R.v[2]) * sqrt_rho_inv;
            /* compute velocity along the line connecting the nodes, and max/min wave speeds */
            double v_line_roe = vx_roe*n_unit[0] + vy_roe*n_unit[1] + vz_roe*n_unit[2];
#ifndef EOS_GENERAL
            double h_roe  = (sqrt_rho_L*h_L  + sqrt_rho_R*h_R) * sqrt_rho_inv;
            double cs_roe = sqrt(DMAX(1.e-30, GAMMA_G9*(h_roe - 0.5*(vx_roe*vx_roe+vy_roe*vy_roe+vz_roe*vz_roe))));
#else
            double cs_roe = (sqrt_rho_L*cs_L  + sqrt_rho_R*cs_R) * sqrt_rho_inv;
#endif
            S_R = DMAX(v_line_R + cs_R , v_line_roe + cs_roe);
            S_L = DMIN(v_line_L - cs_L , v_line_roe - cs_roe);
            rho_wt_R =  Riemann_vec.R.rho * (S_R - v_line_R);
            rho_wt_L = -Riemann_vec.L.rho * (S_L - v_line_L); /* note the sign */
            /* contact wave speed (speed at contact surface): */
            Riemann_out->S_M = ((rho_wt_R*v_line_R + rho_wt_L*v_line_L) + (PT_L - PT_R)) / (rho_wt_R + rho_wt_L);
            /* S_M = v_line_L* = v_line_R* = v_line_M --- this is the speed at interface */
            /* contact pressure (pressure at contact surface): */
            Riemann_out->P_M = Riemann_vec.L.rho * (v_line_L-S_L)*(v_line_L-Riemann_out->S_M) + PT_L;
            if(Riemann_out->P_M <= MIN_REAL_NUMBER) {Riemann_out->P_M = MIN_REAL_NUMBER; Riemann_out->S_M = S_L = S_R = 0;}
            /* p_M = p_L* = p_R*  */
            
            if((Riemann_out->P_M <= 0)||(isnan(Riemann_out->P_M))||(Riemann_out->P_M>press_tot_limiter))
            {
                /* failed again! try the simple primitive-variable estimate (as we would for Rusanov) */
                Riemann_out->P_M = 0.5*((PT_L + PT_R) + (v_line_L-v_line_R)*0.25*(Riemann_vec.L.rho+Riemann_vec.R.rho)*(cs_L+cs_R));
                /* compute the new wave speeds from it */
                Riemann_out->S_M = 0.5*(v_line_R+v_line_L) + 2.0*(PT_L-PT_R)/((Riemann_vec.L.rho+Riemann_vec.R.rho)*(cs_L+cs_R));
                double S_plus = DMAX(DMAX(fabs(v_line_L - cs_L), fabs(v_line_R - cs_R)), DMAX(fabs(v_line_L + cs_L), fabs(v_line_R + cs_R)));
                S_L=-S_plus; S_R=S_plus; if(Riemann_out->S_M<S_L) Riemann_out->S_M=S_L; if(Riemann_out->S_M>S_R) Riemann_out->S_M=S_R;
                if(Riemann_out->P_M <= MIN_REAL_NUMBER) {Riemann_out->P_M = MIN_REAL_NUMBER; Riemann_out->S_M = S_L = S_R = 0;}
            }
        }
    }
    *S_L_out = S_L; *S_R_out = S_R;
    return;
} // if((Riemann_out->P_M <= 0)||(isnan(Riemann_out->P_M))) //





/* --------------------------------------------------------------------------------- */
/* ... fluxes from HLLC four-wave sampling ... */
/*  (written by P. Hopkins) */
/* --------------------------------------------------------------------------------- */
void HLLC_fluxes(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                 double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R, double S_L, double S_R)
{
    /* note that we are solving everything in the REST FRAME of the interface, then de-boosting the Riemann solution in the main loop */
    double P_M = Riemann_out->P_M; double S_M = Riemann_out->S_M;
    if((P_M <= 0)||(isnan(P_M))) return;
    
    double nfac,eK,dv2=0,v_line_frame=0; int k;
    if((S_M==S_L)||(S_M==S_R)||(S_M==v_line_frame))
    {
        /* trap for this case, which gives NAN below but is actually very simple */
        Riemann_out->Fluxes.rho = 0;
        Riemann_out->Fluxes.p = Riemann_out->P_M * Riemann_out->S_M;
        for(k=0;k<3;k++)
            Riemann_out->Fluxes.v[k] = Riemann_out->P_M * n_unit[k];
        return;
    }
    if(v_line_frame < S_L)
    {
        Riemann_out->Fluxes.rho = Riemann_vec.L.rho * (v_line_L - v_line_frame);
        Riemann_out->Fluxes.p = Riemann_vec.L.rho * h_L * (v_line_L - v_line_frame) + Riemann_vec.L.p * v_line_frame;
        for(k=0;k<3;k++)
            Riemann_out->Fluxes.v[k] = Riemann_out->Fluxes.rho * Riemann_vec.L.v[k] + Riemann_vec.L.p * n_unit[k];
#ifdef SAVE_FACE_DENSITY
        Riemann_out->Face_Density = Riemann_vec.L.rho;
#endif
    } else {
        if((S_L <= v_line_frame)&&(v_line_frame <= S_M))
        {
            nfac = Riemann_vec.L.rho * (S_L-v_line_L)/(S_L-S_M);
            if(nfac < 0) nfac=0; /* protect against too large expansion estimate */
            Riemann_out->Fluxes.rho = Riemann_vec.L.rho * (v_line_L - S_L) + nfac * (S_L - v_line_frame);
            
            eK = Riemann_vec.L.rho * h_L - Riemann_vec.L.p;
            Riemann_out->Fluxes.p = (Riemann_vec.L.rho * h_L * v_line_L - eK * S_L) + (S_L - v_line_frame) * nfac *
                (eK/Riemann_vec.L.rho + (S_M-v_line_L) * (S_M + Riemann_vec.L.p/(Riemann_vec.L.rho * (S_L-v_line_L))));
            
            dv2 = nfac * (S_L - v_line_frame) * (S_M - v_line_L) + Riemann_vec.L.p;
            for(k=0;k<3;k++)
                Riemann_out->Fluxes.v[k] = Riemann_out->Fluxes.rho * Riemann_vec.L.v[k] + dv2 * n_unit[k];
#ifdef SAVE_FACE_DENSITY
            Riemann_out->Face_Density = nfac;
#endif
        } else {
            if((S_M <= v_line_frame)&&(v_line_frame <= S_R))
            {
                nfac = Riemann_vec.R.rho * (S_R-v_line_R)/(S_R-S_M);
                if(nfac < 0) nfac=0; /* protect against too large expansion estimate */
                Riemann_out->Fluxes.rho = Riemann_vec.R.rho * (v_line_R - S_R) + nfac * (S_R - v_line_frame);
                
                eK = Riemann_vec.R.rho * h_R - Riemann_vec.R.p;
                Riemann_out->Fluxes.p = (Riemann_vec.R.rho * h_R * v_line_R - eK * S_R) + (S_R - v_line_frame) * nfac *
                    (eK/Riemann_vec.R.rho + (S_M-v_line_R) * (S_M + Riemann_vec.R.p/(Riemann_vec.R.rho * (S_R-v_line_R))));
                
                dv2 = nfac * (S_R - v_line_frame) * (S_M - v_line_R) + Riemann_vec.R.p;
                for(k=0;k<3;k++)
                    Riemann_out->Fluxes.v[k] = Riemann_out->Fluxes.rho * Riemann_vec.R.v[k] + dv2 * n_unit[k];
#ifdef SAVE_FACE_DENSITY
                Riemann_out->Face_Density = nfac;
#endif
            } else {
                Riemann_out->Fluxes.rho = Riemann_vec.R.rho * (v_line_R - v_line_frame);
                Riemann_out->Fluxes.p = Riemann_vec.R.rho * h_R * (v_line_R - v_line_frame) + Riemann_vec.R.p * v_line_frame;
                for(k=0;k<3;k++)
                    Riemann_out->Fluxes.v[k] = Riemann_out->Fluxes.rho * Riemann_vec.R.v[k] + Riemann_vec.R.p * n_unit[k];
#ifdef SAVE_FACE_DENSITY
                Riemann_out->Face_Density = Riemann_vec.R.rho;
#endif
            }
        }
    }
    return;
}





/* --------------------------------------------------------------------------------- */
/*  exact Riemann solver here -- deals with all the problematic states! */
/*  (written by V. Springel for AREPO; as are the extensions to the exact solver below) */
/* --------------------------------------------------------------------------------- */
void Riemann_solver_exact(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double n_unit[3],
                       double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R)
{
    /* first, we need to check for all the special/exceptional cases that will cause things to go haywire */
    if((Riemann_vec.L.p == 0 && Riemann_vec.R.p == 0) || (Riemann_vec.L.rho==0 && Riemann_vec.R.rho==0))
    {
        /* we're in a Vaccuum! */
        Riemann_out->P_M = Riemann_out->S_M = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        Riemann_out->Fluxes.rho = Riemann_out->Fluxes.p = Riemann_out->Fluxes.v[0] = Riemann_out->Fluxes.v[1] = Riemann_out->Fluxes.v[2] = 0;
#endif
        return;
    }
    /* the usual situation is here:: */
    if((Riemann_vec.L.rho > 0) && (Riemann_vec.R.rho > 0))
    {
        if(iterative_Riemann_solver(Riemann_vec, Riemann_out, v_line_L, v_line_R, cs_L, cs_R))
        {
            /* this is the 'normal' Reimann solution */
            sample_reimann_standard(0.0,Riemann_vec,Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
        }
        else
        {
            /* ICs lead to vacuum, need to sample vacuum solution */
            sample_reimann_vaccum_internal(0.0,Riemann_vec,Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
        }
    } else {
        /* one of the densities is zero or negative */
        if((Riemann_vec.L.rho<0)||(Riemann_vec.R.rho<0))
            endrun(1234545);
        if(Riemann_vec.L.rho>0)
            sample_reimann_vaccum_right(0.0,Riemann_vec,Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
        if(Riemann_vec.R.rho>0)
            sample_reimann_vaccum_left(0.0,Riemann_vec,Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
    }
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    /* if we got a valid solution, this solver returns face states: need to convert these to fluxes */
    convert_face_to_flux(Riemann_out, n_unit);
#endif
}


/* --------------------------------------------------------------------------------- */
/* part of exact Riemann solver: */
/* left state is a vacuum, but right state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_vaccum_left(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                                double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    double S_R = v_line_R - GAMMA_G4 * cs_R;
#ifndef HYDRO_MESHLESS_FINITE_VOLUME
    /* in this code mode, we are -always- moving with the contact discontinuity so density flux = 0; 
     this constrains where we reside in the solution fan */
    Riemann_out->P_M = 0;
    Riemann_out->S_M = S_R;
    return;
#endif

    if(S_R > S)
    {
        /* vacuum */
        Riemann_out->P_M = 0;
        Riemann_out->S_M = S_R;
        Riemann_out->Fluxes.rho = 0;
    } else {
        /* right fan */
        double S_R_check = v_line_R + cs_R;
        if(S_R_check > S)
        {
            /* rarefaction fan right state */
            double C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
            Riemann_out->P_M = Riemann_vec.R.p * pow(C_eff / cs_R, GAMMA_G3);
            Riemann_out->S_M = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
            Riemann_out->Fluxes.rho = Riemann_vec.R.rho * pow(C_eff / cs_R, GAMMA_G4);
        } else {
            /* right data state */
            Riemann_out->P_M = Riemann_vec.R.p;
            Riemann_out->S_M = v_line_R;
            Riemann_out->Fluxes.rho = Riemann_vec.R.rho;
        }
    }
    Riemann_out->Fluxes.p = Riemann_out->P_M;
    int k;
    for(k=0;k<3;k++)
        Riemann_out->Fluxes.v[k] = Riemann_vec.R.v[k] + (Riemann_out->S_M - v_line_R) * n_unit[k];
    return;
}


/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/* right state is a vacuum, but left state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_vaccum_right(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                                double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    //double S_L = v_line_L - GAMMA_G4 * cs_L;
    double S_L = v_line_L + GAMMA_G4 * cs_L; // above line was a sign error, caught by Bert Vandenbroucke
#ifndef HYDRO_MESHLESS_FINITE_VOLUME
    /* in this code mode, we are -always- moving with the contact discontinuity so density flux = 0;
     this constrains where we reside in the solution fan */
    Riemann_out->P_M = 0;
    Riemann_out->S_M = S_L;
    return;
#endif
    
    if(S_L < S)
    {
        /* vacuum */
        Riemann_out->P_M = 0;
        Riemann_out->S_M = S_L;
        Riemann_out->Fluxes.rho = 0;
    } else {
        /* left fan */
        double S_L_check = v_line_L - cs_L;
        if(S_L_check < S)
        {
            /* rarefaction fan left state */
            double C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
            Riemann_out->P_M = Riemann_vec.L.p * pow(C_eff / cs_L, GAMMA_G3);
            Riemann_out->S_M = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
            Riemann_out->Fluxes.rho = Riemann_vec.L.rho * pow(C_eff / cs_L, GAMMA_G4);
        } else {
            /* left data state */
            Riemann_out->P_M = Riemann_vec.L.p;
            Riemann_out->S_M = v_line_L;
            Riemann_out->Fluxes.rho = Riemann_vec.L.rho;
        }
    }
    Riemann_out->Fluxes.p = Riemann_out->P_M;
    int k;
    for(k=0;k<3;k++)
        Riemann_out->Fluxes.v[k] = Riemann_vec.L.v[k] + (Riemann_out->S_M - v_line_L) * n_unit[k];
    return;
}



/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/* solution generations a vacuum inside the fan: sample the vacuum appropriately */
/*   (note that these solutions are identical to the left/right solutions above) */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_vaccum_internal(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                                    double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    double S_L = v_line_L + GAMMA_G4 * cs_L;
    double S_R = v_line_R - GAMMA_G4 * cs_R;
    if(S <= S_L)
    {
        /* left fan */
        sample_reimann_vaccum_right(S,Riemann_vec,Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
    }
    else if(S >= S_R)
    {
        /* right fan */
        sample_reimann_vaccum_left(S,Riemann_vec,Riemann_out,n_unit,v_line_L,v_line_R,cs_L,cs_R);
    }
    else
    {
        /* vacuum in between */
        Riemann_out->P_M = 0;
        Riemann_out->S_M = S;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        Riemann_out->Fluxes.rho = 0;
        Riemann_out->Fluxes.p = Riemann_out->P_M;
        int k;
        for(k=0;k<3;k++)
            Riemann_out->Fluxes.v[k] = (Riemann_vec.L.v[k] + (Riemann_vec.R.v[k]-Riemann_vec.L.v[k]) * (S-S_L)/(S_R-S_L)) *
            (1-n_unit[k]) + S * n_unit[k];
#endif
    }
}




/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/*  This is the "normal" Riemann fan, with no vacuum on L or R state! */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_standard(double S, struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out,
                             double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
#ifndef HYDRO_MESHLESS_FINITE_VOLUME
    /* we don't actually need to evaluate the fluxes, and we already have P_M and S_M, which define the 
     contact discontinuity where the rho flux = 0; so can simply exit this routine */
    return;
#endif
    int k; double C_eff,S_eff;
    if(S <= Riemann_out->S_M)  /* sample point is left of contact discontinuity */
    {
        if(Riemann_out->P_M <= Riemann_vec.L.p)	/* left fan (rarefaction) */
        {
            double S_check_L = v_line_L - cs_L;
            if(S <= S_check_L) /* left data state */
            {
                Riemann_out->Fluxes.p = Riemann_vec.L.p;
                Riemann_out->Fluxes.rho = Riemann_vec.L.rho;
                for(k=0;k<3;k++)
                    Riemann_out->Fluxes.v[k] = Riemann_vec.L.v[k];
                return;
            }
            else
            {
                double C_eff_L = cs_L * pow(Riemann_out->P_M / Riemann_vec.L.p, GAMMA_G1);
                double S_tmp_L = Riemann_out->S_M - C_eff_L;
                
                if(S > S_tmp_L)	/* middle left state */
                {
                    Riemann_out->Fluxes.rho = Riemann_vec.L.rho * pow(Riemann_out->P_M / Riemann_vec.L.p, GAMMA_G8);
                    Riemann_out->Fluxes.p = Riemann_out->P_M;
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.L.v[k] + (Riemann_out->S_M-v_line_L)*n_unit[k];
                    return;
                }
                else		/* left state inside fan */
                {
                    S_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
                    C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
                    Riemann_out->Fluxes.rho = Riemann_vec.L.rho * pow(C_eff / cs_L, GAMMA_G4);
                    Riemann_out->Fluxes.p = Riemann_vec.L.p * pow(C_eff / cs_L, GAMMA_G3);
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.L.v[k] + (S_eff-v_line_L)*n_unit[k];
                    return;
                }
            }
        }
        else			/* left shock */
        {
            if(Riemann_vec.L.p > 0)
            {
                double pml = Riemann_out->P_M / Riemann_vec.L.p;
                double S_L = v_line_L - cs_L * sqrt(GAMMA_G2 * pml + GAMMA_G1);
                
                if(S <= S_L)	/* left data state */
                {
                    Riemann_out->Fluxes.p = Riemann_vec.L.p;
                    Riemann_out->Fluxes.rho = Riemann_vec.L.rho;
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.L.v[k];
                    return;
                }
                else		/* middle left state behind shock */
                {
                    Riemann_out->Fluxes.rho = Riemann_vec.L.rho * (pml + GAMMA_G6) / (pml * GAMMA_G6 + 1.0);
                    Riemann_out->Fluxes.p = Riemann_out->P_M;
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.L.v[k] + (Riemann_out->S_M-v_line_L)*n_unit[k];
                    return;
                }
            }
            else
            {
                Riemann_out->Fluxes.rho = Riemann_vec.L.rho / GAMMA_G6;
                Riemann_out->Fluxes.p = Riemann_out->P_M;
                for(k=0;k<3;k++)
                    Riemann_out->Fluxes.v[k] = Riemann_vec.L.v[k] + (Riemann_out->S_M-v_line_L)*n_unit[k];
                return;
            }
        }
    }
    else    /* sample point is right of contact discontinuity */
    {
        if(Riemann_out->P_M > Riemann_vec.R.p)	/* right shock */
        {
            if(Riemann_vec.R.p > 0)
            {
                double pmr = Riemann_out->P_M / Riemann_vec.R.p;
                double S_R = v_line_R + cs_R * sqrt(GAMMA_G2 * pmr + GAMMA_G1);
                
                if(S >= S_R)	/* right data state */
                {
                    Riemann_out->Fluxes.p = Riemann_vec.R.p;
                    Riemann_out->Fluxes.rho = Riemann_vec.R.rho;
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.R.v[k];
                    return;
                }
                else		/* middle right state behind shock */
                {
                    Riemann_out->Fluxes.rho = Riemann_vec.R.rho * (pmr + GAMMA_G6) / (pmr * GAMMA_G6 + 1.0);
                    Riemann_out->Fluxes.p = Riemann_out->P_M;
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.R.v[k] + (Riemann_out->S_M-v_line_R)*n_unit[k];
                    return;
                }
            }
            else
            {
                Riemann_out->Fluxes.rho = Riemann_vec.R.rho / GAMMA_G6;
                Riemann_out->Fluxes.p = Riemann_out->P_M;
                for(k=0;k<3;k++)
                    Riemann_out->Fluxes.v[k] = Riemann_vec.R.v[k] + (Riemann_out->S_M-v_line_R)*n_unit[k];
                return;
            }
        }
        else			/* right fan */
        {
            double S_check_R = v_line_R + cs_R;
            if(S >= S_check_R)		/* right data state */
            {
                Riemann_out->Fluxes.p = Riemann_vec.R.p;
                Riemann_out->Fluxes.rho = Riemann_vec.R.rho;
                for(k=0;k<3;k++)
                    Riemann_out->Fluxes.v[k] = Riemann_vec.R.v[k];
                return;
            }
            else
            {
                double C_eff_R = cs_R * pow(Riemann_out->P_M / Riemann_vec.R.p, GAMMA_G1);
                double S_tmp_R = Riemann_out->S_M + C_eff_R;

                if(S <= S_tmp_R)	/* middle right state */
                {
                    Riemann_out->Fluxes.rho = Riemann_vec.R.rho * pow(Riemann_out->P_M / Riemann_vec.R.p, GAMMA_G8);
                    Riemann_out->Fluxes.p = Riemann_out->P_M;
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.R.v[k] + (Riemann_out->S_M-v_line_R)*n_unit[k];
                    return;
                }
                else		/* fan right state */
                {
                    S_eff = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
                    C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
                    Riemann_out->Fluxes.rho = Riemann_vec.R.rho * pow(C_eff / cs_R, GAMMA_G4);
                    Riemann_out->Fluxes.p = Riemann_vec.R.p * pow(C_eff / cs_R, GAMMA_G3);
                    for(k=0;k<3;k++)
                        Riemann_out->Fluxes.v[k] = Riemann_vec.R.v[k] + (S_eff-v_line_R)*n_unit[k];
                    return;
                }
            }
        }
    }
}



/* --------------------------------------------------------------------------------- */
/* the exact (iterative) Riemann solver: this is slower, but exact. 
 however there is a small chance of the iteration diverging,
 so we still cannot completely gaurantee a valid solution */
/*  (written by P. Hopkins; however this is adapted from the iterative solver in 
        ATHENA by J. Stone) */
/* --------------------------------------------------------------------------------- */
int iterative_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    /* before going on, let's compare this to an exact Riemann solution calculated iteratively */
    double Pg,Pg_prev,W_L,W_R,Z_L,Z_R,tol,pratio; int niter_Riemann=0;
    double a0,a1,a2,dvel,check_vel;
    dvel = v_line_R - v_line_L;
    check_vel = GAMMA_G4 * (cs_R + cs_L) - dvel;
    /* if check_vel<0, this will produce a vacuum: need to use vacuum-specific subroutine */
    if(check_vel < 0) return 0;
    
    tol=100.0;
    Pg = guess_for_pressure(Riemann_vec, Riemann_out, v_line_L, v_line_R, cs_L, cs_R);
    while((tol>TOL_ITER)&&(niter_Riemann<NMAX_ITER))
    {
        Pg_prev=Pg;
        if(Pg>Riemann_vec.L.p)
        {
            /* shock wave */
            a0 = GAMMA_G5 / Riemann_vec.L.rho;
            a1 = GAMMA_G6 * Riemann_vec.L.p;
            a2 = sqrt(a0 / (Pg+a1));
            W_L = (Pg-Riemann_vec.L.p) * a2;
            Z_L = a2 * (1.0 - 0.5*(Pg-Riemann_vec.L.p)/(a1+Pg));
        } else {
            /* rarefaction wave */
            pratio = Pg / Riemann_vec.L.p;
            W_L = GAMMA_G4 * cs_L * (pow(pratio, GAMMA_G1)-1);
            Z_L = 1 / (Riemann_vec.L.rho*cs_L) * pow(Pg/Riemann_vec.L.p, -GAMMA_G2);
        }
        if(Pg>Riemann_vec.R.p)
        {
            /* shock wave */
            a0 = GAMMA_G5 / Riemann_vec.R.rho;
            a1 = GAMMA_G6 * Riemann_vec.R.p;
            a2 = sqrt(a0 / (Pg+a1));
            W_R = (Pg-Riemann_vec.R.p) * a2;
            Z_R = a2 * (1.0 - 0.5*(Pg-Riemann_vec.R.p)/(a1+Pg));
        } else {
            /* rarefaction wave */
            pratio = Pg / Riemann_vec.R.p;
            W_R = GAMMA_G4 * cs_R * (pow(pratio, GAMMA_G1)-1);
            Z_R = 1 / (Riemann_vec.R.rho*cs_R) * pow(pratio, -GAMMA_G2);
        }
        if(niter_Riemann < NMAX_ITER / 2)
            Pg -= (W_L + W_R + dvel) / (Z_L + Z_R);
        else
            Pg -= 0.5 * (W_L + W_R + dvel) / (Z_L + Z_R);
        
        if(Pg < 0.1 * Pg_prev)
            Pg = 0.1 * Pg_prev;
        
        tol = 2.0 * fabs((Pg-Pg_prev)/(Pg+Pg_prev));
        niter_Riemann++;
    }
    if(niter_Riemann<NMAX_ITER)
    {
        Riemann_out->P_M = Pg;
        Riemann_out->S_M = 0.5*(v_line_L+v_line_R) + 0.5*(W_R-W_L);
        return 1;
    } else {
        return 0;
    }
}



/* --------------------------------------------------------------------------------- */
/* get a pressure guess to begin iteration, for the iterative exact Riemann solver(s) */
/*   (written by V. Springel for AREPO, with minor modifications) */
/* --------------------------------------------------------------------------------- */
double guess_for_pressure(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    double pmin, pmax;
    /* start with the usual lowest-order guess for the contact wave pressure */
    double pv = 0.5*(Riemann_vec.L.p+Riemann_vec.R.p) - 0.125*(v_line_R-v_line_L)*(Riemann_vec.L.p+Riemann_vec.R.p)*(cs_L+cs_R);
    pmin = DMIN(Riemann_vec.L.p,Riemann_vec.R.p);
    pmax = DMAX(Riemann_vec.L.p,Riemann_vec.R.p);
    
    /* if one side is vacuum, guess half the mean */
    if(pmin<=0)
        return 0.5*(pmin+pmax);

    /* if the two are sufficiently close, and pv is between both values, return it */
    double qrat = pmax / pmin;
    if(qrat <= 2.0 && (pmin <= pv && pv <= pmax))
        return pv;
    
    if(pv < pmin)
    {
        /* use two-rarefaction solution */
        double pnu = (cs_L+cs_R) - GAMMA_G7 * (v_line_R - v_line_L);
        double pde = cs_L / pow(Riemann_vec.L.p, GAMMA_G1) + cs_R / pow(Riemann_vec.R.p, GAMMA_G1);
        return pow(pnu / pde, GAMMA_G3);
    }
    else
    {
        /* two-shock approximation  */
        double gel = sqrt((GAMMA_G5 / Riemann_vec.L.rho) / (GAMMA_G6 * Riemann_vec.L.p + pv));
        double ger = sqrt((GAMMA_G5 / Riemann_vec.R.rho) / (GAMMA_G6 * Riemann_vec.R.p + pv));
        double x = (gel * Riemann_vec.L.p + ger * Riemann_vec.R.p - (v_line_R - v_line_L)) / (gel + ger);
        if(x < pmin || x > pmax)
            x = pmin;
        return x;
    }
}


/* -------------------------------------------------------------------------------------------------------------- */
/*  Part of exact Riemann solver: */
 /*    take the face state we have calculated from the exact Riemann solution and get the corresponding fluxes */
/*   (written by V. Springel for AREPO, with minor modifications) */
 /* -------------------------------------------------------------------------------------------------------------- */
void convert_face_to_flux(struct Riemann_outputs *Riemann_out, double n_unit[3])
{
    double rho, P, v[3], v_line=0, v_frame=0, h=0; int k;
    rho = Riemann_out->Fluxes.rho;
    P = Riemann_out->Fluxes.p;
    for(k=0;k<3;k++)
    {
        v[k] = Riemann_out->Fluxes.v[k];
        v_line += v[k] * n_unit[k];
        h += v[k] * v[k];
    }
    v_line -= v_frame;
    h *= 0.5 * rho; /* h is the kinetic energy density */
    h += (GAMMA_G0/GAMMA_G9) * P; /* now h is the enthalpy */
    /* now we just compute the standard fluxes for a given face state */
    Riemann_out->Fluxes.p = h * v_line;
    Riemann_out->Fluxes.rho = rho * v_line;
    for(k=0;k<3;k++)
        Riemann_out->Fluxes.v[k] = Riemann_out->Fluxes.rho * v[k] + P * n_unit[k];
    return;
}








#ifdef MAGNETIC
/* -------------------------------------------------------------------------------------------------------------- */
/* the MHD Riemann solvers: includes HLLD, HLL, and Lax_Friedrich, which it will try in that order */
/*   (written by P. Hopkins) */
/* -------------------------------------------------------------------------------------------------------------- */
/* HLLD: (MHD) */
void HLLD_Riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, double press_tot_limiter)
{
    double v_frame = 0; /* frame velocity (in the direction of its normal): this will be used below! */
    double SMALL_NUMBER = 1.0e-11;

    /* first, correct the reconstructed B-field to mean value (required by div.B criterion across shocks) */
    /* begin by computing the fast magnetosonic wave speeds */
    int k;
    double S_L,S_R,PT_L,PT_R,e_L,e_R,S_M,P_M,Bx,vxL,vxR,tmp2,
            vdotB_L,vdotB_R,B2_L,B2_R,cs2rho_L,cs2rho_R,tmp,c_eff,Bx2,rho_wt_L,rho_wt_R;
    struct Conserved_var_Riemann U_s,U_ss,V_s, *Interface_State;
    Interface_State = &Riemann_vec.L;
    B2_L = Riemann_vec.L.B[0]*Riemann_vec.L.B[0]+Riemann_vec.L.B[1]*Riemann_vec.L.B[1]+Riemann_vec.L.B[2]*Riemann_vec.L.B[2];
    B2_R = Riemann_vec.R.B[0]*Riemann_vec.R.B[0]+Riemann_vec.R.B[1]*Riemann_vec.R.B[1]+Riemann_vec.R.B[2]*Riemann_vec.R.B[2];

    /* define rho*cs^2 for wavespeeds below: note that for the HLLD solver, this must be limited at the pressure 
        (minimum soundspeed is the isothermal soundspeed), for good behavior of the wavespeed guesses */
    cs2rho_L = DMAX(Riemann_vec.L.p , Riemann_vec.L.cs*Riemann_vec.L.cs * Riemann_vec.L.rho);
    cs2rho_R = DMAX(Riemann_vec.R.p , Riemann_vec.R.cs*Riemann_vec.R.cs * Riemann_vec.R.rho);
    /* left fast speed (before B-interface correction) */
    tmp = cs2rho_L + B2_L;
    Bx2 = Riemann_vec.L.B[0]*Riemann_vec.L.B[0];
    Riemann_out->cfast_L = sqrt(0.5*(tmp + sqrt(DMAX(0,tmp*tmp - 4*cs2rho_L*Bx2))) / Riemann_vec.L.rho);
    /* right fast speed (before B-interface correction) */
    tmp = cs2rho_R + B2_R;
    Bx2 = Riemann_vec.R.B[0]*Riemann_vec.R.B[0];
    Riemann_out->cfast_R = sqrt(0.5*(tmp + sqrt(DMAX(0,tmp*tmp - 4*cs2rho_R*Bx2))) / Riemann_vec.R.rho);
    /* effective sound speed is the maximum of the two */
    c_eff = DMAX(Riemann_out->cfast_L,Riemann_out->cfast_R);
    
    /* use this to compute the corrected face-centered fields (recall Bx = constant through all states) */
    Bx = 0.5*(Riemann_vec.L.B[0]+Riemann_vec.R.B[0]);
    Riemann_out->B_normal_corrected = Bx;

#ifdef DIVBCLEANING_DEDNER
    /* use the solution for the modified Bx, given the action of the phi-field; 
        however this must be slope-limited to ensure the 'corrected' Bx remains stable */
#if !defined(MHD_CONSTRAINED_GRADIENT) || defined(COOLING) || defined(GALSF)
    /* this is the formulation from E. Gaburov assuming -two- wavespeeds (cL and cR); this down-weights the correction term
     under some circumstances, which appears to increase stability */
    double cL = Riemann_out->cfast_L; // may want to try with correction for approach speed; this helps to stabilize things in the cL=cR case
    double cR = Riemann_out->cfast_R;
    double corr_norm = 1;
    double cinv = 1.0 / (cL+cR);
    Riemann_out->B_normal_corrected = cinv*(Riemann_vec.L.B[0]*cL + Riemann_vec.R.B[0]*cR);
    double corr_p = cinv*(Riemann_vec.L.phi-Riemann_vec.R.phi);
    double corr_p_abs = fabs(corr_p);
    double corr_b_abs = DEDNER_IMPLICIT_LIMITER * fabs(Riemann_out->B_normal_corrected);
    if(corr_p_abs > corr_b_abs) {corr_norm *= corr_b_abs/corr_p_abs;}
    Riemann_out->B_normal_corrected += corr_norm * corr_p;
    Riemann_out->phi_normal_mean = corr_norm * cinv * (cL*Riemann_vec.R.phi + cR*Riemann_vec.L.phi);
    Riemann_out->phi_normal_db = corr_norm * cinv * cL * cR * (Riemann_vec.L.B[0]-Riemann_vec.R.B[0]);
#else
    /* this is the 'default' Dedner formulation, which uses a single (maximum) wavespeed for the interface */
    c_eff += DMAX( 0. , Riemann_vec.L.v[0]-Riemann_vec.R.v[0] );
    // add approach speed, since this is a signal velocity (and that's what we use for growing phi)
    double corr_norm = 1.0;
    double corr_p = 0.5*(Riemann_vec.L.phi-Riemann_vec.R.phi)/c_eff;
    double corr_p_abs = fabs(corr_p);
    double corr_b_abs = DEDNER_IMPLICIT_LIMITER * fabs(Bx);
    if(corr_p_abs > corr_b_abs) {corr_norm *= corr_b_abs/corr_p_abs;}
    Riemann_out->B_normal_corrected += corr_norm * corr_p;
    Riemann_out->phi_normal_mean = 0.5 * corr_norm * (Riemann_vec.R.phi+Riemann_vec.L.phi);
    Riemann_out->phi_normal_db = corr_norm * 0.5*(c_eff*(Riemann_vec.L.B[0]-Riemann_vec.R.B[0]));
#endif
#endif
    
    /* and set the normal component of B to the corrected value */
    Bx = Riemann_out->B_normal_corrected; // need to re-set this using the updated value of Bx //
    Riemann_vec.L.B[0] = Bx;
    Riemann_vec.R.B[0] = Bx;
    Bx2 = Bx*Bx;
    
    /* define useful quantities and obtain the wavespeeds */
    vxL = Riemann_vec.L.v[0];
    vxR = Riemann_vec.R.v[0];
    /* compute the total pressure (use the updated B-field values!) */
    B2_L = Riemann_vec.L.B[0]*Riemann_vec.L.B[0]+Riemann_vec.L.B[1]*Riemann_vec.L.B[1]+Riemann_vec.L.B[2]*Riemann_vec.L.B[2];
    B2_R = Riemann_vec.R.B[0]*Riemann_vec.R.B[0]+Riemann_vec.R.B[1]*Riemann_vec.R.B[1]+Riemann_vec.R.B[2]*Riemann_vec.R.B[2];
    PT_L = Riemann_vec.L.p + 0.5*B2_L;
    PT_R = Riemann_vec.R.p + 0.5*B2_R;

    /* re-compute the fast-magnetosonic wave speeds with the corrected Bx value */
    tmp = cs2rho_L + B2_L;
    Riemann_out->cfast_L = sqrt(0.5*(tmp + sqrt(DMAX(0,tmp*tmp - 4*cs2rho_L*Bx2))) / Riemann_vec.L.rho);
    tmp = cs2rho_R + B2_R;
    Riemann_out->cfast_R = sqrt(0.5*(tmp + sqrt(DMAX(0,tmp*tmp - 4*cs2rho_R*Bx2))) / Riemann_vec.R.rho);
    c_eff = DMAX(Riemann_out->cfast_L, Riemann_out->cfast_R);

    /* now make guesses for S_L and S_R */
    S_L = DMIN(vxL,vxR) - c_eff;
    S_R = DMAX(vxL,vxR) + c_eff;
    /* ok, now we can approximate the contact wavespeed and pressure */
    rho_wt_L = Riemann_vec.L.rho*(S_L-vxL);
    rho_wt_R = Riemann_vec.R.rho*(S_R-vxR);
    S_M = ((PT_R-PT_L) + rho_wt_L*vxL - rho_wt_R*vxR) / (rho_wt_L - rho_wt_R);
    P_M = PT_L + rho_wt_L*(S_M-vxL);
    /* P_M = (PT_L*rho_wt_R - PT_R*rho_wt_L + rho_wt_L*rho_wt_R*(vxR - vxL)) / (rho_wt_R - rho_wt_L); */
    
    /* trap for vacuum solution, important! */
    int within_sampled_vacuum = 0; // no vacuum (default)
    if((vxR - vxL) > c_eff) {P_M = MIN_REAL_NUMBER; within_sampled_vacuum = 1;}
    
    /* now we need to check if these guesses are reasonable; otherwise, use a different wavespeed estimate */
    if((P_M<=0)||(isnan(P_M))||(P_M>press_tot_limiter))
    {
        double sqrt_rho_L = sqrt(Riemann_vec.L.rho);
        double sqrt_rho_R = sqrt(Riemann_vec.R.rho);
        double sqrt_rho_inv = 1 / (sqrt_rho_L + sqrt_rho_R);
        double vx_roe = (sqrt_rho_L*vxL + sqrt_rho_R*vxR) * sqrt_rho_inv;
        /* compute velocity along the line connecting the nodes, and max/min wave speeds */
        double cs_roe = (sqrt_rho_L*Riemann_out->cfast_L  + sqrt_rho_R*Riemann_out->cfast_R) * sqrt_rho_inv;
        S_R = DMAX(vxR + Riemann_out->cfast_R , vx_roe + cs_roe);
        S_L = DMIN(vxL - Riemann_out->cfast_L , vx_roe - cs_roe);
        
        rho_wt_L = Riemann_vec.L.rho*(S_L-vxL);
        rho_wt_R = Riemann_vec.R.rho*(S_R-vxR);
        S_M = ((PT_R-PT_L) + rho_wt_L*vxL - rho_wt_R*vxR) / (rho_wt_L - rho_wt_R);
        P_M = PT_L + rho_wt_L*(S_M-vxL);
        
        /* ok, are the roe-averaged wavespeed guesses reasonable? */
        if((P_M<=0)||(isnan(P_M))||(P_M>press_tot_limiter))
        {
            double S_plus = DMAX(fabs(vxL),fabs(vxR)) + c_eff;
            S_L=-S_plus; S_R=S_plus;
            
            rho_wt_L = Riemann_vec.L.rho*(S_L-vxL);
            rho_wt_R = Riemann_vec.R.rho*(S_R-vxR);
            S_M = ((PT_R-PT_L) + rho_wt_L*vxL - rho_wt_R*vxR) / (rho_wt_L - rho_wt_R);
            P_M = PT_L + rho_wt_L*(S_M-vxL);
        }
    }
    if(P_M <= MIN_REAL_NUMBER) {P_M = MIN_REAL_NUMBER; within_sampled_vacuum = 1;}
    
    /* done trying different wavespeed estimates. we'll enter the flux computation if we have a valid pressure. 
        if not, this will return to the main hydro routine, and a lower-order reconstruction will be attempted */
    
    if(within_sampled_vacuum == 1)
    {
        
        /* ok, we are supersonically separating, the only safe thing to do is to assign no fluxes of conserved quantities */
        P_M = MIN_REAL_NUMBER; S_M = v_frame = Riemann_out->B_normal_corrected = Riemann_out->phi_normal_mean = Riemann_out->phi_normal_db = 0; // vanishing //
        memset(&Riemann_out->Fluxes, 0, sizeof(struct Conserved_var_Riemann)); // set all fluxes to vanish //
        memset(Interface_State, 0, sizeof(struct Conserved_var_Riemann)); // interface state is also vanishing (for subsequent fluxes) //
        
    } else {
        
#if defined(HYDRO_MESHLESS_FINITE_MASS)
        v_frame = S_M; /* in the lagrangian scheme, we must calculate fluxes consistent with the assumption
                        that there is zero mass flux. this gives that result for the HLLC/D flux */
#endif
        if((P_M > 0)&&(!isnan(P_M)))
        {
            /* alright, we have a valid solution! we can now compute the HLLD fluxes */
            if(v_frame <= S_M)
            {
                /* v_frame <= S_M : we're in the left fan */
                tmp = Riemann_vec.L.v[0]*Riemann_vec.L.v[0]+Riemann_vec.L.v[1]*Riemann_vec.L.v[1]+Riemann_vec.L.v[2]*Riemann_vec.L.v[2];
                e_L = Riemann_vec.L.u*Riemann_vec.L.rho + 0.5*B2_L + 0.5*Riemann_vec.L.rho*tmp;
                vdotB_L = Riemann_vec.L.v[0]*Riemann_vec.L.B[0]+Riemann_vec.L.v[1]*Riemann_vec.L.B[1]+Riemann_vec.L.v[2]*Riemann_vec.L.B[2];
                Riemann_out->Fluxes.rho = Riemann_vec.L.rho * vxL;
                Riemann_out->Fluxes.p = (e_L + PT_L) * vxL - vdotB_L * Bx;
                Riemann_out->Fluxes.v[0] = Riemann_out->Fluxes.rho * vxL - Bx2 + PT_L;
                Riemann_out->Fluxes.v[1] = Riemann_out->Fluxes.rho * Riemann_vec.L.v[1] - Riemann_vec.L.B[1]*Bx;
                Riemann_out->Fluxes.v[2] = Riemann_out->Fluxes.rho * Riemann_vec.L.v[2] - Riemann_vec.L.B[2]*Bx;
                Riemann_out->Fluxes.B[0] = 0;
                Riemann_out->Fluxes.B[1] = Riemann_vec.L.B[1] * vxL - Bx * Riemann_vec.L.v[1];
                Riemann_out->Fluxes.B[2] = Riemann_vec.L.B[2] * vxL - Bx * Riemann_vec.L.v[2];
                if(v_frame < S_L)
                {
                    /* left state */
                    if(v_frame != 0)
                    {
                        /* correct for frame velocity if it is non-zero */
                        Riemann_out->Fluxes.rho -= v_frame * Riemann_vec.L.rho;
                        Riemann_out->Fluxes.p -= v_frame * e_L;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] -= v_frame * Riemann_vec.L.rho * Riemann_vec.L.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] -= v_frame * Riemann_vec.L.B[k];
                        }
                        Interface_State = &Riemann_vec.L;
                    }
                } else {
                    U_s.rho = rho_wt_L / (S_L - S_M);
                    double sqrt_rho_star = sqrt(U_s.rho);
                    double S_star = S_M - fabs(Bx) / sqrt_rho_star;
                    tmp = rho_wt_L * (S_L - S_M) - Bx2;
                    U_s.v[0] = S_M;
                    U_s.B[0] = Bx;
                    if(fabs(tmp) < SMALL_NUMBER * P_M)
                    {
                        U_s.v[1] = Riemann_vec.L.v[1];
                        U_s.v[2] = Riemann_vec.L.v[2];
                        U_s.B[1] = Riemann_vec.L.B[1];
                        U_s.B[2] = Riemann_vec.L.B[2];
                    } else {
                        tmp = 1/tmp;
                        tmp2 = tmp * (S_M - vxL) * Bx;
                        U_s.v[1] = Riemann_vec.L.v[1] - Riemann_vec.L.B[1] * tmp2;
                        U_s.v[2] = Riemann_vec.L.v[2] - Riemann_vec.L.B[2] * tmp2;
                        tmp2 = tmp * (rho_wt_L*(S_L-vxL) - Bx2);
                        U_s.B[1] = Riemann_vec.L.B[1] * tmp2;
                        U_s.B[2] = Riemann_vec.L.B[2] * tmp2;
                    }
                    double vdotB_star = U_s.v[0]*U_s.B[0] + U_s.v[1]*U_s.B[1] + U_s.v[2]*U_s.B[2];
                    U_s.p = (e_L*(S_L-vxL) - PT_L*vxL + P_M*S_M + Bx*(vdotB_L-vdotB_star)) / (S_L-S_M);
                    /* F_L - S_L*U_L :: */
                    Riemann_out->Fluxes.rho -= S_L * Riemann_vec.L.rho;
                    Riemann_out->Fluxes.p -= S_L * e_L;
                    for(k=0;k<3;k++)
                    {
                        Riemann_out->Fluxes.v[k] -= S_L * Riemann_vec.L.rho * Riemann_vec.L.v[k];
                        if(k>0) Riemann_out->Fluxes.B[k] -= S_L * Riemann_vec.L.B[k];
                    }
                    
                    if((v_frame <= S_star) || (0.5*Bx2 < SMALL_NUMBER * P_M))
                    {
                        /* left Alfven wave :: F_L - S_L*U_L + (S_L-a_x)*U_star */
                        /* note also: if Bx=0, the star-star state is the same as the star-state, so can just finish here! */
                        double dS = S_L - v_frame;
                        Riemann_out->Fluxes.rho += dS * U_s.rho;
                        Riemann_out->Fluxes.p += dS * U_s.p;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] += dS * U_s.rho * U_s.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] += dS * U_s.B[k];
                        }
                        Interface_State = &U_s;
                        
                    } else {
                        /* ok the star and star-star states are different, need to compute them separately */
                        /* left contact wave :: F_L - S_L*U_L - (S_Lstar-S_L)*U_star + (S_Lstar-a_x)*U_starstar */
                        double dS = S_L - S_star;
                        Riemann_out->Fluxes.rho += dS * U_s.rho;
                        Riemann_out->Fluxes.p += dS * U_s.p;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] += dS * U_s.rho * U_s.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] += dS * U_s.B[k];
                        }
                        /* now we just need the star-star state to finish */
                        U_ss.rho = U_s.rho;
                        U_ss.p = U_s.p;
                        for(k=0;k<3;k++) {U_ss.v[k]=U_s.v[k]; U_ss.B[k]=U_s.B[k];}
                        U_ss.rho = U_s.rho;
                        U_ss.v[0] = S_M;
                        U_ss.B[0] = Bx;
                        double sign_Bx = 1.0;
                        if(Bx < 0) sign_Bx = -1.0;
                        
                        /* getting the middle states requires both sides: get the necessary R quantities */
                        V_s.rho = rho_wt_R / (S_R - S_M);
                        double sqrt_rho_star_alt = sqrt(V_s.rho);
                        double irhowt = 1/(sqrt_rho_star + sqrt_rho_star_alt);
                        tmp = rho_wt_R * (S_R - S_M) - Bx2;
                        if(fabs(tmp) < SMALL_NUMBER * P_M)
                        {
                            V_s.v[1] = Riemann_vec.R.v[1];
                            V_s.v[2] = Riemann_vec.R.v[2];
                            V_s.B[1] = Riemann_vec.R.B[1];
                            V_s.B[2] = Riemann_vec.R.B[2];
                        } else {
                            tmp = 1/tmp;
                            tmp2 = tmp * (S_M - vxR) * Bx;
                            V_s.v[1] = Riemann_vec.R.v[1] - Riemann_vec.R.B[1] * tmp2;
                            V_s.v[2] = Riemann_vec.R.v[2] - Riemann_vec.R.B[2] * tmp2;
                            tmp2 = tmp * (rho_wt_R*(S_R-vxR) - Bx2);
                            V_s.B[1] = Riemann_vec.R.B[1] * tmp2;
                            V_s.B[2] = Riemann_vec.R.B[2] * tmp2;
                        }
                        for(k=1;k<3;k++)
                        {
                            U_ss.v[k] = irhowt * (sqrt_rho_star*U_s.v[k] + sqrt_rho_star_alt*V_s.v[k] +
                                                  (V_s.B[k]-U_s.B[k])*sign_Bx);
                            U_ss.B[k] = irhowt * (sqrt_rho_star*V_s.B[k] + sqrt_rho_star_alt*U_s.B[k] +
                                                  sqrt_rho_star*sqrt_rho_star_alt*(V_s.v[k]-U_s.v[k])*sign_Bx);
                        }
                        double vdotB_ss = U_ss.v[0]*U_ss.B[0] + U_ss.v[1]*U_ss.B[1] + U_ss.v[2]*U_ss.B[2];
                        U_ss.p = U_s.p - sqrt_rho_star * (vdotB_star - vdotB_ss) * sign_Bx;
                        dS = S_star - v_frame;
                        Riemann_out->Fluxes.rho += dS * U_ss.rho;
                        Riemann_out->Fluxes.p += dS * U_ss.p;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] += dS * U_ss.rho * U_ss.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] += dS * U_ss.B[k];
                        }
                        Interface_State = &U_ss;
                    }
                }
                
            } else {
                /* v_frame > S_M : we're in the right fan */
                tmp = Riemann_vec.R.v[0]*Riemann_vec.R.v[0]+Riemann_vec.R.v[1]*Riemann_vec.R.v[1]+Riemann_vec.R.v[2]*Riemann_vec.R.v[2];
                e_R = Riemann_vec.R.u*Riemann_vec.R.rho + 0.5*B2_R + 0.5*Riemann_vec.R.rho*tmp;
                vdotB_R = Riemann_vec.R.v[0]*Riemann_vec.R.B[0]+Riemann_vec.R.v[1]*Riemann_vec.R.B[1]+Riemann_vec.R.v[2]*Riemann_vec.R.B[2];
                Riemann_out->Fluxes.rho = Riemann_vec.R.rho * vxR;
                Riemann_out->Fluxes.p = (e_R + PT_R) * vxR - vdotB_R * Bx;
                Riemann_out->Fluxes.v[0] = Riemann_out->Fluxes.rho * vxR - Bx2 + PT_R;
                Riemann_out->Fluxes.v[1] = Riemann_out->Fluxes.rho * Riemann_vec.R.v[1] - Riemann_vec.R.B[1]*Bx;
                Riemann_out->Fluxes.v[2] = Riemann_out->Fluxes.rho * Riemann_vec.R.v[2] - Riemann_vec.R.B[2]*Bx;
                Riemann_out->Fluxes.B[0] = 0;
                Riemann_out->Fluxes.B[1] = Riemann_vec.R.B[1] * vxR - Bx * Riemann_vec.R.v[1];
                Riemann_out->Fluxes.B[2] = Riemann_vec.R.B[2] * vxR - Bx * Riemann_vec.R.v[2];
                if(v_frame >= S_R)
                {
                    /* right state */
                    if(v_frame != 0)
                    {
                        /* correct for frame velocity if it is non-zero */
                        Riemann_out->Fluxes.rho -= v_frame * Riemann_vec.R.rho;
                        Riemann_out->Fluxes.p -= v_frame * e_R;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] -= v_frame * Riemann_vec.R.rho * Riemann_vec.R.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] -= v_frame * Riemann_vec.R.B[k];
                        }
                        Interface_State = &Riemann_vec.R;
                    }
                } else {
                    U_s.rho = rho_wt_R / (S_R - S_M);
                    double sqrt_rho_star = sqrt(U_s.rho);
                    double inv_sqrt_rho_star = 1 / sqrt_rho_star;
                    double S_star = S_M + fabs(Bx) * inv_sqrt_rho_star;
                    U_s.v[0] = S_M;
                    U_s.B[0] = Bx;
                    tmp = rho_wt_R * (S_R - S_M) - Bx2;
                    if(fabs(tmp) < SMALL_NUMBER * P_M)
                    {
                        U_s.v[1] = Riemann_vec.R.v[1];
                        U_s.v[2] = Riemann_vec.R.v[2];
                        U_s.B[1] = Riemann_vec.R.B[1];
                        U_s.B[2] = Riemann_vec.R.B[2];
                    } else {
                        tmp = 1/tmp;
                        tmp2 = tmp * (S_M - vxR) * Bx;
                        U_s.v[1] = Riemann_vec.R.v[1] - Riemann_vec.R.B[1] * tmp2;
                        U_s.v[2] = Riemann_vec.R.v[2] - Riemann_vec.R.B[2] * tmp2;
                        tmp2 = tmp * (rho_wt_R*(S_R-vxR) - Bx2);
                        U_s.B[1] = Riemann_vec.R.B[1] * tmp2;
                        U_s.B[2] = Riemann_vec.R.B[2] * tmp2;
                    }
                    double vdotB_star = U_s.v[0]*U_s.B[0] + U_s.v[1]*U_s.B[1] + U_s.v[2]*U_s.B[2];
                    U_s.p = (e_R*(S_R-vxR) - PT_R*vxR + P_M*S_M + Bx*(vdotB_R-vdotB_star)) / (S_R-S_M);
                    /* F_R - S_R*U_R :: */
                    Riemann_out->Fluxes.rho -= S_R * Riemann_vec.R.rho;
                    Riemann_out->Fluxes.p -= S_R * e_R;
                    for(k=0;k<3;k++)
                    {
                        Riemann_out->Fluxes.v[k] -= S_R * Riemann_vec.R.rho * Riemann_vec.R.v[k];
                        if(k>0) Riemann_out->Fluxes.B[k] -= S_R * Riemann_vec.R.B[k];
                    }
                    
                    if((v_frame >= S_star) || (0.5*Bx2 < SMALL_NUMBER * P_M))
                    {
                        /* right Alfven wave :: F_R - S_R*U_R + (S_R-a_x)*U_star */
                        /* note also: if Bx=0, the star-star state is the same as the star-state, so can just finish here! */
                        double dS = S_R - v_frame;
                        Riemann_out->Fluxes.rho += dS * U_s.rho;
                        Riemann_out->Fluxes.p += dS * U_s.p;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] += dS * U_s.rho * U_s.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] += dS * U_s.B[k];
                        }
                        Interface_State = &U_s;
                        
                    } else {
                        /* ok the star and star-star states are different, need to compute them separately */
                        /* right contact wave :: F_R - S_R*U_R - (S_Rstar-S_R)*U_star + (S_Rstar-a_x)*U_starstar */
                        double dS = S_R - S_star;
                        Riemann_out->Fluxes.rho += dS * U_s.rho;
                        Riemann_out->Fluxes.p += dS * U_s.p;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] += dS * U_s.rho * U_s.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] += dS * U_s.B[k];
                        }
                        /* now we just need the star-star state to finish */
                        U_ss.rho = U_s.rho;
                        U_ss.v[0] = S_M;
                        U_ss.B[0] = Bx;
                        double sign_Bx = 1.0;
                        if(Bx < 0) sign_Bx = -1.0;
                        
                        /* getting the middle states requires both sides: get the necessary R quantities */
                        V_s.rho = rho_wt_L / (S_L - S_M);
                        double sqrt_rho_star_alt = sqrt(V_s.rho);
                        double irhowt = 1/(sqrt_rho_star + sqrt_rho_star_alt);
                        tmp = rho_wt_L * (S_L - S_M) - Bx2;
                        if(fabs(tmp) < SMALL_NUMBER * P_M)
                        {
                            V_s.v[1] = Riemann_vec.L.v[1];
                            V_s.v[2] = Riemann_vec.L.v[2];
                            V_s.B[1] = Riemann_vec.L.B[1];
                            V_s.B[2] = Riemann_vec.L.B[2];
                        } else {
                            tmp = 1/tmp;
                            tmp2 = tmp * (S_M - vxL) * Bx;
                            V_s.v[1] = Riemann_vec.L.v[1] - Riemann_vec.L.B[1] * tmp2;
                            V_s.v[2] = Riemann_vec.L.v[2] - Riemann_vec.L.B[2] * tmp2;
                            tmp2 = tmp * (rho_wt_L*(S_L-vxL) - Bx2);
                            V_s.B[1] = Riemann_vec.L.B[1] * tmp2;
                            V_s.B[2] = Riemann_vec.L.B[2] * tmp2;
                        }
                        for(k=1;k<3;k++)
                        {
                            U_ss.v[k] = irhowt * (sqrt_rho_star*U_s.v[k] + sqrt_rho_star_alt*V_s.v[k] +
                                                  (U_s.B[k]-V_s.B[k])*sign_Bx);
                            U_ss.B[k] = irhowt * (sqrt_rho_star*V_s.B[k] + sqrt_rho_star_alt*U_s.B[k] +
                                                  sqrt_rho_star*sqrt_rho_star_alt*(U_s.v[k]-V_s.v[k])*sign_Bx);
                        }
                        double vdotB_ss = U_ss.v[0]*U_ss.B[0] + U_ss.v[1]*U_ss.B[1] + U_ss.v[2]*U_ss.B[2];
                        U_ss.p = U_s.p + sqrt_rho_star * (vdotB_star - vdotB_ss) * sign_Bx;
                        dS = S_star - v_frame;
                        Riemann_out->Fluxes.rho += dS * U_ss.rho;
                        Riemann_out->Fluxes.p += dS * U_ss.p;
                        for(k=0;k<3;k++)
                        {
                            Riemann_out->Fluxes.v[k] += dS * U_ss.rho * U_ss.v[k];
                            if(k>0) Riemann_out->Fluxes.B[k] += dS * U_ss.B[k];
                        }
                        Interface_State = &U_ss;
                    }
                }
            }
        } // if((P_M > 0)&&(!isnan(P_M))) //
    } // within_sampled_vacuum check //

    /* alright, we've gotten successful HLLD fluxes! */
    Riemann_out->Fluxes.B[0] = -v_frame * Bx;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME)
    //Riemann_out->Fluxes.phi = -Interface_State->v[0] * Interface_State->rho * Riemann_out->phi_normal_mean; // potentially improved phi-flux for MFV with mass-based fluxes
    //Riemann_out->Fluxes.phi = -v_frame * Riemann_out->phi_normal_mean; // need to use the proper phi from the updated problem // mass-based phi-fluxes don't require this 
    //Riemann_out->Fluxes.phi -= All.DivBcleanHyperbolicSigma * c_eff*c_eff * Bx;
#endif
    Riemann_out->S_M=v_frame;
    Riemann_out->P_M=P_M;
#ifdef SAVE_FACE_DENSITY
    Riemann_out->Face_Density = Interface_State->rho;
#endif
#ifdef SAVE_FACE_BFIELD
    for(k=0;k<3;k++) {Riemann_out->Face_B[k] = Interface_State->B[k];}
#endif
#ifdef SAVE_FACE_VFIELD
    for(k=0;k<3;k++) {Riemann_out->Face_Vel[k] = Interface_State->v[k];}
#endif
    return;

    /* originally we included HLL and LF solvers here, but those are way too prone to give unphysical
        results, and lead to crashing (plus they are more diffusive than a low-order 
        reconstruction); therefore if HLLD fails, we prefer to simply re-calculate with a lower-order 
        reconstruction at the face */
    
} /* yay! we're done writing our HLLD solver! */







/* --------------------------------------------------------------------------------------------------------------
 *  obtain the rotation matrix to rotate into the frame, and do the rotation (as well as routine to rotate back)
 *    (based on the code by V. Springel in AREPO; a nearly identical implementation was also developed by E. Gaburov
 *      based on the Weighted-Particle MHD code)
 * -------------------------------------------------------------------------------------------------------------- */
void rotate_states_to_face(struct Input_vec_Riemann *Riemann_vec, double n_unit[3], struct rotation_matrix *rot_matrix)
{
    rot_matrix->n[0] = n_unit[0];
    rot_matrix->n[1] = n_unit[1];
    rot_matrix->n[2] = n_unit[2];
    /* now we can construct a basis orthonormal to this */
    if((rot_matrix->n[0]==0)&&(rot_matrix->n[1]==0))
    {
        /* trap for the pathological case */
        rot_matrix->m[0] = 1;
        rot_matrix->m[1] = rot_matrix->m[2] = 0;
        rot_matrix->p[1] = 1;
        rot_matrix->p[0] = rot_matrix->p[2] = 0;
    } else {
        rot_matrix->m[0] = -rot_matrix->n[1];
        rot_matrix->m[1] = rot_matrix->n[0];
        rot_matrix->m[2] = 0;
        double mm = sqrt(rot_matrix->m[0]*rot_matrix->m[0] + rot_matrix->m[1]*rot_matrix->m[1]);
        rot_matrix->m[0] /= mm;
        rot_matrix->m[1] /= mm;
        
        rot_matrix->p[0] = rot_matrix->n[1] * rot_matrix->m[2] - rot_matrix->n[2] * rot_matrix->m[1];
        rot_matrix->p[1] = rot_matrix->n[2] * rot_matrix->m[0] - rot_matrix->n[0] * rot_matrix->m[2];
        rot_matrix->p[2] = rot_matrix->n[0] * rot_matrix->m[1] - rot_matrix->n[1] * rot_matrix->m[0];
    }
    
    /* now we have an orthonormal rotation matrix -- we can rotate the states */
    int k; double v[3];
    for(k=0;k<3;k++) {v[k]=Riemann_vec->L.v[k];}
    Riemann_vec->L.v[0] = v[0]*rot_matrix->n[0] + v[1]*rot_matrix->n[1] + v[2]*rot_matrix->n[2];
    Riemann_vec->L.v[1] = v[0]*rot_matrix->m[0] + v[1]*rot_matrix->m[1] + v[2]*rot_matrix->m[2];
    Riemann_vec->L.v[2] = v[0]*rot_matrix->p[0] + v[1]*rot_matrix->p[1] + v[2]*rot_matrix->p[2];
    for(k=0;k<3;k++) {v[k]=Riemann_vec->R.v[k];}
    Riemann_vec->R.v[0] = v[0]*rot_matrix->n[0] + v[1]*rot_matrix->n[1] + v[2]*rot_matrix->n[2];
    Riemann_vec->R.v[1] = v[0]*rot_matrix->m[0] + v[1]*rot_matrix->m[1] + v[2]*rot_matrix->m[2];
    Riemann_vec->R.v[2] = v[0]*rot_matrix->p[0] + v[1]*rot_matrix->p[1] + v[2]*rot_matrix->p[2];
#ifdef MAGNETIC
    for(k=0;k<3;k++) {v[k]=Riemann_vec->L.B[k];}
    Riemann_vec->L.B[0] = v[0]*rot_matrix->n[0] + v[1]*rot_matrix->n[1] + v[2]*rot_matrix->n[2];
    Riemann_vec->L.B[1] = v[0]*rot_matrix->m[0] + v[1]*rot_matrix->m[1] + v[2]*rot_matrix->m[2];
    Riemann_vec->L.B[2] = v[0]*rot_matrix->p[0] + v[1]*rot_matrix->p[1] + v[2]*rot_matrix->p[2];
    for(k=0;k<3;k++) {v[k]=Riemann_vec->R.B[k];}
    Riemann_vec->R.B[0] = v[0]*rot_matrix->n[0] + v[1]*rot_matrix->n[1] + v[2]*rot_matrix->n[2];
    Riemann_vec->R.B[1] = v[0]*rot_matrix->m[0] + v[1]*rot_matrix->m[1] + v[2]*rot_matrix->m[2];
    Riemann_vec->R.B[2] = v[0]*rot_matrix->p[0] + v[1]*rot_matrix->p[1] + v[2]*rot_matrix->p[2];
#endif
}
void rotate_fluxes_back_to_lab(struct Riemann_outputs *Riemann_out, struct rotation_matrix rot_matrix)
{
    /* for an orthonormal rotation matrix A, we have A_transpose = A_inverse, so this is easy */
    int k; double v[3];
    for(k=0;k<3;k++) {v[k]=Riemann_out->Fluxes.v[k];}
    Riemann_out->Fluxes.v[0] = v[0]*rot_matrix.n[0] + v[1]*rot_matrix.m[0] + v[2]*rot_matrix.p[0];
    Riemann_out->Fluxes.v[1] = v[0]*rot_matrix.n[1] + v[1]*rot_matrix.m[1] + v[2]*rot_matrix.p[1];
    Riemann_out->Fluxes.v[2] = v[0]*rot_matrix.n[2] + v[1]*rot_matrix.m[2] + v[2]*rot_matrix.p[2];
#ifdef MAGNETIC
    for(k=0;k<3;k++) {v[k]=Riemann_out->Fluxes.B[k];}
    Riemann_out->Fluxes.B[0] = v[0]*rot_matrix.n[0] + v[1]*rot_matrix.m[0] + v[2]*rot_matrix.p[0];
    Riemann_out->Fluxes.B[1] = v[0]*rot_matrix.n[1] + v[1]*rot_matrix.m[1] + v[2]*rot_matrix.p[1];
    Riemann_out->Fluxes.B[2] = v[0]*rot_matrix.n[2] + v[1]*rot_matrix.m[2] + v[2]*rot_matrix.p[2];
#endif
#ifdef SAVE_FACE_BFIELD
    for(k=0;k<3;k++) {v[k]=Riemann_out->Face_B[k];}
    Riemann_out->Face_B[0] = v[0]*rot_matrix.n[0] + v[1]*rot_matrix.m[0] + v[2]*rot_matrix.p[0];
    Riemann_out->Face_B[1] = v[0]*rot_matrix.n[1] + v[1]*rot_matrix.m[1] + v[2]*rot_matrix.p[1];
    Riemann_out->Face_B[2] = v[0]*rot_matrix.n[2] + v[1]*rot_matrix.m[2] + v[2]*rot_matrix.p[2];
#endif
#ifdef SAVE_FACE_VFIELD
    for(k=0;k<3;k++) {v[k]=Riemann_out->Face_Vel[k];}
    Riemann_out->Face_Vel[0] = v[0]*rot_matrix.n[0] + v[1]*rot_matrix.m[0] + v[2]*rot_matrix.p[0];
    Riemann_out->Face_Vel[1] = v[0]*rot_matrix.n[1] + v[1]*rot_matrix.m[1] + v[2]*rot_matrix.p[1];
    Riemann_out->Face_Vel[2] = v[0]*rot_matrix.n[2] + v[1]*rot_matrix.m[2] + v[2]*rot_matrix.p[2];
#endif
}


#endif // MAGNETIC //

