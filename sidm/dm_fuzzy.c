#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file dm_fuzzy.c
 *  \brief routines needed for fuzzy-DM implementation
 *         This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef DM_FUZZY



/* --------------------------------------------------------------------------
 Actual evaluation of fluxes from the quantum pressure tensor
 -------------------------------------------------------------------------- */



void do_dm_fuzzy_flux_computation(double HLLwt, double dt, double prev_a, double dv[3],
                                  double GradRho_L[3], double GradRho_R[3],
                                  double GradRho2_L[3][3], double GradRho2_R[3][3],
                                  double rho_L, double rho_R, double dv_Right_minus_Left,
                                  double Area[3], double fluxes[3], double AGS_Numerical_QuantumPotential_L, double AGS_Numerical_QuantumPotential_R, double *dt_egy_Numerical_QuantumPotential)
{
    double f00 = 0.5 * All.ScalarField_hbar_over_mass; f00*=f00; // this encodes the coefficient with the mass of the particle: units vel*L = hbar / particle_mass
    double gamma_eff=5./3., PL_dot_A[3]={0}, PR_dot_A[3]={0}, PL_dot_AA=0, PR_dot_AA=0, Face_Area_Norm=0, rSi=1./(rho_L+rho_R), QL=0, QR=0; int m,k;
    for(m=0;m<3;m++) {QL+=GradRho2_L[m][m] - 0.5*GradRho_L[m]*GradRho_L[m]/rho_L; QR+=GradRho2_R[m][m] - 0.5*GradRho_R[m]*GradRho_R[m]/rho_R;} // compute the quantum potential (multiplied by rho to be an energy density to match units of PN)
    double PN_L=(gamma_eff-1.)*(AGS_Numerical_QuantumPotential_L-f00*QL), PN_R=(gamma_eff-1.)*(AGS_Numerical_QuantumPotential_R-f00*QR); PN_L=DMAX(PN_L,0); PN_R=DMAX(PN_R,0); // compute actual scalar-pressure terms and limit to be positive definite
    for(m=0;m<3;m++)
    {
        for(k=0;k<3;k++)
        {
            PL_dot_A[m] += Area[k] * f00*(GradRho_L[m]*GradRho_L[k]/rho_L - GradRho2_L[m][k]); // convert grad^2_rho ~ rho/L^2 from code units to physical [should already all be physical here]
            PR_dot_A[m] += Area[k] * f00*(GradRho_R[m]*GradRho_R[k]/rho_R - GradRho2_R[m][k]); // convert grad^2_rho ~ rho/L^2 from code units to physical [should already all be physical here]
        }
        Face_Area_Norm += Area[m]*Area[m]; PL_dot_AA += PL_dot_A[m]*Area[m]; PR_dot_AA += PR_dot_A[m]*Area[m]; // compute area and normal scalar component of pressure tensor
        fluxes[m] = (rho_L*(PR_dot_A[m]+PN_R*Area[m]) + rho_R*(PL_dot_A[m]+PN_L*Area[m])) * rSi; // term needed for Pstar-dot-A
    }
    PL_dot_AA /= Face_Area_Norm; PR_dot_AA /= Face_Area_Norm; Face_Area_Norm=sqrt(Face_Area_Norm); // define projected normal pressure components
    double PL_norm = PL_dot_AA + PN_L, PR_norm = PR_dot_AA + PN_L, cs_L = gamma_eff*PL_norm/rho_L, cs_R = gamma_eff*PR_norm/rho_R, cs=DMAX(cs_L,cs_R), ceff=cs; // get projected components, plus pure pressure, to sound speed as defined for adiabatic perturb here
    if(dv_Right_minus_Left < 0) {ceff+=fabs(dv_Right_minus_Left);} // this is the necessary upwind step for approaching cells
    for(m=0;m<3;m++) {fluxes[m]-=rho_L*rho_R*rSi*ceff*dv_Right_minus_Left*Area[m];} // numerical HLLC flux in star frame from upwinding appropriately
    double S_M = ((PL_norm-PR_norm)/ceff + 0.5*(rho_R-rho_L)*dv_Right_minus_Left)*rSi; // contact wave speed (oriented in Area[m] direction)
    *dt_egy_Numerical_QuantumPotential=0; for(m=0;m<3;m++) {*dt_egy_Numerical_QuantumPotential += fluxes[m]*(S_M*Area[m]/Face_Area_Norm - 0.5*dv[m]);} // numerical flux for energy
    return;
}



void do_dm_fuzzy_flux_computation_old(double HLLwt, double dt, double m0, double prev_a, double dp[3], double dv[3],
                                  double GradRho_L[3], double GradRho_R[3],
                                  double GradRho2_L[3][3], double GradRho2_R[3][3],
                                  double rho_L, double rho_R, double dv_Right_minus_Left,
                                  double Area[3], double fluxes[3], double AGS_Numerical_QuantumPotential, double *dt_egy_Numerical_QuantumPotential)
{
    if(dt <= 0) return; // no timestep, no flux
    int m,n;
    double f00 = 0.5 * All.ScalarField_hbar_over_mass; // this encodes the coefficient with the mass of the particle: units vel*L = hbar / particle_mass
    // (0.5=1/2, pre-factor from dimensionless equations; 591569 = hbar/eV in cm^2/s; add mass in eV, and put in code units
    double f2 = f00*f00, rhoL_i=1./rho_L, rhoR_i=1./rho_R, r2=0, rSi=1./(rho_L+rho_R); *dt_egy_Numerical_QuantumPotential=0;
    for(m=0;m<3;m++) {r2+=dp[m]*dp[m]; fluxes[m]=0;} /* zero fluxes and calculate separation */
    if(r2 <= 0) return; // same element
    double r=sqrt(r2), wavespeed=2.*f00*(M_PI/r); // approximate k = 2pi/lambda = 2pi/(2*dr) as the maximum k the code will allow locally */
    /* note that the QPT admits waves parallel to k, with wavespeed omega = pm 2*f00*k, so include these for HLLC solution */
    if(dv_Right_minus_Left > wavespeed) return; // elements are receding super-sonically, no way to communicate pressure //

    double fluxmag=0, Face_Area_Norm=0; for(m=0;m<3;m++) {Face_Area_Norm+=Area[m]*Area[m];}
    Face_Area_Norm=sqrt(Face_Area_Norm);

    for(m=0;m<3;m++)
    {
        for(n=0;n<3;n++)
        {
            double QPT_L = f2*(rhoL_i*GradRho_L[m]*GradRho_L[n] - GradRho2_L[m][n]); // convert grad^2_rho ~ rho/L^2 from code units to physical [should already all be physical here]
            double QPT_R = f2*(rhoR_i*GradRho_R[m]*GradRho_R[n] - GradRho2_R[m][n]); // convert grad^2_rho ~ rho/L^2 from code units to physical [should already all be physical here]
            /* calculate 'star' solution (interface moving with contact wave, since we have a Lagrangian code)
             for HLLC reimann problem based on these pressure tensors */
            double P_star = (QPT_L*rho_R + QPT_R*rho_L) * rSi; // if there were no waves and all at rest
            fluxes[m] += Area[n] * P_star; /* momentum flux into direction 'm' given by Area.Pressure */
            // sign convention here: -Area[m] * (positive definite) => repulsive force, +Area[m] => attractive
        }
        fluxmag += fluxes[m]*fluxes[m];
    }
    fluxmag = sqrt(fluxmag);
    double fluxmax = 100. * Face_Area_Norm * f2 * 0.5*(rho_L+rho_R) / (r*r); // limiter to prevent crazy values where derivatives are ill-posed (e.g. discontinuities)
    if(fluxmag > fluxmax) {for(m=0;m<3;m++) {fluxes[m] *= fluxmax/fluxmag;}}

    for(m=0;m<3;m++)
    {
        double ftmp = (2./3.)*AGS_Numerical_QuantumPotential*Area[m]; // 2/3 b/c the equation-of-state of the 'quantum pressure tensor' is gamma=5/3 under isotropic compression/expansion //
        double fmax = 0.5 * m0 * fabs(dv[m]) / dt; fmax = DMAX(fmax, 100.*fluxmag); fmax = DMAX(DMIN(fmax , 100.*prev_a), fluxmag); if(fabs(ftmp) > fmax) {ftmp *= fmax/fabs(ftmp);} // limit pressure-induced acceleration to prevent unphysical cases
        *dt_egy_Numerical_QuantumPotential -= 0.5*ftmp*dv[m]; // PdV work from this pressure term //
        fluxes[m] += ftmp; // add numerical 'pressure' stored from previous timesteps //
    }
    fluxmag=0; for(m=0;m<3;m++) {fluxmag += fluxes[m]*fluxes[m];} if(fluxmag > 0) {fluxmag = sqrt(fluxmag);} else {fluxmag = 0;}

    /* now we have to introduce the numerical diffusivity (the up-wind mixing part from the Reimann problem);
     this can have one of a couple forms, but the most accurate and stable appears to be the traditional HLLC form which we use by default below */
    if(dv_Right_minus_Left < 0) // converging flow, upwind dissipation terms appear //
    {
        // estimate wavenumber needed to calculate wavespeed below, for dissipation terms //
        double g_kL = sqrt(GradRho_L[0]*GradRho_L[0]+GradRho_L[1]*GradRho_L[1]+GradRho_L[2]*GradRho_L[2]); // gradient magnitude
        double g_kR = sqrt(GradRho_R[0]*GradRho_R[0]+GradRho_R[1]*GradRho_R[1]+GradRho_R[2]*GradRho_R[2]); // gradient magnitude
        double k_g1 = 0.5*(rhoL_i*g_kL + rhoR_i*g_kR); // gradient scale-length as proxy for wavenumber [crude] //
        double k_eff=0, k_g2=0, k_g3=0, k_d1 = fabs(rho_L-rho_R)/(0.5*(rho_L+rho_R)) / r; // even cruder delta-based k-estimate
        if(isnan(k_d1)) {k_d1=0;} // trap for division by zero above (should give zero wavespeed)
        if(isnan(k_g1)) {k_g1=0;} // trap
        double k2L = GradRho2_L[0][0]+GradRho2_L[1][1]+GradRho2_L[2][2]; // laplacian
        double k2R = GradRho2_R[0][0]+GradRho2_R[1][1]+GradRho2_R[2][2]; // laplacian
        k_g2 = sqrt( (fabs(k2L)+fabs(k2R)) * rSi ); // sqrt of second derivative -- again crude, doesn't necessarily recover full k but lower value
        k_g3 = 0.5 * sqrt(fabs(k2L-k2R) / (0.5*(g_kL+g_kR)*r)); // third-derivative to first-derivative ratio: exact for resolved wave
        double k_gtan = (fabs(k2L)+fabs(k2R)) / (MIN_REAL_NUMBER + fabs(g_kL) + fabs(g_kR));
        if(isnan(k_g2)) {k_g2=0;} // trap
        if(isnan(k_g3)) {k_g3=0;} // trap
        if(isnan(k_gtan)) {k_gtan=0;} // trap
        k_eff = DMIN(DMAX(k_gtan , DMAX(DMAX((3.+HLLwt)*DMAX(k_d1,k_g1), (2.+HLLwt)*k_g2)  , (1.5+HLLwt)*k_g3)) , 1./r);
        if(isnan(k_eff)) {k_eff=0;} // trap
        double Pstar = (-dv_Right_minus_Left)*(f00*k_eff + (-dv_Right_minus_Left))*rho_L*rho_R*rSi; // HLLC diffusive term
        fluxmag = 0; for(m=0;m<3;m++) {fluxmag += fluxes[m]*fluxes[m];}
        fluxmag = sqrt(fluxmag);
        for(m=0;m<3;m++)
        {
            double f_dir = Area[m]*Pstar, fmax = 0.5 * m0 * fabs(dv[m]) / dt; // assume the face points along the line between particles (very similar, but slightly more stable/diffusive if faces are highly-irregular)
            fmax = DMAX(fmax, 10.*fluxmag); fmax = DMAX(DMIN(fmax , 40.*prev_a), fluxmag); // limit diffusive flux to multiplier of physical flux
            if(fabs(f_dir) > fmax) {f_dir *= fmax/fabs(f_dir);} // limit diffusive flux to avoid overshoot (numerical stability of the diffusion terms) //
            fluxes[m] += f_dir; /* momentum flux into direction 'm' given by Area.Pressure */
            *dt_egy_Numerical_QuantumPotential -= 0.5 * f_dir * dv[m];
        }
    } // approach velocities lead to up-wind mixing
    return;
}




/* kicks for fuzzy-dm integration: just put relevant drift-kick operators here to keep the code clean
 mode=0 -> 'kick', mode=1 -> 'drift' */
void do_dm_fuzzy_drift_kick(int i, double dt, int mode)
{
    if(mode==0)
    {
        // calculate various energies: quantum potential QP0, 'stored' numerical pressure NQ0, kinetic energy KE0
        double dNQ=P[i].AGS_Dt_Numerical_QuantumPotential*dt, NQ0=P[i].AGS_Numerical_QuantumPotential, NQ1=NQ0+dNQ, KE0=0.5*P[i].Mass*(P[i].Vel[0]*P[i].Vel[0]+P[i].Vel[1]*P[i].Vel[1]+P[i].Vel[2]*P[i].Vel[2])*All.cf_a2inv;
        double f00 = 0.5 * All.ScalarField_hbar_over_mass; // this encodes the coefficient with the mass of the particle: units vel*L = hbar / particle_mass
        double d2rho = P[i].AGS_Gradients2_Density[0][0] + P[i].AGS_Gradients2_Density[1][1] + P[i].AGS_Gradients2_Density[2][2]; // laplacian
        double drho2 = P[i].AGS_Gradients_Density[0]*P[i].AGS_Gradients_Density[0] + P[i].AGS_Gradients_Density[1]*P[i].AGS_Gradients_Density[1] + P[i].AGS_Gradients_Density[2]*P[i].AGS_Gradients_Density[2];
        double QP0 = (f00*f00 / P[i].AGS_Density) * (d2rho - 0.5*drho2/P[i].AGS_Density); // quantum 'potential'
        NQ1 = DMAX(0,DMAX(NQ1,0.1*NQ0)); NQ1 = DMIN(NQ1,1.1*DMAX(DMAX(KE0+NQ0,fabs(QP0)),KE0+NQ0+QP0)); // limit kick to not produce unphysical energy over-or-under-shoot
        P[i].AGS_Numerical_QuantumPotential = NQ1;
    }

#if (DM_FUZZY > 0) /* if using direct-wavefunction integration methods */
    double vol_inv = P[i].AGS_Density / P[i].Mass;
    if(mode == 0)
    {
        //double psimag_mass_old = (P[i].AGS_Psi_Re*P[i].AGS_Psi_Re + P[i].AGS_Psi_Im*P[i].AGS_Psi_Im) * vol_inv;
        P[i].AGS_Psi_Re += P[i].AGS_Dt_Psi_Re * dt;
        P[i].AGS_Psi_Im += P[i].AGS_Dt_Psi_Im * dt;
        double mass_old = P[i].Mass, dmass = P[i].AGS_Dt_Psi_Mass * dt, mass_new = mass_old + dmass;
        dmass = DMIN(DMAX(dmass,-0.5*mass_old),0.5*mass_old);
        mass_new = mass_old + dmass;
        double psimag_mass_new = (P[i].AGS_Psi_Re*P[i].AGS_Psi_Re + P[i].AGS_Psi_Im*P[i].AGS_Psi_Im) * vol_inv;
#if (DM_FUZZY == 2)
        mass_new = psimag_mass_new; /* uses direct [NON-MASS-CONSERVING] integration of psi field */
#endif
        double psi_corr_fac = sqrt(mass_new / (MIN_REAL_NUMBER + psimag_mass_new));
        P[i].Mass = mass_new; P[i].AGS_Psi_Re *= psi_corr_fac; P[i].AGS_Psi_Im *= psi_corr_fac;

        P[i].AGS_Density = P[i].Mass * vol_inv;
        P[i].AGS_Psi_Re_Pred = P[i].AGS_Psi_Re;
        P[i].AGS_Psi_Im_Pred = P[i].AGS_Psi_Im;
    } else {
        /* in drift mode, AGS_Density should automatically be drifted already by the predictor step, but not the other quantities here */
        P[i].AGS_Psi_Re_Pred += P[i].AGS_Dt_Psi_Re * dt;
        P[i].AGS_Psi_Im_Pred += P[i].AGS_Dt_Psi_Im * dt;
        P[i].AGS_Density *= 1. + DMIN(DMAX(P[i].AGS_Dt_Psi_Mass*dt/P[i].Mass,-0.5),0.5);
    }
#endif
}


/* initialize wavefunction values in the code ICs */
void do_dm_fuzzy_initialization(void)
{
#if (DM_FUZZY > 0)
    int i;
    for(i = 0; i < NumPart; i++)
    {
        double volume = P[i].AGS_Density / P[i].Mass, psimag = sqrt(P[i].AGS_Density), phase = 0; int k=0;
        /* approximation for initial phase below is fine for slowly-varying k, otherwise not ideal */
        for(k=0;k<3;k++) {phase += P[i].Pos[k] * P[i].Vel[k] / All.ScalarField_hbar_over_mass;}

        P[i].AGS_Psi_Re = psimag * volume * cos(phase); /* remember, we evolve the volume-integrated value of psi */
        P[i].AGS_Psi_Im = psimag * volume * sin(phase);

        P[i].AGS_Dt_Psi_Mass = 0; P[i].AGS_Dt_Psi_Re = 0; P[i].AGS_Dt_Psi_Im = 0;
        P[i].AGS_Psi_Re_Pred = P[i].AGS_Psi_Re; P[i].AGS_Psi_Im_Pred = P[i].AGS_Psi_Im;
    }
#endif
}



void dm_fuzzy_reconstruct_and_slopelimit(double *u_R, double du_R[3], double *u_L, double du_L[3],
                                         double q_R, double dq_R[3], double d2q_R[3][3],
                                         double q_L, double dq_L[3], double d2q_L[3][3],
                                         double dx[3])
{
    double t_L,t_R; int k;
    dm_fuzzy_reconstruct_and_slopelimit_sub(&t_R,&t_L,q_R,dq_R,q_L,dq_L,dx);
    *u_R=t_R; *u_L=t_L;
    for(k=0;k<3;k++)
    {
        dm_fuzzy_reconstruct_and_slopelimit_sub(&t_R,&t_L,dq_R[k],d2q_R[k],dq_L[k],d2q_L[k],dx);
        dq_R[k]=t_R; dq_L[k]=t_L;
    }
    return;
}


void dm_fuzzy_reconstruct_and_slopelimit_sub(double *u_R_f, double *u_L_f, double q_R, double dq_R_0[3], double q_L, double dq_L_0[3], double dx[3])
{
    int k;
    double dq_L=0; for(k=0;k<3;k++) {dq_L += 0.5*dx[k]*dq_L_0[k];}
    double dq_R=0; for(k=0;k<3;k++) {dq_R -= 0.5*dx[k]*dq_R_0[k];}
    double q0=q_L, u_L=0, u_R=0; q_L-=q0; q_R-=q0;
    //double qmid = 0.5*q_R;

    if(dq_L*q_R<0) {dq_L=0;}
    if(dq_R*q_R<0) {dq_R=0;}
    u_L = dq_L; u_R = q_R + dq_R;
    if(q_R > 0)
    {
        if(u_L > u_R) {double tmp=0.5*(u_L+u_R); u_L=tmp; u_R=tmp;}
        if(u_L > q_R) {u_L=q_R;}
        if(u_R > q_R) {u_R=q_R;}
        if(u_L < 0) {u_L=0;}
        if(u_R < 0) {u_R=0;}
    } else {
        if(u_L < u_R) {double tmp=0.5*(u_L+u_R); u_L=tmp; u_R=tmp;}
        if(u_L > 0) {u_L=0;}
        if(u_R > 0) {u_R=0;}
        if(u_L < q_R) {u_L=q_R;}
        if(u_R < q_R) {u_R=q_R;}
    }
    if(q_R==0) {u_L=u_R=0;}
    u_L += q0; u_R += q0;
    *u_L_f = u_L; *u_R_f = u_R;
    return;
}




/* --------------------------------------------------------------------------
 Everything below here is a giant block to define the sub-routines needed
 to calculate the higher-order matrix gradient estimators for the density
 field, around each DM element (based on its interacting neighbor set,
 within the AGS_Hsml volume). This will give the density gradients
 AGS_Gradients_Density needed to actually compute the quantum pressure tensor
 -------------------------------------------------------------------------- */


/* define a common 'gradients' structure to hold everything we're going to take derivatives of */
struct Quantities_for_Gradients_DM
{
    MyDouble AGS_Density, AGS_Gradients_Density[3];
#if (DM_FUZZY > 0)
    MyDouble AGS_Psi_Re, AGS_Gradients_Psi_Re[3], AGS_Psi_Im, AGS_Gradients_Psi_Im[3];
#endif
};

/* this is a temporary structure for quantities used ONLY in the loop below, for example for computing the slope-limiters (for the Reimann problem) */
static struct temporary_dmgradients_data_topass
{
    struct Quantities_for_Gradients_DM Maxima;
    struct Quantities_for_Gradients_DM Minima;
}
*DMGradDataPasser;

struct kernel_DMGrad {double dp[3],r,wk_i, wk_j, dwk_i, dwk_j,h_i;};



#define CORE_FUNCTION_NAME DMGrad_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_DMGrad    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_DMGrad  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(ags_density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */


/* this structure defines the variables that need to be sent -from- the 'searching' element */
struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3], AGS_Hsml;
    struct Quantities_for_Gradients_DM GQuant;
    int NodeList[NODELISTLENGTH], Type;
}
*DATAIN_NAME, *DATAGET_NAME;

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
static inline void particle2in_DMGrad(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    in->AGS_Hsml = PPP[i].AGS_Hsml;
    in->Type = P[i].Type;
    in->GQuant.AGS_Density = P[i].AGS_Density;
    for(k=0;k<3;k++) {in->GQuant.AGS_Gradients_Density[k] = P[i].AGS_Gradients_Density[k];}
#if (DM_FUZZY > 0)
    in->GQuant.AGS_Psi_Re = P[i].AGS_Psi_Re_Pred * P[i].AGS_Density / P[i].Mass;
    for(k=0;k<3;k++) {in->GQuant.AGS_Gradients_Psi_Re[k] = P[i].AGS_Gradients_Psi_Re[k];}
    in->GQuant.AGS_Psi_Im = P[i].AGS_Psi_Im_Pred * P[i].AGS_Density / P[i].Mass;
    for(k=0;k<3;k++) {in->GQuant.AGS_Gradients_Psi_Im[k] = P[i].AGS_Gradients_Psi_Im[k];}
#endif
}


/* this structure defines the variables that need to be sent -back to- the 'searching' element */
struct OUTPUT_STRUCT_NAME
{
    struct Quantities_for_Gradients_DM Gradients[3];
    struct Quantities_for_Gradients_DM Maxima;
    struct Quantities_for_Gradients_DM Minima;
}
*DATARESULT_NAME, *DATAOUT_NAME;

#define ASSIGN_ADD_PRESET(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))

/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
static inline void out2particle_DMGrad(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    if(loop_iteration <= 0)
    {
        int k;
        MAX_ADD(DMGradDataPasser[i].Maxima.AGS_Density,out->Maxima.AGS_Density,mode);
        MIN_ADD(DMGradDataPasser[i].Minima.AGS_Density,out->Minima.AGS_Density,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(P[i].AGS_Gradients_Density[k],out->Gradients[k].AGS_Density,mode);}
#if (DM_FUZZY > 0)
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(P[i].AGS_Gradients_Psi_Re[k],out->Gradients[k].AGS_Psi_Re,mode);}
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(P[i].AGS_Gradients_Psi_Im[k],out->Gradients[k].AGS_Psi_Im,mode);}
#endif
    } else {
        int k,k2;
        for(k=0;k<3;k++)
        {
            MAX_ADD(DMGradDataPasser[i].Maxima.AGS_Gradients_Density[k],out->Maxima.AGS_Gradients_Density[k],mode);
            MIN_ADD(DMGradDataPasser[i].Minima.AGS_Gradients_Density[k],out->Minima.AGS_Gradients_Density[k],mode);
            for(k2=0;k2<3;k2++) {ASSIGN_ADD_PRESET(P[i].AGS_Gradients2_Density[k2][k],out->Gradients[k].AGS_Gradients_Density[k2],mode);}
#if (DM_FUZZY > 0)
            for(k2=0;k2<3;k2++) {ASSIGN_ADD_PRESET(P[i].AGS_Gradients2_Psi_Re[k2][k],out->Gradients[k].AGS_Gradients_Psi_Re[k2],mode);}
            for(k2=0;k2<3;k2++) {ASSIGN_ADD_PRESET(P[i].AGS_Gradients2_Psi_Im[k2][k],out->Gradients[k].AGS_Gradients_Psi_Im[k2],mode);}
#endif
        }
        // do we need limiters here for the density gradients? Not clear if this all needs computing
    }
}

/* this actually builds the final gradients out of the data passed */
void construct_gradient_DMGrad(double *grad, int i);
void construct_gradient_DMGrad(double *grad, int i)
{
    /* use the NV_T matrix-based gradient estimator */
    int k; double v_tmp[3];
    for(k=0;k<3;k++) {v_tmp[k] = grad[k];}
    for(k=0;k<3;k++) {grad[k] = P[i].NV_T[k][0]*v_tmp[0] + P[i].NV_T[k][1]*v_tmp[1] + P[i].NV_T[k][2]*v_tmp[2];}
}


/* this subroutine does the actual neighbor-element calculations (this is the 'core' of the loop, essentially) */
/*!   -- this subroutine contains no writes to shared memory -- */
int DMGrad_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    /* define variables */
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double hinv, hinv3, hinv4, r2, u;
    struct kernel_DMGrad kernel;
    struct INPUT_STRUCT_NAME local;
    struct OUTPUT_STRUCT_NAME out;
    /* zero memory and import data for local target */
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    memset(&kernel, 0, sizeof(struct kernel_DMGrad));
    if(mode == 0) {particle2in_DMGrad(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    /* check if we should bother doing a neighbor loop */
    if(local.AGS_Hsml <= 0) return 0;
    if(local.GQuant.AGS_Density <= 0) return 0;
    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.AGS_Hsml;
    double h2_i = kernel.h_i*kernel.h_i;
    kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
    int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(local.Type); // determine allowed particle types for search for adaptive gravitational softening terms

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.AGS_Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, AGS_kernel_shared_BITFLAG);
            if(numngb_inbox < 0) {return -2;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if((P[j].Mass <= 0)||(P[j].AGS_Density <= 0)) {continue;} /* make sure neighbor is valid */
                /* calculate position relative to target */
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0]; kernel.dp[1] = local.Pos[1] - P[j].Pos[1]; kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); /*  now find the closest image in the given box size  */
                r2 = kernel.dp[0]*kernel.dp[0] + kernel.dp[1]*kernel.dp[1] + kernel.dp[2]*kernel.dp[2];
                if((r2 <= 0) || (r2 >= h2_i)) continue;
                /* calculate kernel quantities needed below */
                kernel.r = sqrt(r2); u = kernel.r * hinv;
                kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, -1);
                /* DIFFERENCE & SLOPE LIMITING: need to check maxima and minima of particle values in the kernel, to avoid 'overshoot' with our gradient estimators. this check should be among all interacting pairs */
                if(loop_iteration <= 0)
                {
                    double d_rho = P[j].AGS_Density - local.GQuant.AGS_Density;
                    MINMAX_CHECK(d_rho,out.Minima.AGS_Density,out.Maxima.AGS_Density);
                    for(k=0;k<3;k++) {out.Gradients[k].AGS_Density += -kernel.wk_i * kernel.dp[k] * d_rho;} /* sign is important here! */
#if (DM_FUZZY > 0)
                    d_rho = P[j].AGS_Psi_Re_Pred * P[j].AGS_Density / P[j].Mass - local.GQuant.AGS_Psi_Re;
                    for(k=0;k<3;k++) {out.Gradients[k].AGS_Psi_Re += -kernel.wk_i * kernel.dp[k] * d_rho;}
                    d_rho = P[j].AGS_Psi_Im_Pred * P[j].AGS_Density / P[j].Mass - local.GQuant.AGS_Psi_Im;
                    for(k=0;k<3;k++) {out.Gradients[k].AGS_Psi_Im += -kernel.wk_i * kernel.dp[k] * d_rho;}
#endif
                } else {
                    int k2; double d_grad_rho;
                    for(k=0;k<3;k++)
                    {
                        d_grad_rho = P[j].AGS_Gradients_Density[k] - local.GQuant.AGS_Gradients_Density[k];
                        MINMAX_CHECK(d_grad_rho,out.Minima.AGS_Gradients_Density[k],out.Maxima.AGS_Gradients_Density[k]);
                        for(k2=0;k2<3;k2++) {out.Gradients[k2].AGS_Gradients_Density[k] += -kernel.wk_i * kernel.dp[k2] * d_grad_rho;}
#if (DM_FUZZY > 0)
                        d_grad_rho = P[j].AGS_Gradients_Psi_Re[k] - local.GQuant.AGS_Gradients_Psi_Re[k];
                        for(k2=0;k2<3;k2++) {out.Gradients[k2].AGS_Gradients_Psi_Re[k] += -kernel.wk_i * kernel.dp[k2] * d_grad_rho;}
                        d_grad_rho = P[j].AGS_Gradients_Psi_Im[k] - local.GQuant.AGS_Gradients_Psi_Im[k];
                        for(k2=0;k2<3;k2++) {out.Gradients[k2].AGS_Gradients_Psi_Im[k] += -kernel.wk_i * kernel.dp[k2] * d_grad_rho;}
#endif
                    }
                } // loop_iteration
            } // numngb_inbox loop
        } // while(startnode)
        /* continue to open leaves if needed */
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DATAGET_NAME[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
            }
        }
    }
    /* Collect the result at the right place */
    if(mode == 0) {out2particle_DMGrad(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
    return 0;
}




void DMGrad_gradient_calc(void)
{
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    PRINT_STATUS(" ..calculating higher-order gradients for DM density field\n");
    /* initialize data, if needed */
    if(All.Time==All.TimeBegin) {int i; for(i=FirstActiveParticle; i>=0; i=NextActiveParticle[i]) {P[i].AGS_Numerical_QuantumPotential=0;}}

    /* allocate memory shared across all loops */
    DMGradDataPasser = (struct temporary_dmgradients_data_topass *) mymalloc("DMGradDataPasser",NumPart * sizeof(struct temporary_dmgradients_data_topass));
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */

    /* loop over the number of iterations needed to actually compute the gradients fully */
    for(loop_iteration=0; loop_iteration<2; loop_iteration++) // need 2 iterations to compute gradients-of-gradients
    {
        #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */

        /* do post-loop operations on the results */
        int i;
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        {
            if(loop_iteration <= 0)
            {
                /* now we can properly calculate (second-order accurate) gradients of hydrodynamic quantities from this loop */
                construct_gradient_DMGrad(P[i].AGS_Gradients_Density,i);
#if (DM_FUZZY > 0)
                construct_gradient_DMGrad(P[i].AGS_Gradients_Psi_Re,i);
                construct_gradient_DMGrad(P[i].AGS_Gradients_Psi_Im,i);
#endif
                /* finally, we need to apply a sensible slope limiter to the gradients, to prevent overshooting */
                /* (actually not clear that we need to slope-limit these, because we are not using the gradients for reconstruction.
                    testing this now. if not, we can remove the limiter information entirely and save some time in these computations) */
            } else {
                int k;
                for(k=0;k<3;k++)
                {
                    /* construct the gradient-of-gradient */
                    construct_gradient_DMGrad(P[i].AGS_Gradients2_Density[k],i);
#if (DM_FUZZY > 0)
                    construct_gradient_DMGrad(P[i].AGS_Gradients2_Psi_Re[k],i);
                    construct_gradient_DMGrad(P[i].AGS_Gradients2_Psi_Im[k],i);
#endif
                }
                /* symmetrize the gradients */
                int k0[3]={0,0,1},k1[3]={1,2,2}; double tmp;
                for(k=0;k<3;k++)
                {
                    tmp = 0.5 * (P[i].AGS_Gradients2_Density[k0[k]][k1[k]] + P[i].AGS_Gradients2_Density[k1[k]][k0[k]]);
                    P[i].AGS_Gradients2_Density[k0[k]][k1[k]] = P[i].AGS_Gradients2_Density[k1[k]][k0[k]] = tmp;
#if (DM_FUZZY > 0)
                    tmp = 0.5 * (P[i].AGS_Gradients2_Psi_Re[k0[k]][k1[k]] + P[i].AGS_Gradients2_Psi_Re[k1[k]][k0[k]]);
                    P[i].AGS_Gradients2_Psi_Re[k0[k]][k1[k]] = P[i].AGS_Gradients2_Psi_Re[k1[k]][k0[k]] = tmp;
                    tmp = 0.5 * (P[i].AGS_Gradients2_Psi_Im[k0[k]][k1[k]] + P[i].AGS_Gradients2_Psi_Im[k1[k]][k0[k]]);
                    P[i].AGS_Gradients2_Psi_Im[k0[k]][k1[k]] = P[i].AGS_Gradients2_Psi_Im[k1[k]][k0[k]] = tmp;
#endif
                }
            }
        }
    } // end of loop_iteration

    /* de-allocate memory and collect timing information */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    myfree(DMGradDataPasser); /* free the temporary structure we created for the MinMax and additional data passing */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait; CPU_Step[CPU_AGSDENSCOMM] += timecomm;
    CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */




#endif // DM_FUZZY
