#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*
 
 This module contains the relevant module physics for various elastic,
   visco-elastic, and plastic body simulations.
 
 This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 
 */



#ifdef EOS_TILLOTSON

/* routine that defines the Tillotson EOS parameters for various solid and liquid materials */
void tillotson_eos_init(void)
{   // order: parameter a,b,u0,rho0,A,B,u_s,u_s_prime,alpha,beta,elastic shear modulus, hugoniot elastic limit [all cgs] //

    // initialize pre-computed material properties //
    double qtmp[6][12]={
        {0.5,1.3,1.60e11,2.700,1.80e11,1.80e11,3.50e10,1.80e11,5.0,5.0,2.17e11,3.8e10},  // granite
        {0.5,1.5,4.87e12,2.700,2.67e11,2.67e11,4.72e10,1.82e11,5.0,5.0,2.27e11,3.5e10},  // basalt
        {0.5,1.5,9.50e10,7.860,1.28e12,1.05e12,1.42e10,8.45e10,5.0,5.0,7.75e11,8.5e10},  // iron
        {0.3,0.1,1.00e11,0.917,9.47e10,9.47e10,7.730e9,3.04e10,10.,5.0,2.80e10,1.0e10},  // ice
        {0.5,1.4,5.50e12,3.500,1.31e12,4.90e11,4.50e10,1.50e11,5.0,5.0,8.13e11,9.0e10},  // olivine/dunite
        {0.5,0.9,2.00e10,1.000,2.00e11,1.00e11,4.000e9,2.00e10,5.0,5.0,1.000e0,1.00e0}}; // water
    int j_t,k_t;
    for(j_t=1;j_t<7;j_t++)
    {
        for(k_t=0;k_t<12;k_t++)
        {
            All.Tillotson_EOS_params[j_t][k_t] = qtmp[j_t-1][k_t];
            if((k_t==2)||(k_t==6)||(k_t==7)) {All.Tillotson_EOS_params[j_t][k_t] /= UNIT_SPECEGY_IN_CGS;}
            if(k_t==3) {All.Tillotson_EOS_params[j_t][k_t] /= UNIT_DENSITY_IN_CGS;}
            if((k_t==4)||(k_t==5)||(k_t==10)||(k_t==11)) {All.Tillotson_EOS_params[j_t][k_t] /= UNIT_PRESSURE_IN_CGS;}
        }
    }
    return;
}


/* routine to calculate the pressure and sound speed from the Tillotson equation-of-state for solids */
double calculate_eos_tillotson(int i)
{
    int type = SphP[i].CompositionType; /* determine material, which determines relevant coefficients */
    double a=All.Tillotson_EOS_params[type][0], b=All.Tillotson_EOS_params[type][1],
    u0=All.Tillotson_EOS_params[type][2], rho0=All.Tillotson_EOS_params[type][3],
    A0=All.Tillotson_EOS_params[type][4], B0=All.Tillotson_EOS_params[type][5],
    u_s=All.Tillotson_EOS_params[type][6], u_s_prime=All.Tillotson_EOS_params[type][7],
    alpha=All.Tillotson_EOS_params[type][8], beta=All.Tillotson_EOS_params[type][9]; /* load all the Tillotson EOS parameters for this composition */
    double rho=SphP[i].Density, u=SphP[i].InternalEnergyPred; /* define gas quantities */
    double eta=rho/rho0, mu=eta-1, u_u0eta2=1+u/(u0*eta*eta), p0=u*rho, z=1/eta-1, press=0, cs=0, press_min=1.e-10*u0*rho0; /* useful variables below */
    double Pc = (a + b/u_u0eta2)*p0 + A0*mu + B0*mu*mu; /* pressure in region I,II (fully-condensed states) */
    double c2c_rho = (1+a+b/u_u0eta2)*Pc + A0+B0*(eta*eta-1) + b*(u_u0eta2-1)*(2*p0-Pc)/(u_u0eta2*u_u0eta2); /* sound speed squared (times density) in this region */
    double Pe = a*p0 + (b/u_u0eta2*p0 + A0*mu * exp(-beta*z)) * exp(-alpha*z*z); /* pressure in region IV (fully-vapor states) */
    double c2e_rho = (1+a+b/u_u0eta2*exp(-alpha*z*z))*Pe + A0*eta*(1+mu*(beta+2*alpha*z-eta)/(eta*eta))*exp(-beta*z-alpha*z*z)
    + b*p0/(u_u0eta2*u_u0eta2*eta*eta)*(2*alpha*z*u_u0eta2*eta + (Pe/(u0*rho)-2*u/u0))*exp(-alpha*z*z); /* sound speed squared (times density) in this region */
    double Px = (Pe*(u-u_s) + Pc*(u_s_prime - u)) / (u_s_prime - u_s); /* pressure in region III (interpolation) */
    double c2x_rho = (c2c_rho*(u-u_s) + c2e_rho*(u_s_prime - u)) / (u_s_prime - u_s); /* sound speed squared (times density) in region III (interpolation) */
    if(u <= u_s) {press=Pc; cs=c2c_rho;} else {if(u >= u_s_prime) {press=Pe; cs=c2e_rho;} else {press=Px; cs=c2x_rho;}} /* check which regime we are in */
    //if(press < press_min) {press=cs=press_min;} /* enforce minimum pressure */
    if(cs <= 0) {cs=press_min;} /* enforce non-zero sound speed */
    SphP[i].SoundSpeed = cs/rho; /* save sound speed for later use */
#ifdef EOS_ELASTIC
    SphP[i].SoundSpeed += All.Tillotson_EOS_params[SphP[i].CompositionType][10] / rho; /* add elastic component to soundspeed */
#endif
    SphP[i].SoundSpeed = sqrt(SphP[i].SoundSpeed);
    return press; /* return pressure */
}

#endif



#ifdef EOS_ELASTIC
/* routine to update the deviatoric stress tensor */
void elastic_body_update_driftkick(int i, double dt_entr, int mode)
{
    int j,k,l,NDim=NUMDIMS;
    double dv0[3][3], R[3][3], S[3][3], S_new[3][3], dS=0, mu, Y0, J2=0, I1=0;
#ifdef EOS_TILLOTSON
    mu = All.Tillotson_EOS_params[SphP[i].CompositionType][10]; Y0 = All.Tillotson_EOS_params[SphP[i].CompositionType][11]; // set for composition
#else
    mu = All.Tillotson_EOS_params[0][10]; Y0 = All.Tillotson_EOS_params[0][11];  // set to universal constants
#endif

    if(mode < 2) // drift or kick operation
    {
        for(j=0;j<NDim;j++) {
            for(k=0;k<NDim;k++) {
                // determine which variable we are updating (mode=0/1 is kick/drift)
                if(mode==0) {S_new[j][k]=SphP[i].Elastic_Stress_Tensor[j][k];} else {S_new[j][k]=SphP[i].Elastic_Stress_Tensor_Pred[j][k];}
                S_new[j][k] += dt_entr * SphP[i].Dt_Elastic_Stress_Tensor[j][k]; // apply time evolution
                if(k==j) {I1 += S_new[j][k];} // first invariant of the tensor
                J2 += 0.5 * S_new[j][k]*S_new[j][k]; // second invariant of the tensor
            }}
        // now apply the von Mises yield criterion //
        if(J2 > 0)
        {
            double f_Y = Y0*Y0/(NDim*J2);
            if(f_Y < 1) {for(j=0;j<NDim;j++) {for(k=0;k<NDim;k++) {S_new[j][k] *= f_Y;}}}
        }
        // write out to variable //
        for(j=0;j<NDim;j++) {
            for(k=0;k<NDim;k++) {
                if(mode==0) {SphP[i].Elastic_Stress_Tensor[j][k]=S_new[j][k];} else {SphP[i].Elastic_Stress_Tensor_Pred[j][k]=S_new[j][k];}
            }}

    } else {

        // ok all below is for mode = 2, which is the actual calculation of the time derivative of the stress tensor
        for(j=0;j<NDim;j++) {for(k=0;k<NDim;k++) {dv0[j][k] = SphP[i].Gradients.Velocity[j][k];}}
        for(j=0;j<NDim;j++) {for(k=0;k<NDim;k++) {S[j][k] = SphP[i].Elastic_Stress_Tensor_Pred[j][k];}}
        for(j=0;j<NDim;j++) {for(k=0;k<NDim;k++) {R[j][k] = 0.5*(dv0[j][k] - dv0[k][j]);}}
        double trace_vel=0; for(j=0;j<NDim;j++) {trace_vel += dv0[j][j];}
        for(j=0;j<NDim;j++) // velocity index
        {
            for(k=0;k<NDim;k++) // gradient index
            {
                dS = mu * (dv0[j][k] + dv0[k][j]); // symmetric strain component
                if(k==j) {dS -= 2.*mu*trace_vel/NDim;} // trace component
                for(l=0;l<NDim;l++) {dS += S[j][l]*R[l][k] - R[j][l]*S[l][k];} // rotation components
                SphP[i].Dt_Elastic_Stress_Tensor[j][k] = dS; // save it to variable
            }
        }

    }
    
    return; // all done here
}
#endif



#if defined(EOS_ELASTIC) || defined(EOS_TILLOTSON)
/* routine to get and define the correction factor needed to prevent tensile instability for negative pressures, for arbitrary kernels & dimensions */
double get_negative_pressure_tensilecorrfac(double r, double h_i, double h_j)
{
    double dx_ips=0, wk_0=0, dwk_tmp=0, wk_r=0, r_over_heff=0;
#if (NUMDIMS==1)
    dx_ips = 2. / All.DesNumNgb; // 1D inter-node separation for desired NNgb, relative to radius of compact support
#elif(NUMDIMS==2)
    dx_ips = sqrt(M_PI / All.DesNumNgb); // 2D inter-node separation for desired NNgb, relative to radius of compact support
#else
    dx_ips = pow(4.*M_PI/3. / All.DesNumNgb, 1./3.); // 3D inter-node separation for desired NNgb, relative to radius of compact support
#endif
    kernel_main(dx_ips, 1., 1., &wk_0, &dwk_tmp, -1); // use kernels because of their stability properties: here weight for 'mean separation'
    r_over_heff = r / DMAX(h_i, h_j);
    kernel_main(r_over_heff, 1., 1., &wk_r, &dwk_tmp, -1); // here weight for actual half-separation
    return 0.2 * pow(wk_r / wk_0, 4); // correction factor for n=4 from Monaghan et al. 2000, Gray et al. 2001
}
#endif
