#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif



/*! \file gradients.c
 *  \brief calculate gradients of hydro quantities
 *
 *  This file contains the "second hydro loop", where the gas hydro quantity
 *   gradients are calculated. All gradients now use the second-order accurate
 *   moving-least-squares formulation, and are calculated here consistently.
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#define ASSIGN_ADD_PRESET(x,y,mode) (x+=y)
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#if defined(COOLING)
#define SHOULD_I_USE_SPH_GRADIENTS(condition_number) ((condition_number > CONDITION_NUMBER_DANGER) ? (1):(0))
#else
#define SHOULD_I_USE_SPH_GRADIENTS(condition_number) ((condition_number > CONDITION_NUMBER_DANGER) ? (0):(0))
#endif


#if defined(MHD_CONSTRAINED_GRADIENT)
#if (MHD_CONSTRAINED_GRADIENT > 1)
#define NUMBER_OF_GRADIENT_ITERATIONS 3
#else
#define NUMBER_OF_GRADIENT_ITERATIONS 2
#endif
#else
#define NUMBER_OF_GRADIENT_ITERATIONS 1
#endif



#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

#define NV_MYSIGN(x) (( x > 0 ) - ( x < 0 ))

void *GasGrad_evaluate_primary(void *p, int gradient_iteration);
void *GasGrad_evaluate_secondary(void *p, int gradient_iteration);


/* function that tells us whether a given element should be active for gradient calculation*/
int GasGrad_isactive(int i)
{
    if(P[i].Type != 0) return 0;
    if(P[i].Mass <= 0) return 0;
    if(SphP[i].Density <= 0 || PPP[i].Hsml <= 0) return 0;
#if defined(GALSF_SUBGRID_WINDS) && !defined(TURB_DIFF_DYNAMIC)
    if(SphP[i].DelayTime > 0) return 0;
#endif
    return 1;
}

/* define a common 'gradients' structure to hold
 everything we're going to take derivatives of */
struct Quantities_for_Gradients
{
    MyDouble Density;
    MyDouble Pressure;
    MyDouble Velocity[3];
#ifdef MAGNETIC
    MyDouble B[3];
#ifdef DIVBCLEANING_DEDNER
    MyDouble Phi;
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
    MyDouble Metallicity[NUM_METAL_SPECIES];
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
    MyFloat Rad_E_gamma[N_RT_FREQ_BINS];
    MyFloat Rad_E_gamma_ET[N_RT_FREQ_BINS][6];
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
    MyDouble InternalEnergy;
#endif
#ifdef COSMIC_RAY_FLUID
    MyDouble CosmicRayPressure[N_CR_PARTICLE_BINS];
#endif
#ifdef DOGRAD_SOUNDSPEED
    MyDouble SoundSpeed;
#endif
#ifdef TURB_DIFF_DYNAMIC
    MyDouble Velocity_bar[3];
#endif
};

struct kernel_GasGrad
{
    double dp[3],r,wk_i, wk_j, dwk_i, dwk_j,h_i;
};

struct GasGraddata_in
{
    MyDouble Pos[3];
    MyFloat Mass;
    MyFloat Hsml;
    integertime Timestep;
#ifdef MHD_CONSTRAINED_GRADIENT
    MyDouble ConditionNumber;
    MyLongDouble NV_T[3][3];
    MyFloat BGrad[3][3];
#ifdef MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV
    MyFloat PhiGrad[3];
#endif
#endif
    int NodeList[NODELISTLENGTH];
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
    MyFloat NV_DivVel;
#endif
    struct Quantities_for_Gradients GQuant;

#ifdef TURB_DIFF_DYNAMIC
    MyDouble Norm_hat;
#ifdef GALSF_SUBGRID_WINDS
    MyFloat DelayTime;
#endif
#endif
}
*GasGradDataIn, *GasGradDataGet;


struct GasGraddata_out
{
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS)
    MyDouble third_moment_weights[10];
    MyDouble first_second_moment_weights[HYDRO_TENSOR_FACE_CORRECTIONS_NUMBER_MOMWTS];
#endif
#if defined(KERNEL_CRK_FACES)
    MyDouble m0;
    MyDouble m1[3];
    MyDouble m2[6];
    MyDouble dm0[3];
    MyDouble dm1[3][3];
    MyDouble dm2[6][3];
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
    MyFloat GlassAcc[3];
#endif
#ifdef HYDRO_SPH
    MyFloat alpha_limiter;
#ifdef MAGNETIC
#ifdef DIVBCLEANING_DEDNER
    MyFloat divB;
#endif
    MyFloat DtB[3];
#endif
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
    MyFloat Face_Area[3];
    MyFloat FaceDotB;
    MyFloat FaceCrossX[3][3];
#endif
    struct Quantities_for_Gradients Gradients[3];
    struct Quantities_for_Gradients Maxima;
    struct Quantities_for_Gradients Minima;
    MyFloat MaxDistance;
#ifdef TURB_DIFF_DYNAMIC
    MyDouble Velocity_hat[3];
#endif
}
*GasGradDataResult, *GasGradDataOut;



struct GasGraddata_out_iter
{
#ifdef MHD_CONSTRAINED_GRADIENT
    MyFloat FaceDotB;
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
    MyDouble PhiGrad[3];
#endif
#else
    MyFloat dummy;
#endif
}
*GasGradDataResult_iter, *GasGradDataOut_iter;



/* this is a temporary structure for quantities used ONLY in the loop below,
 for example for computing the slope-limiters (for the Reimann problem) */
static struct temporary_data_topass
{
    struct Quantities_for_Gradients Maxima;
    struct Quantities_for_Gradients Minima;
    MyFloat MaxDistance;
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS)
    MyDouble third_moment_weights[10];
    MyDouble first_second_moment_weights[HYDRO_TENSOR_FACE_CORRECTIONS_NUMBER_MOMWTS];
#endif
#if defined(KERNEL_CRK_FACES)
    MyDouble m0;
    MyDouble m1[3];
    MyDouble m2[6];
    MyDouble dm0[3];
    MyDouble dm1[3][3];
    MyDouble dm2[6][3];
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
    MyFloat GlassAcc[3];
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
    MyDouble FaceDotB;
    MyDouble FaceCrossX[3][3];
    MyDouble BGrad[3][3];
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
    MyDouble PhiGrad[3];
#endif
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
    MyFloat Gradients_Rad_E_gamma[N_RT_FREQ_BINS][3];
#endif
#ifdef TURB_DIFF_DYNAMIC
    MyDouble GradVelocity_bar[3][3];
#endif
}
*GasGradDataPasser;



static inline void particle2in_GasGrad(struct GasGraddata_in *in, int i, int gradient_iteration);
static inline void out2particle_GasGrad(struct GasGraddata_out *out, int i, int mode, int gradient_iteration);
static inline void out2particle_GasGrad_iter(struct GasGraddata_out_iter *out, int i, int mode, int gradient_iteration);



static inline void particle2in_GasGrad(struct GasGraddata_in *in, int i, int gradient_iteration)
{
    int k;
    for(k = 0; k < 3; k++) {in->Pos[k] = P[i].Pos[k];}
    in->Hsml = PPP[i].Hsml;
    in->Mass = P[i].Mass;
    if(in->Mass < 0) {in->Mass = 0;}

#ifdef TURB_DIFF_DYNAMIC
    in->Norm_hat = SphP[i].Norm_hat;
#ifdef GALSF_SUBGRID_WINDS
    in->DelayTime = SphP[i].DelayTime;
#endif
#endif

    if(SHOULD_I_USE_SPH_GRADIENTS(SphP[i].ConditionNumber)) {in->Mass *= -1;}
    in->Timestep = GET_PARTICLE_INTEGERTIME(i);
#ifdef MHD_CONSTRAINED_GRADIENT
    in->ConditionNumber = SphP[i].ConditionNumber;
    if(gradient_iteration > 0) {if(SphP[i].FlagForConstrainedGradients <= 0) {in->Mass = 0;}}
    int j;
    for(j=0;j<3;j++)
    {
        for(k=0;k<3;k++)
        {
            in->BGrad[j][k] = SphP[i].Gradients.B[j][k];
            in->NV_T[j][k] = SphP[i].NV_T[j][k];
        }
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
        in->PhiGrad[j] = SphP[i].Gradients.Phi[j];
#endif
    }
    if(gradient_iteration > 0)
    {
        for(k = 0; k < 3; k++) {in->GQuant.B[k] = Get_Gas_BField(i,k);}
        in->GQuant.Density = SphP[i].Density;
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
        in->GQuant.Phi = Get_Gas_PhiField(i);
#endif
    }
#endif

    if(gradient_iteration == 0)
    {
        in->GQuant.Density = SphP[i].Density;
        in->GQuant.Pressure = SphP[i].Pressure;
        for(k = 0; k < 3; k++) {in->GQuant.Velocity[k] = SphP[i].VelPred[k];}
#ifdef TURB_DIFF_DYNAMIC
        for (k = 0; k < 3; k++) {in->GQuant.Velocity_bar[k] = SphP[i].Velocity_bar[k];}
#endif
#ifdef MAGNETIC
        for(k = 0; k < 3; k++) {in->GQuant.B[k] = Get_Gas_BField(i,k);}
#ifdef DIVBCLEANING_DEDNER
        in->GQuant.Phi = Get_Gas_PhiField(i);
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        for(k = 0; k < NUM_METAL_SPECIES; k++) {in->GQuant.Metallicity[k] = P[i].Metallicity[k];}
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
        for(k = 0; k < N_RT_FREQ_BINS; k++)
        {
        	in->GQuant.Rad_E_gamma[k] = SphP[i].Rad_E_gamma_Pred[k];
        	int k_et; for(k_et = 0; k_et < 6; k_et++) {in->GQuant.Rad_E_gamma_ET[k][k_et] = SphP[i].ET[k][k_et];}
        }
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
        in->GQuant.InternalEnergy = SphP[i].InternalEnergyPred;
#endif
#ifdef COSMIC_RAY_FLUID
        for(k=0;k<N_CR_PARTICLE_BINS;k++) {in->GQuant.CosmicRayPressure[k] = Get_Gas_CosmicRayPressure(i,k);}
#endif
#ifdef DOGRAD_SOUNDSPEED
        in->GQuant.SoundSpeed = Get_Gas_effective_soundspeed_i(i);
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
       in->NV_DivVel = SphP[i].NV_DivVel;
#endif
    } // gradient_iteration == 0
}



//#define MAX_ADD(x,y,mode) (mode == 0 ? (x=y) : (((x)<(y)) ? (x=y) : (x))) // these definitions applied before the symmetric re-formulation of this routine
//#define MIN_ADD(x,y,mode) (mode == 0 ? (x=y) : (((x)>(y)) ? (x=y) : (x)))
#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))


static inline void out2particle_GasGrad_iter(struct GasGraddata_out_iter *out, int i, int mode, int gradient_iteration)
{
#ifdef MHD_CONSTRAINED_GRADIENT
    {
        ASSIGN_ADD_PRESET(GasGradDataPasser[i].FaceDotB,out->FaceDotB,mode);
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
        int k;
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].PhiGrad[k],out->PhiGrad[k],mode);}
#endif
    }
#endif
}



static inline void out2particle_GasGrad(struct GasGraddata_out *out, int i, int mode, int gradient_iteration)
{
#ifdef MHD_CONSTRAINED_GRADIENT
    {
        ASSIGN_ADD_PRESET(GasGradDataPasser[i].FaceDotB,out->FaceDotB,mode);
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
        int k;
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].PhiGrad[k],out->Gradients[k].Phi,mode);}
#endif
    }
#endif

    if(gradient_iteration == 0)
    {
        int j,k;
        MAX_ADD(GasGradDataPasser[i].MaxDistance,out->MaxDistance,mode);
#ifdef TURB_DIFF_DYNAMIC
        for (j = 0; j < 3; j++) {
            MAX_ADD(GasGradDataPasser[i].Maxima.Velocity_bar[j], out->Maxima.Velocity_bar[j], mode);
            MIN_ADD(GasGradDataPasser[i].Minima.Velocity_bar[j], out->Minima.Velocity_bar[j], mode);
            ASSIGN_ADD_PRESET(SphP[i].Velocity_hat[j], out->Velocity_hat[j], mode);
            for (k = 0; k < 3; k++) {
                ASSIGN_ADD_PRESET(GasGradDataPasser[i].GradVelocity_bar[j][k], out->Gradients[k].Velocity_bar[j], mode);
            }
        }
#endif

#if defined(HYDRO_TENSOR_FACE_CORRECTIONS)
        for(k=0;k<10;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].third_moment_weights[k],out->third_moment_weights[k],mode);}
        for(k=0;k<HYDRO_TENSOR_FACE_CORRECTIONS_NUMBER_MOMWTS;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].first_second_moment_weights[k],out->first_second_moment_weights[k],mode);}
#endif
#if defined(KERNEL_CRK_FACES)
        ASSIGN_ADD_PRESET(GasGradDataPasser[i].m0,out->m0,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].dm0[k],out->dm0[k],mode);}
        for(j=0;j<3;j++)
        {
            ASSIGN_ADD_PRESET(GasGradDataPasser[i].m1[j],out->m1[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].dm1[j][k],out->dm1[j][k],mode);}
        }
        for(j=0;j<6;j++)
        {
            ASSIGN_ADD_PRESET(GasGradDataPasser[i].m2[j],out->m2[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].dm2[j][k],out->dm2[j][k],mode);}
        }
#endif

#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].GlassAcc[k],out->GlassAcc[k],mode);}
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        ASSIGN_ADD_PRESET(SphP[i].alpha_limiter, out->alpha_limiter, mode);
#endif

        MAX_ADD(GasGradDataPasser[i].Maxima.Density,out->Maxima.Density,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.Density,out->Minima.Density,mode);
        MAX_ADD(GasGradDataPasser[i].Maxima.Pressure,out->Maxima.Pressure,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.Pressure,out->Minima.Pressure,mode);
        for(k=0;k<3;k++)
        {
            ASSIGN_ADD_PRESET(SphP[i].Gradients.Density[k],out->Gradients[k].Density,mode);
            ASSIGN_ADD_PRESET(SphP[i].Gradients.Pressure[k],out->Gradients[k].Pressure,mode);
        }
#ifdef DOGRAD_INTERNAL_ENERGY
        MAX_ADD(GasGradDataPasser[i].Maxima.InternalEnergy,out->Maxima.InternalEnergy,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.InternalEnergy,out->Minima.InternalEnergy,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(SphP[i].Gradients.InternalEnergy[k],out->Gradients[k].InternalEnergy,mode);}
#endif
#ifdef COSMIC_RAY_FLUID
        for(j=0;j<N_CR_PARTICLE_BINS;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.CosmicRayPressure[j],out->Maxima.CosmicRayPressure[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.CosmicRayPressure[j],out->Minima.CosmicRayPressure[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(SphP[i].Gradients.CosmicRayPressure[j][k],out->Gradients[k].CosmicRayPressure[j],mode);}
        }
#endif
#ifdef DOGRAD_SOUNDSPEED
        MAX_ADD(GasGradDataPasser[i].Maxima.SoundSpeed,out->Maxima.SoundSpeed,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.SoundSpeed,out->Minima.SoundSpeed,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(SphP[i].Gradients.SoundSpeed[k],out->Gradients[k].SoundSpeed,mode);}
#endif

        for(j=0;j<3;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.Velocity[j],out->Maxima.Velocity[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.Velocity[j],out->Minima.Velocity[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(SphP[i].Gradients.Velocity[j][k],out->Gradients[k].Velocity[j],mode);}
        }

#ifdef MAGNETIC

#ifdef HYDRO_SPH
#ifdef DIVBCLEANING_DEDNER
        ASSIGN_ADD_PRESET(SphP[i].divB,out->divB, mode);
#endif
        for(k = 0; k < 3; k++) {ASSIGN_ADD_PRESET(SphP[i].DtB[k],out->DtB[k], mode);}
#endif


#ifdef MHD_CONSTRAINED_GRADIENT
        for(j=0;j<3;j++)
        {
            ASSIGN_ADD_PRESET(SphP[i].Face_Area[j],out->Face_Area[j],mode);
            for(k=0;k<3;k++)
            {
                ASSIGN_ADD_PRESET(GasGradDataPasser[i].BGrad[j][k],out->Gradients[k].B[j],mode);
                ASSIGN_ADD_PRESET(GasGradDataPasser[i].FaceCrossX[j][k],out->FaceCrossX[j][k],mode);
            }
        }
#endif

        for(j=0;j<3;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.B[j],out->Maxima.B[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.B[j],out->Minima.B[j],mode);
            for(k=0;k<3;k++)
            {
#ifndef MHD_CONSTRAINED_GRADIENT
                ASSIGN_ADD_PRESET(SphP[i].Gradients.B[j][k],out->Gradients[k].B[j],mode);
#endif
            }
        }

#ifdef DIVBCLEANING_DEDNER
        MAX_ADD(GasGradDataPasser[i].Maxima.Phi,out->Maxima.Phi,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.Phi,out->Minima.Phi,mode);
#ifndef MHD_CONSTRAINED_GRADIENT_MIDPOINT
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(SphP[i].Gradients.Phi[k],out->Gradients[k].Phi,mode);}
#endif
#endif
#endif // closes MAGNETIC

#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        for(j=0;j<NUM_METAL_SPECIES;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.Metallicity[j],out->Maxima.Metallicity[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.Metallicity[j],out->Minima.Metallicity[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(SphP[i].Gradients.Metallicity[j][k],out->Gradients[k].Metallicity[j],mode);}
        }
#endif

#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
        for(j=0;j<N_RT_FREQ_BINS;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.Rad_E_gamma[j],out->Maxima.Rad_E_gamma[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.Rad_E_gamma[j],out->Minima.Rad_E_gamma[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].Gradients_Rad_E_gamma[j][k],out->Gradients[k].Rad_E_gamma[j],mode);}
        }
		/* the gradient dotted into the Eddington tensor is more complicated: let's handle this below */
        {
        	int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
        	{
        		int k_xyz; for(k_xyz=0;k_xyz<3;k_xyz++)
        		{
        			int j_xyz,i_xyz,k_et_loop[3]; // recall, for ET: 0=xx,1=yy,2=zz,3=xy,4=yz,5=xz
					if(k_xyz==0) {k_et_loop[0]=0; k_et_loop[1]=3; k_et_loop[2]=5;}
					if(k_xyz==1) {k_et_loop[0]=3; k_et_loop[1]=1; k_et_loop[2]=4;}
					if(k_xyz==2) {k_et_loop[0]=5; k_et_loop[1]=4; k_et_loop[2]=2;}
        			for(j_xyz=0;j_xyz<3;j_xyz++)
        			{
        				for(i_xyz=0;i_xyz<3;i_xyz++)
        				{
        					SphP[i].Gradients.Rad_E_gamma_ET[k_freq][k_xyz] += SphP[i].NV_T[j_xyz][i_xyz] * out->Gradients[i_xyz].Rad_E_gamma_ET[k_freq][k_et_loop[j_xyz]];
						}
        			}
        		}
        	}
        }
#endif
    } // gradient_iteration == 0
}




void local_slopelimiter(double *grad, double valmax, double valmin, double alim, double h, double shoot_tol, int pos_preserve, double d_max, double val_cen)
{
    int k;
    double d_abs = 0.0;
    for(k=0;k<3;k++) {d_abs += grad[k]*grad[k];}
    if(d_abs > 0)
    {
        d_abs=sqrt(d_abs); double cfac = 1 / (alim * h * d_abs); /* inverse change over distance for limiter */
        double fabs_max = fabs(valmax), fabs_min = fabs(valmin), abs_max=fabs_max, abs_min=fabs_min, f_corr_overshoot;
        if(abs_max<abs_min) {abs_max=fabs_min; abs_min=fabs_max;} /* get largest positive/negative deviations, determine smaller in absolute value */
        f_corr_overshoot = DMIN(abs_min + shoot_tol*abs_max, abs_max); /* = abs_min for shoot_tol = 0; don't let gradient deviate by more than this in size, slightly larger if 'shoot_tol' allows some overshoot tolerance */
        cfac *= f_corr_overshoot; /* multiply by the correction factor of interest */
        if(pos_preserve == 1) /* demand that the limited slope be strictly positivity-preserving over the maximal range to any neighbors */
        {
            double fmin = DMIN(val_cen, DMAX(0, DMAX(MIN_REAL_NUMBER*val_cen, DMIN(0.5*(val_cen+valmin), val_cen-f_corr_overshoot)))); /* minimum value: smaller of overshoot target or half positive-definite value, but cannot go negative in larger range */
            cfac = DMIN( (((val_cen-fmin) / d_max) / d_abs) , cfac ); /* use more conservative limiter, of cfac above or this, over longer range d_max, to restrict here */
        }
        if(cfac < 1) {for(k=0;k<3;k++) {grad[k] *= cfac;}} /* scalar gradient correction */
    }
}


void construct_gradient(double *grad, int i)
{
    /* check if the matrix is well-conditioned: otherwise we will use the 'standard SPH-like' derivative estimation */
    if(SHOULD_I_USE_SPH_GRADIENTS(SphP[i].ConditionNumber))
    {
        /* the condition number was bad, so we used SPH-like gradients */
        if(SphP[i].Density > 0) {int k; for(k=0;k<3;k++) {grad[k] *= PPP[i].DhsmlNgbFactor / SphP[i].Density;}}
    } else {
        /* ok, the condition number was good so we used the matrix-like gradient estimator */
        int k; double v_tmp[3];
        for(k=0;k<3;k++) {v_tmp[k] = grad[k];}
        for(k=0;k<3;k++) {grad[k] = SphP[i].NV_T[k][0]*v_tmp[0] + SphP[i].NV_T[k][1]*v_tmp[1] + SphP[i].NV_T[k][2]*v_tmp[2];}
    }
}




void hydro_gradient_calc(void)
{
    CPU_Step[CPU_DENSMISC] += measure_time(); double t0 = my_second();
    int i, j, k, k1, ndone, ndone_flag, recvTask, place, save_NextParticle;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
    double timecomp, timecomm, timewait, tstart, tend, t1;
    long long n_exported = 0;
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
    double NV_dt,NV_dummy,NV_limiter,NV_A,divVel_physical,h_eff,alphaloc,cs_nv;
#endif
#ifdef TURB_DIFF_DYNAMIC
    double smoothInv = 1.0 / All.TurbDynamicDiffSmoothing;
#endif

    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    GasGradDataPasser = (struct temporary_data_topass *) mymalloc("GasGradDataPasser",N_gas * sizeof(struct temporary_data_topass));
    NTaskTimesNumPart = maxThreads * NumPart; size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct GasGraddata_in) + sizeof(struct GasGraddata_out) +
                                                             sizemax(sizeof(struct GasGraddata_in),sizeof(struct GasGraddata_out))));
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    /* before doing any operations, need to zero the appropriate memory so we can correctly do pair-wise operations */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type==0)
        {
            int k2;
            memset(&GasGradDataPasser[i], 0, sizeof(struct temporary_data_topass));
#ifdef HYDRO_SPH
#ifdef MAGNETIC
            for(k = 0; k < 3; k++) {SphP[i].DtB[k] = 0;}
#endif
#ifdef DIVBCLEANING_DEDNER
            SphP[i].divB = 0;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
            SphP[i].alpha_limiter = 0;
#endif
#endif
#ifdef TURB_DIFF_DYNAMIC
            /* NEED Velocity_bar CORRECT HERE */
            for (k = 0; k < 3; k++) {
                SphP[i].Velocity_bar[k] *= All.TurbDynamicDiffSmoothing;
                SphP[i].Velocity_hat[k] = SphP[i].Velocity_bar[k] * smoothInv;
            }
#endif

            /* and zero out the gradients structure itself */
            for(k=0;k<3;k++)
            {
                SphP[i].Gradients.Density[k] = 0;
                SphP[i].Gradients.Pressure[k] = 0;
                for(k2=0;k2<3;k2++) {SphP[i].Gradients.Velocity[k2][k] = 0;}
#ifdef DOGRAD_INTERNAL_ENERGY
                SphP[i].Gradients.InternalEnergy[k] = 0;
#endif
#ifdef COSMIC_RAY_FLUID
                for(k2=0;k2<N_CR_PARTICLE_BINS;k2++) {SphP[i].Gradients.CosmicRayPressure[k2][k] = 0;}
#endif
#ifdef DOGRAD_SOUNDSPEED
                SphP[i].Gradients.SoundSpeed[k] = 0;
#endif
#ifdef MAGNETIC
#ifndef MHD_CONSTRAINED_GRADIENT
                for(k2=0;k2<3;k2++) {SphP[i].Gradients.B[k2][k] = 0;}
#else
                SphP[i].Face_Area[k] = 0;
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
                SphP[i].Gradients.Phi[k] = 0;
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
                for(k2=0;k2<NUM_METAL_SPECIES;k2++) {SphP[i].Gradients.Metallicity[k2][k] = 0;}
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
                for(k2=0;k2<N_RT_FREQ_BINS;k2++) {SphP[i].Gradients.Rad_E_gamma_ET[k2][k] = 0;}
#endif
            }
        }



    /* prepare to do the requisite number of sweeps over the particle distribution */
    int gradient_iteration;
    for(gradient_iteration = 0; gradient_iteration < NUMBER_OF_GRADIENT_ITERATIONS; gradient_iteration++)
    {
        // need to zero things used in the iteration (anything appearing in out2particle_GasGrad_iter)
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type==0)
            {
#ifdef MHD_CONSTRAINED_GRADIENT
                GasGradDataPasser[i].FaceDotB = 0;
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
                for(k=0;k<3;k++) {GasGradDataPasser[i].PhiGrad[k] = 0;}
#endif
#endif
            }

        // now we actually begin the main gradient loop //
        NextParticle = FirstActiveParticle;	/* begin with this index */
        do
        {
            BufferFullFlag = 0; Nexport = 0; save_NextParticle = NextParticle; tstart = my_second();
            for(j = 0; j < NTask; j++) {Send_count[j] = 0; Exportflag[j] = -1;} /* do local particles and prepare export list */
#ifdef PTHREADS_NUM_THREADS
            pthread_t mythreads[PTHREADS_NUM_THREADS - 1]; int threadid[PTHREADS_NUM_THREADS - 1]; pthread_attr_t attr;
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            pthread_mutex_init(&mutex_nexport, NULL); pthread_mutex_init(&mutex_partnodedrift, NULL);
            TimerFlag = 0; for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {threadid[j] = j + 1; pthread_create(&mythreads[j], &attr, GasGrad_evaluate_primary, &threadid[j]);}
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
#ifdef _OPENMP
                int mainthreadid = omp_get_thread_num();
#else
                int mainthreadid = 0;
#endif
                GasGrad_evaluate_primary(&mainthreadid, gradient_iteration);	/* do local particles and prepare export list */
            }
#ifdef PTHREADS_NUM_THREADS
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
#endif
            tend = my_second(); timecomp1 += timediff(tstart, tend);

            if(BufferFullFlag) /* we've filled the buffer or reached the end of the list, prepare for communications */
            {
                int last_nextparticle = NextParticle; NextParticle = save_NextParticle; /* figure out where we are */
                while(NextParticle >= 0)
                {
                    if(NextParticle == last_nextparticle) {break;}
                    if(ProcessedFlag[NextParticle] != 1) {break;}
                    ProcessedFlag[NextParticle] = 2; NextParticle = NextActiveParticle[NextParticle];
                }
                if(NextParticle == save_NextParticle) {endrun(113308);} /* in this case, the buffer is too small to process even a single particle */

                int new_export = 0; /* actually calculate exports [so we can tell other tasks] */
                for(j = 0, k = 0; j < Nexport; j++)
                {
                    if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                    {
                        if(k < j + 1) {k = j + 1;}
                        for(; k < Nexport; k++)
                            if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                            {
                                int old_index = DataIndexTable[j].Index;
                                DataIndexTable[j] = DataIndexTable[k]; DataNodeList[j] = DataNodeList[k]; DataIndexTable[j].IndexGet = j; new_export++;
                                DataIndexTable[k].Index = old_index; k++;
                                break;
                            }
                    }
                    else {new_export++;}
                }
                Nexport = new_export; /* counting exports... */
            }
            n_exported += Nexport;
            for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
            for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
            MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare); /* construct export count tables */
            tstart = my_second();
            MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD); /* broadcast import/export counts */
            tend = my_second(); timewait1 += timediff(tstart, tend);

            for(j = 0, Send_offset[0] = 0; j < NTask; j++) {if(j > 0) {Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];}} /* calculate export table offsets */
            GasGradDataIn = (struct GasGraddata_in *) mymalloc("GasGradDataIn", Nexport * sizeof(struct GasGraddata_in)); /* allocate memory for exports */
            if(gradient_iteration==0) /* allocate memory for exports: here we have a different structure for the different iterations, which makes this especially complicated */
            {
                GasGradDataOut = (struct GasGraddata_out *) mymalloc("GasGradDataOut", Nexport * sizeof(struct GasGraddata_out));
            } else {
                GasGradDataOut_iter = (struct GasGraddata_out_iter *) mymalloc("GasGradDataOut_iter", Nexport * sizeof(struct GasGraddata_out_iter));
            }
            for(j = 0; j < Nexport; j++) /* prepare particle data for export [fill in the structures to be passed] */
            {
                place = DataIndexTable[j].Index;
                particle2in_GasGrad(&GasGradDataIn[j], place, gradient_iteration);
                memcpy(GasGradDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
            }

            /* ok now we have to figure out if there is enough memory to handle all the tasks sending us their data, and if not, break it into sub-chunks */
            int N_chunks_for_import, ngrp_initial, ngrp;
            for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) /* sub-chunking loop opener */
            {
                int flagall;
                N_chunks_for_import = (1 << PTask) - ngrp_initial;
                do {
                    int flag = 0; Nimport = 0;
                    for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++)
                    {
                        recvTask = ThisTask ^ ngrp;
                        if(recvTask < NTask) {if(Recv_count[recvTask] > 0) {Nimport += Recv_count[recvTask];}}
                    }
                    size_t space_needed = Nimport * sizeof(struct GasGraddata_in) + Nimport * sizeof(struct GasGraddata_out) + 16384; /* extra bitflag is a padding, to avoid overflows */
                    if(space_needed > FreeBytes) {flag = 1;}

                    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                    if(flagall) {N_chunks_for_import /= 2;} else {break;}
                } while(N_chunks_for_import > 0);
                if(N_chunks_for_import == 0) {printf("Memory is insufficient for even one import-chunk: N_chunks_for_import=%d  ngrp_initial=%d  Nimport=%ld  FreeBytes=%lld , but we need to allocate=%lld \n",N_chunks_for_import, ngrp_initial, Nimport, (long long)FreeBytes,(long long)(Nimport * sizeof(struct GasGraddata_in) + Nimport * sizeof(struct GasGraddata_out) + 16384)); endrun(9999);}
                if(flagall) {if(ThisTask==0) PRINT_WARNING("Splitting import operation into sub-chunks as we are hitting memory limits (check this isn't imposing large communication cost)");}

                /* now allocated the import and results buffers */
                GasGradDataGet = (struct GasGraddata_in *) mymalloc("GasGradDataGet", Nimport * sizeof(struct GasGraddata_in));
                if(gradient_iteration==0)
                {
                    GasGradDataResult = (struct GasGraddata_out *) mymalloc("GasGradDataResult", Nimport * sizeof(struct GasGraddata_out));
                } else {
                    GasGradDataResult_iter = (struct GasGraddata_out_iter *) mymalloc("GasGradDataResult_iter", Nimport * sizeof(struct GasGraddata_out_iter));
                }


                tstart = my_second(); Nimport = 0; /* reset because this will be cycled below to calculate the recieve offsets (Recv_offset) */
                for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) /* exchange particle data */
                {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask)
                    {
                        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) /* get the particles */
                        {
                            MPI_Sendrecv(&GasGradDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct GasGraddata_in), MPI_BYTE, recvTask, TAG_GRADLOOP_A,
                                         &GasGradDataGet[Nimport], Recv_count[recvTask] * sizeof(struct GasGraddata_in), MPI_BYTE, recvTask, TAG_GRADLOOP_A,
                                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            Nimport += Recv_count[recvTask];
                        }
                    }
                }
                tend = my_second(); timecommsumm1 += timediff(tstart, tend);



                /* now do the particles that were sent to us */
                tstart = my_second(); NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
                for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_create(&mythreads[j], &attr, GasGrad_evaluate_secondary, &threadid[j]);}
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
                {
#ifdef _OPENMP
                    int mainthreadid = omp_get_thread_num();
#else
                    int mainthreadid = 0;
#endif
                    GasGrad_evaluate_secondary(&mainthreadid, gradient_iteration);
                }
#ifdef PTHREADS_NUM_THREADS
                for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
                pthread_mutex_destroy(&mutex_partnodedrift);
                pthread_mutex_destroy(&mutex_nexport);
                pthread_attr_destroy(&attr);
#endif
                tend = my_second(); timecomp2 += timediff(tstart, tend); tstart = my_second();
                MPI_Barrier(MPI_COMM_WORLD); /* insert MPI Barrier here - will be forced by comms below anyways but this allows for clean timing measurements */
                tend = my_second(); timewait2 += timediff(tstart, tend);

                tstart = my_second(); Nimport = 0;
                for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) /* send the results for imported elements back to their host tasks */
                {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask)
                    {
                        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                        {
                            if(gradient_iteration==0) /* send the results */
                            {
                                MPI_Sendrecv(&GasGradDataResult[Nimport], Recv_count[recvTask] * sizeof(struct GasGraddata_out), MPI_BYTE, recvTask, TAG_GRADLOOP_B,
                                             &GasGradDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct GasGraddata_out), MPI_BYTE, recvTask, TAG_GRADLOOP_B,
                                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            } else {
                                MPI_Sendrecv(&GasGradDataResult_iter[Nimport], Recv_count[recvTask] * sizeof(struct GasGraddata_out_iter), MPI_BYTE, recvTask, TAG_GRADLOOP_C,
                                             &GasGradDataOut_iter[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct GasGraddata_out_iter), MPI_BYTE, recvTask, TAG_GRADLOOP_C,
                                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            }
                            Nimport += Recv_count[recvTask];
                        }
                    }
                }
                tend = my_second(); timecommsumm2 += timediff(tstart, tend);
                if(gradient_iteration==0) {myfree(GasGradDataResult);} else {myfree(GasGradDataResult_iter);} /* free the structures used to send data back to tasks, its sent */
                myfree(GasGradDataGet); /* free the structures used to send data back to tasks, its sent */

            } /* close the sub-chunking loop: for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) */


            /* we have all our results back from the elements we exported: add the result to the local elements */
            tstart = my_second();
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index;
                if(gradient_iteration==0)
                {
                    out2particle_GasGrad(&GasGradDataOut[j], place, 1, gradient_iteration);
                } else {
                    out2particle_GasGrad_iter(&GasGradDataOut_iter[j], place, 1, gradient_iteration);
                }
            }
            tend = my_second(); timecomp1 += timediff(tstart, tend);
            if(gradient_iteration==0) {myfree(GasGradDataOut);} else {myfree(GasGradDataOut_iter);} /* free the structures used to receive results, weve used it */
            myfree(GasGradDataIn); /* free the structures used to prepare our initial export data, we're done here! */

            if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;} /* figure out if we are done with the particular active set here */
            tstart = my_second();
            MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /* call an allreduce to figure out if all tasks are also done here, otherwise we need to iterate */
            tend = my_second(); timewait2 += timediff(tstart, tend);
        }
        while(ndone < NTask);


        /* here, we insert intermediate operations on the results, from the iterations we have completed */
#ifdef MHD_CONSTRAINED_GRADIENT
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 0)
            {
                SphP[i].FlagForConstrainedGradients = 1;
                /* copy everything from the structure holding B-gradients (needed so they dont change mid-loop) */
                for(k=0;k<3;k++)
                {
                    for(k1=0;k1<3;k1++)
                    {
                        SphP[i].Gradients.B[k][k1] = GasGradDataPasser[i].BGrad[k][k1];
                    }
                }
                /* build the gradient */
                for(k=0;k<3;k++) {construct_gradient(SphP[i].Gradients.B[k],i);}
                /* slope limit it */
                double v_tmp = P[i].Mass / SphP[i].Density;
                double tmp_d = sqrt(1.0e-37 + (2. * All.cf_atime/All.cf_afac1 * SphP[i].Pressure*v_tmp*v_tmp) +
                                    SphP[i].BPred[0]*SphP[i].BPred[0]+SphP[i].BPred[1]*SphP[i].BPred[1]+SphP[i].BPred[2]*SphP[i].BPred[2]);
                double tmp = 3.0e3 * fabs(SphP[i].divB) * PPP[i].Hsml / tmp_d;
                double alim; alim = 1. + DMIN(1.,tmp*tmp);
#if (MHD_CONSTRAINED_GRADIENT <= 1)
                double dbmax=0, dbgrad=0;
                double dh=0.25*PPP[i].Hsml; // need to be more aggressive with new wt_i,wt_j formalism
                for(k=0;k<3;k++)
                {
                    double b0 = Get_Gas_BField(i,k);
                    double dd = 2. * fabs(b0) * DMIN(fabs(GasGradDataPasser[i].Minima.B[k]) , fabs(GasGradDataPasser[i].Maxima.B[k]));
                    dbmax = DMIN(fabs(dbmax+dd),fabs(dbmax-dd));
                    for(k1=0;k1<3;k1++) {dbgrad += 2.*dh * fabs(b0*SphP[i].Gradients.B[k][k1]);}
                }
                dbmax /= dbgrad;
                for(k1=0;k1<3;k1++)
                {
                    double d_abs=0; for(k=0;k<3;k++) {d_abs += SphP[i].Gradients.B[k1][k]*SphP[i].Gradients.B[k1][k];}
                    if(d_abs > 0)
                    {
                        double cfac = 1 / (0.25 * PPP[i].Hsml * sqrt(d_abs));
                        cfac *= DMIN(fabs(GasGradDataPasser[i].Maxima.B[k1]) , fabs(GasGradDataPasser[i].Minima.B[k1]));
                        double c_eff = DMIN( cfac , DMAX(cfac/alim , dbmax) );
                        if(c_eff < 1) {for(k=0;k<3;k++) {SphP[i].Gradients.B[k1][k] *= c_eff;}}
                    } else {
                        for(k=0;k<3;k++) {SphP[i].Gradients.B[k1][k]=0;}
                    }
                }
#endif
                /* check the particle area closure, which will inform whether it is safe to use the constrained gradients */
                double area = fabs(SphP[i].Face_Area[0]) + fabs(SphP[i].Face_Area[1]) + fabs(SphP[i].Face_Area[2]);
                area /= Get_Particle_Expected_Area(PPP[i].Hsml);
                /* set the relevant flags to decide whether or not we use the constrained gradients */
                if(area > 0.5) {SphP[i].FlagForConstrainedGradients = 0;}
                if(SphP[i].ConditionNumber > 1000.) {SphP[i].FlagForConstrainedGradients = 0;}
                if(SHOULD_I_USE_SPH_GRADIENTS(SphP[i].ConditionNumber)) {SphP[i].FlagForConstrainedGradients = 0;} /* this must be here, since in this case the SPH gradients are used, which will not work with this method */

                /* now check, and if ok, enter the gradient re-calculation */
                if(SphP[i].FlagForConstrainedGradients == 1)
                {
                    double GB0[3][3];
                    double fsum = 0.0, dmag = 0.0;
                    double h_eff = Get_Particle_Size(i);
                    for(k=0;k<3;k++)
                    {
                        double grad_limiter_mag = Get_Gas_BField(i,k) / h_eff;
                        dmag += grad_limiter_mag * grad_limiter_mag;
                        for(k1=0;k1<3;k1++)
                        {
                            GB0[k][k1] = SphP[i].Gradients.B[k][k1];
                            dmag += GB0[k][k1] * GB0[k][k1];
                            fsum += GasGradDataPasser[i].FaceCrossX[k][k1] * GasGradDataPasser[i].FaceCrossX[k][k1];
                        }
                    }
                    if((fsum <= 0) || (dmag <= 0))
                    {
                        SphP[i].FlagForConstrainedGradients = 0;
                    } else {
                        dmag = 2.0 * sqrt(dmag); // limits the maximum magnitude of the correction term we will allow //
                        fsum = -1 / fsum;
                        int j_gloop;
                        for(j_gloop = 0; j_gloop < 5; j_gloop++)
                        {
                            /* calculate the correction terms */
                            double asum=GasGradDataPasser[i].FaceDotB;
                            for(k=0;k<3;k++)
                            {
                                for(k1=0;k1<3;k1++)
                                {
                                    asum += SphP[i].Gradients.B[k][k1] * GasGradDataPasser[i].FaceCrossX[k][k1];
                                }
                            }
                            double prefac = 1.0 * asum * fsum;
                            double ecorr[3][3];
                            double cmag=0;
                            for(k=0;k<3;k++)
                            {
                                for(k1=0;k1<3;k1++)
                                {
                                    ecorr[k][k1] = prefac * GasGradDataPasser[i].FaceCrossX[k][k1];
                                    double grad_limiter_mag = (SphP[i].Gradients.B[k][k1] + ecorr[k][k1]) - GB0[k][k1];
                                    cmag += grad_limiter_mag * grad_limiter_mag;
                                }
                            }
                            cmag = sqrt(cmag);
                            /* limit the correction term, based on the maximum calculated above */
                            double nnorm = 1.0;
                            if(cmag > dmag) nnorm *= dmag / cmag;
                            /* finally, we can apply the correction */
                            for(k=0;k<3;k++)
                            {
                                for(k1=0;k1<3;k1++)
                                {
                                    SphP[i].Gradients.B[k][k1] = GB0[k][k1] + nnorm*(SphP[i].Gradients.B[k][k1]+ecorr[k][k1] - GB0[k][k1]);
                                }
                                /* slope-limit the corrected gradients again, but with a more tolerant slope-limiter */
#if (MHD_CONSTRAINED_GRADIENT <= 1)
                                local_slopelimiter(SphP[i].Gradients.B[k],GasGradDataPasser[i].Maxima.B[k],GasGradDataPasser[i].Minima.B[k],0.25, PPP[i].Hsml, 0.25, 0, 0, 0);
#endif
                            }
                        } // closes j_gloop loop
                    } // closes fsum/dmag check
                } // closes FlagForConstrainedGradients check
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
                double a_limiter = 0.25; if(SphP[i].ConditionNumber>100) a_limiter=DMIN(0.5, 0.25 + 0.25 * (SphP[i].ConditionNumber-100)/100);
                /* copy everything from the structure holding phi-gradients (needed so they dont change mid-loop) */
                for(k=0;k<3;k++) {SphP[i].Gradients.Phi[k] = GasGradDataPasser[i].PhiGrad[k];}
                /* build and limit the gradient */
                construct_gradient(SphP[i].Gradients.Phi,i);
                local_slopelimiter(SphP[i].Gradients.Phi,GasGradDataPasser[i].Maxima.Phi,GasGradDataPasser[i].Minima.Phi,a_limiter,PPP[i].Hsml,0.0, 0, 0, 0);
#endif
            } // closes Ptype == 0 check
#endif
    } // closes gradient_iteration

    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);


    /* do final operations on results: these are operations that can be done after the complete set of iterations */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type == 0)
        {
            /* now we can properly calculate (second-order accurate) gradients of hydrodynamic quantities from this loop */
            construct_gradient(SphP[i].Gradients.Density,i);
            construct_gradient(SphP[i].Gradients.Pressure,i);
            for(k=0;k<3;k++) {construct_gradient(SphP[i].Gradients.Velocity[k],i);}
#ifdef TURB_DIFF_DYNAMIC
            for (k = 0; k < 3; k++) {
                construct_gradient(GasGradDataPasser[i].GradVelocity_bar[k], i);
            }
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
            construct_gradient(SphP[i].Gradients.InternalEnergy,i);
#endif
#ifdef COSMIC_RAY_FLUID
            for(k=0;k<N_CR_PARTICLE_BINS;k++) {construct_gradient(SphP[i].Gradients.CosmicRayPressure[k],i);}
            int is_particle_local_extremum[N_CR_PARTICLE_BINS]={0}; is_particle_local_extremum[0]=0; // test for local extremum to revert to lower-order reconstruction if necessary
#endif
#ifdef DOGRAD_SOUNDSPEED
            construct_gradient(SphP[i].Gradients.SoundSpeed,i);
#endif
#ifdef MAGNETIC
#ifndef MHD_CONSTRAINED_GRADIENT
            for(k=0;k<3;k++) {construct_gradient(SphP[i].Gradients.B[k],i);}
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
            construct_gradient(SphP[i].Gradients.Phi,i);
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
            for(k=0;k<NUM_METAL_SPECIES;k++) {construct_gradient(SphP[i].Gradients.Metallicity[k],i);}
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
            for(k=0;k<N_RT_FREQ_BINS;k++) {construct_gradient(GasGradDataPasser[i].Gradients_Rad_E_gamma[k],i);}
#endif

            /* now the gradients are calculated: below are simply useful operations on the results */
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
            /* this is here because for the models of BH growth and self-shielding of stars, we need to calculate GradRho: we don't bother doing it in density.c if we're already calculating it here! but note, this is the -un-limited- gradient here */
            for(k=0;k<3;k++) {P[i].GradRho[k] = SphP[i].Gradients.Density[k];}
#endif

#if defined(TURB_DRIVING) || defined(OUTPUT_VORTICITY)
            SphP[i].Vorticity[0] = SphP[i].Gradients.Velocity[1][2] - SphP[i].Gradients.Velocity[2][1];
            SphP[i].Vorticity[1] = SphP[i].Gradients.Velocity[2][0] - SphP[i].Gradients.Velocity[0][2];
            SphP[i].Vorticity[2] = SphP[i].Gradients.Velocity[0][1] - SphP[i].Gradients.Velocity[1][0];
#endif

#ifdef SPH_TP12_ARTIFICIAL_RESISTIVITY
            /* use the magnitude of the B-field gradients relative to kernel length to calculate artificial resistivity */
            double GradBMag=0.0;
            double BMag=0.0;
            for(k=0;k<3;k++)
            {
                for(j=0;j<3;j++)
                {
                    GradBMag += SphP[i].Gradients.B[k][j]*SphP[i].Gradients.B[k][j];
                }
                BMag += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);
            }
            SphP[i].Balpha = PPP[i].Hsml * sqrt(GradBMag/(BMag+1.0e-33));
            SphP[i].Balpha = DMIN(SphP[i].Balpha, 0.1 * All.ArtMagDispConst);
            SphP[i].Balpha = DMAX(SphP[i].Balpha, 0.005);
#endif


#ifdef HYDRO_SPH

#ifdef MAGNETIC
            if(SphP[i].Density > 0)
            {
                for(k=0;k<3;k++) SphP[i].DtB[k] *= PPP[i].DhsmlNgbFactor * P[i].Mass / (SphP[i].Density * SphP[i].Density) / All.cf_atime; // induction equation (convert from Bcode*vcode/rcode to Bphy/tphys) //
#ifdef DIVBCLEANING_DEDNER
                /* full correct form of D(phi)/Dt = -ch*ch*div.dot.B - phi/tau - (1/2)*phi*div.dot.v */
                /* PFH: here's the div.dot.B term: make sure div.dot.B def'n matches appropriate grad_phi conjugate pair: recommend direct diff div.dot.B */
                SphP[i].divB *= PPP[i].DhsmlNgbFactor * P[i].Mass / (SphP[i].Density * SphP[i].Density);
                if((!isnan(SphP[i].divB))&&(PPP[i].Hsml>0)&&(SphP[i].divB!=0)&&(SphP[i].Density>0))
                {
                    double tmp_ded = 0.5 * SphP[i].MaxSignalVel * All.cf_afac3; // has units of v_physical now
                    /* do a check to make sure divB isn't something wildly divergent (owing to particles being too close) */
                    double b2_max = 0.0;
                    for(k=0;k<3;k++) {b2_max += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);}
                    b2_max = 100.0 * fabs( sqrt(b2_max) * All.cf_a2inv * P[i].Mass / (SphP[i].Density*All.cf_a3inv) * 1.0 / (PPP[i].Hsml*All.cf_atime) );
                    if(fabs(SphP[i].divB) > b2_max) {SphP[i].divB *= b2_max / fabs(SphP[i].divB);}
                    /* ok now can apply this to get the growth rate of phi */
                    // SphP[i].DtPhi = -tmp_ded * tmp_ded * All.DivBcleanHyperbolicSigma * SphP[i].divB;
                    SphP[i].DtPhi = -tmp_ded * tmp_ded * All.DivBcleanHyperbolicSigma * SphP[i].divB * SphP[i].Density*All.cf_a3inv; // mass-based phi-flux
                    // phiphi above now has units of [Bcode]*[vcode]^2/[rcode]=(Bcode*vcode)*vcode/rcode; needs to have units of [Phicode]*[vcode]/[rcode]
                    // [PhiGrad]=[Phicode]/[rcode] = [DtB] = [Bcode]*[vcode]/[rcode] IFF [Phicode]=[Bcode]*[vcode]; this also makes the above self-consistent //
                    // (implicitly, this gives the correct evolution in comoving, adiabatic coordinates where the sound speed is the relevant speed at which
                    //   the 'damping wave' propagates. another choice (provided everything else is self-consistent) is fine, it just makes different assumptions
                    //   about the relevant 'desired' timescale for damping wave propagation in the expanding box) //
                } else {
                    SphP[i].DtPhi=0; SphP[i].divB=0; for(k=0;k<3;k++) {SphP[i].DtB[k] = 0;}
                }
                SphP[i].divB = 0.0; // now we re-zero it, since a -different- divB definition must be used in hydro to subtract the tensile terms */
#endif
            } else {
                for(k=0;k<3;k++) SphP[i].DtB[k] = 0;
#ifdef DIVBCLEANING_DEDNER
                SphP[i].divB = 0; SphP[i].DtPhi = 0;
#endif
            }
#endif


#ifdef SPHAV_CD10_VISCOSITY_SWITCH
            SphP[i].alpha_limiter /= SphP[i].Density;
            NV_dt =  GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); // physical
            NV_dummy = fabs(1.0 * pow(1.0 - SphP[i].alpha_limiter,4.0) * SphP[i].NV_DivVel); // NV_ quantities are in physical units
            NV_limiter = NV_dummy*NV_dummy / (NV_dummy*NV_dummy + SphP[i].NV_trSSt);
            NV_A = DMAX(-SphP[i].NV_dt_DivVel, 0.0);
            divVel_physical = SphP[i].NV_DivVel;

            // add a simple limiter here: alpha_loc is 'prepped' but only switches on when the divergence goes negative: want to add hubble flow here //
            if(All.ComovingIntegrationOn) {divVel_physical += 3*All.cf_hubble_a;} // hubble-flow correction added
            if(divVel_physical>=0.0) {NV_A = 0.0;}

            h_eff = Get_Particle_Size(i) * All.cf_atime / 0.5; // 'default' parameter choices are scaled for a cubic spline //
            cs_nv = Get_Gas_effective_soundspeed_i(i) * All.cf_afac3; // converts to physical velocity units //
            alphaloc = All.ViscosityAMax * h_eff*h_eff*NV_A / (0.36*cs_nv*cs_nv*(0.05/SPHAV_CD10_VISCOSITY_SWITCH) + h_eff*h_eff*NV_A);
            // 0.25 in front of vsig is the 'noise parameter' that determines the relative amplitude which will trigger the switch:
            //    that choice was quite large (requires approach velocity rate-of-change is super-sonic); better to use c_s (above), and 0.05-0.25 //
            // NV_A is physical 1/(time*time), but Hsml and vsig can be comoving, so need appropriate correction terms above //

            if(SphP[i].alpha < alphaloc) {SphP[i].alpha = alphaloc;}
                else if (SphP[i].alpha > alphaloc) {SphP[i].alpha = alphaloc + (SphP[i].alpha - alphaloc) * exp(-NV_dt * (0.5*fabs(SphP[i].MaxSignalVel)*All.cf_afac3)/(0.5*h_eff) * SPHAV_CD10_VISCOSITY_SWITCH);}

            if(SphP[i].alpha < All.ViscosityAMin) {SphP[i].alpha = All.ViscosityAMin;}
            SphP[i].alpha_limiter = DMAX(NV_limiter,All.ViscosityAMin/SphP[i].alpha);
#else
            /* compute the traditional Balsara limiter (now that we have velocity gradients) */
            double divVel = All.cf_a2inv * fabs(SphP[i].Gradients.Velocity[0][0] + SphP[i].Gradients.Velocity[1][1] + SphP[i].Gradients.Velocity[2][2]);
            if(All.ComovingIntegrationOn) {divVel += 3*All.cf_hubble_a;} // hubble-flow correction added (physical units)
            double CurlVel[3], MagCurl;
            CurlVel[0] = SphP[i].Gradients.Velocity[1][2] - SphP[i].Gradients.Velocity[2][1];
            CurlVel[1] = SphP[i].Gradients.Velocity[2][0] - SphP[i].Gradients.Velocity[0][2];
            CurlVel[2] = SphP[i].Gradients.Velocity[0][1] - SphP[i].Gradients.Velocity[1][0];
            MagCurl = All.cf_a2inv * sqrt(CurlVel[0]*CurlVel[0] + CurlVel[1]*CurlVel[1] + CurlVel[2]*CurlVel[2]);
            double fac_mu = 1 / (All.cf_afac3 * All.cf_atime);
            SphP[i].alpha_limiter = divVel / (divVel + MagCurl + 0.0001 * Get_Gas_effective_soundspeed_i(i) / (Get_Particle_Size(i)) / fac_mu);
#endif
#endif



#if (defined(CONDUCTION_SPITZER) || defined(VISCOSITY_BRAGINSKII) || defined(MHD_NON_IDEAL))
            double ion_frac; ion_frac=1;
#if defined(COOLING) /* get the ionized fraction. NOTE we CANNOT call 'ThermalProperties' or functions like 'Get_Ionized_Fraction' here in gradients.c, as we have not done self-shielding steps yet and most modules will yield unphysical answers! */
            ion_frac = SphP[i].Ne / (1. + 2.*yhelium(i)); /* quick estimator. this is actually what we need for conduction since its the free electrons conducting, and we want number relative to fully-ionized gas */
#endif
#endif


#ifdef CONDUCTION
            {
                SphP[i].Kappa_Conduction = All.ConductionCoeff;
#ifdef CONDUCTION_SPITZER
                /* calculate the thermal conductivities: use the Spitzer formula */
                SphP[i].Kappa_Conduction *= ion_frac * pow(SphP[i].InternalEnergyPred, 2.5);

                /* account for saturation (when the mean free path of electrons is large): estimate whether we're in that limit with the gradients */
                double electron_free_path = All.ElectronFreePathFactor * SphP[i].InternalEnergyPred * SphP[i].InternalEnergyPred / (SphP[i].Density * All.cf_a3inv);
                double du_conduction=0;
                for(k=0;k<3;k++) {du_conduction += SphP[i].Gradients.InternalEnergy[k] * SphP[i].Gradients.InternalEnergy[k];}
                double temp_scale_length = SphP[i].InternalEnergyPred / sqrt(du_conduction) * All.cf_atime;
#ifdef MAGNETIC
                // following Jono Squire's notes, the 'Whistler instability' limits the heat flux at high-beta; Komarov et al., arXiv:1711.11462 (2017) //
                double beta_i=0; for(k=0;k<3;k++) {beta_i += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);}
                beta_i *= All.cf_afac1 / (All.cf_atime * SphP[i].Density * Get_Gas_thermal_soundspeed_i(i)*Get_Gas_thermal_soundspeed_i(i));
                SphP[i].Kappa_Conduction /= (1 + (4.2 + 1./(3.*beta_i)) * electron_free_path / temp_scale_length); // should be in physical units //
#else
                SphP[i].Kappa_Conduction /= (1 + 4.2 * electron_free_path / temp_scale_length); // should be in physical units //
#endif

#ifdef DIFFUSION_OPTIMIZERS
                double cs = Get_Gas_effective_soundspeed_i(i);
#ifdef MAGNETIC
                double vA_2 = 0.0; for(k=0;k<3;k++) {vA_2 += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);}
                vA_2 *= All.cf_afac1 / (All.cf_atime * SphP[i].Density);
                cs = DMIN(1.e4*cs , sqrt(cs*cs+vA_2));
#endif
                cs *= All.cf_afac3;
                SphP[i].Kappa_Conduction = DMIN(SphP[i].Kappa_Conduction , 42.85 * SphP[i].Density*All.cf_a3inv * cs * DMIN(20.*Get_Particle_Size(i)*All.cf_atime , temp_scale_length));
#endif
#endif
            }
#endif



#ifdef VISCOSITY
            {
                SphP[i].Eta_ShearViscosity = All.ShearViscosityCoeff;
                SphP[i].Zeta_BulkViscosity = All.BulkViscosityCoeff;
#ifdef VISCOSITY_BRAGINSKII
                /* calculate the viscosity coefficients: use the Braginskii shear tensor formulation expanded to first order */
                SphP[i].Eta_ShearViscosity *= ion_frac * pow(SphP[i].InternalEnergyPred, 2.5);
                SphP[i].Zeta_BulkViscosity = 0;

                /* again need to account for possible saturation (when the mean free path of ions is large): estimate whether we're in that limit with the gradients */
                double ion_free_path = All.ElectronFreePathFactor * SphP[i].InternalEnergyPred * SphP[i].InternalEnergyPred / (SphP[i].Density * All.cf_a3inv);
                /* need an estimate of the internal energy gradient scale length, which we get by d(P/rho) = P/rho * (dP/P - drho/rho) */
                double dv_magnitude=0, v_magnitude=0;
#ifdef MAGNETIC
                double bhat[3]={0},beta_i=0,bmag=0; for(k=0;k<3;k++) {bhat[k]=Get_Gas_BField(i,k); bmag+=bhat[k]*bhat[k];}
                double double_dot_dv=0; if(bmag>0) {bmag = sqrt(bmag); for(k=0;k<3;k++) {bhat[k]/=bmag;}}
                beta_i = bmag*bmag * All.cf_afac1 / (All.cf_atime * SphP[i].Density * Get_Gas_thermal_soundspeed_i(i)*Get_Gas_thermal_soundspeed_i(i));
#endif
                for(k=0;k<3;k++)
                {
                    for(k1=0;k1<3;k1++)
                    {
                        dv_magnitude += SphP[i].Gradients.Velocity[k][k1]*SphP[i].Gradients.Velocity[k][k1];
#ifdef MAGNETIC
                        double_dot_dv += SphP[i].Gradients.Velocity[k][k1] * bhat[k]*bhat[k1] * All.cf_a2inv; // physical units
#endif
                    }
                    v_magnitude += SphP[i].VelPred[k]*SphP[i].VelPred[k];
                }
                double vel_scale_length = sqrt( v_magnitude / dv_magnitude ) * All.cf_atime;
                SphP[i].Eta_ShearViscosity /= (1 + 4.2 * ion_free_path / vel_scale_length); // should be in physical units //
                /* also limit to saturation magnitude ~ signal_speed / lambda_MFP^2 */
                double cs = Get_Gas_effective_soundspeed_i(i);
#ifdef MAGNETIC
                // following Jono Squire's notes, the mirror and firehose instabilities limit pressure anisotropies [which scale as the viscous term inside the gradient: nu_braginskii*(bhat.bhat:grad.v)] to >-2*P_magnetic and <1*P_magnetic
                double P_effective_visc = SphP[i].Eta_ShearViscosity * double_dot_dv;
                double P_magnetic = 0.5 * (bmag*All.cf_a2inv) * (bmag*All.cf_a2inv);
                if(P_effective_visc < -2.*P_magnetic) {SphP[i].Eta_ShearViscosity = 2.*P_magnetic / fabs(double_dot_dv);}
                if(P_effective_visc > P_magnetic) {SphP[i].Eta_ShearViscosity = P_magnetic / fabs(double_dot_dv);}
                double vA_2 = 0.0; for(k=0;k<3;k++) {vA_2 += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);}
                vA_2 *= All.cf_afac1 / (All.cf_atime * SphP[i].Density);
                cs = DMIN(1.e4*cs , sqrt(cs*cs+vA_2));
#endif
                cs *= All.cf_afac3;
#ifdef DIFFUSION_OPTIMIZERS
                double eta_sat = (SphP[i].Density*All.cf_a3inv) * cs / (ion_free_path * (1 + 4.2 * ion_free_path / vel_scale_length));
                if(eta_sat <= 0) SphP[i].Eta_ShearViscosity=0;
                if(SphP[i].Eta_ShearViscosity>0) {SphP[i].Eta_ShearViscosity = 1. / (1./SphP[i].Eta_ShearViscosity + 1./eta_sat);} // again, all physical units //
                //SphP[i].Eta_ShearViscosity = DMIN(SphP[i].Eta_ShearViscosity , SphP[i].Density*All.cf_a3inv * cs * DMAX(Get_Particle_Size(i)*All.cf_atime , vel_scale_length));
#endif
#endif
            }
#endif



#ifdef MHD_NON_IDEAL
#ifdef COOLING	    
	    if(All.Time > All.TimeBegin) // only try to get self-consistent resistivities after the first timestep, when we have calculated the self-consistent ionization state - on the first timestep default to 0 (=ideal MHD)
#endif	      
            {
                /* calculations below follow Wardle 2007 and Keith & Wardle 2014, for the equation sets */
                double mean_molecular_weight = 2.38; // molecular H2, +He with solar mass fractions and metals
#ifdef COOLING
                double T_eff_atomic = 1.23 * (5./3.-1.) * U_TO_TEMP_UNITS * SphP[i].InternalEnergyPred; /* we'll use this to make a quick approximation to the actual mean molecular weight here */
                double nH_cgs = SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, T_transition=DMIN(8000.,nH_cgs), f_mol=1./(1. + T_eff_atomic*T_eff_atomic/(T_transition*T_transition));
                mean_molecular_weight = 4. / (1. + (3. + 4.*SphP[i].Ne - 2.*f_mol) * HYDROGEN_MASSFRAC);
#endif
                double temperature = mean_molecular_weight * (GAMMA(i)-1.) * U_TO_TEMP_UNITS * SphP[i].InternalEnergyPred; // will use appropriate EOS to estimate temperature
                double zeta_cr = 1.0e-17; // cosmic ray ionization rate (fixed as constant for non-CR runs)
#ifdef COSMIC_RAY_FLUID
                double u_cr=0; for(k=0;k<N_CR_PARTICLE_BINS;k++) {u_cr += SphP[i].CosmicRayEnergyPred[k];}
                zeta_cr = u_cr * 2.2e-6 * ((1. / P[i].Mass * SphP[i].Density * All.cf_a3inv) * (UNIT_PRESSURE_IN_CGS)); // convert to ionization rate
#endif
                double a_grain_micron = 0.1, f_dustgas = 0.01; // effective size of grains that matter at these densities
                double m_ion = 24.3; // Mg dominates ions in dense gas [where this is relevant]; this is ion mass in units of proton mass
#ifdef METALS
		        f_dustgas = 0.5 * P[i].Metallicity[0] * return_dust_to_metals_ratio_vs_solar(i); // appropriate dust-to-metals ratio
#endif
		        // now everything should be fully-determined (given the inputs above and the known properties of the gas) //
                double m_neutral = mean_molecular_weight; // in units of the proton mass
		double ag01 = a_grain_micron/0.1, m_grain = 7.51e9 * ag01*ag01*ag01; // grain mass [internal density =3 g/cm^3]
		double rho = SphP[i].Density*All.cf_a3inv * UNIT_DENSITY_IN_CGS, n_eff = rho / PROTONMASS_CGS; // density in cgs
                // calculate ionization fraction in dense gas; use rate coefficients k to estimate grain charge
                double k0 = 1.95e-4 * ag01*ag01 * sqrt(temperature); // prefactor for rate coefficient for electron-grain collisions
                double ngr_ngas = (m_neutral/m_grain) * f_dustgas; // number of grains per neutral
                double psi_prefac = 167.1 / (ag01 * temperature); // e*e/(a_grain*k_boltzmann*T): Z_grain = psi/psi_prefac where psi is constant determines charge
                double alpha = zeta_cr * psi_prefac / (ngr_ngas*ngr_ngas * k0 * (n_eff/m_neutral)); // coefficient for equation that determines Z_grain
                // psi solves the equation: psi = alpha * (exp[psi] - y/(1+psi)) where y=sqrt(m_ion/m_electron); note the solution for small alpha is independent of m_ion, only large alpha
                //   (where the non-ideal effects are weak, generally) produces a difference: at very high-T, appropriate m_ion should be hydrogen+helium, but in this limit our cooling
                //    routines will already correctly determine the ionization states. so we can safely adopt Mg as our ion of consideration
                double y=sqrt(m_ion*PROTONMASS_CGS/ELECTRONMASS_CGS), psi_0 = 0.5188025-0.804386*log(y), psi=psi_0; // solution for large alpha [>~10]
                if(alpha<0.002) {psi=alpha*(1.-y)/(1.+alpha*(1.+y));} else if(alpha<10.) {psi=psi_0/(1.+0.027/alpha);} // accurate approximation for intermediate values we can use here
                double k_e = k0 * exp(psi); // e-grain collision rate coefficient
                double k_i = k0 * sqrt(ELECTRONMASS_CGS / (m_ion*PROTONMASS_CGS)) * (1 - psi); // i-grain collision rate coefficient
                double n_elec = zeta_cr / (ngr_ngas * k_e); // electron number density
                double n_ion = zeta_cr / (ngr_ngas * k_i); // ion number density
                double Z_grain = psi / psi_prefac; // mean grain charge (note this is signed, will be negative)
#ifdef COOLING
                double mu_eff=2.38, x_elec=DMAX(1e-16, SphP[i].Ne*HYDROGEN_MASSFRAC*mu_eff), R=x_elec*psi_prefac/ngr_ngas; psi_0=-3.787124454911839; n_elec=x_elec*n_eff/mu_eff; // R is essentially the ratio of negative charge in e- to dust: determines which regime we're in to set quantities below
                if(R > 100.) {psi=psi_0;} else if(R < 0.002) {psi=R*(1.-y)/(1.+2.*y*R);} else {psi=psi_0/(1.+pow(R/0.18967,-0.5646));} // simple set of functions to solve for psi, given R above, using the same equations used to determine low-temp ion fractions
                n_ion = n_elec * y * exp(psi)/(1.-psi); Z_grain = psi / psi_prefac; // we can immediately now calculate these from the above
#endif
                // now define more variables we will need below //
                double gizmo2gauss = UNIT_B_IN_GAUSS; // convert to B-field to gauss (units)
                double B_Gauss = 0; for(k=0;k<3;k++) {B_Gauss += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);} // get magnitude of B //
                if(B_Gauss<=0) {B_Gauss=0;} else {B_Gauss = sqrt(B_Gauss) * All.cf_a2inv * gizmo2gauss;} // B-field magnitude in Gauss
                double xe = n_elec / n_eff;
                double xi = n_ion / n_eff;
                double xg = ngr_ngas;
                // get collision rates/cross sections for different species //
                double nu_g = 7.90e-6 * ag01*ag01 * sqrt(temperature/m_neutral) / (m_neutral+m_grain); // Pinto & Galli 2008
                double nu_ei = 51.*xe*pow(temperature,-1.5); // Pandey & Wardle 2008 (e-ion)
                double nu_e = nu_ei + 6.21e-9*pow(temperature/100.,0.65)/m_neutral; // Pinto & Galli 2008 for latter (e-neutral)
                double nu_i = (xe/xi)*nu_ei + 1.57e-9/(m_neutral+m_ion); // // Pandey & Wardle 2008 for former (e-ion), Pinto & Galli 2008 for latter (i-neutral)
                // use the cross sections to determine the hall parameters and conductivities //
                double beta_prefac = ELECTRONCHARGE_CGS * B_Gauss / (PROTONMASS_CGS * C_LIGHT_CGS * n_eff);
                double beta_i = beta_prefac / (m_ion * nu_i); // standard beta factors (Hall parameters)
                double beta_e = beta_prefac / (ELECTRONMASS_CGS/PROTONMASS_CGS * nu_e);
                double beta_g = beta_prefac / (m_grain * nu_g) * fabs(Z_grain);
                double be_inv = 1/(1 + beta_e*beta_e), bi_inv = 1/(1 + beta_i*beta_i), bg_inv = 1/(1 + beta_g*beta_g);
                double sigma_O = xe*beta_e + xi*beta_i + xg*fabs(Z_grain)*beta_g; // ohmic conductivity
                double sigma_H = -xe*be_inv + xi*bi_inv - xg*Z_grain*bg_inv; // hall conductivity
                double sigma_P = xe*beta_e*be_inv + xi*beta_i*bi_inv + xg*fabs(Z_grain)*beta_g*bg_inv; // pedersen conductivity
                // now we can finally calculate the diffusivities //
                double eta_prefac = B_Gauss * C_LIGHT_CGS / (4 * M_PI * ELECTRONCHARGE_CGS * n_eff );
                double eta_O = eta_prefac / sigma_O;
                double sigma_perp2 = sigma_H*sigma_H + sigma_P*sigma_P;
                double eta_H = eta_prefac * sigma_H / sigma_perp2;
                double eta_A = eta_prefac * (sigma_P/sigma_perp2 - 1/sigma_O);
                eta_O = DMAX(0,eta_O); eta_H = DMAX(0,eta_H); eta_A = DMAX(0,eta_A); // check against unphysical negative diffusivities
                // convert units to code units
                double units_cgs_to_code = UNIT_TIME_IN_CGS / (UNIT_LENGTH_IN_CGS * UNIT_LENGTH_IN_CGS); // convert coefficients (L^2/t) to code units [physical]
                double eta_ohmic = eta_O*units_cgs_to_code, eta_hall = eta_H*units_cgs_to_code, eta_ad = eta_A*units_cgs_to_code;

                SphP[i].Eta_MHD_OhmicResistivity_Coeff = eta_ohmic;     /*!< Ohmic resistivity coefficient [physical units of L^2/t] */
                SphP[i].Eta_MHD_HallEffect_Coeff = eta_hall;            /*!< Hall effect coefficient [physical units of L^2/t] */
                SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = eta_ad;      /*!< Hall effect coefficient [physical units of L^2/t] */
            }
#ifdef COOLING
	    else {SphP[i].Eta_MHD_OhmicResistivity_Coeff = SphP[i].Eta_MHD_HallEffect_Coeff = SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;} // =0 on the first timestep, since we don't know the ionization yet
#endif	    
#endif


#ifdef RADTRANSFER
            {
                int k_freq; for(k_freq = 0; k_freq < N_RT_FREQ_BINS; k_freq++)
                {
                    /* calculate the opacity */
                    SphP[i].Rad_Kappa[k_freq] = rt_kappa(i,k_freq); // physical units //
#if defined(RT_FLUXLIMITER) && defined(RT_COMPGRAD_EDDINGTON_TENSOR)
                    /* compute the flux-limiter for radiation transport: also convenient here to compute the relevant opacities for all particles */
                    double lambda = 1;
                    if(SphP[i].Rad_E_gamma_Pred[k_freq] > 0) /* can compute gradient length scale */
                    {
                        double R_ET = sqrt(SphP[i].Gradients.Rad_E_gamma_ET[k_freq][0] * SphP[i].Gradients.Rad_E_gamma_ET[k_freq][0] + SphP[i].Gradients.Rad_E_gamma_ET[k_freq][1] * SphP[i].Gradients.Rad_E_gamma_ET[k_freq][1] + SphP[i].Gradients.Rad_E_gamma_ET[k_freq][2] * SphP[i].Gradients.Rad_E_gamma_ET[k_freq][2]) / (MIN_REAL_NUMBER + SphP[i].Rad_E_gamma_Pred[k_freq] * SphP[i].Density/(MIN_REAL_NUMBER+P[i].Mass));
                        R_ET = 3.*DMAX(R_ET , 1.e-6/Get_Particle_Size(i)) / (1.e-55 + All.cf_atime*SphP[i].Rad_Kappa[k_freq]*(SphP[i].Density*All.cf_a3inv)); // limit to be > 0, divide by kappa-rho to get desired dimensionless ratio
                        lambda = DMIN(1., DMAX( 3.*(2. + R_ET) / (6. + 3.*R_ET + R_ET*R_ET), MIN_REAL_NUMBER )); // slope-limiter
#ifdef RT_OTVET         /* note that the OTVET eddington tensor is close to the correct value for the optically-thin limit. for the diffusion limit
                            it may be incorrect. we can therefore interpolate using an M1-like relation below, based on the gradients above (used
                            to determine which limit we are actually in: ratio f=|flux|/(c_eff*Energy_density_rad): f<<1 = diffusion limit, f~1 = free-streaming limit: this is our slope-limiter above */
                        double chi=DMAX(1./3.,DMIN(1.,(3.+4.*lambda*lambda)/(5.+2.*sqrt(4.-3.*lambda*lambda)))), chifac_iso=3.*(1-chi)/2., chifac_ot=(3.*chi-1.)/2.;
                        for(k=0;k<3;k++) {SphP[i].Gradients.Rad_E_gamma_ET[k_freq][k] = chifac_ot*SphP[i].Gradients.Rad_E_gamma_ET[k_freq][k] + chifac_iso/3.*GasGradDataPasser[i].Gradients_Rad_E_gamma[k_freq][k];}
#ifdef RT_DIFFUSION_CG
                        for(k=0;k<6;k++) {SphP[i].ET[k_freq][k] *= chifac_ot; if(k<3) {SphP[i].ET[k_freq][k] += chifac_iso/3.;}} // diagonal components // (this only makes sense if ET is freq-dependent) [note this will cause instability in the explicit methods; only use for CG where ET is explicitly called and this is done only on global timesteps]
#endif
#endif // ifdef otvet
                    }
                    SphP[i].Rad_Flux_Limiter[k_freq] = lambda;
#endif // ifdef fluxlimiter

#if defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_OTVET)
                    /* set the output gradient grad.(D*Prad) = D.(grad Prad), i.e. move the tensor outside the gradient. while not strictly self-consistent, this is more stable, and correct at the level of the ad-hoc M1 or OTVET closure, b/c otherwise M1 introduces unphysical behaviors from the gradients of the tensor where the ad-hoc closure relation causes changes to D [see e.g. Hopkins 'Anisotropic Diffusion in Mesh-Free Numerical Magnetohydrodynamics' Fig 8 and associated discussion]. Also works for FLD trivially. */
                    eddington_tensor_dot_vector(SphP[i].ET[k_freq],GasGradDataPasser[i].Gradients_Rad_E_gamma[k_freq],SphP[i].Gradients.Rad_E_gamma_ET[k_freq]);
#endif
#if defined(GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION) /* yet another hack for this problem to get the boundaries to play nicely once dust evacuated -- this is a bit redundant with other hacks, but here for safety */
                    if(SphP[i].Interpolated_Opacity[0] < 1.e-3 * All.Dust_to_Gas_Mass_Ratio*0.75*All.Grain_Q_at_MaxGrainSize/(All.Grain_Internal_Density*All.Grain_Size_Max)) {double gmax=-1; if(P[i].GravAccel[GRAV_DIRECTION_RDI] < gmax) {P[i].GravAccel[GRAV_DIRECTION_RDI]=gmax;}}
#endif
                }
            }
#endif // ifdef radtransfer

#if defined(EOS_ELASTIC) // update time-derivative of stress tensor (needs to be done before slope-limiting to use full velocity gradient information) //
            elastic_body_update_driftkick(i,1.,2);
#endif

            /* finally, we need to apply a sensible slope limiter to the gradients, to prevent overshooting */
            double stol = 0.0, stol_tmp, stol_diffusion; stol_diffusion = 0.1; stol_tmp = stol;
            double h_lim = PPP[i].Hsml, d_max = DMAX(PPP[i].Hsml,GasGradDataPasser[i].MaxDistance); h_lim = d_max;
            /* fraction of H at which maximum reconstruction is allowed (=0.5 for 'standard'); for pure hydro we can
             be a little more aggresive and the equations are still stable (but this is as far as you want to push it) */
            double a_limiter = 0.25; if(SphP[i].ConditionNumber>100) a_limiter=DMIN(0.5, 0.25 + 0.25 * (SphP[i].ConditionNumber-100)/100);
#if defined(SELFGRAVITY_OFF) && (!defined(MAGNETIC) && !defined(GALSF))
            h_lim=PPP[i].Hsml; stol=0.1;
#endif
#if (SLOPE_LIMITER_TOLERANCE == 2)
            h_lim = PPP[i].Hsml; a_limiter *= 0.5; stol = 0.125;
#endif
#if (SLOPE_LIMITER_TOLERANCE == 0)
            a_limiter *= 2.0; stol = 0.0;
#endif

#if (SINGLE_STAR_SINK_FORMATION & 4)
            SphP[i].Density_Relative_Maximum_in_Kernel = GasGradDataPasser[i].Maxima.Density;
#endif
            local_slopelimiter(SphP[i].Gradients.Density,GasGradDataPasser[i].Maxima.Density,GasGradDataPasser[i].Minima.Density,a_limiter,h_lim,0, 1,d_max,SphP[i].Density);
            int pressure_is_positive_definite = 1;
#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC)
            pressure_is_positive_definite = 0; /* some physics allow negative pressures - account for that here */
#endif
            local_slopelimiter(SphP[i].Gradients.Pressure,GasGradDataPasser[i].Maxima.Pressure,GasGradDataPasser[i].Minima.Pressure,a_limiter,h_lim,stol, pressure_is_positive_definite,d_max,SphP[i].Pressure);
            stol_tmp = stol;
#if defined(VISCOSITY)
            stol_tmp = DMAX(stol,stol_diffusion);
#endif
#ifdef TURB_DIFF_DYNAMIC
            for (k1=0;k1<3;k1++) {local_slopelimiter(GasGradDataPasser[i].GradVelocity_bar[k1], GasGradDataPasser[i].Maxima.Velocity_bar[k1], GasGradDataPasser[i].Minima.Velocity_bar[k1], a_limiter, h_lim, stol, 0,0,0);}
#endif
            for(k1=0;k1<3;k1++) {local_slopelimiter(SphP[i].Gradients.Velocity[k1],GasGradDataPasser[i].Maxima.Velocity[k1],GasGradDataPasser[i].Minima.Velocity[k1],a_limiter,h_lim,stol_tmp, 0,0,0);}
#ifdef DOGRAD_INTERNAL_ENERGY
            stol_tmp = stol;
#if defined(CONDUCTION)
            stol_tmp = DMAX(stol,stol_diffusion);
#endif
            local_slopelimiter(SphP[i].Gradients.InternalEnergy,GasGradDataPasser[i].Maxima.InternalEnergy,GasGradDataPasser[i].Minima.InternalEnergy,a_limiter,h_lim,stol_tmp, 1,d_max,SphP[i].InternalEnergyPred);
#endif
#ifdef COSMIC_RAY_FLUID
            stol_tmp = stol;
#ifndef CRFLUID_M1
            for(k1=0;k1<N_CR_PARTICLE_BINS;k1++)
            {
                local_slopelimiter(SphP[i].Gradients.CosmicRayPressure[k1],GasGradDataPasser[i].Maxima.CosmicRayPressure[k1],GasGradDataPasser[i].Minima.CosmicRayPressure[k1],DMAX(1.,a_limiter),h_lim,0., 1,d_max,Get_Gas_CosmicRayPressure(i,k1));
                if((GasGradDataPasser[i].Maxima.CosmicRayPressure[k1]==0)||(GasGradDataPasser[i].Minima.CosmicRayPressure[k1]==0)) {is_particle_local_extremum[k1] = 1;}
            }
#endif
#endif
#ifdef DOGRAD_SOUNDSPEED
            local_slopelimiter(SphP[i].Gradients.SoundSpeed,GasGradDataPasser[i].Maxima.SoundSpeed,GasGradDataPasser[i].Minima.SoundSpeed,a_limiter,h_lim,stol, 1,d_max,Get_Gas_effective_soundspeed_i(i));
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
            for(k1=0;k1<NUM_METAL_SPECIES;k1++) {local_slopelimiter(SphP[i].Gradients.Metallicity[k1],GasGradDataPasser[i].Maxima.Metallicity[k1],GasGradDataPasser[i].Minima.Metallicity[k1],a_limiter,h_lim,DMAX(stol,stol_diffusion), 1,d_max,P[i].Metallicity[k1]);}
#endif
#if defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_EVOLVE_FLUX)
            for(k1=0;k1<N_RT_FREQ_BINS;k1++)
            {
#if 0 /* stricter limiter here */
                local_slopelimiter(SphP[i].Gradients.Rad_E_gamma_ET[k1],GasGradDataPasser[i].Maxima.Rad_E_gamma[k1],GasGradDataPasser[i].Minima.Rad_E_gamma[k1],DMAX(1.,a_limiter),DMAX(h_lim,d_max),0, 1,d_max,SphP[i].Rad_E_gamma_Pred[k1]*SphP[i].Density/P[i].Mass);
                local_slopelimiter(GasGradDataPasser[i].Gradients_Rad_E_gamma[k1],GasGradDataPasser[i].Maxima.Rad_E_gamma[k1],GasGradDataPasser[i].Minima.Rad_E_gamma[k1],DMAX(1.,a_limiter),DMAX(h_lim,d_max),0, 1,d_max,SphP[i].Rad_E_gamma_Pred[k1]*SphP[i].Density/P[i].Mass);
#else /* weaker limiter here [choose with current tests] */
                local_slopelimiter(SphP[i].Gradients.Rad_E_gamma_ET[k1],GasGradDataPasser[i].Maxima.Rad_E_gamma[k1],GasGradDataPasser[i].Minima.Rad_E_gamma[k1],a_limiter,h_lim,stol, 1,d_max,SphP[i].Rad_E_gamma_Pred[k1]*SphP[i].Density/P[i].Mass);
                local_slopelimiter(GasGradDataPasser[i].Gradients_Rad_E_gamma[k1],GasGradDataPasser[i].Maxima.Rad_E_gamma[k1],GasGradDataPasser[i].Minima.Rad_E_gamma[k1],a_limiter,h_lim,DMAX(stol,stol_diffusion), 1,d_max,SphP[i].Rad_E_gamma_Pred[k1]*SphP[i].Density/P[i].Mass);
#endif
            }
#endif
#ifdef MAGNETIC
#ifndef MHD_CONSTRAINED_GRADIENT
            double v_tmp = P[i].Mass / SphP[i].Density;
            double tmp_d = sqrt(1.0e-37 + (2. * All.cf_atime/All.cf_afac1 * SphP[i].Pressure*v_tmp*v_tmp) +
                                SphP[i].BPred[0]*SphP[i].BPred[0]+SphP[i].BPred[1]*SphP[i].BPred[1]+SphP[i].BPred[2]*SphP[i].BPred[2]);
            double q = fabs(SphP[i].divB) * PPP[i].Hsml / tmp_d, alim2 = a_limiter * (1. + q*q); if(alim2 > 0.5) alim2=0.5;
            stol_tmp = stol;
#ifdef MHD_NON_IDEAL
            stol_tmp = DMAX(stol,stol_diffusion);
#endif
            for(k1=0;k1<3;k1++) {local_slopelimiter(SphP[i].Gradients.B[k1],GasGradDataPasser[i].Maxima.B[k1],GasGradDataPasser[i].Minima.B[k1],alim2,h_lim,stol_tmp, 0,0,0);}
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
            local_slopelimiter(SphP[i].Gradients.Phi,GasGradDataPasser[i].Maxima.Phi,GasGradDataPasser[i].Minima.Phi,a_limiter,h_lim,stol, 0,0,0);
#endif
#endif




#ifdef TURB_DIFFUSION
            {
                /* estimate local turbulent diffusion coefficient from velocity gradients using Smagorinsky mixing model:
                    we do this after slope-limiting to prevent the estimated velocity gradients from being unphysically large */
                double h_turb = Get_Particle_Size(i) * All.cf_atime; // physical
                if(h_turb > 0)
                {
                    // overall normalization //
                    double C_Smagorinsky_Lilly = 0.15; // this is the standard Smagorinsky-Lilly constant, calculated from Kolmogorov theory: should be 0.1-0.2 //
                    double turb_prefactor = 0.25 * All.TurbDiffusion_Coefficient * C_Smagorinsky_Lilly*C_Smagorinsky_Lilly * sqrt(2.0);
                    // then scale with inter-particle spacing //
                    turb_prefactor *= h_turb*h_turb;
                    // calculate frobenius norm of symmetric shear velocity gradient tensor //
                    double shear_factor = sqrt((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) *
                                                        (SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) +
                                                        (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) *
                                                        (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) +
                                                        (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2]) *
                                                        (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
                                               (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] +
                                                         SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] +
                                                         SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) -
                                                        (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] +
                                                         SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] +
                                                         SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2])));
                    // slope-limit and convert to physical units //
                    double shearfac_max = 0.5 * sqrt(SphP[i].VelPred[0]*SphP[i].VelPred[0]+SphP[i].VelPred[1]*SphP[i].VelPred[1]+SphP[i].VelPred[2]*SphP[i].VelPred[2]) / h_turb;
                    shear_factor = DMIN(shear_factor , shearfac_max * All.cf_atime) * All.cf_a2inv; // physical
#ifdef TURB_DIFF_DYNAMIC
                    int u, v; double trace = 0;
                    shearfac_max = 0.5 * sqrt(SphP[i].Velocity_bar[0] * SphP[i].Velocity_bar[0] + SphP[i].Velocity_bar[1] * SphP[i].Velocity_bar[1]+SphP[i].Velocity_bar[2] * SphP[i].Velocity_bar[2]) * All.cf_atime / h_turb;
                    for (u = 0; u < 3; u++) {
                        for (v = 0; v < 3; v++) {
                            SphP[i].VelShear_bar[u][v] = 0.5 * (GasGradDataPasser[i].GradVelocity_bar[u][v] + GasGradDataPasser[i].GradVelocity_bar[v][u]);
                            if (SphP[i].VelShear_bar[u][v] < 0) {SphP[i].VelShear_bar[u][v] = DMAX(SphP[i].VelShear_bar[u][v], -shearfac_max);}
                            else {SphP[i].VelShear_bar[u][v] = DMIN(SphP[i].VelShear_bar[u][v], shearfac_max);}
                            if (u == v) {trace += SphP[i].VelShear_bar[u][u];}}}
                    /* If it was already trace-free, don't zero out the diagonal components */
                    if (trace != 0 && NUMDIMS > 1) {for (u = 0; u < NUMDIMS; u++) {SphP[i].VelShear_bar[u][u] -= 1.0 / NUMDIMS * trace;}}
                    for (u = 0; u < 3; u++) { /* Don't want to recalculate these a bunch later on, so save them */
                        SphP[i].Velocity_hat[u] *= All.TurbDynamicDiffSmoothing;
                        for (v = 0; v < 3; v++) {SphP[i].MagShear_bar += SphP[i].VelShear_bar[u][v] * SphP[i].VelShear_bar[u][v];}}
                    SphP[i].MagShear = sqrt(2.0) * shear_factor / All.cf_a2inv; // Don't want this physical
                    SphP[i].MagShear_bar = DMIN(sqrt(2.0 * SphP[i].MagShear_bar), shearfac_max); turb_prefactor /= 0.25;
#endif

                    // ok, combine to get the diffusion coefficient //
                    SphP[i].TD_DiffCoeff = turb_prefactor * shear_factor; // physical
                } else {
                    SphP[i].TD_DiffCoeff = 0;
                }
#ifdef TURB_DIFF_ENERGY
                SphP[i].Kappa_Conduction = All.ConductionCoeff * SphP[i].TD_DiffCoeff * SphP[i].Density * All.cf_a3inv; // physical
#endif
#ifdef TURB_DIFF_VELOCITY
                SphP[i].Eta_ShearViscosity = All.ShearViscosityCoeff * SphP[i].TD_DiffCoeff * SphP[i].Density * All.cf_a3inv; // physical
                SphP[i].Zeta_BulkViscosity = All.BulkViscosityCoeff * SphP[i].TD_DiffCoeff * SphP[i].Density * All.cf_a3inv; // physical
#endif
            }
#endif



#if defined(COSMIC_RAY_FLUID) && !defined(CRFLUID_EVOLVE_SCATTERINGWAVES) /* note that because of the way this depends on the gradient scale-length, we should calculate it -after- the slope-limiters are applied */
            for(k=0;k<N_CR_PARTICLE_BINS;k++) {SphP[i].CosmicRayDiffusionCoeff[k]=0;}
            if(SphP[i].Density > 0 && P[i].Mass > 0) {CalculateAndAssign_CosmicRay_DiffusionAndStreamingCoefficients(i);}/* only assign diffusivities to 'valid' gas particles */
#ifndef CRFLUID_M1
            for(k=0;k<N_CR_PARTICLE_BINS;k++) {if(is_particle_local_extremum[k]==1) {SphP[i].CosmicRayDiffusionCoeff[k] *= -1;}} // negative here codes for local extrema
#endif
#endif


#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
            /* if the mesh motion is specified to be glass-generating, this is where we apply the appropriate mesh velocity */
            if(All.Time > 0)
            {
                double cs_invelunits = Get_Gas_effective_soundspeed_i(i) * All.cf_afac3 * All.cf_atime; // soundspeed, converted to units of code velocity
                double L_i_code = Get_Particle_Size(i); // particle effective size (in code units)
                double dvel[3]={0}, velnorm=0; for(k=0;k<3;k++) {dvel[k] = L_i_code*L_i_code*GasGradDataPasser[i].GlassAcc[k]; velnorm += dvel[k]*dvel[k];} // calculate quantities to use for glass
                double dtx = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); // need timestep for limiter below
                if(velnorm > 0 && dtx > 0)
                {
                    velnorm = sqrt(velnorm); // normalization for glass 'force'
                    double v00 = 0.5 * DMIN(cs_invelunits*(0.5*velnorm) , All.CourantFac*(L_i_code/dtx)/All.cf_a2inv); // limit added velocity of mesh-generating point to Courant factor
                    for(k=0;k<3;k++) {SphP[i].ParticleVel[k] += v00 * (dvel[k]/velnorm);} // actually add the correction velocity to the mesh velocity
                }
            }
#endif

            
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS)
            {
                double S_in[9], S1[3], Bvec[3]={0}, S2[3][3], S2_inv[3][3], S3_proj[3][3], S3_proj_inv[3][3], CNum_S2, CNum_S3; int k1,k2;
                for(k1=0;k1<HYDRO_TENSOR_FACE_CORRECTIONS_NUMBER_MOMWTS;k1++) {S_in[k1]=GasGradDataPasser[i].first_second_moment_weights[k1];} // we have this saved from the step below
                S1[0]=S_in[0]; S1[1]=S_in[1]; S1[2]=S_in[2]; // these are the vector 1st-moment sums SUM[wt*dp]
                // now assign the matrix based on the 6 independent 2nd-moment sums SUM[wt*dp*dp]
                S2[0][0]=S_in[3]; S2[1][0]=S2[0][1]=S_in[4]; S2[2][0]=S2[0][2]=S_in[5];
                S2[1][1]=S_in[6]; S2[2][1]=S2[1][2]=S_in[7]; S2[2][2]=S_in[8];
                CNum_S2 = matrix_invert_ndims(S2, S2_inv); // invert the S2 matrix [momment-2 matrix]
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS_NGBITER)
                double S2B[3][3], S2B_inv[3][3];
                S2B[0][0]=S_in[3+6]; S2B[1][0]=S2B[0][1]=S_in[4+6]; S2B[2][0]=S2B[0][2]=S_in[5+6];
                S2B[1][1]=S_in[6+6]; S2B[2][1]=S2B[1][2]=S_in[7+6]; S2B[2][2]=S_in[8+6];
                CNum_S2 = matrix_invert_ndims(S2B, S2B_inv); // invert the S2 matrix [moment-2 matrix]
                for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Bvec[k1] += -S2B_inv[k1][k2] * S1[k2];}} // calculate the B-coefficient matrix
#else
                for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Bvec[k1] += -S2_inv[k1][k2] * S1[k2];}} // calculate the B-coefficient matrix
#endif
                double m3wt[10]; for(k=0;k<10;k++) {m3wt[k]=GasGradDataPasser[i].third_moment_weights[k];} // assign these to a shorter variable
                double xxx=m3wt[0], xxy=m3wt[1], xxz=m3wt[2], yyy=m3wt[3], xyy=m3wt[4], yyz=m3wt[5], zzz=m3wt[6], xzz=m3wt[7], yzz=m3wt[8], xyz=m3wt[9]; // assign them intuitive names here, that's all
                S3_proj[0][0]= xxx*Bvec[0] + xxy*Bvec[1] + xxz*Bvec[2]; // now construct the components of this projection matrix
                S3_proj[1][0]=S3_proj[0][1]= xxy*Bvec[0] + xyy*Bvec[1] + xyz*Bvec[2]; //=xyx*Bvec[0] + xyy*Bvec[1] + xyz*Bvec[2];
                S3_proj[2][0]=S3_proj[0][2]= xxz*Bvec[0] + xyz*Bvec[1] + xzz*Bvec[2]; //=xzx*Bvec[0] + xzy*Bvec[1] + xzz*Bvec[2];
                S3_proj[1][1]= xyy*Bvec[0] + yyy*Bvec[1] + yyz*Bvec[2]; //=yyx*Bvec[0] + yyy*Bvec[1] + yyz*Bvec[2];
                S3_proj[2][1]=S3_proj[1][2]= xyz*Bvec[0] + yyz*Bvec[1] + yzz*Bvec[2]; //=yzx*Bvec[0] + yzy*Bvec[1] + yzz*Bvec[2];
                S3_proj[2][2]= xzz*Bvec[0] + yzz*Bvec[1] + zzz*Bvec[2]; //=zzx*Bvec[0] + zzy*Bvec[1] + zzz*Bvec[2];
                for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {S3_proj[k1][k2] = -S3_proj[k1][k2] + S2[k1][k2];}} // add the S2 term that the term above corrects
                CNum_S3 = matrix_invert_ndims(S3_proj, S3_proj_inv); // invert the corrected S2 matrix after correction from the projected S3 matrix down to a modified moment-2 matrix
                for(k1=0;k1<3;k1++) {SphP[i].Tensor_MFM_Face_Corrections[k1] = Bvec[k1];} // load everything into the saved variables we'll use in the hydro step
                SphP[i].Tensor_MFM_Face_Corrections[3] = S3_proj_inv[0][0]; SphP[i].Tensor_MFM_Face_Corrections[4] = S3_proj_inv[0][1]; SphP[i].Tensor_MFM_Face_Corrections[5] = S3_proj_inv[0][2];
                SphP[i].Tensor_MFM_Face_Corrections[6] = S3_proj_inv[1][1]; SphP[i].Tensor_MFM_Face_Corrections[7] = S3_proj_inv[1][2]; SphP[i].Tensor_MFM_Face_Corrections[8] = S3_proj_inv[2][2];
            }
#endif

            
#if defined(KERNEL_CRK_FACES)
            {
                // ok first, load the data from the passer structure into more convenient form //
                double m0, dm0[3], m1[3], dm1[3][3], m2[3][3], m2i[3][3], dm2[3][3][3], Cnum_m2;
                m0 = GasGradDataPasser[i].m0;
                int k_x, k_y;
                for(k=0;k<3;k++)
                {
                    dm0[k] = GasGradDataPasser[i].dm0[k];
                    m1[k] = GasGradDataPasser[i].m1[k];
                    for(k_x=0;k_x<3;k_x++)
                    {
                        dm1[k][k_x] = GasGradDataPasser[i].dm1[k][k_x];
                        int k_tmp;
                        if((k==0)&&(k_x==0)) {k_tmp=0;}
                        if((k==1)&&(k_x==1)) {k_tmp=1;}
                        if((k==2)&&(k_x==2)) {k_tmp=2;}
                        if((k==0)&&(k_x==1)) {k_tmp=3;}
                        if((k==1)&&(k_x==0)) {k_tmp=3;}
                        if((k==0)&&(k_x==2)) {k_tmp=4;}
                        if((k==2)&&(k_x==0)) {k_tmp=4;}
                        if((k==1)&&(k_x==2)) {k_tmp=5;}
                        if((k==2)&&(k_x==1)) {k_tmp=5;}
                        m2[k][k_x] = GasGradDataPasser[i].m2[k_tmp]; m2i[k][k_x] = 0;
                        for(k_y=0;k_y<3;k_y++) {dm2[k][k_x][k_y] = GasGradDataPasser[i].dm2[k_tmp][k_y];}
                    }
                }
                // transform from 'mu' variables to 'm' variables for derivatives:
                for(k=0;k<3;k++) {dm1[k][k] += m0;}
                for(k=0;k<3;k++) {for(k_x=0;k_x<3;k_x++) {dm2[k][k_x][k_x] += m1[k]; dm2[k_x][k][k_x] += m1[k];}}
                Cnum_m2 = matrix_invert_ndims(m2, m2i); // now, invert the m2 matrix into the form we will actually use
                // now start constructing the actual derivatives we need //
                double A = 0, B[3] = {0}, Bdotm1 = 0, dB[3][3]={{0}}, dA[3]={0};
                for(k=0;k<3;k++)
                {
                    for(k_x=0;k_x<3;k_x++) {B[k] += -m2i[k][k_x] * m1[k_x];}
                    Bdotm1 += B[k] * m1[k];
                }
                A = 1. / (m0 + Bdotm1);

                // now the painful part (likely to be errors) -- construct the complicated tensor derivatives contracting all components //
                double minus_m2i_dm1_dotm1[3]={0}, contracted_twotensor[3][3]={{0}}, contracted_twotensor_x[3][3]={{0}}, contracted_twotensor_dotm1[3]={0};
                int k_alpha, k_gamma, k_beta, k_delta;
                for(k_gamma=0; k_gamma<3; k_gamma++)
                {
                    for(k_alpha=0;k_alpha<3;k_alpha++)
                    {
                        for(k_beta=0;k_beta<3;k_beta++)
                        {
                            contracted_twotensor[k_beta][k_gamma] = 0;
                            for(k_delta=0;k_delta<3;k_delta++) {contracted_twotensor[k_beta][k_gamma] += dm2[k_beta][k_delta][k_gamma] * B[k_delta];}
                            contracted_twotensor_x[k_alpha][k_gamma] += dm2[k_alpha][k_beta][k_gamma] * B[k_beta];
                            dB[k_alpha][k_gamma] += -m2i[k_alpha][k_beta] * (dm1[k_beta][k_gamma] + contracted_twotensor[k_beta][k_gamma]);
                        }
                        minus_m2i_dm1_dotm1[k_gamma] += 2.*B[k_alpha]*dm1[k_alpha][k_gamma];
                        contracted_twotensor_dotm1[k_gamma] += B[k_alpha]*contracted_twotensor_x[k_alpha][k_gamma];
                    }
                    dA[k_gamma] = -A*A * (dm0[k_gamma] + minus_m2i_dm1_dotm1[k_gamma] + contracted_twotensor_dotm1[k_gamma]);
                }

                // collect the final vector and tensor terms actually needed for the face construction
                double vector_corr[3] = {0}, tensor_corr[3][3] = {{0}};
                for(k=0;k<3;k++)
                {
                    vector_corr[k] = dA[k] + A*B[k];
                    for(k_x=0;k_x<3;k_x++) {tensor_corr[k][k_x] = B[k]*dA[k_x] + A*dB[k][k_x];}
                }
                // assign these to an ordered list (for ease of reference) and to particle. order: A, B[3], (dA+A*B)[3], (dA.B+A.dB)[3][3]
                SphP[i].Tensor_CRK_Face_Corrections[0] = A;
                for(k=0;k<3;k++) {SphP[i].Tensor_CRK_Face_Corrections[1+k] = B[k];}
                for(k=0;k<3;k++) {SphP[i].Tensor_CRK_Face_Corrections[1+3+k] = vector_corr[k];}
                for(k=0;k<3;k++) {for(k_x=0;k_x<3;k_x++) {SphP[i].Tensor_CRK_Face_Corrections[1+3+3+3*k+k_x] = tensor_corr[k][k_x];}}
            }
#endif

        }


    /* free the temporary structure we created for the MinMax and additional data passing */
    myfree(GasGradDataPasser);

    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall = timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;

    CPU_Step[CPU_DENSCOMPUTE] += timecomp;
    CPU_Step[CPU_DENSWAIT] += timewait;
    CPU_Step[CPU_DENSCOMM] += timecomm;
    CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}


/* this is the main work routine for the gradients calculations */
/*!   -- this subroutine ONLY should write to shared memory when the local pairwise 'swap_to_j' flag is set. that flag should never be active in OPENMP runs, by the definitions below, for thread safety.
    comparing this to multithreaded code using thread locks or atomic for safety shows the latter provides no performance gain and often a loss, so this is better for safety and speed. if you
    are adding to this routine, you must protect any writes to shared memory appropriately [including j variables and any global variables as well] -- */

int GasGrad_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
                     int *ngblist, int gradient_iteration)
{
    int startnode, numngb, listindex = 0;
    int j, k, k2, n, swap_to_j;
    double hinv, hinv3, hinv4, r2, u, hinv_j, hinv3_j, hinv4_j;
    struct kernel_GasGrad kernel;
    struct GasGraddata_in local;
    struct GasGraddata_out out;
    struct GasGraddata_out_iter out_iter;
    if(gradient_iteration==0)
    {
        memset(&out, 0, sizeof(struct GasGraddata_out));
    } else {
        memset(&out_iter, 0, sizeof(struct GasGraddata_out_iter));
    }
    memset(&kernel, 0, sizeof(struct kernel_GasGrad));

    if(mode == 0)
        particle2in_GasGrad(&local, target, gradient_iteration);
    else
        local = GasGradDataGet[target];

    /* check if we should bother doing a neighbor loop */
    if(gradient_iteration == 0)
        if(local.GQuant.Density <= 0) return 0;

    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.Hsml;
    double h2_i = kernel.h_i*kernel.h_i;
    kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
    int sph_gradients_flag_i = 0;
    int sph_gradients_flag_j = 0;
    if(local.Mass < 0) {sph_gradients_flag_i=1; local.Mass*=-1;}
    double V_i;
    V_i = local.Mass / local.GQuant.Density;

    int kernel_mode_i = -1; // only need to calculate wk, by default
    if(sph_gradients_flag_i) kernel_mode_i = 0; // for sph, only need dwk
#if defined(HYDRO_SPH) || defined(KERNEL_CRK_FACES)
    kernel_mode_i = 0; // for some circumstances, we require both wk and dwk //
#endif


    /* Now start the actual neighbor computation for this particle */

    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = GasGradDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
#ifdef TURB_DIFF_DYNAMIC
            if (gradient_iteration == 0) {
                numngb = ngb_treefind_pairs_threads(local.Pos, All.TurbDynamicDiffFac * kernel.h_i, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            }
            else
#endif
            {
                numngb = ngb_treefind_pairs_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            }
            if(numngb < 0) {return -2;}

            for(n = 0; n < numngb; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(GasGrad_isactive(j)==0) continue;

                integertime TimeStep_J; TimeStep_J = GET_PARTICLE_INTEGERTIME(j);
#if 0 //!defined(BOX_SHEARING) && !defined(_OPENMP) // (shearing box means the fluxes at the boundaries are not actually symmetric, so can't do this; OpenMP on some new compilers goes bad here because pointers [e.g. P...] are not thread-safe shared with predictive operations, and vectorization means no gain here with OMP anyways) //
                if(local.Timestep > TimeStep_J) continue; /* compute from particle with smaller timestep */
                /* use relative positions to break degeneracy */
                if(local.Timestep == TimeStep_J)
                {
                    int n0=0; if(local.Pos[n0] == P[j].Pos[n0]) {n0++; if(local.Pos[n0] == P[j].Pos[n0]) n0++;}
                    if(local.Pos[n0] < P[j].Pos[n0]) continue;
                }
                swap_to_j = TimeBinActive[P[j].TimeBin];
#else
                swap_to_j = 0;
#endif

                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); /*  now find the closest image in the given box size  */
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                double h_j = PPP[j].Hsml;
#if !defined(HYDRO_SPH) && !defined(KERNEL_CRK_FACES)
                if(r2 <= 0) continue;
#else
                if(r2 <= 0) {swap_to_j = 0;}
#endif
#ifdef TURB_DIFF_DYNAMIC
#ifdef GALSF_SUBGRID_WINDS
                if (gradient_iteration == 0 && ((SphP[j].DelayTime == 0 && local.DelayTime == 0) || (SphP[j].DelayTime > 0 && local.DelayTime > 0))) {
#else
                if (gradient_iteration == 0) {
#endif
                    double hhat_i, hhat_j, hhatinv_i, hhatinv3_i, hhatinv4_i, wkhat_i, dwkhat_i;
                    hhat_i = All.TurbDynamicDiffFac * kernel.h_i; hhat_j = All.TurbDynamicDiffFac * h_j;
                    if((r2 >= (hhat_i * hhat_i)) && (r2 >= (hhat_j * hhat_j))) continue;
                    double h_avg = 0.5 * (hhat_i + hhat_j), particle_distance = sqrt(r2);
                    kernel_hinv(h_avg, &hhatinv_i, &hhatinv3_i, &hhatinv4_i); u = DMIN(particle_distance * hhatinv_i, 1.0);
                    if(u<1) {kernel_main(u, hhatinv3_i, hhatinv4_i, &wkhat_i, &dwkhat_i, 0);} else {wkhat_i=dwkhat_i=0;} /* wkhat is symmetric in this case W_{ij} = W_{ji} */
                    double mean_weight = wkhat_i * 0.5 * (SphP[j].Norm_hat + local.Norm_hat) / (local.Norm_hat * SphP[j].Norm_hat);
                    double weight_i = P[j].Mass * mean_weight, weight_j = local.Mass * mean_weight, Velocity_bar_diff[3];
                    if(particle_distance < h_avg) {
                        for(k=0;k<3;k++) {Velocity_bar_diff[k] = SphP[j].Velocity_bar[k] - local.GQuant.Velocity_bar[k]; out.Velocity_hat[k] += Velocity_bar_diff[k] * weight_i;}
                        if(swap_to_j) {for(k=0;k<3;k++) {SphP[j].Velocity_hat[k] -= Velocity_bar_diff[k] * weight_j;}}
                    }
                } /* closes gradient_iteration == 0 */
#endif
                if((r2 >= h2_i) && (r2 >= h_j * h_j)) continue;

                kernel.r = sqrt(r2);
                if(kernel.r < kernel.h_i)
                {
                    u = kernel.r * hinv;
                    kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, kernel_mode_i);
                }
                else
                {
                    kernel.dwk_i = kernel.wk_i = 0;
                }
#if defined(MHD_CONSTRAINED_GRADIENT) || defined(KERNEL_CRK_FACES) || defined(HYDRO_TENSOR_FACE_CORRECTIONS_NGBITER)
                if(kernel.r < h_j)
#else
                if((kernel.r < h_j) && (swap_to_j))
#endif
                {
                    /* ok, we need the j-particle weights, but first check what kind of gradient we are calculating */
                    sph_gradients_flag_j = SHOULD_I_USE_SPH_GRADIENTS(SphP[j].ConditionNumber);
                    int kernel_mode_j;
#if defined(HYDRO_SPH) || defined(KERNEL_CRK_FACES)
                    kernel_mode_j = 0; // for some circumstances, we require both wk and dwk //
#else
                    if(sph_gradients_flag_j) {kernel_mode_j=0;} else {kernel_mode_j=-1;}
#endif
                    kernel_hinv(h_j, &hinv_j, &hinv3_j, &hinv4_j);
                    u = kernel.r * hinv_j;
                    kernel_main(u, hinv3_j, hinv4_j, &kernel.wk_j, &kernel.dwk_j, kernel_mode_j);
                }
                else
                {
                    kernel.dwk_j = kernel.wk_j = 0;
                }
                double Particle_Size_j, Particle_Size_i;  Particle_Size_j=Get_Particle_Size(j); Particle_Size_i=pow(local.Mass/local.GQuant.Density, 1./NUMDIMS);

#if defined(MHD_CONSTRAINED_GRADIENT)
                double V_j = P[j].Mass / SphP[j].Density, Face_Area_Vec[3], Face_Area_Norm, cnumcrit2 = ((double)CONDITION_NUMBER_DANGER)*((double)CONDITION_NUMBER_DANGER) - local.ConditionNumber*local.ConditionNumber;

#include "compute_finitevol_faces.h" /* insert code block for computing Face_Area_Vec, Face_Area_Norm, n_unit, etc. */
                
                for(k=0;k<3;k++){
                    if(gradient_iteration == 0)
                    {
                        out.Face_Area[k] += Face_Area_Vec[k];
                        if(swap_to_j) SphP[j].Face_Area[k] -= Face_Area_Vec[k];

                        for(k2=0;k2<3;k2++)
                        {
                            double q = -0.5 * Face_Area_Vec[k] * kernel.dp[k2];
                            out.FaceCrossX[k][k2] += q;
                            if(swap_to_j) GasGradDataPasser[j].FaceCrossX[k][k2] += q;
                        }
                    }

                    /* now use the gradients to construct the B_L,R states */
                    double Bjk = Get_Gas_BField(j,k);
                    double db_c=0, db_cR=0;
                    for(k2=0;k2<3;k2++)
                    {
                        db_c += 0.5 * SphP[j].Gradients.B[k][k2] * kernel.dp[k2];
                        db_cR -= 0.5 * local.BGrad[k][k2]  * kernel.dp[k2];
                    }

                    /* now we apply our slope-limiter to the B_L,R reconstruction */
                    double Q_L, Q_R;
                    if(Bjk == local.GQuant.B[k])
                    {
                        Q_L = Q_R = Bjk;
                    } else {
                        Q_L = Bjk + db_c;
                        Q_R = local.GQuant.B[k] + db_cR;
                        double Qmax, Qmin, Qmed = 0.5*(local.GQuant.B[k] + Bjk);
                        if(local.GQuant.B[k] < Bjk) {Qmax=Bjk; Qmin=local.GQuant.B[k];} else {Qmax=local.GQuant.B[k]; Qmin=Bjk;}
                        double fac = MHD_CONSTRAINED_GRADIENT_FAC_MINMAX * (Qmax-Qmin);
                        fac += MHD_CONSTRAINED_GRADIENT_FAC_MAX_PM * fabs(Qmed);
                        double Qmax_eff = Qmax + fac;
                        double Qmin_eff = Qmin - fac;
                        fac = MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV * (Qmax-Qmin);
                        fac += MHD_CONSTRAINED_GRADIENT_FAC_MED_PM * fabs(Qmed);
                        double Qmed_max = Qmed + fac;
                        double Qmed_min = Qmed - fac;
                        if(Qmed_max>Qmax_eff) Qmed_max=Qmax_eff;
                        if(Qmed_min<Qmin_eff) Qmed_min=Qmin_eff;
                        if(local.GQuant.B[k] < Bjk)
                        {
                            if(Q_L>Qmax_eff) Q_L=Qmax_eff;
                            if(Q_L<Qmed_min) Q_L=Qmed_min;
                            if(Q_R<Qmin_eff) Q_R=Qmin_eff;
                            if(Q_R>Qmed_max) Q_R=Qmed_max;
                        } else {
                            if(Q_L<Qmin_eff) Q_L=Qmin_eff;
                            if(Q_L>Qmed_max) Q_L=Qmed_max;
                            if(Q_R>Qmax_eff) Q_R=Qmax_eff;
                            if(Q_R<Qmed_min) Q_R=Qmed_min;
                        }
                    }

                    if(gradient_iteration==0)
                    {
                        out.FaceDotB += Face_Area_Vec[k] * (local.GQuant.B[k] + Q_L);
                    } else {
                        out_iter.FaceDotB += Face_Area_Vec[k] * (local.GQuant.B[k] + Q_L);
                    }
                    if(swap_to_j) GasGradDataPasser[j].FaceDotB -= Face_Area_Vec[k] * (Bjk + Q_R);
                }

#if defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
                /* this will fit the gradient at the -midpoint- as opposed to at the j locations, i.e.
                 attempting to minimize the quantity phi_L - phi_R, at face locations */
                double dphi = Get_Gas_PhiField(j) - local.GQuant.Phi;
                if(gradient_iteration == 0)
                {
                    MINMAX_CHECK(dphi,out.Minima.Phi,out.Maxima.Phi);
                    if(swap_to_j) {MINMAX_CHECK(-dphi,GasGradDataPasser[j].Minima.Phi,GasGradDataPasser[j].Maxima.Phi);}
                }

                // dphi = phi_j - phi_i :: if phi_i = 0, dphi = phi_j //
                double dphi_grad_j = 0, dphi_grad_i = 0;
                for(k=0;k<3;k++)
                {
                    dphi_grad_j += 0.5 * kernel.dp[k] * SphP[j].Gradients.Phi[k];
                    dphi_grad_i -= 0.5 * kernel.dp[k] * local.PhiGrad[k];
                }
                if(dphi > 0)
                {
                    if(dphi_grad_j>0) {dphi_grad_j=0;} else {if(dphi_grad_j<0.5*dphi) dphi_grad_j=0.5*dphi;}
                    if(dphi_grad_i<0) {dphi_grad_i=0;} else {if(dphi_grad_i>0.5*dphi) dphi_grad_i=0.5*dphi;}
                } else {
                    if(dphi_grad_j<0) {dphi_grad_j=0;} else {if(dphi_grad_j>0.5*dphi) dphi_grad_j=0.5*dphi;}
                    if(dphi_grad_i>0) {dphi_grad_i=0;} else {if(dphi_grad_i<0.5*dphi) dphi_grad_i=0.5*dphi;}
                }
                double dphi_j = dphi + dphi_grad_j;
                double dphi_i = dphi - dphi_grad_i;
                if(sph_gradients_flag_i) {dphi_j *= -2*kernel.wk_i;} else {dphi_j *= kernel.dwk_i/kernel.r * P[j].Mass;}
                if(sph_gradients_flag_j) {dphi_i *= -2*kernel.wk_j;} else {dphi_i *= kernel.dwk_j/kernel.r * local.Mass;}
                if(gradient_iteration == 0) {for(k=0;k<3;k++) {out.Gradients[k].Phi += dphi_j * kernel.dp[k];}} else {for(k=0;k<3;k++) {out_iter.PhiGrad[k] += dphi_j * kernel.dp[k];}}
                if(swap_to_j) {for(k=0;k<3;k++) {GasGradDataPasser[j].PhiGrad[k] += dphi_i * kernel.dp[k];}}
#endif
#endif // MHD_CONSTRAINED_GRADIENT

                if(gradient_iteration == 0)
                {
                    /* ------------------------------------------------------------------------------------------------ */
                    /* DIFFERENCE & SLOPE LIMITING: need to check maxima and minima of particle values in the kernel, to avoid
                     'overshoot' with our gradient estimators. this check should be among all interacting pairs */

                    if(kernel.r > out.MaxDistance) {out.MaxDistance = kernel.r;}
                    if(swap_to_j) {if(kernel.r > GasGradDataPasser[j].MaxDistance) {GasGradDataPasser[j].MaxDistance = kernel.r;}}

                    double d_rho = SphP[j].Density - local.GQuant.Density;
                    MINMAX_CHECK(d_rho,out.Minima.Density,out.Maxima.Density);
                    if(swap_to_j) {MINMAX_CHECK(-d_rho,GasGradDataPasser[j].Minima.Density,GasGradDataPasser[j].Maxima.Density);}

                    double dp = SphP[j].Pressure - local.GQuant.Pressure;
                    MINMAX_CHECK(dp,out.Minima.Pressure,out.Maxima.Pressure);
                    if(swap_to_j) {MINMAX_CHECK(-dp,GasGradDataPasser[j].Minima.Pressure,GasGradDataPasser[j].Maxima.Pressure);}

#ifdef TURB_DIFF_DYNAMIC
                    double dv_bar[3]; for(k=0;k<3;k++) {dv_bar[k] = SphP[j].Velocity_bar[k] - local.GQuant.Velocity_bar[k];} /* Need to calculate the filtered velocity gradient for the filtered shear */
                    NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,dv_bar,-1); /* wrap velocities for shearing boxes if needed */
                    for(k=0;k<3;k++) {MINMAX_CHECK(dv_bar[k], out.Minima.Velocity_bar[k], out.Maxima.Velocity_bar[k]);
                        if(swap_to_j) {MINMAX_CHECK(-dv_bar[k], GasGradDataPasser[j].Minima.Velocity_bar[k], GasGradDataPasser[j].Maxima.Velocity_bar[k]);}}
#endif

                    
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS)
                    double m3wt[10], m2wt[6], wt_i=kernel.wk_i, wt_j=kernel.wk_j, wt_i_A=wt_i, wt_i_B=wt_i, wt_j_A=wt_j, wt_j_B=wt_j;
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS_NGBITER)
                    {
                        double V_i = local.Mass/local.GQuant.Density, V_j = P[j].Mass/SphP[j].Density;
                        wt_i_A = (V_i*kernel.wk_i + V_j*kernel.wk_j) / (2.*V_i);
                        wt_i_B = (V_i*kernel.wk_i - V_j*kernel.wk_j) / (2.*V_i);
                        wt_j_A = wt_i_A * V_i/V_j; wt_j_B = wt_i_B * V_i/V_j;
                    }
#endif
                    m3wt[0]=kernel.dp[0]*kernel.dp[0]*kernel.dp[0]; m3wt[1]=kernel.dp[0]*kernel.dp[0]*kernel.dp[1]; m3wt[2]=kernel.dp[0]*kernel.dp[0]*kernel.dp[2];
                    m3wt[3]=kernel.dp[1]*kernel.dp[1]*kernel.dp[1]; m3wt[4]=kernel.dp[1]*kernel.dp[1]*kernel.dp[0]; m3wt[5]=kernel.dp[1]*kernel.dp[1]*kernel.dp[2];
                    m3wt[6]=kernel.dp[2]*kernel.dp[2]*kernel.dp[2]; m3wt[7]=kernel.dp[2]*kernel.dp[2]*kernel.dp[0]; m3wt[8]=kernel.dp[2]*kernel.dp[2]*kernel.dp[1];
                    m3wt[9]=kernel.dp[0]*kernel.dp[1]*kernel.dp[2];
                    for(k=0;k<10;k++) {out.third_moment_weights[k] += wt_i_B * m3wt[k];}
                    for(k=0;k<3;k++) {out.first_second_moment_weights[k] += wt_i_A * kernel.dp[k];}
                    m2wt[0]=kernel.dp[0]*kernel.dp[0]; m2wt[1]=kernel.dp[0]*kernel.dp[1]; m2wt[2]=kernel.dp[0]*kernel.dp[2];
                    m2wt[3]=kernel.dp[1]*kernel.dp[1]; m2wt[4]=kernel.dp[1]*kernel.dp[2]; m2wt[5]=kernel.dp[2]*kernel.dp[2];
                    for(k=0;k<6;k++) {out.first_second_moment_weights[k+3] += wt_i_A * m2wt[k];}
                    if(swap_to_j) {
                        for(k=0;k<10;k++) {GasGradDataPasser[j].third_moment_weights[k] -= wt_j_B * m3wt[k];}
                        for(k=0;k<3;k++) {GasGradDataPasser[j].first_second_moment_weights[k] -= wt_j_A * kernel.dp[k];}
                        for(k=0;k<6;k++) {GasGradDataPasser[j].first_second_moment_weights[k+3] += wt_j_A * m2wt[k];}
                    }
#if defined(HYDRO_TENSOR_FACE_CORRECTIONS_NGBITER)
                    for(k=0;k<6;k++) {out.first_second_moment_weights[k+3+6] += wt_i_B * m2wt[k];}
                    if(swap_to_j) {for(k=0;k<6;k++) {out.first_second_moment_weights[k+3+6] -= wt_i_B * m2wt[k];}}
#endif
#endif

#if defined(KERNEL_CRK_FACES)
                    {
                        double V_i = local.Mass/local.GQuant.Density, V_j = P[j].Mass/SphP[j].Density;
                        double wk_ij = 0.5*(kernel.wk_i + kernel.wk_j), dwk_ij = 0.5*(kernel.dwk_i + kernel.dwk_j), rinv = 1./(MIN_REAL_NUMBER + kernel.r);
                        double Vj_wki = V_j*wk_ij, Vj_dwki = V_j*dwk_ij*rinv, Vi_wkj = V_i*wk_ij, Vi_dwkj = V_i*dwk_ij*rinv;
                        out.m0 += Vj_wki;
                        for(k=0;k<3;k++) {out.dm0[k] += Vj_dwki*kernel.dp[k];}
                        for(k2=0;k2<3;k2++)
                        {
                            out.m1[k2] += Vj_wki*kernel.dp[k2];
                            for(k=0;k<3;k++) {out.dm1[k2][k] += Vj_dwki*kernel.dp[k2]*kernel.dp[k];}
                        }
                        for(k2=0;k2<6;k2++)
                        {
                            int kk0[6]={0,1,2,0,0,1};
                            int kk1[6]={0,1,2,1,2,2};
                            out.m2[k2] += Vj_wki*kernel.dp[kk0[k2]]*kernel.dp[kk1[k2]];
                            for(k=0;k<3;k++) {out.dm2[k2][k] += Vj_dwki*kernel.dp[kk0[k2]]*kernel.dp[kk1[k2]]*kernel.dp[k];}
                        }
                        if(swap_to_j)
                        {
                            GasGradDataPasser[j].m0 += Vi_wkj;
                            for(k=0;k<3;k++) {GasGradDataPasser[j].dm0[k] -= Vi_dwkj*kernel.dp[k];}
                            for(k2=0;k2<3;k2++)
                            {
                                GasGradDataPasser[j].m1[k2] -= Vi_wkj*kernel.dp[k2];
                                for(k=0;k<3;k++) {GasGradDataPasser[j].dm1[k2][k] += Vi_dwkj*kernel.dp[k2]*kernel.dp[k];}
                            }
                            for(k2=0;k2<6;k2++)
                            {
                                int kk0[6]={0,1,2,0,0,1};
                                int kk1[6]={0,1,2,1,2,2};
                                GasGradDataPasser[j].m2[k2] += Vi_wkj*kernel.dp[kk0[k2]]*kernel.dp[kk1[k2]];
                                for(k=0;k<3;k++) {GasGradDataPasser[j].dm2[k2][k] -= Vi_dwkj*kernel.dp[kk0[k2]]*kernel.dp[kk1[k2]]*kernel.dp[k];}
                            }
                        }
                    }
#endif

                    double dv[3]; for(k=0;k<3;k++) {dv[k] = SphP[j].VelPred[k] - local.GQuant.Velocity[k];}
                    NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,dv,-1); /* wrap velocities for shearing boxes if needed */
                    for(k=0;k<3;k++) {
                        MINMAX_CHECK(dv[k],out.Minima.Velocity[k],out.Maxima.Velocity[k]);
                        if(swap_to_j) {MINMAX_CHECK(-dv[k],GasGradDataPasser[j].Minima.Velocity[k],GasGradDataPasser[j].Maxima.Velocity[k]);}
                    }

#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
                    for(k=0;k<3;k++)
                    {
                        double GlassAcc = kernel.dp[k] / (kernel.r*kernel.r*kernel.r); // acceleration to apply to force cells into a glass
                        out.GlassAcc[k] += GlassAcc;
                        if(swap_to_j) {GasGradDataPasser[j].GlassAcc[k] -= GlassAcc;}
                    }
#endif

#ifdef DOGRAD_INTERNAL_ENERGY
                    double du = SphP[j].InternalEnergyPred - local.GQuant.InternalEnergy;
                    MINMAX_CHECK(du,out.Minima.InternalEnergy,out.Maxima.InternalEnergy);
                    if(swap_to_j) {MINMAX_CHECK(-du,GasGradDataPasser[j].Minima.InternalEnergy,GasGradDataPasser[j].Maxima.InternalEnergy);}
#endif
#ifdef COSMIC_RAY_FLUID
                    double dpCR[N_CR_PARTICLE_BINS];
                    for(k=0;k<N_CR_PARTICLE_BINS;k++)
                    {
                        dpCR[k] = Get_Gas_CosmicRayPressure(j,k) - local.GQuant.CosmicRayPressure[k];
                        MINMAX_CHECK(dpCR[k],out.Minima.CosmicRayPressure[k],out.Maxima.CosmicRayPressure[k]);
                        if(swap_to_j) {MINMAX_CHECK(-dpCR[k],GasGradDataPasser[j].Minima.CosmicRayPressure[k],GasGradDataPasser[j].Maxima.CosmicRayPressure[k]);}
                    }
#endif
#ifdef DOGRAD_SOUNDSPEED
                    double dc = Get_Gas_effective_soundspeed_i(j) - local.GQuant.SoundSpeed;
                    MINMAX_CHECK(dc,out.Minima.SoundSpeed,out.Maxima.SoundSpeed);
                    if(swap_to_j) {MINMAX_CHECK(-dc,GasGradDataPasser[j].Minima.SoundSpeed,GasGradDataPasser[j].Maxima.SoundSpeed);}
#endif
#ifdef MAGNETIC
                    double Bj[3],dB[3];
                    for(k=0;k<3;k++)
                    {
                        Bj[k] = Get_Gas_BField(j,k);
                        dB[k] = Bj[k] - local.GQuant.B[k];
                        MINMAX_CHECK(dB[k],out.Minima.B[k],out.Maxima.B[k]);
                        if(swap_to_j) {MINMAX_CHECK(-dB[k],GasGradDataPasser[j].Minima.B[k],GasGradDataPasser[j].Maxima.B[k]);}
                    }
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
                    double dphi = Get_Gas_PhiField(j) - local.GQuant.Phi;
                    MINMAX_CHECK(dphi,out.Minima.Phi,out.Maxima.Phi);
                    if(swap_to_j) {MINMAX_CHECK(-dphi,GasGradDataPasser[j].Minima.Phi,GasGradDataPasser[j].Maxima.Phi);}
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
                    double dmetal[NUM_METAL_SPECIES];
                    for(k = 0; k < NUM_METAL_SPECIES; k++)
                    {
                        dmetal[k] = P[j].Metallicity[k] - local.GQuant.Metallicity[k];
                        MINMAX_CHECK(dmetal[k],out.Minima.Metallicity[k],out.Maxima.Metallicity[k]);
                        if(swap_to_j) {MINMAX_CHECK(-dmetal[k],GasGradDataPasser[j].Minima.Metallicity[k],GasGradDataPasser[j].Maxima.Metallicity[k]);}
                    }
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
                    double dnET[N_RT_FREQ_BINS][6];
                    double dn[N_RT_FREQ_BINS];
                    double V_i_inv = 1/V_i, V_j_inv = SphP[j].Density/P[j].Mass;
                    for(k = 0; k < N_RT_FREQ_BINS; k++)
                    {
                        int k_dE; for(k_dE=0;k_dE<6;k_dE++) {dnET[k][k_dE] = SphP[j].Rad_E_gamma_Pred[k]*SphP[j].ET[k][k_dE]*V_j_inv - local.GQuant.Rad_E_gamma[k]*local.GQuant.Rad_E_gamma_ET[k][k_dE]*V_i_inv;}
                        dn[k] = SphP[j].Rad_E_gamma_Pred[k]*V_j_inv - local.GQuant.Rad_E_gamma[k]*V_i_inv;
                        MINMAX_CHECK(dn[k],out.Minima.Rad_E_gamma[k],out.Maxima.Rad_E_gamma[k]);
                        if(swap_to_j) {MINMAX_CHECK(-dn[k],GasGradDataPasser[j].Minima.Rad_E_gamma[k],GasGradDataPasser[j].Maxima.Rad_E_gamma[k]);}
                    }
#endif
                    /* end of difference and slope-limiter (min/max) block */
                    /* ------------------------------------------------------------------------------------------------ */


                    /* ------------------------------------------------------------------------------------------------ */
                    /*  Here we insert additional operations we want to fit into the gradients loop. at the moment, all of these are SPH-specific */
#ifdef HYDRO_SPH
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
                    out.alpha_limiter += NV_MYSIGN(SphP[j].NV_DivVel) * P[j].Mass * kernel.wk_i;
                    if(swap_to_j) SphP[j].alpha_limiter += NV_MYSIGN(local.NV_DivVel) * local.Mass * kernel.wk_j;
#endif
#ifdef MAGNETIC
                    double mji_dwk_r = P[j].Mass * kernel.dwk_i / kernel.r;
                    double mij_dwk_r = local.Mass * kernel.dwk_j / kernel.r;
                    for(k=0;k<3;k++)
                    {
                        for(k2=0;k2<3;k2++)
                        {
                            out.DtB[k] += local.GQuant.B[k2] * mji_dwk_r * kernel.dp[k2] * dv[k];
                            if(swap_to_j) SphP[j].DtB[k] += Bj[k2] * mij_dwk_r * kernel.dp[k2] * dv[k];
                        }
#ifdef DIVBCLEANING_DEDNER
                        out.divB += dB[k] * kernel.dp[k] * mji_dwk_r;
                        if(swap_to_j) SphP[j].divB += dB[k] * kernel.dp[k] * mij_dwk_r;
#endif
                    }
#endif
#endif
                    /* end of additional/miscellaneous operators block */
                    /* ------------------------------------------------------------------------------------------------ */


                    /* ------------------------------------------------------------------------------------------------ */
                    /* Finally, save actual output for GRADIENTS */

                    /* first do particle i */
                    if(kernel.r < kernel.h_i)
                    {
                        if(sph_gradients_flag_i && kernel.r > 0) {kernel.wk_i = -kernel.dwk_i/kernel.r * P[j].Mass;} // sph-like weights for gradients //
                        for(k=0;k<3;k++)
                        {
                            double wk_xyz_i = -kernel.wk_i * kernel.dp[k]; /* sign is important here! */
                            out.Gradients[k].Density += wk_xyz_i * d_rho;
                            out.Gradients[k].Pressure += wk_xyz_i * dp;
                            for(k2=0;k2<3;k2++) {out.Gradients[k].Velocity[k2] += wk_xyz_i * dv[k2];}
#ifdef TURB_DIFF_DYNAMIC
                            for (k2 = 0; k2 < 3; k2++) {out.Gradients[k].Velocity_bar[k2] += wk_xyz_i * dv_bar[k2];}
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
                            out.Gradients[k].InternalEnergy += wk_xyz_i * du;
#endif
#ifdef COSMIC_RAY_FLUID
                            for(k2=0;k2<N_CR_PARTICLE_BINS;k2++) {out.Gradients[k].CosmicRayPressure[k2] += wk_xyz_i * dpCR[k2];}
#endif
#ifdef DOGRAD_SOUNDSPEED
                            out.Gradients[k].SoundSpeed += wk_xyz_i * dc;
#endif
#ifdef MAGNETIC
                            for(k2=0;k2<3;k2++) {out.Gradients[k].B[k2] += wk_xyz_i * dB[k2];}
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
                            out.Gradients[k].Phi += wk_xyz_i * dphi;
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
                            for(k2=0;k2<NUM_METAL_SPECIES;k2++) {out.Gradients[k].Metallicity[k2] += wk_xyz_i * dmetal[k2];}
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
                            for(k2=0;k2<N_RT_FREQ_BINS;k2++)
                            {
                            	out.Gradients[k].Rad_E_gamma[k2] += wk_xyz_i * dn[k2];
                            	int k_et; for(k_et=0;k_et<6;k_et++) out.Gradients[k].Rad_E_gamma_ET[k2][k_et] += wk_xyz_i * dnET[k2][k_et];
                            }
#endif
                        }
                    }

                    /* next do particle j */
                    if((kernel.r < h_j) && (swap_to_j))
                    {
                        if(sph_gradients_flag_j && kernel.r > 0) {kernel.wk_j = -kernel.dwk_j/kernel.r * local.Mass;} // sph-like weights for gradients //
                        for(k=0;k<3;k++)
                        {
                            double wk_xyz_j = -kernel.wk_j * kernel.dp[k]; /* sign is important here! (note dp-dd signs cancel) */
                            SphP[j].Gradients.Density[k] += wk_xyz_j * d_rho;
                            SphP[j].Gradients.Pressure[k] += wk_xyz_j * dp;
                            for(k2=0;k2<3;k2++) {SphP[j].Gradients.Velocity[k2][k] += wk_xyz_j * dv[k2];}
#ifdef TURB_DIFF_DYNAMIC
                            for (k2 = 0; k2 < 3; k2++) {GasGradDataPasser[j].GradVelocity_bar[k2][k] += wk_xyz_j * dv_bar[k2];}
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
                            SphP[j].Gradients.InternalEnergy[k] += wk_xyz_j * du;
#endif
#ifdef COSMIC_RAY_FLUID
                            for(k2=0;k2<N_CR_PARTICLE_BINS;k2++) {SphP[j].Gradients.CosmicRayPressure[k2][k] += wk_xyz_j * dpCR[k2];}
#endif
#ifdef DOGRAD_SOUNDSPEED
                            SphP[j].Gradients.SoundSpeed[k] += wk_xyz_j * dc;
#endif
#ifdef MAGNETIC
#ifdef MHD_CONSTRAINED_GRADIENT
                            for(k2=0;k2<3;k2++) {GasGradDataPasser[j].BGrad[k2][k] += wk_xyz_j * dB[k2];}
#else
                            for(k2=0;k2<3;k2++) {SphP[j].Gradients.B[k2][k] += wk_xyz_j * dB[k2];}
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
                            SphP[j].Gradients.Phi[k] += wk_xyz_j * dphi;
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
                            for(k2=0;k2<NUM_METAL_SPECIES;k2++) {SphP[j].Gradients.Metallicity[k2][k] += wk_xyz_j * dmetal[k2];}
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
                            for(k2=0;k2<N_RT_FREQ_BINS;k2++)
                            {
                            	GasGradDataPasser[j].Gradients_Rad_E_gamma[k2][k] += wk_xyz_j * dn[k2];
								/* below we have the gradient dotted into the Eddington tensor (more complicated than a scalar gradient, but should recover full anisotropy */
								int k_freq=k2,k_xyz,j_xyz,i_xyz=k,k_et_loop[3]; // recall, for ET: 0=xx,1=yy,2=zz,3=xy,4=yz,5=xz
								for(k_xyz=0;k_xyz<3;k_xyz++)
								{
									if(k_xyz==0) {k_et_loop[0]=0; k_et_loop[1]=3; k_et_loop[2]=5;}
									if(k_xyz==1) {k_et_loop[0]=3; k_et_loop[1]=1; k_et_loop[2]=4;}
									if(k_xyz==2) {k_et_loop[0]=5; k_et_loop[1]=4; k_et_loop[2]=2;}
									for(j_xyz=0;j_xyz<3;j_xyz++)
									{
										SphP[j].Gradients.Rad_E_gamma_ET[k_freq][k_xyz] += SphP[j].NV_T[j_xyz][i_xyz] * wk_xyz_j * dnET[k_freq][k_et_loop[j_xyz]];
									}
								}
                            }
#endif
                        }
                    }

                    /* end of GRADIENTS calculation block */
                    /* ------------------------------------------------------------------------------------------------ */


                } // (r2 < h2i || r2 < h2j) && gradient_iteration==0
            } // numngb loop
        } // while(startnode)

        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = GasGradDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }


    /* ------------------------------------------------------------------------------------------------ */
    /* Now collect the result at the right place */
    if(gradient_iteration==0)
    {
        if(mode == 0)
            out2particle_GasGrad(&out, target, 0, gradient_iteration);
        else
            GasGradDataResult[target] = out;
    } else {
        if(mode == 0)
            out2particle_GasGrad_iter(&out_iter, target, 0, gradient_iteration);
        else
            GasGradDataResult_iter[target] = out_iter;
    }
    /* ------------------------------------------------------------------------------------------------ */

    return 0;
}





void *GasGrad_evaluate_primary(void *p, int gradient_iteration)
{
#define CONDITION_FOR_EVALUATION if(GasGrad_isactive(i))
#define EVALUATION_CALL GasGrad_evaluate(i,0,exportflag,exportnodecount,exportindex,ngblist,gradient_iteration)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *GasGrad_evaluate_secondary(void *p, int gradient_iteration)
{
#define EVALUATION_CALL GasGrad_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist, gradient_iteration);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}
