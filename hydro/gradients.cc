#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../core/step_phases.h"
#include "../system/gpu_particles_arena.h"
#include "../mesh/kernel.h"
#include "../mesh/neighbor_list.h"
/* gradient_evaluate_gpu writes results as GasGraddata_out_ structs (defined in
   gradient_functions.h, layout-identical to GasGraddata_out defined below).
   We pass a void* buffer to avoid including gradient_functions.h which would
   conflict with the local struct definitions. The cast to GasGraddata_out*
   is safe because the two structs have identical field layout. */
extern void gradient_evaluate_gpu(struct particle_data *, struct gas_cell_data *,
                                  int, int *, int, int64_t *, int *, int64_t, void *, int);
#include "../mesh/ghost_symlist_lifecycle.h"
#include "compute_finitevol_faces_functions.h"



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


#define NV_MYSIGN(x) (( x > 0 ) - ( x < 0 ))

/* function that tells us whether a given element should be active for gradient calculation*/
int GasGrad_isactive(int i, struct particle_data *pp, struct gas_cell_data *cell)
{
    if(pp[i].Type != 0) return 0;
    if(pp[i].Mass <= 0) return 0;
    if(cell[i].Density <= 0 || pp[i].KernelRadius <= 0) return 0;
#if defined(GALSF_SUBGRID_WINDS) && !defined(TURB_DIFF_DYNAMIC)
    if(cell[i].DelayTime > 0) return 0;
#endif
    return 1;
}

/* define a common 'gradients' structure to hold
 everything we're going to take derivatives of */
struct Quantities_for_Gradients
{
    MyDouble Density;
    MyDouble Pressure;
    Vec3<MyDouble> Velocity;
#ifdef MAGNETIC
    Vec3<MyDouble> B;
#ifdef DIVBCLEANING_DEDNER
    MyDouble Phi;
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
    MyDouble Metallicity[NUM_METAL_SPECIES];
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
    MyFloat Rad_E_gamma[N_RT_FREQ_BINS];
    SymmetricTensor2<MyFloat> Rad_E_gamma_ET[N_RT_FREQ_BINS];
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
    MyFloat Rad_Flux[N_RT_FREQ_BINS][3];
#endif
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
    MyDouble InternalEnergy;
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
    MyDouble ElectronNumberDensity;  /* n_e for Biermann battery: provides grad(n_e) used in dB/dt|_Bier ~ grad(T_e) x grad(n_e). */
    MyDouble ElectronTemperature;    /* T_e for Biermann battery (currently == T_gas). */
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
    Vec3<MyDouble> E_battery_T2;     /* Tier-2 battery EMF (sum of radiative + dust). Gradient pass produces grad(E_battery_T2); curl is taken in hydro_toplevel. */
#endif
#ifdef COSMIC_RAY_FLUID
    MyDouble CosmicRayPressure[N_CR_PARTICLE_BINS];
#endif
#ifdef DOGRAD_SOUNDSPEED
    MyDouble SoundSpeed;
#endif
#ifdef TURB_DIFF_DYNAMIC
    Vec3<MyDouble> Velocity_bar;
#endif
};

struct kernel_GasGrad
{
    Vec3<double> dp; double r,wk_i, wk_j, dwk_i, dwk_j,h_i;
};

struct GasGraddata_in
{
    Vec3<MyDouble> Pos;
    MyFloat Mass;
    MyFloat KernelRadius;
#ifdef MHD_CONSTRAINED_GRADIENT
    MyDouble ConditionNumber;
    SymmetricTensor2<MyDouble> NV_T;
    MyFloat BGrad[3][3];
#ifdef MHD_MODIFIED_GRADIENT
    MyFloat MG_cgcoeff;
#endif
#ifdef MHD_CONSTRAINED_GRADIENT_FAC_MEDDEV
    Vec3<MyFloat> PhiGrad;
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
};


struct GasGraddata_out
{
#if defined(KERNEL_CRK_FACES)
    MyDouble m0;
    MyDouble m1[3];
    MyDouble m2[6];
    MyDouble dm0[3];
    MyDouble dm1[3][3];
    MyDouble dm2[6][3];
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
    Vec3<MyFloat> GlassAcc;
#endif
#ifdef HYDRO_SPH
    MyFloat alpha_limiter;
#ifdef MAGNETIC
#ifdef DIVBCLEANING_DEDNER
    MyFloat divB;
#endif
    Vec3<MyFloat> DtB;
#endif
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
    Vec3<MyFloat> Face_Area;
    MyFloat FaceDotB;
    MyFloat FaceCrossX[3][3];
#endif
    struct Quantities_for_Gradients Gradients[3];
    struct Quantities_for_Gradients Maxima;
    struct Quantities_for_Gradients Minima;
    MyFloat MaxDistance;
#ifdef TURB_DIFF_DYNAMIC
    Vec3<MyDouble> Velocity_hat;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
    MyFloat AGS_zeta;
#endif
};



struct GasGraddata_out_iter
{
#ifdef MHD_CONSTRAINED_GRADIENT
    MyFloat FaceDotB;
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
    Vec3<MyDouble> PhiGrad;
#endif
#else
    MyFloat dummy;
#endif
};



/* this is a temporary structure for quantities used ONLY in the loop below,
 for example for computing the slope-limiters (for the Reimann problem) */
static struct temporary_data_topass
{
    struct Quantities_for_Gradients Maxima;
    struct Quantities_for_Gradients Minima;
    MyFloat MaxDistance;
#if defined(KERNEL_CRK_FACES)
    MyDouble m0;
    MyDouble m1[3];
    MyDouble m2[6];
    MyDouble dm0[3];
    MyDouble dm1[3][3];
    MyDouble dm2[6][3];
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
    Vec3<MyFloat> GlassAcc;
#endif
#ifdef MHD_CONSTRAINED_GRADIENT
    MyDouble FaceDotB;
    MyDouble FaceCrossX[3][3];
    MyDouble BGrad[3][3];
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
    Vec3<MyDouble> PhiGrad;
#endif
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
    Vec3<MyFloat> Gradients_Rad_E_gamma[N_RT_FREQ_BINS];
#endif
#ifdef TURB_DIFF_DYNAMIC
    Vec3<MyDouble> GradVelocity_bar[3];
#endif
}
*GasGradDataPasser;



static inline void out2particle_GasGrad(struct GasGraddata_out *out, int i, int mode, int gradient_iteration);
static inline void out2particle_GasGrad_iter(struct GasGraddata_out_iter *out, int i, int mode, int gradient_iteration);

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
        ASSIGN_ADD_PRESET(GasGradDataPasser[i].PhiGrad,out->PhiGrad,mode);
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
            ASSIGN_ADD_PRESET(CellP[i].Velocity_hat[j], out->Velocity_hat[j], mode);
            for (k = 0; k < 3; k++) {
                ASSIGN_ADD_PRESET(GasGradDataPasser[i].GradVelocity_bar[j][k], out->Gradients[k].Velocity_bar[j], mode);
            }
        }
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
        ASSIGN_ADD_PRESET(GasGradDataPasser[i].GlassAcc,out->GlassAcc,mode);
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        ASSIGN_ADD_PRESET(CellP[i].alpha_limiter, out->alpha_limiter, mode);
#endif

        MAX_ADD(GasGradDataPasser[i].Maxima.Density,out->Maxima.Density,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.Density,out->Minima.Density,mode);
        MAX_ADD(GasGradDataPasser[i].Maxima.Pressure,out->Maxima.Pressure,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.Pressure,out->Minima.Pressure,mode);
        for(k=0;k<3;k++)
        {
            ASSIGN_ADD_PRESET(CellP[i].Gradients.Density[k],out->Gradients[k].Density,mode);
            ASSIGN_ADD_PRESET(CellP[i].Gradients.Pressure[k],out->Gradients[k].Pressure,mode);
        }
#ifdef DOGRAD_INTERNAL_ENERGY
        MAX_ADD(GasGradDataPasser[i].Maxima.InternalEnergy,out->Maxima.InternalEnergy,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.InternalEnergy,out->Minima.InternalEnergy,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.InternalEnergy[k],out->Gradients[k].InternalEnergy,mode);}
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
        MAX_ADD(GasGradDataPasser[i].Maxima.ElectronNumberDensity,out->Maxima.ElectronNumberDensity,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.ElectronNumberDensity,out->Minima.ElectronNumberDensity,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.ElectronNumberDensity[k],out->Gradients[k].ElectronNumberDensity,mode);}
        MAX_ADD(GasGradDataPasser[i].Maxima.ElectronTemperature,out->Maxima.ElectronTemperature,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.ElectronTemperature,out->Minima.ElectronTemperature,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.ElectronTemperature[k],out->Gradients[k].ElectronTemperature,mode);}
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
        for(j=0;j<3;j++) {
            MAX_ADD(GasGradDataPasser[i].Maxima.E_battery_T2[j],out->Maxima.E_battery_T2[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.E_battery_T2[j],out->Minima.E_battery_T2[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.E_battery_T2[j][k],out->Gradients[k].E_battery_T2[j],mode);}
        }
#endif
#ifdef COSMIC_RAY_FLUID
        for(j=0;j<N_CR_PARTICLE_BINS;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.CosmicRayPressure[j],out->Maxima.CosmicRayPressure[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.CosmicRayPressure[j],out->Minima.CosmicRayPressure[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.CosmicRayPressure[j][k],out->Gradients[k].CosmicRayPressure[j],mode);}
        }
#endif
#ifdef DOGRAD_SOUNDSPEED
        MAX_ADD(GasGradDataPasser[i].Maxima.SoundSpeed,out->Maxima.SoundSpeed,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.SoundSpeed,out->Minima.SoundSpeed,mode);
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.SoundSpeed[k],out->Gradients[k].SoundSpeed,mode);}
#endif

        for(j=0;j<3;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.Velocity[j],out->Maxima.Velocity[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.Velocity[j],out->Minima.Velocity[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.Velocity[j][k],out->Gradients[k].Velocity[j],mode);}
        }

#ifdef MAGNETIC

#ifdef HYDRO_SPH
#ifdef DIVBCLEANING_DEDNER
        ASSIGN_ADD_PRESET(CellP[i].divB,out->divB, mode);
#endif
        ASSIGN_ADD_PRESET(CellP[i].DtB,out->DtB, mode);
#endif


#ifdef MHD_CONSTRAINED_GRADIENT
        for(j=0;j<3;j++)
        {
            ASSIGN_ADD_PRESET(CellP[i].Face_Area[j],out->Face_Area[j],mode);
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
                ASSIGN_ADD_PRESET(CellP[i].Gradients.B[j][k],out->Gradients[k].B[j],mode);
#endif
            }
        }

#ifdef DIVBCLEANING_DEDNER
        MAX_ADD(GasGradDataPasser[i].Maxima.Phi,out->Maxima.Phi,mode);
        MIN_ADD(GasGradDataPasser[i].Minima.Phi,out->Minima.Phi,mode);
#ifndef MHD_CONSTRAINED_GRADIENT_MIDPOINT
        for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.Phi[k],out->Gradients[k].Phi,mode);}
#endif
#endif
#endif // closes MAGNETIC

#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        for(j=0;j<NUM_METAL_SPECIES;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.Metallicity[j],out->Maxima.Metallicity[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.Metallicity[j],out->Minima.Metallicity[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.Metallicity[j][k],out->Gradients[k].Metallicity[j],mode);}
        }
#endif

#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
        for(j=0;j<N_RT_FREQ_BINS;j++)
        {
            MAX_ADD(GasGradDataPasser[i].Maxima.Rad_E_gamma[j],out->Maxima.Rad_E_gamma[j],mode);
            MIN_ADD(GasGradDataPasser[i].Minima.Rad_E_gamma[j],out->Minima.Rad_E_gamma[j],mode);
            for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(GasGradDataPasser[i].Gradients_Rad_E_gamma[j][k],out->Gradients[k].Rad_E_gamma[j],mode);}
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
            {int k_d; for(k_d=0;k_d<3;k_d++) {
                MAX_ADD(GasGradDataPasser[i].Maxima.Rad_Flux[j][k_d], out->Maxima.Rad_Flux[j][k_d], mode);
                MIN_ADD(GasGradDataPasser[i].Minima.Rad_Flux[j][k_d], out->Minima.Rad_Flux[j][k_d], mode);
                for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(CellP[i].Gradients.Rad_Flux_Grad[j][k_d][k], out->Gradients[k].Rad_Flux[j][k_d], mode);}
            }}
#endif
        }
		/* the gradient dotted into the Eddington tensor is more complicated: let's handle this below */
        {
        	int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
        	{
        		int k_xyz; for(k_xyz=0;k_xyz<3;k_xyz++)
        		{
        			int j_xyz,i_xyz;
        			for(j_xyz=0;j_xyz<3;j_xyz++)
        			{
        				for(i_xyz=0;i_xyz<3;i_xyz++)
        				{
        					CellP[i].Gradients.Rad_E_gamma_ET[k_freq][k_xyz] += CellP[i].NV_T[j_xyz][i_xyz] * out->Gradients[i_xyz].Rad_E_gamma_ET[k_freq][k_xyz][j_xyz];
						}
        			}
        		}
        	}
        }
#endif
        
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
        ASSIGN_ADD_PRESET(P[i].AGS_zeta, out->AGS_zeta,   mode);
#endif

    } // gradient_iteration == 0
}




void local_slopelimiter(double *grad, double valmax, double valmin, double alim, double h, double shoot_tol, int pos_preserve, double d_max, double val_cen)
{
    Vec3<double>& g = *reinterpret_cast<Vec3<double>*>(grad);
    double d_abs = g.norm_sq();
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
        if(cfac < 1) {g *= cfac;} /* scalar gradient correction */
    }
}


void construct_gradient(Vec3<MyDouble>& grad, int i)
{
    /* check if the matrix is well-conditioned: otherwise we will use the 'standard SPH-like' derivative estimation */
    if(SHOULD_I_USE_SPH_GRADIENTS(CellP[i].ConditionNumber))
    {
        /* the condition number was bad, so we used SPH-like gradients */
        if(CellP[i].Density > 0) {grad *= P[i].DrkernNgbFactor / CellP[i].Density;}
    } else {
        /* ok, the condition number was good so we used the matrix-like gradient estimator */
        grad = CellP[i].NV_T.matvec(grad);
    }
}




void hydro_gradient_calc(void)
{
    CPU_Step[CPU_DENSMISC] += measure_time(); double t0 = my_second();
    /* Phase 7+ outer-wrapper sub-bucket timing — env-gated; no-op when off. */
    double t_grad_outer_start = my_second();
    /* Neighbor-list path: allocate active-index array, refresh ghosts, build symmetric CSR.
       The CSR is reused by hydro_force. Helper is a no-op on the tree-walk build. */
    double gsl_safety = gizmo_ghost_safety_factor();
    double t_diag_symlist_start = my_second(); /* DIAG */
    gizmo_gradients_prep_symlist(gsl_safety, gsl_safety);
    double t_grad_after_symlist = my_second();
    gizmo_step_phase_record("gradient_prep_symlist", timediff(t_diag_symlist_start, t_grad_after_symlist));
    if(ThisTask == 0 && gizmo_verbose_diag()) {
        printf("[DIAG_SYMNL step=%d N=%d pairs=%lld] symlist_build=%.3f\n",
               (int)All.NumCurrentTiStep, gizmo_sym_num_active, (long long)(gizmo_sym_neighbor_list.total_pairs),
               timediff(t_diag_symlist_start, my_second()));
        fflush(stdout);
    }
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
    All.BunchSize = (long) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct GasGraddata_in) + sizeof(struct GasGraddata_out) +
                                                             sizemax(sizeof(struct GasGraddata_in),sizeof(struct GasGraddata_out))));
    Ngblist.resize(NTaskTimesNumPart);
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    /* before doing any operations, need to zero the appropriate memory so we can correctly do pair-wise operations */
    for (int i : ActiveParticleList)
        if(P[i].Type==0)
        {
            int k2;
            memset(&GasGradDataPasser[i], 0, sizeof(struct temporary_data_topass));
#ifdef HYDRO_SPH
#ifdef MAGNETIC
            CellP[i].DtB = {};
#endif
#ifdef DIVBCLEANING_DEDNER
            CellP[i].divB = 0;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
            CellP[i].alpha_limiter = 0;
#endif
#endif
#ifdef TURB_DIFF_DYNAMIC
            /* NEED Velocity_bar CORRECT HERE */
            CellP[i].Velocity_bar *= All.TurbDynamicDiffSmoothing;
            CellP[i].Velocity_hat = CellP[i].Velocity_bar * smoothInv;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
            P[i].AGS_zeta = 0;
#endif

            /* and zero out the gradients structure itself */
            CellP[i].Gradients.Density = {};
            CellP[i].Gradients.Pressure = {};
            CellP[i].Gradients.Velocity = {};
#ifdef DOGRAD_INTERNAL_ENERGY
            CellP[i].Gradients.InternalEnergy = {};
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
            CellP[i].Gradients.ElectronNumberDensity = {};
            CellP[i].Gradients.ElectronTemperature = {};
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
            CellP[i].Gradients.E_battery_T2 = {};
#endif
#ifdef COSMIC_RAY_FLUID
            for(k2=0;k2<N_CR_PARTICLE_BINS;k2++) {CellP[i].Gradients.CosmicRayPressure[k2] = {};}
#endif
#ifdef DOGRAD_SOUNDSPEED
            CellP[i].Gradients.SoundSpeed = {};
#endif
#ifdef MAGNETIC
#ifndef MHD_CONSTRAINED_GRADIENT
            CellP[i].Gradients.B = {};
#else
            CellP[i].Face_Area = {};
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
            CellP[i].Gradients.Phi = {};
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
            for(k2=0;k2<NUM_METAL_SPECIES;k2++) {CellP[i].Gradients.Metallicity[k2] = {};}
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
            for(k2=0;k2<N_RT_FREQ_BINS;k2++) {CellP[i].Gradients.Rad_E_gamma_ET[k2] = {};}
#endif
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
            for(k2=0;k2<N_RT_FREQ_BINS;k2++) {
                int k3; for(k3=0;k3<3;k3++) {CellP[i].Gradients.Rad_E_gamma_Grad[k2][k3] = 0;}
                int k4; for(k4=0;k4<3;k4++) {int k5; for(k5=0;k5<3;k5++) {CellP[i].Gradients.Rad_Flux_Grad[k2][k4][k5] = 0;}}
            }
#endif
        }

    /* CORRECTNESS FIX (caught by GIZMO_GPU_ARENA_DEBUG=1, 2026-05-03):
     * The zero-out loop above mutates host P[i].AGS_zeta and CellP[i].Gradients.*
     * for active gas particles. The GPU arena (populated by the symlist memcpy
     * earlier) is NOT updated, so a subsequent fast-path acquire by
     * gradient_evaluate_gpu would read stale values. Invalidate arena here so
     * the next acquire forces a slow-path memcpy and re-seeds from host.
     * Phase 8a Round 2 will replace this with a mirror-update + mark_clean
     * (write zeros to BOTH host and arena) to recover the slow-path memcpy. */
    gpu_particles_arena_invalidate();

    /* prepare to do the requisite number of sweeps over the particle distribution */
    int gradient_iteration;
    for(gradient_iteration = 0; gradient_iteration < NUMBER_OF_GRADIENT_ITERATIONS; gradient_iteration++)
    {
        // need to zero things used in the iteration (anything appearing in out2particle_GasGrad_iter)
        for (int i : ActiveParticleList)
            if(P[i].Type==0)
            {
#ifdef MHD_CONSTRAINED_GRADIENT
                GasGradDataPasser[i].FaceDotB = 0;
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
                GasGradDataPasser[i].PhiGrad = {};
#endif
#endif
            }

        /* Neighbor-list path: use cached symmetric CSR list with GPU/Kokkos dispatch.
           Works for all gradient iterations (0 and >0 for MHD_CONSTRAINED_GRADIENT). */
        {
            struct GasGraddata_out *grad_out = (struct GasGraddata_out *) mymalloc("grad_out",
                (gizmo_sym_num_active > 0 ? gizmo_sym_num_active : 1) * sizeof(struct GasGraddata_out));

            gradient_evaluate_gpu(P, CellP, NumPart,
                                  gizmo_sym_active_indices, gizmo_sym_num_active,
                                  gizmo_sym_neighbor_list.offsets,
                                  gizmo_sym_neighbor_list.neighbors,
                                  gizmo_sym_neighbor_list.total_pairs,
                                  (void *)grad_out, gradient_iteration);

            if(gradient_iteration == 0) {
                for(int aa = 0; aa < gizmo_sym_num_active; aa++) {
                    int ii = gizmo_sym_active_indices[aa];
                    out2particle_GasGrad(&grad_out[aa], ii, 0, gradient_iteration);
                }
            } else {
                /* For iteration >0: extract only FaceDotB (and PhiGrad if MIDPOINT) from the
                   full output struct into the compact iter struct for scatter. */
                for(int aa = 0; aa < gizmo_sym_num_active; aa++) {
                    int ii = gizmo_sym_active_indices[aa];
                    struct GasGraddata_out_iter out_iter;
                    memset(&out_iter, 0, sizeof(out_iter));
#ifdef MHD_CONSTRAINED_GRADIENT
                    out_iter.FaceDotB = grad_out[aa].FaceDotB;
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
                    out_iter.PhiGrad[0] = grad_out[aa].Gradients[0].Phi;
                    out_iter.PhiGrad[1] = grad_out[aa].Gradients[1].Phi;
                    out_iter.PhiGrad[2] = grad_out[aa].Gradients[2].Phi;
#endif
#endif
                    out2particle_GasGrad_iter(&out_iter, ii, 0, gradient_iteration);
                }
            }
            myfree(grad_out);
        }


        /* here, we insert intermediate operations on the results, from the iterations we have completed */
#ifdef MHD_CONSTRAINED_GRADIENT
        for (int i : ActiveParticleList)
            if(P[i].Type == 0)
            {
                CellP[i].FlagForConstrainedGradients = 1;
                /* copy everything from the structure holding B-gradients (needed so they dont change mid-loop) */
                for(k=0;k<3;k++)
                {
                    for(k1=0;k1<3;k1++)
                    {
                        CellP[i].Gradients.B[k][k1] = GasGradDataPasser[i].BGrad[k][k1];
                    }
                }
                /* build the gradient */
                for(k=0;k<3;k++) {construct_gradient(CellP[i].Gradients.B[k],i);}
                /* slope limit it */
                double v_tmp = P[i].Mass / CellP[i].Density;
                double tmp_d = sqrt(1.0e-37 + (2. * All.cf_atime/ CellP[i].Pressure*v_tmp*v_tmp) + CellP[i].BPred.norm_sq());
                double tmp = 3.0e3 * fabs(CellP[i].divB) * P[i].KernelRadius / tmp_d;
                double alim; alim = 1. + DMIN(1.,tmp*tmp);
#if (MHD_CONSTRAINED_GRADIENT <= 1)
                double dbmax=0, dbgrad=0;
                double dh=0.25*P[i].KernelRadius; // need to be more aggressive with new wt_i,wt_j formalism
                for(k=0;k<3;k++)
                {
                    double b0 = CellP[i].Bfield_component(k);
                    double dd = 2. * fabs(b0) * DMIN(fabs(GasGradDataPasser[i].Minima.B[k]) , fabs(GasGradDataPasser[i].Maxima.B[k]));
                    dbmax = DMIN(fabs(dbmax+dd),fabs(dbmax-dd));
                    for(k1=0;k1<3;k1++) {dbgrad += 2.*dh * fabs(b0*CellP[i].Gradients.B[k][k1]);}
                }
                dbmax /= dbgrad;
                for(k1=0;k1<3;k1++)
                {
                    double d_abs = CellP[i].Gradients.B[k1].norm_sq();
                    if(d_abs > 0)
                    {
                        double cfac = 1 / (0.25 * P[i].KernelRadius * sqrt(d_abs));
                        cfac *= DMIN(fabs(GasGradDataPasser[i].Maxima.B[k1]) , fabs(GasGradDataPasser[i].Minima.B[k1]));
                        double c_eff = DMIN( cfac , DMAX(cfac/alim , dbmax) );
                        if(c_eff < 1) {CellP[i].Gradients.B[k1] *= c_eff;}
                    } else {
                        CellP[i].Gradients.B[k1] = {};
                    }
                }
#endif
                /* check the particle area closure, which will inform whether it is safe to use the constrained gradients */
                double area = fabs(CellP[i].Face_Area[0]) + fabs(CellP[i].Face_Area[1]) + fabs(CellP[i].Face_Area[2]);
                area /= Get_Particle_Expected_Area(P[i].KernelRadius);
                /* set the relevant flags to decide whether or not we use the constrained gradients */
                if(area > 0.5) {CellP[i].FlagForConstrainedGradients = 0;}
                if(CellP[i].ConditionNumber > 1000.) {CellP[i].FlagForConstrainedGradients = 0;}
                if(SHOULD_I_USE_SPH_GRADIENTS(CellP[i].ConditionNumber)) {CellP[i].FlagForConstrainedGradients = 0;} /* this must be here, since in this case the SPH gradients are used, which will not work with this method */

                /* now check, and if ok, enter the gradient re-calculation */
                if(CellP[i].FlagForConstrainedGradients == 1)
                {
                {
                    /* When MHD_MODIFIED_GRADIENT is enabled but the MG global solve was skipped this
                       timestep (small active fraction), fall back to the standard CG iterative correction
                       for active cells. When MG did run, skip CG (the global solve handles div(B)). */
                    int do_cg_correction = 1; /* default: always do CG when MHD_MODIFIED_GRADIENT is not defined */
#ifdef MHD_MODIFIED_GRADIENT
                    do_cg_correction = All.Flag_SkipMGSolve; /* 1 = MG was skipped, so do CG; 0 = MG ran, skip CG */
#endif
                    if(do_cg_correction)
                    {
                        double GB0[3][3];
                        double fsum = 0.0, dmag = 0.0;
                        double h_eff = P[i].Get_Particle_Size();
                        for(k=0;k<3;k++)
                        {
                            double grad_limiter_mag = CellP[i].Bfield_component(k) / h_eff;
                            dmag += grad_limiter_mag * grad_limiter_mag;
                            for(k1=0;k1<3;k1++)
                            {
                                GB0[k][k1] = CellP[i].Gradients.B[k][k1];
                                dmag += GB0[k][k1] * GB0[k][k1];
                                fsum += GasGradDataPasser[i].FaceCrossX[k][k1] * GasGradDataPasser[i].FaceCrossX[k][k1];
                            }
                        }
                        if((fsum <= 0) || (dmag <= 0))
                        {
                            CellP[i].FlagForConstrainedGradients = 0;
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
                                        asum += CellP[i].Gradients.B[k][k1] * GasGradDataPasser[i].FaceCrossX[k][k1];
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
                                        double grad_limiter_mag = (CellP[i].Gradients.B[k][k1] + ecorr[k][k1]) - GB0[k][k1];
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
                                        CellP[i].Gradients.B[k][k1] = GB0[k][k1] + nnorm*(CellP[i].Gradients.B[k][k1]+ecorr[k][k1] - GB0[k][k1]);
                                    }
                                    /* slope-limit the corrected gradients again, but with a more tolerant slope-limiter */
#if (MHD_CONSTRAINED_GRADIENT <= 1)
                                    local_slopelimiter(CellP[i].Gradients.B[k],GasGradDataPasser[i].Maxima.B[k],GasGradDataPasser[i].Minima.B[k],0.25, P[i].KernelRadius, 0.25, 0, 0, 0);
#endif
                                }
                            } // closes j_gloop loop
                        } // closes fsum/dmag check
                    } // closes do_cg_correction check
                } // closes inner block
                } // closes FlagForConstrainedGradients check
#ifdef MHD_CONSTRAINED_GRADIENT_MIDPOINT
                double a_limiter = 0.25; if(CellP[i].ConditionNumber>100) a_limiter=DMIN(0.5, 0.25 + 0.25 * (CellP[i].ConditionNumber-100)/100);
                /* copy everything from the structure holding phi-gradients (needed so they dont change mid-loop) */
                CellP[i].Gradients.Phi = GasGradDataPasser[i].PhiGrad;
                /* build and limit the gradient */
                construct_gradient(CellP[i].Gradients.Phi,i);
                local_slopelimiter(CellP[i].Gradients.Phi,GasGradDataPasser[i].Maxima.Phi,GasGradDataPasser[i].Minima.Phi,a_limiter,P[i].KernelRadius,0.0, 0, 0, 0);
#endif
            } // closes Ptype == 0 check
#endif
    } // closes gradient_iteration

    myfree(DataNodeList);
    myfree(DataIndexTable);
    


    /* do final operations on results: these are operations that can be done after the complete set of iterations */
    for (int i : ActiveParticleList)
        if(P[i].Type == 0)
        {
            /* now we can properly calculate (second-order accurate) gradients of hydrodynamic quantities from this loop */
            construct_gradient(CellP[i].Gradients.Density,i);
            construct_gradient(CellP[i].Gradients.Pressure,i);
            for(k=0;k<3;k++) {construct_gradient(CellP[i].Gradients.Velocity[k],i);}
#ifdef TURB_DIFF_DYNAMIC
            for(k=0;k<3;k++) {construct_gradient(GasGradDataPasser[i].GradVelocity_bar[k], i);}
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
            construct_gradient(CellP[i].Gradients.InternalEnergy,i);
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
            construct_gradient(CellP[i].Gradients.ElectronNumberDensity,i);
            construct_gradient(CellP[i].Gradients.ElectronTemperature,i);
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
            for(k=0;k<3;k++) {construct_gradient(CellP[i].Gradients.E_battery_T2[k],i);}
#endif
#ifdef COSMIC_RAY_FLUID
            for(k=0;k<N_CR_PARTICLE_BINS;k++) {construct_gradient(CellP[i].Gradients.CosmicRayPressure[k],i);}
            int is_particle_local_extremum[N_CR_PARTICLE_BINS]={0}; is_particle_local_extremum[0]=0; // test for local extremum to revert to lower-order reconstruction if necessary
#endif
#ifdef DOGRAD_SOUNDSPEED
            construct_gradient(CellP[i].Gradients.SoundSpeed,i);
#endif
#ifdef MAGNETIC
#ifndef MHD_CONSTRAINED_GRADIENT
            for(k=0;k<3;k++) {construct_gradient(CellP[i].Gradients.B[k],i);}
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
            construct_gradient(CellP[i].Gradients.Phi,i);
#endif
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
            for(k=0;k<NUM_METAL_SPECIES;k++) {construct_gradient(CellP[i].Gradients.Metallicity[k],i);}
#endif
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
            for(k=0;k<N_RT_FREQ_BINS;k++) {construct_gradient(GasGradDataPasser[i].Gradients_Rad_E_gamma[k],i);}
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
            {int k_f; for(k_f=0;k_f<N_RT_FREQ_BINS;k_f++) {int k_d; for(k_d=0;k_d<3;k_d++) {construct_gradient(CellP[i].Gradients.Rad_Flux_Grad[k_f][k_d],i);}}}
#endif
#endif

            /* now the gradients are calculated: below are simply useful operations on the results */
#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
            /* this is here because for the models of sink growth and self-shielding of stars, we need to calculate GradRho: we don't bother doing it in density.c if we're already calculating it here! but note, this is the -un-limited- gradient here */
            P[i].GradRho = CellP[i].Gradients.Density;
#endif

#if defined(TURB_DRIVING) || defined(OUTPUT_VORTICITY)
            CellP[i].Vorticity = CellP[i].Gradients.Velocity.curl();
#endif

#ifdef SPH_TP12_ARTIFICIAL_RESISTIVITY
            /* use the magnitude of the B-field gradients relative to kernel length to calculate artificial resistivity */
            double GradBMag = CellP[i].Gradients.B.frobenius_norm_sq(); double rho_over_m = CellP[i].Density / P[i].Mass; double BMag = CellP[i].BPred.norm_sq() * rho_over_m * rho_over_m;
            CellP[i].Balpha = DMAX(DMIN(P[i].KernelRadius * sqrt(GradBMag/(BMag+1.0e-33)), 0.1 * All.ArtMagDispConst), 0.005);
#endif
            
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
            /* note non-gas particles are handled separately, in the ags_rkern routine. here the zeta terms ONLY control errors if we maintain the 'correct' neighbor number: for boundary particles, it can actually be worse. so we need to check whether we should use it or not */
            double ngb_eff = pow(P[i].NumNgb, NUMDIMS); // calculate the actual neighbor number, needed here
            if((fabs(ngb_eff-All.DesNumNgb)/All.DesNumNgb < 0.05) && (P[i].KernelRadius > 1.001*All.MinKernelRadius) && (P[i].KernelRadius < 0.999*All.MaxKernelRadius)) {
                double ndenNGB = ngb_eff / ( VOLUME_NORM_COEFF_FOR_NDIMS * pow(P[i].KernelRadius,NUMDIMS) ); P[i].AGS_zeta *= P[i].Mass * P[i].KernelRadius / (NUMDIMS * ndenNGB) * P[i].DrkernNgbFactor;
            } else {P[i].AGS_zeta=0;}
#endif


#ifdef HYDRO_SPH
#ifdef MAGNETIC
            if(CellP[i].Density > 0)
            {
                CellP[i].DtB *= P[i].DrkernNgbFactor * P[i].Mass / (CellP[i].Density * CellP[i].Density) / All.cf_atime; // induction equation (convert from Bcode*vcode/rcode to Bphy/tphys) //
#ifdef DIVBCLEANING_DEDNER
                /* full correct form of D(phi)/Dt = -ch*ch*div.dot.B - phi/tau - (1/2)*phi*div.dot.v */ /* PFH: here's the div.dot.B term: make sure div.dot.B def'n matches appropriate grad_phi conjugate pair: recommend direct diff div.dot.B */
                CellP[i].divB *= P[i].DrkernNgbFactor * P[i].Mass / (CellP[i].Density * CellP[i].Density);
                if((!isnan(CellP[i].divB))&&(P[i].KernelRadius>0)&&(CellP[i].divB!=0)&&(CellP[i].Density>0)) {
                    double tmp_ded = 0.5 * CellP[i].MaxSignalVel ; /* has units of v_physical now *//* do a check to make sure divB isn't something wildly divergent (owing to particles being too close) */
                    double rho_over_m_divb = CellP[i].Density / P[i].Mass; double b2_max = CellP[i].BPred.norm_sq() * rho_over_m_divb * rho_over_m_divb;
                    b2_max = 100.0 * fabs( sqrt(b2_max) * All.cf_a2inv * P[i].Mass / (CellP[i].Density*All.cf_a3inv) * 1.0 / (P[i].KernelRadius*All.cf_atime) );
                    if(fabs(CellP[i].divB) > b2_max) {CellP[i].divB *= b2_max / fabs(CellP[i].divB);} /* ok now can apply this to get the growth rate of phi */
                    CellP[i].DtPhi = -tmp_ded * tmp_ded * All.DivBcleanHyperbolicSigma * CellP[i].divB * CellP[i].Density*All.cf_a3inv; // mass-based phi-flux
                    // phiphi above now has units of [Bcode]*[vcode]^2/[rcode]=(Bcode*vcode)*vcode/rcode; needs to have units of [Phicode]*[vcode]/[rcode] so [PhiGrad]=[Phicode]/[rcode] = [DtB] = [Bcode]*[vcode]/[rcode] IFF [Phicode]=[Bcode]*[vcode]; this also makes the above self-consistent //
                    // (implicitly, this gives the correct evolution in comoving, adiabatic coordinates where the sound speed is the relevant speed at which the 'damping wave' propagates. another choice (provided everything else is self-consistent) is fine, it just makes different assumptions about the relevant 'desired' timescale for damping wave propagation in the expanding box) //
                } else {CellP[i].DtPhi=0; CellP[i].divB=0; CellP[i].DtB = {};}
                CellP[i].divB = 0.0; /* now we re-zero it, since a -different- divB definition must be used in hydro to subtract the tensile terms */
#endif
            } else {
                CellP[i].DtB = {};
#ifdef DIVBCLEANING_DEDNER
                CellP[i].divB = 0; CellP[i].DtPhi = 0;
#endif
            }
#endif

#ifdef SPHAV_CD10_VISCOSITY_SWITCH
            CellP[i].alpha_limiter /= CellP[i].Density;
            NV_dt =  get_particle_timestep_in_physical(i); // physical
            NV_dummy = fabs(1.0 * pow(1.0 - CellP[i].alpha_limiter,4.0) * CellP[i].NV_DivVel); // NV_ quantities are in physical units
            NV_limiter = NV_dummy*NV_dummy / (NV_dummy*NV_dummy + CellP[i].NV_trSSt);
            NV_A = DMAX(-CellP[i].NV_dt_DivVel, 0.0);
            divVel_physical = CellP[i].NV_DivVel;
            // add a simple limiter here: alpha_loc is 'prepped' but only switches on when the divergence goes negative: want to add hubble flow here //
            if(All.ComovingIntegrationOn) {divVel_physical += 3*All.cf_hubble_a;} // hubble-flow correction added
            if(divVel_physical>=0.0) {NV_A = 0.0;}
            h_eff = P[i].Get_Particle_Size() * All.cf_atime / 0.5; // 'default' parameter choices are scaled for a cubic spline, but code will attempt to scale appropriately to other kernel choices //
            cs_nv = CellP[i].effective_soundspeed() ; // converts to physical velocity units //
            alphaloc = All.ViscosityAMax * h_eff*h_eff*NV_A / (0.36*cs_nv*cs_nv + h_eff*h_eff*NV_A);
            // 0.25 in front of vsig is the 'noise parameter' that determines the relative amplitude which will trigger the switch: that choice was quite large (requires approach velocity rate-of-change is super-sonic); better to use c_s (above), and 0.05-0.25 //
            // NV_A is physical 1/(time*time), but KernelRadius and vsig can be comoving, so need appropriate correction terms above //
            if(CellP[i].alpha < alphaloc) {CellP[i].alpha = alphaloc;}
                else if (CellP[i].alpha > alphaloc) {CellP[i].alpha = alphaloc + (CellP[i].alpha - alphaloc) * exp(-NV_dt * (0.5*fabs(CellP[i].MaxSignalVel))/(0.5*h_eff) * 0.05);}
            if(CellP[i].alpha < All.ViscosityAMin) {CellP[i].alpha = All.ViscosityAMin;}
            CellP[i].alpha_limiter = DMAX(NV_limiter,All.ViscosityAMin/CellP[i].alpha);
#else
            /* compute the traditional Balsara limiter (now that we have velocity gradients) */
            double divVel = All.cf_a2inv * fabs(CellP[i].Gradients.Velocity.trace());
            if(All.ComovingIntegrationOn) {divVel += 3*All.cf_hubble_a;} // hubble-flow correction added (physical units)
            Vec3<double> CurlVel; double MagCurl;
            CurlVel = CellP[i].Gradients.Velocity.curl();
            MagCurl = All.cf_a2inv * CurlVel.norm();
            double fac_mu = 1 / ( All.cf_atime);
            CellP[i].alpha_limiter = divVel / (divVel + MagCurl + 0.0001 * CellP[i].effective_soundspeed() / (P[i].Get_Particle_Size()) / fac_mu);
#endif
#endif

            calculate_and_assign_conduction_and_viscosity_coefficients(i, P, CellP);
            calculate_and_assign_nonideal_mhd_coefficients(i, P, CellP);
            
#ifdef RADTRANSFER
            {
                int k_freq; for(k_freq = 0; k_freq < N_RT_FREQ_BINS; k_freq++)
                {
                    /* calculate the opacity */
                    CellP[i].Rad_Kappa[k_freq] = rt_kappa(i,k_freq, P, CellP); // physical units //
#if defined(RT_FLUXLIMITER) && defined(RT_COMPGRAD_EDDINGTON_TENSOR)
                    /* compute the flux-limiter for radiation transport: also convenient here to compute the relevant opacities for all particles */
                    double lambda = 1;
                    if(CellP[i].Rad_E_gamma_Pred[k_freq] > 0) /* can compute gradient length scale */
                    {
                        double R_ET = CellP[i].Gradients.Rad_E_gamma_ET[k_freq].norm() / (MIN_REAL_NUMBER + CellP[i].Rad_E_gamma_Pred[k_freq] * CellP[i].Density/(MIN_REAL_NUMBER+P[i].Mass));
                        R_ET = 3.*DMAX(R_ET , 1.e-6/P[i].Get_Particle_Size()) / (1.e-55 + All.cf_atime*CellP[i].Rad_Kappa[k_freq]*(CellP[i].Density*All.cf_a3inv)); // limit to be > 0, divide by kappa-rho to get desired dimensionless ratio
                        lambda = DMIN(1., DMAX( 3.*(2. + R_ET) / (6. + 3.*R_ET + R_ET*R_ET), MIN_REAL_NUMBER )); // slope-limiter
#ifdef RT_OTVET         /* note that the OTVET eddington tensor is close to the correct value for the optically-thin limit. for the diffusion limit
                            it may be incorrect. we can therefore interpolate using an M1-like relation below, based on the gradients above (used
                            to determine which limit we are actually in: ratio f=|flux|/(c_eff*Energy_density_rad): f<<1 = diffusion limit, f~1 = free-streaming limit: this is our slope-limiter above */
                        double chi=DMAX(1./3.,DMIN(1.,(3.+4.*lambda*lambda)/(5.+2.*sqrt(4.-3.*lambda*lambda)))), chifac_iso=3.*(1-chi)/2., chifac_ot=(3.*chi-1.)/2.;
                        CellP[i].Gradients.Rad_E_gamma_ET[k_freq] = chifac_ot*CellP[i].Gradients.Rad_E_gamma_ET[k_freq] + (chifac_iso/3.)*GasGradDataPasser[i].Gradients_Rad_E_gamma[k_freq];
#endif // ifdef otvet
                    }
                    CellP[i].Rad_Flux_Limiter[k_freq] = lambda;
#endif // ifdef fluxlimiter

#if defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_OTVET)
                    {Vec3<MyDouble> g{GasGradDataPasser[i].Gradients_Rad_E_gamma[k_freq][0], GasGradDataPasser[i].Gradients_Rad_E_gamma[k_freq][1], GasGradDataPasser[i].Gradients_Rad_E_gamma[k_freq][2]};
                    CellP[i].Gradients.Rad_E_gamma_ET[k_freq] = CellP[i].ET[k_freq].matvec(g);} /* set the output gradient grad.(D*Prad) = D.(grad Prad) */
#endif
#if defined(GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION) /* yet another hack for this problem to get the boundaries to play nicely once dust evacuated -- this is a bit redundant with other hacks, but here for safety */
                    if(CellP[i].Interpolated_Opacity[0] < 1.e-3 * All.Dust_to_Gas_Mass_Ratio*0.75*All.Grain_Q_at_MaxGrainSize/((All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS)*(All.Grain_Size_Max/UNIT_LENGTH_IN_CGS))) {double gmax=-1; if(P[i].GravAccel[GRAV_DIRECTION_RDI] < gmax) {P[i].GravAccel[GRAV_DIRECTION_RDI]=gmax;}} // the interpolated opacity here should be in code units by default
#endif
                }
            }
#endif // ifdef radtransfer

            
#if defined(EOS_ELASTIC) // update time-derivative of stress tensor (needs to be done before slope-limiting to use full velocity gradient information) //
            elastic_body_update_driftkick(i,1.,2);
#endif

            
            /* finally, we need to apply a sensible slope limiter to the gradients, to prevent overshooting */
            double stol = 0.0, stol_tmp, stol_diffusion; stol_diffusion = 0.1; stol_tmp = stol;
            double h_lim = P[i].KernelRadius, d_max = DMAX(P[i].KernelRadius,GasGradDataPasser[i].MaxDistance); h_lim = d_max;
            /* fraction of H at which maximum reconstruction is allowed (=0.5 for 'standard'); for pure hydro we can
             be a little more aggresive and the equations are still stable (but this is as far as you want to push it) */
            double a_limiter = 0.25; if(CellP[i].ConditionNumber>100) a_limiter=DMIN(0.5, 0.25 + 0.25 * (CellP[i].ConditionNumber-100)/100);
#if defined(SELFGRAVITY_OFF) && (!defined(MAGNETIC) && !defined(GALSF))
            h_lim=P[i].KernelRadius; stol=0.1;
#endif
#if (SLOPE_LIMITER_TOLERANCE == 2)
            h_lim = P[i].KernelRadius; a_limiter *= 0.5; stol = 0.125;
#endif
#if (SLOPE_LIMITER_TOLERANCE == 0)
            a_limiter *= 2.0; stol = 0.0;
#endif

#if (SINGLE_STAR_SINK_FORMATION & 4)
            CellP[i].Density_Relative_Maximum_in_Kernel = GasGradDataPasser[i].Maxima.Density;
#endif
            local_slopelimiter(CellP[i].Gradients.Density,GasGradDataPasser[i].Maxima.Density,GasGradDataPasser[i].Minima.Density,a_limiter,h_lim,0, 1,d_max,CellP[i].Density);
            int pressure_is_positive_definite = 1;
#if defined(EOS_TILLOTSON) || defined(EOS_ELASTIC) || defined(EOS_ANEOS)
            pressure_is_positive_definite = 0; /* some physics allow negative pressures - account for that here */
#endif
            local_slopelimiter(CellP[i].Gradients.Pressure,GasGradDataPasser[i].Maxima.Pressure,GasGradDataPasser[i].Minima.Pressure,a_limiter,h_lim,stol, pressure_is_positive_definite,d_max,CellP[i].Pressure);
            stol_tmp = stol;
#if defined(VISCOSITY)
            stol_tmp = DMAX(stol,stol_diffusion);
#endif
#ifdef TURB_DIFF_DYNAMIC
            for (k1=0;k1<3;k1++) {local_slopelimiter(GasGradDataPasser[i].GradVelocity_bar[k1], GasGradDataPasser[i].Maxima.Velocity_bar[k1], GasGradDataPasser[i].Minima.Velocity_bar[k1], a_limiter, h_lim, stol, 0,0,0);}
#endif
            for(k1=0;k1<3;k1++) {local_slopelimiter(CellP[i].Gradients.Velocity[k1],GasGradDataPasser[i].Maxima.Velocity[k1],GasGradDataPasser[i].Minima.Velocity[k1],a_limiter,h_lim,stol_tmp, 0,0,0);}
#ifdef DOGRAD_INTERNAL_ENERGY
            stol_tmp = stol;
#if defined(CONDUCTION)
            stol_tmp = DMAX(stol,stol_diffusion);
#endif
            local_slopelimiter(CellP[i].Gradients.InternalEnergy,GasGradDataPasser[i].Maxima.InternalEnergy,GasGradDataPasser[i].Minima.InternalEnergy,a_limiter,h_lim,stol_tmp, 1,d_max,CellP[i].InternalEnergyPred);
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 1)
            local_slopelimiter(CellP[i].Gradients.ElectronNumberDensity,GasGradDataPasser[i].Maxima.ElectronNumberDensity,GasGradDataPasser[i].Minima.ElectronNumberDensity,a_limiter,h_lim,stol, 1,d_max,CellP[i].n_e());
            local_slopelimiter(CellP[i].Gradients.ElectronTemperature,GasGradDataPasser[i].Maxima.ElectronTemperature,GasGradDataPasser[i].Minima.ElectronTemperature,a_limiter,h_lim,stol, 1,d_max,CellP[i].T_e());
#endif
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & (2|4|8))
            for(k1=0;k1<3;k1++) {local_slopelimiter(CellP[i].Gradients.E_battery_T2[k1],GasGradDataPasser[i].Maxima.E_battery_T2[k1],GasGradDataPasser[i].Minima.E_battery_T2[k1],a_limiter,h_lim,stol, 0,0,0);}
#endif
#ifdef DOGRAD_SOUNDSPEED
            local_slopelimiter(CellP[i].Gradients.SoundSpeed,GasGradDataPasser[i].Maxima.SoundSpeed,GasGradDataPasser[i].Minima.SoundSpeed,a_limiter,h_lim,stol, 1,d_max,CellP[i].effective_soundspeed());
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
            for(k1=0;k1<NUM_METAL_SPECIES;k1++) {local_slopelimiter(CellP[i].Gradients.Metallicity[k1],GasGradDataPasser[i].Maxima.Metallicity[k1],GasGradDataPasser[i].Minima.Metallicity[k1],a_limiter,h_lim,DMAX(stol,stol_diffusion), 1,d_max,P[i].Metallicity[k1]);}
#endif
#if defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_EVOLVE_FLUX)
            for(k1=0;k1<N_RT_FREQ_BINS;k1++)
            {
                local_slopelimiter(CellP[i].Gradients.Rad_E_gamma_ET[k1],GasGradDataPasser[i].Maxima.Rad_E_gamma[k1],GasGradDataPasser[i].Minima.Rad_E_gamma[k1],a_limiter,h_lim,stol, 1,d_max,CellP[i].Rad_E_gamma_Pred[k1]*CellP[i].Density/P[i].Mass);
                local_slopelimiter(GasGradDataPasser[i].Gradients_Rad_E_gamma[k1],GasGradDataPasser[i].Maxima.Rad_E_gamma[k1],GasGradDataPasser[i].Minima.Rad_E_gamma[k1],a_limiter,h_lim,DMAX(stol,stol_diffusion), 1,d_max,CellP[i].Rad_E_gamma_Pred[k1]*CellP[i].Density/P[i].Mass);
            }
#endif
#if defined(RT_M1_SECONDORDER) && defined(RT_EVOLVE_FLUX)
            {
                double V_i_inv_rt = CellP[i].Density / P[i].Mass;
                for(k1=0;k1<N_RT_FREQ_BINS;k1++)
                {
                    double val_cen_e = CellP[i].Rad_E_gamma_Pred[k1] * V_i_inv_rt;
                    local_slopelimiter(GasGradDataPasser[i].Gradients_Rad_E_gamma[k1],GasGradDataPasser[i].Maxima.Rad_E_gamma[k1],GasGradDataPasser[i].Minima.Rad_E_gamma[k1],a_limiter,h_lim,stol, 1,d_max,val_cen_e);
                    for(k=0;k<3;k++) {CellP[i].Gradients.Rad_E_gamma_Grad[k1][k] = GasGradDataPasser[i].Gradients_Rad_E_gamma[k1][k];}
                    int k_d; for(k_d=0;k_d<3;k_d++) {
                        double val_cen_f = CellP[i].Rad_Flux_Pred[k1][k_d] * V_i_inv_rt;
                        local_slopelimiter(CellP[i].Gradients.Rad_Flux_Grad[k1][k_d],GasGradDataPasser[i].Maxima.Rad_Flux[k1][k_d],GasGradDataPasser[i].Minima.Rad_Flux[k1][k_d],a_limiter,h_lim,stol, 0,d_max,val_cen_f);
                    }
                }
            }
#endif
#ifdef MAGNETIC
#ifndef MHD_CONSTRAINED_GRADIENT
            double v_tmp = P[i].Mass / CellP[i].Density;
            double tmp_d = sqrt(1.0e-37 + (2. * All.cf_atime/ CellP[i].Pressure*v_tmp*v_tmp) +
                                CellP[i].BPred.norm_sq());
            double q = fabs(CellP[i].divB) * P[i].KernelRadius / tmp_d, alim2 = a_limiter * (1. + q*q); if(alim2 > 0.5) alim2=0.5;
            stol_tmp = stol;
#ifdef MHD_NON_IDEAL
            stol_tmp = DMAX(stol,stol_diffusion);
#endif
            for(k1=0;k1<3;k1++) {local_slopelimiter(CellP[i].Gradients.B[k1],GasGradDataPasser[i].Maxima.B[k1],GasGradDataPasser[i].Minima.B[k1],alim2,h_lim,stol_tmp, 0,0,0);}
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(MHD_CONSTRAINED_GRADIENT_MIDPOINT)
            local_slopelimiter(CellP[i].Gradients.Phi,GasGradDataPasser[i].Maxima.Phi,GasGradDataPasser[i].Minima.Phi,a_limiter,h_lim,stol, 0,0,0);
#endif
#endif



#ifdef TURB_DIFFUSION
#ifdef TURB_DIFF_DYNAMIC
            {int k1,k2; for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {CellP[i].VelShear_bar[k1][k2] = 0.5 * (GasGradDataPasser[i].GradVelocity_bar[k1][k2] + GasGradDataPasser[i].GradVelocity_bar[k2][k1]);}}} // need to initialize this before sending to routine below
#endif
            calculate_and_assign_turbulent_diffusion_coefficients(i, P, CellP);
#endif


#if defined(COSMIC_RAY_FLUID) && !defined(CRFLUID_EVOLVE_SCATTERINGWAVES) /* note that because of the way this depends on the gradient scale-length, we should calculate it -after- the slope-limiters are applied */
            for(k=0;k<N_CR_PARTICLE_BINS;k++) {CellP[i].CosmicRayDiffusionCoeff[k]=0;}
            if(CellP[i].Density > 0 && P[i].Mass > 0) {CalculateAndAssign_CosmicRay_DiffusionAndStreamingCoefficients(i, P, CellP);}/* only assign diffusivities to 'valid' gas particles */
#endif


#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==6)
            /* if the mesh motion is specified to be glass-generating, this is where we apply the appropriate mesh velocity */
            if(All.Time > 0)
            {
                double cs_invelunits = CellP[i].effective_soundspeed()  * All.cf_atime; // soundspeed, converted to units of code velocity
                double L_i_code = P[i].Get_Particle_Size(); // particle effective size (in code units)
                Vec3<double> dvel = L_i_code*L_i_code*GasGradDataPasser[i].GlassAcc; double velnorm = dvel.norm(); // calculate quantities to use for glass
                double dtx = get_particle_timestep_in_physical(i); // need timestep for limiter below
                if(velnorm > 0 && dtx > 0)
                {
                    double v00 = 0.5 * DMIN(cs_invelunits*(0.5*velnorm) , All.CourantFac*(L_i_code/dtx)/All.cf_a2inv); // limit added velocity of mesh-generating point to Courant factor
                    CellP[i].ParticleVel += v00 * (dvel/velnorm); // actually add the correction velocity to the mesh velocity
                }
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
                CellP[i].Tensor_CRK_Face_Corrections[0] = A;
                for(k=0;k<3;k++) {CellP[i].Tensor_CRK_Face_Corrections[1+k] = B[k];}
                for(k=0;k<3;k++) {CellP[i].Tensor_CRK_Face_Corrections[1+3+k] = vector_corr[k];}
                for(k=0;k<3;k++) {for(k_x=0;k_x<3;k_x++) {CellP[i].Tensor_CRK_Face_Corrections[1+3+3+3*k+k_x] = tensor_corr[k][k_x];}}
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

    double t_grad_before_refresh = my_second();
    /* Neighbor-list path: refresh ghosts so hydro_force sees converged gradients on both
       sides of each pair, and rebuild the symmetric CSR. No-op on tree-walk build. */
    gizmo_gradients_refresh_symlist(gsl_safety, gsl_safety);
    double t_grad_outer_end = my_second();
    /* Phase 7+ outer-wrapper sub-buckets — env-gated. The bucket
     * "hydro_gradient_calc" is the timed top-level (recorded in core/accel.cc).
     * These split out the major phases of THIS function, complementing the
     * inner gradient_evaluate_gpu sub-buckets (gradient_arena, gradient_kernel,
     * etc.) that are recorded inside the GPU evaluator. */
    gizmo_step_phase_record("gradient_zero_iter_loops", timediff(t_grad_after_symlist, t_grad_before_refresh));
    gizmo_step_phase_record("gradient_refresh_symlist", timediff(t_grad_before_refresh, t_grad_outer_end));
    gizmo_step_phase_record("gradient_outer_total",     timediff(t_grad_outer_start,    t_grad_outer_end));
}

