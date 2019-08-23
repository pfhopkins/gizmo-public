/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/*! This function is the 'core' of the hydro force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*   In this routine, we find the gas particle neighbors, and do the loop over 
*  neighbors to calculate the hydro fluxes. The actual flux calculation, 
*  and the returned values, should be in PHYSICAL (not comoving) units */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
int hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, k, n, startnode, numngb, kernel_mode, listindex;
    double hinv_i,hinv3_i,hinv4_i,hinv_j,hinv3_j,hinv4_j,V_i,V_j,dt_hydrostep,r2,rinv,rinv_soft,u,Particle_Size_i;
    double v_hll,k_hll,b_hll; v_hll=k_hll=0,b_hll=1;
    struct kernel_hydra kernel;
    struct hydrodata_in local;
    struct hydrodata_out out;
    struct Conserved_var_Riemann Fluxes;
    listindex = 0;
    memset(&out, 0, sizeof(struct hydrodata_out));
    memset(&kernel, 0, sizeof(struct kernel_hydra));
    memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifndef HYDRO_SPH
    struct Input_vec_Riemann Riemann_vec;
    struct Riemann_outputs Riemann_out;
    double face_area_dot_vel;
    face_area_dot_vel = 0;
#endif
    double face_vel_i=0, face_vel_j=0, Face_Area_Norm=0, Face_Area_Vec[3];

#ifdef HYDRO_MESHLESS_FINITE_MASS
    double epsilon_entropic_eos_big = 0.5; // can be anything from (small number=more diffusive, less accurate entropy conservation) to ~1.1-1.3 (least diffusive, most noisy)
    double epsilon_entropic_eos_small = 1.e-3; // should be << epsilon_entropic_eos_big
#if defined(FORCE_ENTROPIC_EOS_BELOW)
    epsilon_entropic_eos_small = FORCE_ENTROPIC_EOS_BELOW; // if set manually
#elif !defined(SELFGRAVITY_OFF)
    epsilon_entropic_eos_small = 1.e-2; epsilon_entropic_eos_big = 0.6; // with gravity larger tolerance behaves better on hydrostatic equilibrium problems //
#endif
#endif

    if(mode == 0)
    {
        particle2in_hydra(&local, target); // this setup allows for all the fields we need to define (don't hard-code here)
    }
    else
    {
        local = HydroDataGet[target]; // this setup allows for all the fields we need to define (don't hard-code here)
    }
    
    /* certain particles should never enter the loop: check for these */
    if(local.Mass <= 0) return 0;
#ifdef GALSF_SUBGRID_WINDS
    if(local.DelayTime > 0) {return 0;}
#endif
    
    /* --------------------------------------------------------------------------------- */
    /* pre-define Particle-i based variables (so we save time in the loop below) */
    /* --------------------------------------------------------------------------------- */
    kernel.sound_i = local.SoundSpeed;
    kernel.spec_egy_u_i = local.InternalEnergyPred;
    kernel.h_i = local.Hsml;
    kernel_hinv(kernel.h_i, &hinv_i, &hinv3_i, &hinv4_i);
    hinv_j=hinv3_j=hinv4_j=0;
    V_i = local.Mass / local.Density;
    Particle_Size_i = pow(V_i,1./NUMDIMS) * All.cf_atime; // in physical, used below in some routines //
    double Amax_i = MAX_REAL_NUMBER;
#if (NUMDIMS==2)
    Amax_i = 2. * sqrt(V_i/M_PI);
#endif
#if (NUMDIMS==3)
    Amax_i = M_PI * pow((3.*V_i)/(4.*M_PI), 2./3.);
#endif    
    dt_hydrostep = local.Timestep * All.Timebase_interval / All.cf_hubble_a; /* (physical) timestep */
    out.MaxSignalVel = kernel.sound_i;
    kernel_mode = 0; /* need dwk and wk */
    double cnumcrit2 = ((double)CONDITION_NUMBER_DANGER)*((double)CONDITION_NUMBER_DANGER) - local.ConditionNumber*local.ConditionNumber;
#if defined(HYDRO_SPH)
#ifdef HYDRO_PRESSURE_SPH
    kernel.p_over_rho2_i = local.Pressure / (local.EgyWtRho*local.EgyWtRho);
#else 
    kernel.p_over_rho2_i = local.Pressure / (local.Density*local.Density);
#endif
#endif
    
#ifdef MAGNETIC
    kernel.b2_i = local.BPred[0]*local.BPred[0] + local.BPred[1]*local.BPred[1] + local.BPred[2]*local.BPred[2];
#if defined(HYDRO_SPH)
    double magfluxv[3],resistivity_heatflux=0; magfluxv[0]=magfluxv[1]=magfluxv[2]=0;
    kernel.mf_i = local.Mass * fac_magnetic_pressure / (local.Density * local.Density);
    kernel.mf_j = local.Mass * fac_magnetic_pressure;
    // PFH: comoving factors here to convert from B*B/rho to P/rho for accelerations //
    double mm_i[3][3], mm_j[3][3];
    for(k = 0; k < 3; k++)
    {
        for(j = 0; j < 3; j++)
            mm_i[k][j] = local.BPred[k] * local.BPred[j];
    }
    for(k = 0; k < 3; k++)
        mm_i[k][k] -= 0.5 * kernel.b2_i;
#endif
    kernel.alfven2_i = kernel.b2_i * fac_magnetic_pressure / local.Density;
    kernel.alfven2_i = DMIN(kernel.alfven2_i, 1000. * kernel.sound_i*kernel.sound_i);
    double vcsa2_i = kernel.sound_i*kernel.sound_i + kernel.alfven2_i;
#endif // MAGNETIC //

#if defined(RT_EVOLVE_NGAMMA_IN_HYDRO)
    double Fluxes_E_gamma[N_RT_FREQ_BINS];
    double tau_c_i[N_RT_FREQ_BINS];
    for(k=0;k<N_RT_FREQ_BINS;k++) {tau_c_i[k] = Particle_Size_i * local.Kappa_RT[k]*local.Density*All.cf_a3inv;}
#ifdef RT_EVOLVE_FLUX
    double Fluxes_Flux[N_RT_FREQ_BINS][3];
#endif
#endif

    
    /* --------------------------------------------------------------------------------- */
    /* Now start the actual SPH computation for this particle */
    /* --------------------------------------------------------------------------------- */
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
#ifndef DONOTUSENODELIST
        startnode = HydroDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
#else
        startnode = All.MaxPart;	/* root node */
#endif
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            /* --------------------------------------------------------------------------------- */
            /* get the neighbor list */
            /* --------------------------------------------------------------------------------- */
            numngb = ngb_treefind_pairs_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
                                       exportnodecount, exportindex, ngblist);
            if(numngb < 0) return -1;
            
            for(n = 0; n < numngb; n++)
            {
                j = ngblist[n];
                
                /* check if I need to compute this pair-wise interaction from "i" to "j", or skip it and 
                    let it be computed from "j" to "i" */
                integertime TimeStep_J = (P[j].TimeBin ? (((integertime) 1) << P[j].TimeBin) : 0);
                int j_is_active_for_fluxes = 0;
#ifndef BOX_SHEARING // (shearing box means the fluxes at the boundaries are not actually symmetric, so can't do this) //
                if(local.Timestep > TimeStep_J) continue; /* compute from particle with smaller timestep */
                /* use relative positions to break degeneracy */
                if(local.Timestep == TimeStep_J)
                {
                    int n0=0; if(local.Pos[n0] == P[j].Pos[n0]) {n0++; if(local.Pos[n0] == P[j].Pos[n0]) n0++;}
                    if(local.Pos[n0] < P[j].Pos[n0]) continue;
                }
                if(TimeBinActive[P[j].TimeBin]) {j_is_active_for_fluxes = 1;}
#endif
                if(P[j].Mass <= 0) continue;
                if(SphP[j].Density <= 0) continue;
#ifdef GALSF_SUBGRID_WINDS
                if(SphP[j].DelayTime > 0) continue; /* no hydro forces for decoupled wind particles */
#endif
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC  /* find the closest image in the given box size  */
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                kernel.h_j = PPP[j].Hsml;
                
                /* force applied for all particles inside each-others kernels! */
                if((r2 >= kernel.h_i * kernel.h_i) && (r2 >= kernel.h_j * kernel.h_j)) continue;
                if(r2 <= 0) continue;
                
                /* --------------------------------------------------------------------------------- */
                /* ok, now we definitely have two interacting particles */
                /* --------------------------------------------------------------------------------- */
                
                /* --------------------------------------------------------------------------------- */
                /* calculate a couple basic properties needed: separation, velocity difference (needed for timestepping) */
                kernel.r = sqrt(r2);
#ifdef HYDRO_REGULAR_GRID
                if(kernel.r > 1.1 * Particle_Size_i * sqrt(NUMDIMS)) continue; // only do interactions for the immediate neighbors //
#endif
                rinv = 1 / kernel.r;
                /* we require a 'softener' to prevent numerical madness in interpolating functions */
                rinv_soft = 1.0 / sqrt(r2 + 0.0001*kernel.h_i*kernel.h_i);
#ifdef BOX_SHEARING
                /* in a shearing box, need to set dv appropriately for the shearing boundary conditions */
                MyDouble VelPred_j[3]; for(k=0;k<3;k++) {VelPred_j[k]=SphP[j].VelPred[k];}
                if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                MyDouble ParticleVel_j[3]; for(k=0;k<3;k++) {ParticleVel_j[k]=SphP[j].VelPred[k];}
                if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {ParticleVel_j[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {ParticleVel_j[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
#else
                /* faster to just set a pointer directly */
                MyDouble *VelPred_j = SphP[j].VelPred;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                MyDouble *ParticleVel_j = SphP[j].ParticleVel;
#endif
#endif
                kernel.dv[0] = local.Vel[0] - VelPred_j[0];
                kernel.dv[1] = local.Vel[1] - VelPred_j[1];
                kernel.dv[2] = local.Vel[2] - VelPred_j[2];
                kernel.rho_ij_inv = 2.0 / (local.Density + SphP[j].Density);
                
                /* --------------------------------------------------------------------------------- */
                /* sound speed, relative velocity, and signal velocity computation */
                kernel.sound_j = Particle_effective_soundspeed_i(j);
                kernel.vsig = kernel.sound_i + kernel.sound_j;
#ifdef MAGNETIC
                double BPred_j[3];
                for(k=0;k<3;k++) {BPred_j[k]=Get_Particle_BField(j,k);} /* defined j b-field in appropriate units for everything */
#ifdef DIVBCLEANING_DEDNER
                double PhiPred_j = Get_Particle_PhiField(j); /* define j phi-field in appropriate units */
#endif
                kernel.b2_j = BPred_j[0]*BPred_j[0] + BPred_j[1]*BPred_j[1] + BPred_j[2]*BPred_j[2];
                kernel.alfven2_j = kernel.b2_j * fac_magnetic_pressure / SphP[j].Density;
                kernel.alfven2_j = DMIN(kernel.alfven2_j, 1000. * kernel.sound_j*kernel.sound_j);
                double vcsa2_j = kernel.sound_j*kernel.sound_j + kernel.alfven2_j;
                double Bpro2_j = (BPred_j[0]*kernel.dp[0] + BPred_j[1]*kernel.dp[1] + BPred_j[2]*kernel.dp[2]) / kernel.r;
                Bpro2_j *= Bpro2_j;
                double magneticspeed_j = sqrt(0.5 * (vcsa2_j + sqrt(DMAX((vcsa2_j*vcsa2_j -
                        4 * kernel.sound_j*kernel.sound_j * Bpro2_j*fac_magnetic_pressure/SphP[j].Density), 0))));
                double Bpro2_i = (local.BPred[0]*kernel.dp[0] + local.BPred[1]*kernel.dp[1] + local.BPred[2]*kernel.dp[2]) / kernel.r;
                Bpro2_i *= Bpro2_i;
                double magneticspeed_i = sqrt(0.5 * (vcsa2_i + sqrt(DMAX((vcsa2_i*vcsa2_i -
                        4 * kernel.sound_i*kernel.sound_i * Bpro2_i*fac_magnetic_pressure/local.Density), 0))));
                kernel.vsig = magneticspeed_i + magneticspeed_j;
                Bpro2_i /= kernel.b2_i; Bpro2_j /= kernel.b2_j;
#endif
                kernel.vdotr2 = kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2];
                // hubble-flow correction: need in -code- units, hence extra a2 appearing here //
                if(All.ComovingIntegrationOn) kernel.vdotr2 += All.cf_hubble_a2 * r2;
                if(kernel.vdotr2 < 0)
                {
#if defined(HYDRO_SPH) || defined(HYDRO_MESHLESS_FINITE_VOLUME)
                    kernel.vsig -= 3 * fac_mu * kernel.vdotr2 * rinv;
#else
                    kernel.vsig -= fac_mu * kernel.vdotr2 * rinv;
#endif
                }
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
                double KE = kernel.dv[0]*kernel.dv[0] + kernel.dv[1]*kernel.dv[1] + kernel.dv[2]*kernel.dv[2];
                if(KE > out.MaxKineticEnergyNgb) out.MaxKineticEnergyNgb = KE;
                if(j_is_active_for_fluxes) {if(KE > SphP[j].MaxKineticEnergyNgb) SphP[j].MaxKineticEnergyNgb = KE;}
#endif
#ifdef TURB_DIFF_METALS
                double mdot_estimated = 0;
#endif
#if defined(RT_INFRARED)
                double Fluxes_E_gamma_T_weighted_IR = 0;
#endif
                
                /* --------------------------------------------------------------------------------- */
                /* calculate the kernel functions (centered on both 'i' and 'j') */
                if(kernel.r < kernel.h_i)
                {
                    u = kernel.r * hinv_i;
                    kernel_main(u, hinv3_i, hinv4_i, &kernel.wk_i, &kernel.dwk_i, kernel_mode);
                }
                else
                {
                    kernel.dwk_i = 0;
                    kernel.wk_i = 0;
                }
                if(kernel.r < kernel.h_j)
                {
                    kernel_hinv(kernel.h_j, &hinv_j, &hinv3_j, &hinv4_j);
                    u = kernel.r * hinv_j;
                    kernel_main(u, hinv3_j, hinv4_j, &kernel.wk_j, &kernel.dwk_j, kernel_mode);
                }
                else
                {
                    kernel.dwk_j = 0;
                    kernel.wk_j = 0;
                }
                
                /* --------------------------------------------------------------------------------- */
                /* with the overhead numbers above calculated, we now 'feed into' the "core" 
                    hydro computation (SPH, meshless godunov, etc -- doesn't matter, should all take the same inputs) 
                    the core code is -inserted- here from the appropriate .h file, depending on the mode 
                    the code has been compiled in */
                /* --------------------------------------------------------------------------------- */
#ifdef HYDRO_SPH
#include "hydra_core_sph.h"
#else
#include "hydra_core_meshless.h"
#endif
                
#ifdef FREEZE_HYDRO
                memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#endif


                
                
//#ifndef HYDRO_SPH
                /* the following macros are useful for all the diffusion operations below: this is the diffusion term associated
                    with the HLL reimann problem solution. This adds numerical diffusion (albeit limited to the magnitude of the 
                    physical diffusion coefficients), but stabilizes the relevant equations */
#ifdef HYDRO_SPH
        face_vel_i = face_vel_j = 0;
        for(k=0;k<3;k++) 
        {
        face_vel_i += local.Vel[k] * kernel.dp[k] / (kernel.r * All.cf_atime); 
        face_vel_j += SphP[j].VelPred[k] * kernel.dp[k] / (kernel.r * All.cf_atime);
        }
        // SPH: use the sph 'effective areas' oriented along the lines between particles and direct-difference gradients
        Face_Area_Norm = local.Mass * P[j].Mass * fabs(kernel.dwk_i+kernel.dwk_j) / (local.Density * SphP[j].Density);
        for(k=0;k<3;k++) {Face_Area_Vec[k] = Face_Area_Norm * kernel.dp[k]/kernel.r;}
#endif

#ifdef MAGNETIC
                double bhat[3]={0.5*(local.BPred[0]+BPred_j[0])*All.cf_a2inv,0.5*(local.BPred[1]+BPred_j[1])*All.cf_a2inv,0.5*(local.BPred[2]+BPred_j[2])*All.cf_a2inv};
                double bhat_mag=bhat[0]*bhat[0]+bhat[1]*bhat[1]+bhat[2]*bhat[2];
                if(bhat_mag>0) {bhat_mag=sqrt(bhat_mag); bhat[0]/=bhat_mag; bhat[1]/=bhat_mag; bhat[2]/=bhat_mag;}
                v_hll = 0.5*fabs(face_vel_i-face_vel_j) + DMAX(magneticspeed_i,magneticspeed_j);
#define B_dot_grad_weights(grad_i,grad_j) {if(bhat_mag<=0) {b_hll=1;} else {double q_tmp_sum=0,b_tmp_sum=0; for(k=0;k<3;k++) {\
                                           double q_tmp=0.5*(grad_i[k]+grad_j[k]); q_tmp_sum+=q_tmp*q_tmp; b_tmp_sum+=bhat[k]*q_tmp;}\
                                           if((b_tmp_sum!=0)&&(q_tmp_sum>0)) {b_hll=fabs(b_tmp_sum)/sqrt(q_tmp_sum); b_hll*=b_hll;} else {b_hll=0;}}}
#define HLL_DIFFUSION_COMPROMISE_FACTOR 1.1
#else
                v_hll = 0.5*fabs(face_vel_i-face_vel_j) + DMAX(kernel.sound_i,kernel.sound_j);
#define B_dot_grad_weights(grad_i,grad_j) {b_hll=1;}
#define HLL_DIFFUSION_COMPROMISE_FACTOR 1.5
#endif
#define HLL_correction(ui,uj,wt,kappa) (k_hll = v_hll * (wt) * kernel.r * All.cf_atime / fabs(kappa),\
                                        k_hll = (0.2 + k_hll) / (0.2 + k_hll + k_hll*k_hll),\
                                        -1.0*k_hll*Face_Area_Norm*v_hll*((ui)-(uj)))
#if !defined(MAGNETIC) || defined(GALSF) || defined(COOLING) || defined(BLACKHOLES)
#define HLL_DIFFUSION_OVERSHOOT_FACTOR  0.005
#else
#define HLL_DIFFUSION_OVERSHOOT_FACTOR  1.0
#endif
                
#ifdef EOS_ELASTIC
#include "../solids/elastic_stress_tensor_force.h"
#endif

#ifdef MHD_NON_IDEAL
#include "nonideal_mhd.h"
#endif
                
#ifdef CONDUCTION
#include "conduction.h"
#endif

#ifdef VISCOSITY
#include "viscosity.h"
#endif
                
#ifdef TURB_DIFFUSION
#include "../turb/turbulent_diffusion.h"
#endif

#ifdef CHIMES_TURB_DIFF_IONS
#include "chimes_turbulent_ion_diffusion.h" 
#endif 
                
                
#ifdef RT_DIFFUSION_EXPLICIT
#if defined(RT_EVOLVE_INTENSITIES)
#include "../radiation/rt_direct_ray_transport.h"
#else
#include "../radiation/rt_diffusion_explicit.h"
#endif
#endif
                
                
                /* --------------------------------------------------------------------------------- */
                /* now we will actually assign the hydro variables for the evolution step */
                /* --------------------------------------------------------------------------------- */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                double dmass_holder = Fluxes.rho * dt_hydrostep, dmass_limiter;
                if(dmass_holder > 0) {dmass_limiter=P[j].Mass;} else {dmass_limiter=local.Mass;}
                dmass_limiter *= 0.1;
                if(fabs(dmass_holder) > dmass_limiter) {dmass_holder *= dmass_limiter / fabs(dmass_holder);}
                out.dMass += dmass_holder;
                out.DtMass += Fluxes.rho;
#ifndef BOX_SHEARING
                SphP[j].dMass -= dmass_holder;
#endif
                double gravwork[3]; gravwork[0]=Fluxes.rho*kernel.dp[0]; gravwork[1]=Fluxes.rho*kernel.dp[1]; gravwork[2]=Fluxes.rho*kernel.dp[2];
                for(k=0;k<3;k++) {out.GravWorkTerm[k] += gravwork[k];}
#endif
                for(k=0;k<3;k++) {out.Acc[k] += Fluxes.v[k];}
                out.DtInternalEnergy += Fluxes.p;                
#ifdef MAGNETIC
#ifndef HYDRO_SPH
                for(k=0;k<3;k++) {out.Face_Area[k] += Face_Area_Vec[k];}
#endif
#ifndef FREEZE_HYDRO
                for(k=0;k<3;k++) {out.DtB[k]+=Fluxes.B[k];}
                out.divB += Fluxes.B_normal_corrected;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME) // mass-based phi-flux
                out.DtPhi += Fluxes.phi;
#endif
#ifdef HYDRO_SPH
                for(k=0;k<3;k++) {out.DtInternalEnergy+=magfluxv[k]*local.Vel[k]/All.cf_atime;}
                out.DtInternalEnergy += resistivity_heatflux;
#else
                double wt_face_sum = Face_Area_Norm * (-face_area_dot_vel+face_vel_i);
                out.DtInternalEnergy += 0.5 * kernel.b2_i*All.cf_a2inv*All.cf_a2inv * wt_face_sum;
#ifdef DIVBCLEANING_DEDNER
                for(k=0; k<3; k++)
                {
                    out.DtB_PhiCorr[k] += Riemann_out.phi_normal_db * Face_Area_Vec[k];
                    out.DtB[k] += Riemann_out.phi_normal_mean * Face_Area_Vec[k];
                    out.DtInternalEnergy += Riemann_out.phi_normal_mean * Face_Area_Vec[k] * local.BPred[k]*All.cf_a2inv;
                }
#endif
#ifdef MHD_NON_IDEAL
                for(k=0;k<3;k++) {out.DtInternalEnergy += local.BPred[k]*All.cf_a2inv*bflux_from_nonideal_effects[k];}
#endif
#endif
#endif
#endif // magnetic //
                
                /* if this is particle j's active timestep, you should sent them the time-derivative information as well, for their subsequent drift operations */
                if(j_is_active_for_fluxes)
                {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                    SphP[j].DtMass -= Fluxes.rho;
                    for(k=0;k<3;k++) {SphP[j].GravWorkTerm[k] -= gravwork[k];}
#endif
                    for(k=0;k<3;k++) {SphP[j].HydroAccel[k] -= Fluxes.v[k];}
                    SphP[j].DtInternalEnergy -= Fluxes.p;
#ifdef MAGNETIC
#ifndef HYDRO_SPH
                    for(k=0;k<3;k++) {SphP[j].Face_Area[k] -= Face_Area_Vec[k];}
#endif
#ifndef FREEZE_HYDRO
                    for(k=0;k<3;k++) {SphP[j].DtB[k]-=Fluxes.B[k];}
                    SphP[j].divB -= Fluxes.B_normal_corrected;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME) // mass-based phi-flux
                    SphP[j].DtPhi -= Fluxes.phi;
#endif
#ifdef HYDRO_SPH
                    for(k=0;k<3;k++) {SphP[j].DtInternalEnergy-=magfluxv[k]*VelPred_j[k]/All.cf_atime;}
                    SphP[j].DtInternalEnergy += resistivity_heatflux;
#else
                    double wt_face_sum = Face_Area_Norm * (-face_area_dot_vel+face_vel_j);
                    SphP[j].DtInternalEnergy -= 0.5 * kernel.b2_j*All.cf_a2inv*All.cf_a2inv * wt_face_sum;
#ifdef DIVBCLEANING_DEDNER
                    for(k=0; k<3; k++)
                    {
                        SphP[j].DtB_PhiCorr[k] -= Riemann_out.phi_normal_db * Face_Area_Vec[k];
                        SphP[j].DtB[k] -= Riemann_out.phi_normal_mean * Face_Area_Vec[k];
                        SphP[j].DtInternalEnergy -= Riemann_out.phi_normal_mean * Face_Area_Vec[k] * BPred_j[k]*All.cf_a2inv;
                    }
#endif
#ifdef MHD_NON_IDEAL
                    for(k=0;k<3;k++) {SphP[j].DtInternalEnergy -= BPred_j[k]*All.cf_a2inv*bflux_from_nonideal_effects[k];}
#endif
#endif
#endif
#endif // magnetic //

                }

                /* if we have mass fluxes, we need to have metal fluxes if we're using them (or any other passive scalars) */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                if(dmass_holder != 0)
                {
#ifdef METALS
                    if(Fluxes.rho > 0)
                    {
                        /* particle i gains mass from particle j */
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            out.Dyield[k] += (P[j].Metallicity[k] - local.Metallicity[k]) * dmass_holder;
                    } else {
                        /* particle j gains mass from particle i */
                        dmass_holder /= -P[j].Mass;
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            P[j].Metallicity[k] += (local.Metallicity[k] - P[j].Metallicity[k]) * dmass_holder;
                    }
#endif
                }
#endif

                /* --------------------------------------------------------------------------------- */
                /* don't forget to save the signal velocity for time-stepping! */
                /* --------------------------------------------------------------------------------- */
                if(kernel.vsig > out.MaxSignalVel) out.MaxSignalVel = kernel.vsig;
                if(j_is_active_for_fluxes) {if(kernel.vsig > SphP[j].MaxSignalVel) SphP[j].MaxSignalVel = kernel.vsig;}
#ifdef WAKEUP
                if(!(TimeBinActive[P[j].TimeBin]))
                {
                    if(kernel.vsig > WAKEUP*SphP[j].MaxSignalVel) PPPZ[j].wakeup = 1;
#if (SLOPE_LIMITER_TOLERANCE <= 0)
                    if(local.Timestep*WAKEUP < TimeStep_J) PPPZ[j].wakeup = 1;
#endif
                }
#endif
                
                
            } // for(n = 0; n < numngb; n++) //
        } // while(startnode >= 0) //
#ifndef DONOTUSENODELIST
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = HydroDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        } // if(mode == 1) //
#endif
    } // while(startnode >= 0) //
    
    /* Now collect the result at the right place */
    if(mode == 0)
        out2particle_hydra(&out, target, 0);
    else
        HydroDataResult[target] = out;
    
    return 0;
}

