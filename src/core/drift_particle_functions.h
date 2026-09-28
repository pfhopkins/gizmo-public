/* drift_particle_functions.h — the canonical body of the per-particle drift.
 *
 * Single source of truth for both the host drift and, once the drift is offloaded,
 * the device kernel. The body addresses particles through the pointers it is given
 * rather than the P[]/CellP[] globals, so it is correct when it is handed a compacted
 * staging buffer: under staging the index space is the buffer's, and any global reach
 * would silently read a different particle.
 *
 * The drift and gravkick time factors come from a caller-supplied table view for the
 * same reason -- the tables live in host memory and the device sees a mirror of them.
 *
 * Unlike core/predict_functions.h, which is a leaf header included by other physics
 * headers, this one sits at the top of the include order and may pull in the homes of
 * everything the drift calls. That is what lets the body be device-callable here
 * instead of only wherever it happens to be instantiated.
 *
 * Include order: after allvars.h and core/proto.h. */
#pragma once

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

#include "timestep_functions.h"      /* DriftKickTableView, get_{drift,gravkick}_factor_impl, dilation */
#include "../gravity/binary_functions.h" /* odeint_super_timestep */
#include "predict_functions.h"       /* advect_mesh_point_P, apply_special_boundary_conditions_P */
#include "../gravity/ags_functions.h"    /* ags_density_isactive_P, ags_return_{min,max}soft_P, dm_fuzzy */
#include "../eos/eos_functions.h"        /* set_eos_pressure_impl */
#ifdef COSMIC_RAY_FLUID
#include "../eos/cosmic_ray_fluid/cosmic_ray_functions.h"
#endif
#ifdef RADTRANSFER
#include "../radiation/rt_functions.h"
#endif
#ifdef EOS_ELASTIC
#include "../solids/elastic_physics_functions.h"
#endif


KOKKOS_INLINE_FUNCTION
void drift_extra_physics_P(int i, integertime tstart, integertime tend, double dt_entr,
                           struct particle_data *pp, struct gas_cell_data *cell)
{
#ifdef MAGNETIC
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    cell[i].BPred += cell[i].DtB * (dt_entr * BphysVolphys_to_BcodeVolCode); // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1 / All.cf_a3inv; // for mass-based phi fluxes (otherwise coefficient is 1)
    double dtphi_code = (PhiphysVolphys_to_PhicodeVolCode) * cell[i].DtPhi;
    cell[i].PhiPred += dtphi_code  * dt_entr;
    double t_damp = Get_Gas_PhiField_DampingTimeInv_P(i, pp, cell);
    if((t_damp>0) && (!isnan(t_damp)))
    {
        cell[i].PhiPred *= exp( -dt_entr * t_damp );
    }
#endif
#ifdef MHD_ALTERNATIVE_LEAPFROG_SCHEME
    cell[i].B = cell[i].BPred;
#ifdef DIVBCLEANING_DEDNER
    cell[i].Phi=cell[i].PhiPred;
#endif
#endif
#endif
#ifdef COSMIC_RAY_FLUID
    CosmicRay_Update_DriftKick(i, dt_entr, 1, pp, cell);
#endif
#ifdef RADTRANSFER
    rt_update_driftkick(i, dt_entr, 1, pp, cell);
#endif
#ifdef EOS_ELASTIC
    elastic_body_update_driftkick_P(i,dt_entr,1,pp,cell);
#endif
}


KOKKOS_INLINE_FUNCTION
void drift_particle_impl(int i, integertime time1, struct particle_data *pp,
                         struct gas_cell_data *cell, const struct DriftKickTableView *tables,
                         const struct EosTableView *eos_tables)
{
    int j __attribute__((unused)); double dt_drift; integertime time0 = pp[i].Ti_current;
    if(time1 < time0)
    {
        PRINT_WARNING("no prediction into past allowed: i=%d time0=%lld time1=%lld", i, (long long)time0, (long long)time1);
        endrun(90001004);
        return;   /* graceful: skip the drift; bad-stop drains at the next phase poll */
    }
    if(time1 == time0) {return;}
    
    dt_drift = get_drift_factor_impl(time0, time1, timestep_dilation_factor(i, pp), tables);
        
#if (SINGLE_STAR_TIMESTEPPING > 0)
    /* A super-timestepped sink: the drift below moves its binary's centre of mass, and the
       internal orbit is advanced on top of that.  Also consumed by the velocity prediction
       further down, so declared outside the position update. */
    Vec3<double> fewbody_drift_dx = {}, fewbody_kick_dv = {};
    volatile int super_timestepped_sink = ((pp[i].Type == 5) && (pp[i].SuperTimestepFlag >= 2)) ? 1 : 0;
#endif
#if !defined(FREEZE_HYDRO)
    /* A finite-volume gas cell moves with its mesh-generating point, a super-timestepped sink
       with its binary, everything else with its own velocity.  The two special cases are
       different particle types, so both must be live in a build that has both. */
#if (SINGLE_STAR_TIMESTEPPING > 0)
    if(super_timestepped_sink)
    {
        Vec3<double> COM_Vel = pp[i].Vel + pp[i].comp_dv * (pp[i].comp_Mass/(pp[i].Mass+pp[i].comp_Mass)); //center of mass velocity
        pp[i].Pos += COM_Vel * dt_drift; //center of mass drift
        odeint_super_timestep(i, dt_drift, fewbody_kick_dv, fewbody_drift_dx, 1, pp); // do_fewbody_drift
        pp[i].GravAccel = pp[i].COM_GravAccel; //Overwrite the acceleration with center of mass value
        pp[i].Pos += fewbody_drift_dx; //Keplerian evolution
        pp[i].Vel += fewbody_kick_dv; //move on binary.orbit
    }
    else
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    if(pp[i].Type==0) {advect_mesh_point_P(i,dt_drift,pp,cell);}
    else
#endif
    {pp[i].Pos += pp[i].Vel * dt_drift;}
#endif // FREEZE_HYDRO clause
#if (NUMDIMS==1)
    pp[i].Pos[1]=pp[i].Pos[2]=0; // force zero-ing
#endif
#if (NUMDIMS==2)
    pp[i].Pos[2]=0; // force zero-ing
#endif

#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    double dilation = timestep_dilation_factor(i, pp); /* f = 1/a <= 1 */
    if(dilation < 1.) {
        /* the drift above advanced the particle over only the fraction f of the raw interval, since
           dt_drift already carries the f. add back the bulk motion over the remaining (1-f) of the
           raw interval, so that only the motion relative to the surroundings is dilated */
        double cfac = dt_drift * (1./dilation - 1.);
        pp[i].Pos += pp[i].vel_of_nearest_special * cfac;
    }
#endif

    double divv_fac = pp[i].Particle_DivVel * dt_drift;
    /* How far the predicted kernel radius may move in one drift.  The two directions
       are not the same risk: growing follows the flow and simply predicts a larger
       search, while shrinking can empty a kernel of the neighbours the next density
       solve needs.  So growth carries the looser bound, which also tracks the
       divergence it is approximating more closely than the old tight one did.  Adaptive
       softening has used the loose bound on both sides for the same reason. */
    double divv_fac_max_grow   = 4.0; /* was 0.3, more conservative, in older versions of the code, to prevent large jumps, but this leads to lower order and failures of the density guess more often */
    double divv_fac_max_shrink = 1.5; /* was 0.3, more conservative, in older versions of the code, to prevent large jumps, but this leads to lower order and failures of the density guess more often */
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(ags_density_isactive_P(i, pp) && pp[i].Type>0) {divv_fac_max_grow = 4; divv_fac_max_shrink = 4;}
#endif
    if(divv_fac > +divv_fac_max_grow)   divv_fac = +divv_fac_max_grow;
    if(divv_fac < -divv_fac_max_shrink) divv_fac = -divv_fac_max_shrink;
    
#ifdef GRAIN_FLUID
    if((1 << pp[i].Type) & (GRAIN_PTYPES))
    {
        pp[i].KernelRadius *= exp((double)divv_fac / ((double)NUMDIMS));
        if(pp[i].KernelRadius < All.MinKernelRadius) {pp[i].KernelRadius = All.MinKernelRadius;}
        if(pp[i].KernelRadius > All.MaxKernelRadius) {pp[i].KernelRadius = All.MaxKernelRadius;}
    }
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(ags_density_isactive_P(i, pp) && (dt_drift>0)) /* particle is AGS-active */
    {
        double minsoft = ags_return_minsoft_P(i, pp), maxsoft = ags_return_maxsoft_P(i, pp);
        pp[i].AGS_KernelRadius *= exp((double)divv_fac / ((double)NUMDIMS));
        if(pp[i].AGS_KernelRadius < minsoft) {pp[i].AGS_KernelRadius = minsoft;}
        if(pp[i].AGS_KernelRadius > maxsoft) {pp[i].AGS_KernelRadius = maxsoft;}
    } else {pp[i].AGS_KernelRadius = ForceSoftening_KernelRadius_P(i, pp);} /* non-AGS-active particles use fixed softening */
#endif
    
#ifdef DM_FUZZY
    do_dm_fuzzy_drift_kick_P(i, dt_drift, 1, pp);
#endif

#ifdef CBE_INTEGRATOR
    /* CBE-moment predictor SUPPRESSED: the implemented per-basis CBE-moment
     * drift-prediction degraded accuracy (it damped the harmonic breathing test
     * more than leaving it off), so the do_cbe_predict_drift() call is removed
     * pending a corrected revival (separate future work). The general
     * AGS_KernelRadius / gas VelPred prediction is unaffected; do_cbe_predict_drift()
     * remains defined (currently unused) as a placeholder for that revival. */
#endif

    if((pp[i].Type == 0) && (pp[i].Mass > 0))
        {
            double dt_gravkick, dt_gravkick_pm, dt_hydrokick, dt_entr;
            dt_entr = dt_hydrokick = (time1 - time0) * unit_integertime_in_physical(i, pp);
            dt_gravkick = get_gravkick_factor_impl(time0, time1, timestep_dilation_factor(i, pp), tables);
            
#ifdef PMGRID
            dt_gravkick_pm = get_gravkick_factor_impl(time0, time1, timestep_dilation_factor(-1, pp), tables);
            cell[i].VelPred += pp[i].GravAccel*dt_gravkick + pp[i].GravPM*dt_gravkick_pm + cell[i].HydroAccel*(dt_hydrokick*All.cf_atime); /* make sure v is in code units */
#else
            cell[i].VelPred += pp[i].GravAccel * dt_gravkick + cell[i].HydroAccel * (dt_hydrokick*All.cf_atime); /* make sure v is in code units */
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
	        if(super_timestepped_sink) {cell[i].VelPred += fewbody_kick_dv;}
#endif	    
            
#if defined(TURB_DRIVING)
            cell[i].VelPred += cell[i].TurbAccel * dt_gravkick;
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            cell[i].VelPred += cell[i].Rad_Accel * (All.cf_atime * dt_hydrokick);
#endif
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            pp[i].Mass = DMAX(pp[i].Mass + cell[i].DtMass * dt_entr, 0.5 * cell[i].MassTrue); cell[i].Mass = pp[i].Mass;
#endif
            
            cell[i].Density *= exp(-divv_fac);
            double etmp = cell[i].InternalEnergyPred + cell[i].DtInternalEnergy * dt_entr;
#if defined(RADTRANSFER) && defined(RT_EVOLVE_ENERGY) /* block here to deal with tricky cases where radiation energy density is -much- larger than thermal */ 
            int kfreq; double erad_tot=0,tot_e_min=0,enew=0,int_e_min=0,dErad=0,rsol_fac=C_LIGHT_CODE_REDUCED/C_LIGHT_CODE; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad_tot+=cell[i].Rad_E_gamma_Pred[kfreq];}
            if(erad_tot > 0)
            {
                int_e_min=0.025*cell[i].InternalEnergyPred; tot_e_min=0.025*(erad_tot/rsol_fac+cell[i].InternalEnergyPred*pp[i].Mass);
                enew=DMAX(erad_tot/rsol_fac+etmp*pp[i].Mass,tot_e_min); etmp=(enew-erad_tot/rsol_fac)/pp[i].Mass; if(etmp<int_e_min) {dErad=rsol_fac*(etmp-int_e_min); etmp=int_e_min;}
                if(dErad<-0.975*erad_tot) {dErad=-0.975*erad_tot;} cell[i].InternalEnergyPred = etmp; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {cell[i].Rad_E_gamma_Pred[kfreq] *= 1 + dErad/erad_tot;}
            } else {
                if(etmp<0.5*cell[i].InternalEnergyPred) {cell[i].InternalEnergyPred *= 0.5;} else {cell[i].InternalEnergyPred=etmp;}
            }
#else
            if(etmp<0.5*cell[i].InternalEnergyPred) {cell[i].InternalEnergyPred *= 0.5;} else {cell[i].InternalEnergyPred=etmp;}
#endif
            if(cell[i].InternalEnergyPred<All.MinEgySpec) cell[i].InternalEnergyPred=All.MinEgySpec;
            
#ifdef HYDRO_PRESSURE_SPH
            cell[i].EgyWtDensity *= exp(-divv_fac);
#endif
            
#if (HYDRO_FIX_MESH_MOTION > 0)
            pp[i].KernelRadius *= exp((double)divv_fac / ((double)NUMDIMS));
            if(pp[i].KernelRadius < All.MinKernelRadius) {pp[i].KernelRadius = All.MinKernelRadius;}
            if(pp[i].KernelRadius > All.MaxKernelRadius) {pp[i].KernelRadius = All.MaxKernelRadius;}
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            if(1 & ADAPTIVE_GRAVSOFT_FORALL) {pp[i].AGS_KernelRadius = pp[i].KernelRadius;} /* gas is AGS-active, so needs to be set here to match updated KernelRadius */
#endif
#endif
            drift_extra_physics_P(i, time0, time1, dt_entr, pp, cell);

            set_eos_pressure_impl(i, pp, cell, eos_tables);
        }
    
    /* check for reflecting or outflow or otherwise special boundaries: if so, do the reflection/boundary! */
    apply_special_boundary_conditions_P(i,pp[i].Mass,0,pp,cell);

    pp[i].Ti_current = time1;
}
