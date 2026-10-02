#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "timestep_functions.h" /* apply_special_boundary_conditions_P */
#include "../system/gpu_particles_arena.h"
#include "../mesh/gpu_neighbor_list.h" /* gpu_sidx_notify_owned_changed */
#ifdef CBE_INTEGRATOR
#include "../sidm/sidm_gpu_decls.h"
#endif

/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO 
 * (added energy/entropy switch, terms for explicit mass conservation in mass fluxes, 
 *  and updates to additional fluid variables, options for different hydro solvers,
 *  libraries and architectures for different grid/cell structures, rewritten
 *  to allow boundary conditions other than periodic and open, etc.)
 */

void apply_long_range_kick(integertime tstart, integertime tend);

void do_first_halfstep_kick(void)
{
    int i; integertime ti_step, tstart=0, tend=0;
    
#ifdef TURB_DRIVING
    do_turb_driving_step_first_half();
#endif
    
#ifdef PMGRID
    if(All.PM_Ti_begstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
        tstart = All.PM_Ti_begstep;
        tend = tstart + ti_step / 2;
        apply_long_range_kick(tstart, tend);
    }
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(ti_step, tstart, tend)
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    /* as currently written with some revisions to MFV methods, should only update on active timesteps */
    for(i = 0; i < NumPart; i++)
    {
        if((TimeBinActive[P[i].TimeBin]) || (P[i].Type==0)) /* active OR gas, need to check each timestep to ensure manifest conservation */
#else
    for (int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++) /* 'full' kick for active particles */
    { int i = ActiveParticleList[_apl];
#endif
        {
            if(P[i].Mass > 0)
            {
                ti_step = P[i].integertime_step();
                tstart = P[i].Ti_begstep;	/* beginning of step */
                tend = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */
                do_the_kick(i, tstart, tend, P[i].Ti_current, 0);
            }
        }
    } // for(i = 0; i < NumPart; i++) //
#if defined(CBE_INTEGRATOR)
    {
        int n = (int)ActiveParticleList.size();
        int    *cbe_active = (int *)   mymalloc("cbe_kick_idx", (n>0?n:1)*sizeof(int));
        double *cbe_dt     = (double *)mymalloc("cbe_kick_dt",  (n>0?n:1)*sizeof(double));
        int k = 0;
        for(int _a = 0; _a < n; _a++) {
            int ii = ActiveParticleList[_a];
            if(P[ii].Mass > 0 && CBE_INTEGRATOR_DOES_TYPE(P[ii].Type)) {
                cbe_active[k] = ii;
                cbe_dt[k] = (double)(P[ii].integertime_step()/2) * unit_integertime_in_physical(ii, P);
                k++;
            }
        }
        cbe_drift_kick_evaluate_gpu(P, cbe_active, k, cbe_dt);
        myfree(cbe_dt); myfree(cbe_active);
    }
#endif
    gpu_particles_arena_invalidate(); /* host Vel/Entropy kicked; arena stale */
}

void do_second_halfstep_kick(void)
{
    int i; integertime ti_step, tstart=0, tend=0;

#if defined(CBE_INTEGRATOR)
    /* Palindromic operator split: kick1 is gravity->CBE; kick2 is CBE->gravity,
     * so the full step is G C | C G — symmetric in the gravity/CBE coupling
     * (2nd-order split instead of 1st). Run the CBE kick BEFORE BOTH gravity
     * half-kicks here (the PM long-range kick is moved below the CBE block too);
     * the predicted state is resynced to conserved AFTER all gravity below (the
     * pre-gravity CBE-kick reset would otherwise leave it stale). */
    {
        int n = (int)ActiveParticleList.size();
        int    *cbe_active = (int *)   mymalloc("cbe_kick_idx", (n>0?n:1)*sizeof(int));
        double *cbe_dt     = (double *)mymalloc("cbe_kick_dt",  (n>0?n:1)*sizeof(double));
        int k = 0;
        for(int _a = 0; _a < n; _a++) {
            int ii = ActiveParticleList[_a];
            if(P[ii].Mass > 0 && CBE_INTEGRATOR_DOES_TYPE(P[ii].Type)) {
                cbe_active[k] = ii;
                cbe_dt[k] = (double)(P[ii].integertime_step()/2) * unit_integertime_in_physical(ii, P);
                k++;
            }
        }
        cbe_drift_kick_evaluate_gpu(P, cbe_active, k, cbe_dt);
        myfree(cbe_dt); myfree(cbe_active);
    }
#endif
#ifdef PMGRID
    if(All.PM_Ti_endstep == All.Ti_Current)	/* second-half PM long-range kick — AFTER the CBE block so all of gravity (PM long-range + short-range) lands on one side of the CBE kick, keeping kick2 = C->G palindromic with kick1 G->C. */
    {
        ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
        tstart = All.PM_Ti_begstep + ti_step / 2;
        tend = tstart + ti_step / 2;
        apply_long_range_kick(tstart, tend);
    }
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(ti_step, tstart, tend)
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    for(i = 0; i < NumPart; i++)
    {
        if((TimeBinActive[P[i].TimeBin]) || (P[i].Type==0)) /* active OR gas, need to check each timestep to ensure manifest conservation */
#else
    for (int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++) /* 'full' kick for active particles */
    { int i = ActiveParticleList[_apl];
#endif
        {
            if(P[i].Mass > 0)
            {
                ti_step = P[i].integertime_step();
                tstart = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */
                tend = P[i].Ti_begstep + ti_step;	/* end of step */
                do_the_kick(i, tstart, tend, P[i].Ti_current, 1);
                set_predicted_quantities_for_extra_physics(i);
            }
        }
    } // for(i = 0; i < NumPart; i++) //
#if defined(CBE_INTEGRATOR)
    /* Resync predicted CBE state to conserved AFTER the gravity half-kick. The
     * CBE kick (above, before gravity) reset pred=conserved in the pre-gravity
     * frame; gravity has since advanced pi.Vel, so without this the predicted
     * frame velocity would be stale at kick completion. Restores pred==conserved
     * at the end of the step. */
    for(int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++) {
        int i = ActiveParticleList[_apl];
        if(P[i].Mass > 0 && CBE_INTEGRATOR_DOES_TYPE(P[i].Type)) cbe_sync_pred_to_conserved(i);
    }
#endif

#ifdef TURB_DRIVING
    do_turb_driving_step_second_half();
#endif
    gpu_particles_arena_invalidate(); /* host Vel/Entropy kicked; arena stale */
}

#ifdef HERMITE_INTEGRATION
// Initial "prediction" step of Hermite integration, performed after the initial force evaluation 
// Note: the below routines only account for gravitational acceleration - only appropriate for stars or collisionless particles
void do_hermite_prediction(void)
{
    int i; integertime ti_step, tstart=0, tend=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(ti_step, tstart, tend)
#endif
    for (int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++) {
        int i = ActiveParticleList[_apl];
	if(eligible_for_hermite(i, P)) { /* check if we're actually eligible */
	    if(P[i].Mass > 0) { /* skip massless particles scheduled for deletion */
		ti_step = P[i].integertime_step();
		tstart = P[i].Ti_begstep;    /* beginning of step */
		tend = P[i].Ti_begstep + ti_step;    /* end of step */
            double dt_grav = get_gravkick_factor(tstart, tend, i, 0);
#ifdef PMGRID
            //Add the long-range kick from the first half-step, if necessary (since we are overwriting the previous kick operations with the Hermite scheme)
            if(All.PM_Ti_begstep == All.Ti_Current)	/* need to do long-range kick */
            {
                double dt_grav_pm = get_gravkick_factor(All.PM_Ti_begstep, All.PM_Ti_begstep + (All.PM_Ti_endstep - All.PM_Ti_begstep)/2, i, 0);
                P[i].OldVel += P[i].GravPM * dt_grav_pm;
            }
#endif
            P[i].Pos = P[i].OldPos + (P[i].OldVel + (P[i].Hermite_OldAcc + P[i].OldJerk * (dt_grav/3)) * (dt_grav/2)) * dt_grav;
            P[i].Vel = P[i].OldVel + (P[i].Hermite_OldAcc + P[i].OldJerk * (dt_grav/2)) * dt_grav;
		}}} // for (int _apl : ActiveParticleList)
}

void do_hermite_correction(void) // corrector step
{
    int i; integertime ti_step, tstart=0, tend=0;
    int gas_corrected = 0; /* the corrector writes positions, which the kept gas neighbour index must see */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(ti_step, tstart, tend) reduction(+:gas_corrected)
#endif
    for (int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++) {
        int i = ActiveParticleList[_apl];
	if(eligible_for_hermite(i, P)){
                if(P[i].Mass > 0) {
                    if(P[i].Type == 0) {gas_corrected++;}
                    ti_step = P[i].integertime_step();
                    tstart = P[i].Ti_begstep;    /* beginning of step */
                    tend = P[i].Ti_begstep + ti_step;    /* end of step */
                    double dt_grav = get_gravkick_factor(tstart, tend, i, 0);
                    P[i].Vel = P[i].OldVel + (P[i].Hermite_OldAcc + P[i].GravAccel) * (dt_grav * 0.5) + (P[i].OldJerk - P[i].GravJerk) * (dt_grav * dt_grav / 12);
                    P[i].Pos = P[i].OldPos + (P[i].Vel + P[i].OldVel) * (dt_grav * 0.5) + (P[i].Hermite_OldAcc - P[i].GravAccel) * (dt_grav * dt_grav / 12);
#ifdef PMGRID
                    //Add the long-range kick from the second half-step, if necessary (since we are overwriting the previous kick operations with the Hermite scheme)
                    if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
                    {
                        double dt_grav_pm = get_gravkick_factor(All.PM_Ti_begstep + (All.PM_Ti_endstep - All.PM_Ti_begstep)/2, All.PM_Ti_endstep, i, 0);
                        P[i].Vel += P[i].GravPM * dt_grav_pm;
                    }
#endif
		}}} //     for (int _apl : ActiveParticleList)
    if(gas_corrected > 0) {gpu_sidx_notify_owned_changed();}
}
#endif // HERMITE_INTEGRATION


#ifdef PMGRID
void apply_long_range_kick(integertime tstart, integertime tend)
{
    int i;
    double dt_gravkick = get_gravkick_factor(tstart, tend, -1, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Mass > 0)
        {
            Vec3<double> dvel = P[i].GravPM * dt_gravkick; /* do the kick, only collisionless particles */
            P[i].Vel += dvel;
            P[i].dp += dvel * P[i].Mass;
        }
    }
}
#endif


void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent, int mode)
{
    Vec3<double> dp; double dt_entr, dt_gravkick, dt_hydrokick;
    double mass_old, mass_pred, mass_new;
    mass_old = mass_pred = mass_new = P[i].Mass;    
    
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    /* need to do the slightly more complicated update scheme to maintain exact mass conservation */
    if(P[i].Type==0)
    {
        if(CellP[i].dMass != 0) //ent_old = CellP[i].InternalEnergy; for(j=0;j<3;j++) v_old[j] = P[i].Vel[j];
        {
            double dMass=0; // fraction of delta_conserved to couple per kick step (each 'kick' is 1/2-timestep) // double dv[3], v_old[3], dMass, ent_old=0, d_inc = 0.5;
            if(mode != 0) // update the --conserved-- variables of each particle //
            {
                dMass = ((tend - tstart) * unit_integertime_in_physical(i, P)) * CellP[i].DtMass; if(dMass * CellP[i].dMass < 0) {dMass = 0;} // slope-limit: no opposite reconstruction! //
                if((fabs(dMass) > fabs(CellP[i].dMass))) {dMass = CellP[i].dMass;} // try to get close to what the time-integration scheme would give //
                CellP[i].dMass -= dMass;
            } else {dMass = CellP[i].dMass;}
            if(dMass < -0.99*CellP[i].MassTrue) {dMass = -0.99*CellP[i].MassTrue;} // limiter to prevent madness //

            /* load and update the particle masses : particle mass update here, from hydro fluxes */
            mass_old = CellP[i].MassTrue; mass_pred = P[i].Mass; mass_new = mass_old + dMass; CellP[i].MassTrue = mass_new; // UNITS: remember all time derivatives (DtX, dX) are in -physical- units; as are mass, entropy/internal energy, but -not- velocity //
            /* double e_old = mass_old * CellP[i].InternalEnergy; for(j = 0; j< 3; j++) e_old += 0.5*mass_old * (P[i].Vel[j]/All.cf_atime)*(P[i].Vel[j]/All.cf_atime); // physical //
            for(j = 0; j < 3; j++) // momentum-space-kick
            {
                dp[j] = d_inc * CellP[i].dMomentum[j]; // now update the velocity based on the total momentum change
                P[i].Vel[j] = (mass_old*P[i].Vel[j] + dp[j]*All.cf_atime) / mass_new; // call after tabulating dP[j] //
            } // kick for gas internal energy/entropy
            e_old += d_inc * CellP[i].dInternalEnergy; // for(j = 0; j< 3; j++) e_old -= 0.5*mass_new * (P[i].Vel[j]/All.cf_atime)*(P[i].Vel[j]/All.cf_atime); // increment of total (thermal+kinetic) energy; subtract off the new kinetic energy //
            CellP[i].InternalEnergy = e_old / mass_new; CellP[i].enforce_temperature_floor(); // obtain the new internal energy per unit mass, check floor // */
             
            // at the end of this kick, need to re-zero the dInternalEnergy, and other conserved-variable gas/fluid quantities set in the hydro loop, to avoid double-counting them
            if(mode==0) {CellP[i].dMass=0;} /* CellP[i].dInternalEnergy=0; CellP[i].dMomentum[0]=CellP[i].dMomentum[1]=CellP[i].dMomentum[2]=0; */
        }
    } // if(P[i].Type==0) //
#endif
    
    /* only enter the 'normal' kick loop below for genuinely active particles */
    if(TimeBinActive[P[i].TimeBin])
    {
        /* get the timestep (physical units for dt_entr and dt_hydrokick) */
        dt_entr = dt_hydrokick = (tend - tstart) * unit_integertime_in_physical(i, P);
        dt_gravkick = get_gravkick_factor(tstart, tend, i, 0);
        
        if(P[i].Type==0)
        {
            Vec3<double> grav_acc; double dEnt_Gravity = 0;
            grav_acc = P[i].GravAccel * All.cf_a2inv;
#ifdef PMGRID
            grav_acc += P[i].GravPM * All.cf_a2inv;
#endif

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
            dEnt_Gravity += -dot(CellP[i].GravWorkTerm, grav_acc) * All.cf_atime * dt_hydrokick;
#endif
            double du_tot = CellP[i].DtInternalEnergy * dt_hydrokick + dEnt_Gravity;
#if defined(COOLING) && !defined(COOLING_OPERATOR_SPLIT)
            if((mode == 1) && (du_tot != 0) && (dt_hydrokick > 0)) { /* if about to consider second-halfstep kick (just after hydro) */
                CellP[i].CoolingIsOperatorSplitThisTimestep=0; /* default to assume un-split */
#if 0 /* this is a more sophisticated check to see if we need to split the cooling from the hydro work term, but it is not currently used because it is less accurate for complicated non-trivial chemistry, which we now adopt, and the energy conservation issues were more related to older problems in the code that we have since fixed*/
                /* decide if we need to split this particular cell on this particular timestep, since this un-split solver can lead to energy conservation problems if the mechanical heating is much larger than cooling; compare the effective heating rate from the hydro work term to a reference cooling time (which is set to a typical cooling time for the cell, based on its density and the current timestep) */
                CellP[i].CoolingIsOperatorSplitThisTimestep=1; /* default to assume split */
                double DtInternalEnergyEff = du_tot / dt_hydrokick; double DtInternalEnergyReference = 1.e-20*CellP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS / ((UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS_CGS/HYDROGEN_MASSFRAC)); /* define the effective work term and a reference typical cooling time (converted to code units) */
                double DtInternalEnergyDyn = DMAX(2.*CellP[i].InternalEnergy / dt_hydrokick, CellP[i].InternalEnergy * fabs(P[i].Particle_DivVel)*All.cf_a2inv); /* define a dynamical time for the cell, based on either the timestep or the divergence of the velocity field */
                DtInternalEnergyReference = DMAX(DtInternalEnergyReference, DtInternalEnergyDyn); /* use the larger of the two as the reference cooling time */
                if(DtInternalEnergyEff < DtInternalEnergyReference) {CellP[i].CoolingIsOperatorSplitThisTimestep=0;} /* cooling is fast compared to the hydro work term, or the hydro term is negative [cooling], so un-split the operation */
#endif
            }
            if(CellP[i].CoolingIsOperatorSplitThisTimestep==0) {du_tot=0;} /* cooling in unsplit, so zero contribution here */
#endif
            double dEnt = CellP[i].InternalEnergy + du_tot;
            
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            /* if we're using a Riemann solver, we include an energy/entropy-type switch to ensure
                that we don't corrupt the temperature evolution of extremely cold, adiabatic flows */
            /* MHD tests suggest that this switch often does more harm than good: we will
             pay the price of noisier temperature fields (corrupting them when c_s << v_A << v_bulk)
             and they are dynamically irrelevant, in exchange for avoiding potentially much more
             serious errors if this tripped when the B-fields were important */
            double e_thermal,e_kinetic,e_potential;
            e_potential = grav_acc.norm_sq();
            e_potential = P[i].Mass * sqrt(e_potential) * (P[i].Get_Particle_Size()*All.cf_atime); // = M*|a_grav|*h (physical)
            e_kinetic = 0.5 * P[i].Mass * All.cf_a2inv * CellP[i].MaxKineticEnergyNgb;
            e_thermal = DMAX(0.5*CellP[i].InternalEnergy, dEnt) * P[i].Mass;
#ifdef MAGNETIC
            e_thermal += 0.5*CellP[i].B.norm_sq()*CellP[i].Density/(All.cf_atime*P[i].Mass);
#endif
            int do_entropy = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            if(0.01*(e_thermal+e_kinetic) > e_thermal) {do_entropy=1;}
#else
            if(0.005*(e_thermal+e_kinetic) > e_thermal) {do_entropy=1;}
#endif
            do_entropy = 0;
            if(0.01*e_potential > e_thermal) {do_entropy=1;}
            // note that for the Zeldovich problem, either the gravity or kinetic energy switch is sufficient for good resolution;
            //  both are not needed. we find slightly cleaner results on that test keeping the gravity and removing the KE switch
            
            // also check for flows which are totally dominated by the adiabatic component of their temperature evolution //
            // double mach = fabs(CellP[i].MaxSignalVel/CellP[i].effective_soundspeed() - 2.0); //
            // if(mach < 1.1) {do_entropy=1;} // (actually, this switch tends to do more harm than good!) //
            //do_entropy = 0; // seems unstable in tests like interacting blastwaves... //
            if(do_entropy)
            {
                /* use the pure-SPH entropy equation, which is exact up to the mass flux, for adiabatic flows */
                CellP[i].DtInternalEnergy = -(CellP[i].Pressure/CellP[i].Density) * P[i].Particle_DivVel*All.cf_a2inv;
#ifdef MAGNETIC
                CellP[i].DtB = CellP[i].B * ((1./3.) * All.cf_atime * P[i].Particle_DivVel*All.cf_a2inv);
#ifdef DIVBCLEANING_DEDNER
                CellP[i].DtPhi = (1./3.) * (CellP[i].Phi*All.cf_a3inv) * P[i].Particle_DivVel*All.cf_a2inv; // cf_a3inv from mass-based phi-fluxes
#endif
#endif
                if(All.ComovingIntegrationOn) {CellP[i].DtInternalEnergy -= 3*(CellP[i].Gamma-1) * CellP[i].InternalEnergyPred * All.cf_hubble_a;}
                dEnt = CellP[i].InternalEnergy + CellP[i].DtInternalEnergy * dt_hydrokick; /* gravity term not included here, as it makes this unstable */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                CellP[i].dMass = CellP[i].DtMass = 0;
#endif
            }
#endif // closes ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            
#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
            CellP[i].Density_ExplicitInt *= exp(-DMIN(1.5,DMAX(-1.5,P[i].Particle_DivVel*All.cf_a2inv * dt_hydrokick))); /*!< explicitly integrated volume/density variable to be used if integrating the SPH-like form of the continuity directly */
            if(CellP[i].FaceClosureError > 0) {double drho2 = CellP[i].Gradients.Density.norm_sq(); /* the evolved density evolves back to the explicit density on a relaxation time of order the sound-crossing or tension wave-crossing time across the density gradient length */
                if(drho2>0 && CellP[i].Density_ExplicitInt>0 && CellP[i].Density>0) {
                    double Lgrad = CellP[i].Density / sqrt(drho2); Lgrad=DMAX(Lgrad,P[i].KernelRadius); double cs_eff_forrestoringforce=CellP[i].effective_soundspeed(); /* gradient scale length and sound speed */
#if defined(EOS_TILLOTSON)
                    cs_eff_forrestoringforce=DMIN(cs_eff_forrestoringforce , sqrt(All.Tillotson_EOS_params[CellP[i].CompositionType][10] / CellP[i].Density)); /* speed of deviatoric waves, which is most relevant, if defined */
#endif
                    double delta = 0.1 * dt_hydrokick * cs_eff_forrestoringforce / Lgrad, q0=log(CellP[i].Density_ExplicitInt), q1=log(P[i].Mass/CellP[i].FaceClosureError), qn=0; if(delta > 0.005) {qn=q0*exp(-delta) + q1*(1.-exp(-delta));} else {qn=q0 + (q1-q0)*delta*(1.-0.5*delta);} /* evolves in log-space across this span */
                    CellP[i].Density_ExplicitInt = exp(q0); /* set final density */
                }}
#endif

#ifdef RADTRANSFER /* block here to deal with tricky cases where radiation energy density is -much- larger than thermal, re-distribute the energy that would have taken us negative in gas back into radiation */
            int kfreq; double erad_tot=0,emin=0,enew=0,demin=0,dErad=0,rsol_fac=C_LIGHT_CODE_REDUCED/C_LIGHT_CODE;  for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad_tot+=CellP[i].Rad_E_gamma[kfreq];}
            if(erad_tot > 0) // do some checks if this helps or hurts (identical setup in predict)
            {
                double u_before_radblock = CellP[i].InternalEnergy;
                demin=0.025*CellP[i].InternalEnergy; emin=0.025*(erad_tot/rsol_fac + CellP[i].InternalEnergy*P[i].Mass); enew=DMAX(erad_tot/rsol_fac + dEnt*P[i].Mass, emin);
                dEnt=(enew - erad_tot/rsol_fac) / P[i].Mass; if(dEnt < demin) {dErad=rsol_fac*(dEnt-demin); dEnt=demin;}
                if(dErad<-0.975*erad_tot) {dErad=-0.975*erad_tot;} CellP[i].InternalEnergy = dEnt; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {CellP[i].Rad_E_gamma[kfreq] *= 1 + dErad/erad_tot;}
            } else {
                if(dEnt < 0.5*CellP[i].InternalEnergy) {CellP[i].InternalEnergy *= 0.5;} else {CellP[i].InternalEnergy = dEnt;}
            }
#else
            if(dEnt < 0.5*CellP[i].InternalEnergy) {CellP[i].InternalEnergy *= 0.5;} else {CellP[i].InternalEnergy = dEnt;}
#endif
            CellP[i].enforce_temperature_floor(); /* if we've fallen below the minimum temperature, force the 'floor' */
        }
        
        /* now, kick for non-gas/fluid quantities (accounting for momentum conservation if masses are changing) */
        dp = {};
        if(P[i].Type==0)
        {
            dp += P[i].GravAccel * (mass_pred * dt_gravkick);
            dp += CellP[i].HydroAccel * (mass_pred * All.cf_atime * dt_hydrokick); // convert to code units
#ifdef TURB_DRIVING
            dp += CellP[i].TurbAccel * (mass_pred * dt_gravkick);
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            dp += CellP[i].Rad_Accel * (mass_pred * All.cf_atime * dt_hydrokick);
#endif
        } else {
            dp += P[i].GravAccel * (mass_pred * dt_gravkick);
        }
#if (SINGLE_STAR_TIMESTEPPING > 0)  //if we're super-timestepping, the above accounts for the change in COM velocity. Now we do the internal binary velocity change
        if((P[i].Type == 5) && (P[i].SuperTimestepFlag>=2)) {dp += (P[i].COM_GravAccel - P[i].GravAccel) * (mass_pred * dt_gravkick);}
#endif
#ifdef HERMITE_INTEGRATION
        // we augment this to a whole-step kick for the initial Hermite prediction step, which is done alongside the first half-step kick.
        if((1<<P[i].Type) & HERMITE_INTEGRATION)
        {
            if(mode == 0)
            {
                P[i].OldVel = P[i].Vel;
                P[i].OldPos = P[i].Pos;
                P[i].OldJerk = P[i].GravJerk;
                P[i].Hermite_OldAcc = P[i].GravAccel; // this is the value from the first Hermite tree pass for this timestep
            }
        }
#endif
        P[i].Vel += dp * (1.0 / mass_new); /* correctly accounts for mass change if its allowed */

#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
        double dilation = timestep_dilation_factor(i, P); /* f = 1/a <= 1 */
        if(dilation < 1.) {
            /* the kick above spanned only the fraction f of the raw interval, since dt_gravkick
               already carries the f. add back the mean kick of the surroundings over the remaining
               (1-f) of the raw interval, so that only the relative dynamics are dilated */
            double cfac = dt_gravkick * (1./dilation - 1.);
            P[i].Vel += P[i].acc_of_nearest_special * cfac;
        }
#endif

 
        /* check for reflecting or outflow or otherwise special boundaries: if so, do the reflection/boundary! */
        apply_special_boundary_conditions(i,mass_new,1);
        if(P[i].Mass <= 0 || !isfinite(P[i].Mass)) {return;} /* exit if we have zero'd the particle mass, to avoid errors with dividing by zero */

        /* any other gas-specific kicks (e.g. B-fields, radiation) go here */
        if(P[i].Type==0)
        {
            do_kick_for_extra_physics(i, tstart, tend, dt_entr);

            /* after completion of a full step, set the predicted values of gas/fluid quantities
             * to the current values. They will then predicted further along in drift operations */
            if(mode==1)
            {
#ifdef HYDRO_GENERATE_TARGET_MESH // it is often desirable to damp transient velocities when setting up a stable mesh: do so here by un-commenting the line below //
                //for(j=0;j<3;j++) {P[i].Vel[j] *= exp(-0.15);} // coefficient is constant per-timestep: adjust to make as aggressive or weak as desired //
#endif

                CellP[i].VelPred = P[i].Vel; //(mass_old*v_old[j] + dp[j]) / mass_new;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                P[i].Mass = CellP[i].MassTrue; CellP[i].Mass = P[i].Mass; //mass_old + CellP[i].DtMass * dt_hydrokick;
#endif
                CellP[i].InternalEnergyPred = CellP[i].InternalEnergy; //ent_old + CellP[i].DtInternalEnergy * dt_entr;
#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
                CellP[i].Density = CellP[i].Density_ExplicitInt; /*!< explicitly integrated volume/density variable to be used if integrating the SPH-like form of the continuity directly */
#endif
            }
        }
        
        /* set the momentum shift so we know how to move the tree! */
        P[i].dp += dp;
#ifdef DM_FUZZY
        do_dm_fuzzy_drift_kick(i, dt_entr, 0); /* kicks for fuzzy-dm integration */
#endif
        
    } // if(TimeBinActive[P[i].TimeBin]) //
}


void set_predicted_quantities_for_extra_physics(int i)
{
    if(P[i].Type == 0 && P[i].Mass > 0)
    {
        int k, kf; k=0, kf=0;
#if defined(MAGNETIC)
#ifndef MHD_ALTERNATIVE_LEAPFROG_SCHEME
        CellP[i].BPred = CellP[i].B;
#if defined(DIVBCLEANING_DEDNER)
        CellP[i].PhiPred = CellP[i].Phi;
#endif
#endif
#endif
#ifdef COSMIC_RAY_FLUID
        for(kf=0;kf<N_CR_PARTICLE_BINS;kf++)
        {
            CellP[i].CosmicRayEnergyPred[kf] = CellP[i].CosmicRayEnergy[kf];
            CellP[i].CosmicRayFluxPred[kf] = CellP[i].CosmicRayFlux[kf];
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(k=0;k<2;k++) {CellP[i].CosmicRayAlfvenEnergyPred[kf][k] = CellP[i].CosmicRayAlfvenEnergy[kf][k];}
#endif
        }
#endif
        
#if defined(RT_EVOLVE_ENERGY)
        for(kf=0;kf<N_RT_FREQ_BINS;kf++)
        {
            CellP[i].Rad_E_gamma_Pred[kf] = CellP[i].Rad_E_gamma[kf];
#if defined(RT_EVOLVE_FLUX)
            CellP[i].Rad_Flux_Pred[kf] = CellP[i].Rad_Flux[kf];
#endif
        }
        rt_eddington_update_calculation(i, P, CellP);
#endif
#ifdef RT_EVOLVE_INTENSITIES
        for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k=0;k<N_RT_INTENSITY_BINS;k++) {CellP[i].Rad_Intensity_Pred[kf][k] = CellP[i].Rad_Intensity[kf][k];}}
#endif

#ifdef EOS_ELASTIC
        CellP[i].Elastic_Stress_Tensor_Pred = CellP[i].Elastic_Stress_Tensor;
#endif
        
        set_eos_pressure(i, P, CellP);
    }
}



void do_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
    int j; j=0;
#ifdef MAGNETIC
#ifndef MHD_ALTERNATIVE_LEAPFROG_SCHEME
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    CellP[i].B += CellP[i].DtB * (dt_entr * BphysVolphys_to_BcodeVolCode); // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1 / All.cf_a3inv; // for mass-based phi-fluxes (otherwise is just "1")
    /* phi units are [vcode][Bcode]=a^3 * vphys*Bphys */
    if(CellP[i].Density > 0)
    {
        /* now we're going to check for physically reasonable phi values */
        double cs_phys = CellP[i].effective_soundspeed();
        double b_phys = sqrt(CellP[i].Bfield().norm_sq())*All.cf_a2inv;
        double vsig1 = sqrt(cs_phys*cs_phys + b_phys*b_phys/(CellP[i].Density*All.cf_a3inv));
        double vsig2 = 0.5 * fabs(CellP[i].MaxSignalVel);
        double vsig_max = DMAX( DMAX(vsig1,vsig2) , All.FastestWaveSpeed );
        double phi_phys_abs = fabs(Get_Gas_PhiField(i)) * All.cf_a3inv;
        double vb_phy_abs = vsig_max * b_phys;

        if((!isnan(CellP[i].DtPhi))&&(phi_phys_abs>0)&&(vb_phy_abs>0)&&(!isnan(phi_phys_abs))&&(!isnan(vb_phy_abs)))
        {
            double phi_max_tolerance = 10.0;
            if(phi_phys_abs > 1000. * phi_max_tolerance * vb_phy_abs)
            {
                /* this can indicate a problem! issue a warning and zero phi */
                if(phi_phys_abs > 1.0e6 * phi_max_tolerance * vb_phy_abs) {
                    PRINT_WARNING("significant growth detected in phi-field: phi_phys_abs=%g vb_phy_abs=%g vsig_max=%g b_phys=%g particle_id_i=%d dtphi_code=%g Pressure=%g rho=%g x/y/z=%g/%g/%g vx/vy/vz=%g/%g/%g Bx/By/Bz=%g/%g/%g h=%g u=%g m=%g phi=%g bin=%d SigVel=%g a=%g \n",
                       phi_phys_abs,vb_phy_abs,vsig_max,b_phys,i,CellP[i].DtPhi,CellP[i].Pressure,CellP[i].Density,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],
                       P[i].Vel[0],P[i].Vel[1],P[i].Vel[2],CellP[i].B[0],CellP[i].B[1],CellP[i].B[2],
                       P[i].KernelRadius,CellP[i].InternalEnergy,P[i].Mass,CellP[i].Phi,P[i].TimeBin,CellP[i].MaxSignalVel,All.cf_atime);}
                CellP[i].PhiPred = CellP[i].Phi = CellP[i].DtPhi = 0;
            } else {
                if(phi_phys_abs > phi_max_tolerance * vb_phy_abs)
                {
                    /* in this limit, only allow for decay of phi: to avoid over-shooting, we apply the force as damping */
                    if(CellP[i].Phi > 0) {CellP[i].DtPhi=DMIN(CellP[i].DtPhi,0);} else {CellP[i].DtPhi=DMAX(CellP[i].DtPhi,0);}
                    double dtphi_code = dt_entr * PhiphysVolphys_to_PhicodeVolCode * CellP[i].DtPhi;
                    if(CellP[i].Phi != 0) {CellP[i].Phi *= exp( - fabs(dtphi_code) / fabs(CellP[i].Phi) );}
                } else {
                    /* ok, in this regime, we're safe to apply the 'normal' time evolution */
                    double dtphi_code = dt_entr * PhiphysVolphys_to_PhicodeVolCode * CellP[i].DtPhi;
                    CellP[i].Phi += dtphi_code;
                }
            }
        }
    } else {
        CellP[i].Phi = CellP[i].PhiPred = CellP[i].DtPhi = 0;
    }
    /* now apply the usual damping term */
    double t_damp = Get_Gas_PhiField_DampingTimeInv(i);
    if((t_damp>0) && (!isnan(t_damp)) && (dt_entr>0))
    {
        CellP[i].Phi *= exp( -dt_entr * t_damp );
    }
    if(isnan(CellP[i].DtPhi)) {CellP[i].DtPhi=0;}
    if(isnan(CellP[i].Phi)) {CellP[i].Phi=0;}
    if(isnan(CellP[i].PhiPred)) {CellP[i].PhiPred=CellP[i].Phi;}
#endif
#endif
#endif
    
#ifdef COSMIC_RAY_FLUID
#ifndef TRANSPORT_SUBCYCLE
    CosmicRay_Update_DriftKick(i,dt_entr,0,P,CellP);
#endif
#endif

#ifdef RADTRANSFER
#ifndef TRANSPORT_SUBCYCLE  /* when subcycling, RT kicks are handled in the subcycle loop */
    rt_update_driftkick(i,dt_entr,0,P,CellP);
#endif
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
    if(P[i].Pos[2] > DMIN(19., DMAX(1.1*All.Time*C_LIGHT_CODE_REDUCED, DMIN(18.*boxSize_X + (All.Vertical_Grain_Accel*All.Dust_to_Gas_Mass_Ratio - All.Vertical_Gravity_Strength)*All.Time*All.Time/2., 19.)))) {for(j=0;j<N_RT_FREQ_BINS;j++) {CellP[i].Rad_E_gamma[j]*=0.5; CellP[i].Rad_E_gamma_Pred[j]*=0.5;
#ifdef RT_EVOLVE_FLUX
        if(CellP[i].Rad_Flux[j][2] < 0) {CellP[i].Rad_Flux[j][2]=-CellP[i].Rad_Flux[j][2]; CellP[i].Rad_Flux_Pred[j][2]=CellP[i].Rad_Flux[j][2];}
#endif
    }}
#endif
#endif // RADTRANSFER

#ifdef EOS_ELASTIC
    elastic_body_update_driftkick(i,dt_entr,0);
#endif
}

    
    
void apply_special_boundary_conditions(int i, double mass_for_dp, int mode) { apply_special_boundary_conditions_P(i, mass_for_dp, mode, P, CellP); }
