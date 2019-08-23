#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO 
 * (added energy/entropy switch, terms for explicit mass conservation in mass fluxes, 
 *  and updates to additional fluid variables as needed)
 */

void apply_long_range_kick(integertime, integertime);

void do_first_halfstep_kick(void)
{
    int i;
    integertime ti_step, tstart=0, tend=0;
    
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
    
    /* as currently written with some revisions to MFV methods, should only update on active timesteps */
    for(i = 0; i < NumPart; i++)
    {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        if((TimeBinActive[P[i].TimeBin]) || (P[i].Type==0)) /* active OR gas, need to check each timestep to ensure manifest conservation */
#else
        if(TimeBinActive[P[i].TimeBin]) /* 'full' kick for active particles */
#endif
        {
            if(P[i].Mass > 0)
            {
                ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
                tstart = P[i].Ti_begstep;	/* beginning of step */
                tend = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */
                do_the_kick(i, tstart, tend, P[i].Ti_current, 0);
            }
        }
    } // for(i = 0; i < NumPart; i++) //
}

void do_second_halfstep_kick(void)
{
    int i;
    integertime ti_step, tstart=0, tend=0;
    
#ifdef PMGRID
    if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
        tstart = All.PM_Ti_begstep + ti_step / 2;
        tend = tstart + ti_step / 2;
        apply_long_range_kick(tstart, tend);
    }
#endif

    for(i = 0; i < NumPart; i++)
    {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        if((TimeBinActive[P[i].TimeBin]) || (P[i].Type==0)) /* active OR gas, need to check each timestep to ensure manifest conservation */
#else
        if(TimeBinActive[P[i].TimeBin]) /* 'full' kick for active particles */
#endif
        {
            if(P[i].Mass > 0)
            {
                ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
                tstart = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */
                tend = P[i].Ti_begstep + ti_step;	/* end of step */
                do_the_kick(i, tstart, tend, P[i].Ti_current, 1);
                set_predicted_sph_quantities_for_extra_physics(i);
            }
        }
    } // for(i = 0; i < NumPart; i++) //
    
#ifdef TURB_DRIVING
    do_turb_driving_step_second_half();
#endif
}



#ifdef PMGRID
void apply_long_range_kick(integertime tstart, integertime tend)
{
    int i, j;
    double dt_gravkick, dvel[3];
    
    if(All.ComovingIntegrationOn)
        dt_gravkick = get_gravkick_factor(tstart, tend);
    else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;
    
    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Mass > 0)
            for(j = 0; j < 3; j++)	/* do the kick, only collisionless particles */
            {
                dvel[j] = P[i].GravPM[j] * dt_gravkick;
                P[i].Vel[j] += dvel[j];
                P[i].dp[j] += P[i].Mass * dvel[j];
            }
#ifdef GDE_DISTORTIONTENSOR
        do_long_range_phase_space_kick(i, dt_gravkick);
#endif
    }
}
#endif


void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent, int mode)
{
    int j;
    double dp[3], dt_entr, dt_gravkick, dt_hydrokick;
    double mass_old, mass_pred, mass_new;
    mass_old = mass_pred = mass_new = P[i].Mass;    
    
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    /* need to do the slightly more complicated update scheme to maintain exact mass conservation */
    if(P[i].Type==0)
    {
        if(SphP[i].dMass != 0) //ent_old = SphP[i].InternalEnergy; for(j=0;j<3;j++) v_old[j] = P[i].Vel[j];
        {
            double dMass=0; // fraction of delta_conserved to couple per kick step (each 'kick' is 1/2-timestep) // double dv[3], v_old[3], dMass, ent_old=0, d_inc = 0.5;
            if(mode != 0) // update the --conserved-- variables of each particle //
            {
                dMass = (tend - tstart) * All.Timebase_interval / All.cf_hubble_a * SphP[i].DtMass; if(dMass * SphP[i].dMass < 0) {dMass = 0;} // slope-limit: no opposite reconstruction! //
                if((fabs(dMass) > fabs(SphP[i].dMass))) {dMass = SphP[i].dMass;} // try to get close to what the time-integration scheme would give //
                SphP[i].dMass -= dMass;
            } else {dMass = SphP[i].dMass;}
            if(dMass < -0.99*SphP[i].MassTrue) {dMass = -0.99*SphP[i].MassTrue;} // limiter to prevent madness //

            /* load and update the particle masses : particle mass update here, from hydro fluxes */
            mass_old = SphP[i].MassTrue; mass_pred = P[i].Mass; mass_new = mass_old + dMass; SphP[i].MassTrue = mass_new; // UNITS: remember all time derivatives (DtX, dX) are in -physical- units; as are mass, entropy/internal energy, but -not- velocity //
            /* double e_old = mass_old * SphP[i].InternalEnergy; for(j = 0; j< 3; j++) e_old += 0.5*mass_old * (P[i].Vel[j]/All.cf_atime)*(P[i].Vel[j]/All.cf_atime); // physical //
            for(j = 0; j < 3; j++) // momentum-space-kick
            {
                dp[j] = d_inc * SphP[i].dMomentum[j]; // now update the velocity based on the total momentum change
                P[i].Vel[j] = (mass_old*P[i].Vel[j] + dp[j]*All.cf_atime) / mass_new; // call after tabulating dP[j] //
            } // kick for gas internal energy/entropy
            e_old += d_inc * SphP[i].dInternalEnergy; // for(j = 0; j< 3; j++) e_old -= 0.5*mass_new * (P[i].Vel[j]/All.cf_atime)*(P[i].Vel[j]/All.cf_atime); // increment of total (thermal+kinetic) energy; subtract off the new kinetic energy //
            SphP[i].InternalEnergy = e_old / mass_new; check_particle_for_temperature_minimum(i); // obtain the new internal energy per unit mass, check floor // */
             
            // at the end of this kick, need to re-zero the dInternalEnergy, and other conserved-variable SPH quantities set in the hydro loop, to avoid double-counting them
            if(mode==0) {SphP[i].dMass=0;} /* SphP[i].dInternalEnergy=0; SphP[i].dMomentum[0]=SphP[i].dMomentum[1]=SphP[i].dMomentum[2]=0; */
        }
    } // if(P[i].Type==0) //
#endif
    
    /* only enter the 'normal' kick loop below for genuinely active particles */
    if(TimeBinActive[P[i].TimeBin])
    {
        /* get the timestep (physical units for dt_entr and dt_hydrokick) */
        if(All.ComovingIntegrationOn)
        {
            dt_entr = dt_hydrokick = (tend - tstart) * All.Timebase_interval / All.cf_hubble_a;
            dt_gravkick = get_gravkick_factor(tstart, tend);
        }
        else
        {
            dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
        }
        
        if(P[i].Type==0)
        {
            double grav_acc[3], dEnt_Gravity = 0;
            for(j = 0; j < 3; j++)
            {
                grav_acc[j] = All.cf_a2inv * P[i].GravAccel[j];
#ifdef PMGRID
                grav_acc[j] += All.cf_a2inv * P[i].GravPM[j];
#endif
            }
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
            for(j=0;j<3;j++) {dEnt_Gravity += -(SphP[i].GravWorkTerm[j] * All.cf_atime * dt_hydrokick) * grav_acc[j];}
#endif
            double du_tot = SphP[i].DtInternalEnergy * dt_hydrokick + dEnt_Gravity;
#if defined(COOLING) && !defined(COOLING_OPERATOR_SPLIT)
            if(mode == 1) {du_tot = 0;}
#endif
            double dEnt = SphP[i].InternalEnergy + du_tot;
            
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            /* if we're using a Riemann solver, we include an energy/entropy-type switch to ensure
                that we don't corrupt the temperature evolution of extremely cold, adiabatic flows */
            /* MHD tests suggest that this switch often does more harm than good: we will
             pay the price of noisier temperature fields (corrupting them when c_s << v_A << v_bulk)
             and they are dynamically irrelevant, in exchange for avoiding potentially much more
             serious errors if this tripped when the B-fields were important */
            double e_thermal,e_kinetic,e_potential;
            e_potential=0; for(j=0;j<3;j++) {e_potential += grav_acc[j]*grav_acc[j];}
            e_potential = P[i].Mass * sqrt(e_potential) * (Get_Particle_Size(i)*All.cf_atime); // = M*|a_grav|*h (physical)
            e_kinetic = 0.5 * P[i].Mass * All.cf_a2inv * SphP[i].MaxKineticEnergyNgb;
            e_thermal = DMAX(0.5*SphP[i].InternalEnergy, dEnt) * P[i].Mass;
#ifdef MAGNETIC
            for(j=0;j<3;j++) {e_thermal += 0.5*SphP[i].B[j]*SphP[i].B[j]*SphP[i].Density/(All.cf_atime*P[i].Mass);}
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
            // double mach = fabs(SphP[i].MaxSignalVel/Particle_effective_soundspeed_i(i) - 2.0); //
            // if(mach < 1.1) {do_entropy=1;} // (actually, this switch tends to do more harm than good!) //
            //do_entropy = 0; // seems unstable in tests like interacting blastwaves... //
            if(do_entropy)
            {
                /* use the pure-SPH entropy equation, which is exact up to the mass flux, for adiabatic flows */
                SphP[i].DtInternalEnergy = -(SphP[i].Pressure/SphP[i].Density) * P[i].Particle_DivVel*All.cf_a2inv;
#ifdef MAGNETIC
                for(j=0;j<3;j++)
                {
                    SphP[i].DtB[j] = (1./3.) * SphP[i].B[j]*All.cf_atime * P[i].Particle_DivVel*All.cf_a2inv;
                }
#ifdef DIVBCLEANING_DEDNER
                SphP[i].DtPhi = (1./3.) * (SphP[i].Phi*All.cf_a3inv) * P[i].Particle_DivVel*All.cf_a2inv; // cf_a3inv from mass-based phi-fluxes
#endif
#endif
                if(All.ComovingIntegrationOn) SphP[i].DtInternalEnergy -= 3*GAMMA_MINUS1 * SphP[i].InternalEnergyPred * All.cf_hubble_a;
                dEnt = SphP[i].InternalEnergy + SphP[i].DtInternalEnergy * dt_hydrokick; /* gravity term not included here, as it makes this unstable */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[i].dMass = SphP[i].DtMass = 0;
#endif
            }
#endif // closes ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            
            if(dEnt < 0.5*SphP[i].InternalEnergy) {SphP[i].InternalEnergy *= 0.5;} else {SphP[i].InternalEnergy = dEnt;}
            check_particle_for_temperature_minimum(i); /* if we've fallen below the minimum temperature, force the 'floor' */
        }
        
        /* now, kick for non-SPH quantities (accounting for momentum conservation if masses are changing) */
        for(j = 0; j < 3; j++)
        {
            dp[j] = 0;
            if(P[i].Type==0)
            {
                dp[j] += mass_pred * SphP[i].HydroAccel[j] * All.cf_atime * dt_hydrokick; // convert to code units
#ifdef TURB_DRIVING
                dp[j] += mass_pred * SphP[i].TurbAccel[j] * dt_gravkick;
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
                dp[j] += mass_pred * SphP[i].RadAccel[j] * All.cf_atime * dt_hydrokick;
#endif
            }
            dp[j] += mass_pred * P[i].GravAccel[j] * dt_gravkick;
            P[i].Vel[j] += dp[j] / mass_new; /* correctly accounts for mass change if its allowed */
        }

 
        /* check for reflecting boundaries: if so, do the reflection! */
#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z)
        double box_upper[3]; box_upper[0]=box_upper[1]=box_upper[2]=1;
#ifdef BOX_PERIODIC
        box_upper[0]=boxSize_X; box_upper[1]=boxSize_Y; box_upper[2]=boxSize_Z;
#endif
        for(j = 0; j < 3; j++)
        {
            /* skip the non-reflecting boundaries */
#ifndef BOX_REFLECT_X
            if(j==0) continue;
#endif
#ifndef BOX_REFLECT_Y
            if(j==1) continue;
#endif
#ifndef BOX_REFLECT_Z
            if(j==2) continue;
#endif
            if(P[i].Pos[j] <= 0)
            {
                if(P[i].Vel[j]<0) {P[i].Vel[j]=-P[i].Vel[j]; if(P[i].Type==0) {SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0;} dp[j]+=2*P[i].Vel[j]*mass_new;}
                P[i].Pos[j]=(0.+((double)P[i].ID)*1.e-9)*box_upper[j];
            }
            if(P[i].Pos[j] >= box_upper[j])
            {
                if(P[i].Vel[j]>0) {P[i].Vel[j]=-P[i].Vel[j]; if(P[i].Type==0) {SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0;} dp[j]+=2*P[i].Vel[j]*mass_new;}
                P[i].Pos[j]=box_upper[j]*(1.-((double)P[i].ID)*1.e-9);
            }
        }
#endif
        /* any other gas-specific kicks (e.g. B-fields, radiation) go here */
        if(P[i].Type==0)
        {
            do_sph_kick_for_extra_physics(i, tstart, tend, dt_entr);

            /* after completion of a full step, set the predicted values of SPH quantities
             * to the current values. They will then predicted further along in drift operations */
            if(mode==1)
            {
#ifdef HYDRO_GENERATE_TARGET_MESH // it is often desirable to damp transient velocities when setting up a stable mesh: do so here by un-commenting the line below //
                //for(j=0;j<3;j++) {P[i].Vel[j] *= exp(-0.15);} // coefficient is constant per-timestep: adjust to make as aggressive or weak as desired //
#endif

                for(j = 0; j < 3; j++)
                    SphP[i].VelPred[j] = P[i].Vel[j];//(mass_old*v_old[j] + dp[j]) / mass_new;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                P[i].Mass = SphP[i].MassTrue; //mass_old + SphP[i].DtMass * dt_hydrokick;
#endif
                SphP[i].InternalEnergyPred = SphP[i].InternalEnergy; //ent_old + SphP[i].DtInternalEnergy * dt_entr;
            }
        }
        
        /* set the momentum shift so we know how to move the tree! */
        for(j=0;j<3;j++) {P[i].dp[j] += dp[j];}
#ifdef GDE_DISTORTIONTENSOR
        do_the_phase_space_kick(i, dt_gravkick); /* momentum-space correction for following phase-space distribution (call after momentum-space kicks) */
#endif
#ifdef DM_FUZZY
        do_dm_fuzzy_drift_kick(i, dt_entr, 0); /* kicks for fuzzy-dm integration */
#endif
        
    } // if(TimeBinActive[P[i].TimeBin]) //
}


void set_predicted_sph_quantities_for_extra_physics(int i)
{
    if(P[i].Type == 0 && P[i].Mass > 0)
    {
        int k, kf; k=0, kf=0;
#if defined(MAGNETIC)
#ifndef MHD_ALTERNATIVE_LEAPFROG_SCHEME
        for(k=0;k<3;k++) {SphP[i].BPred[k] = SphP[i].B[k];}
#if defined(DIVBCLEANING_DEDNER)
        SphP[i].PhiPred = SphP[i].Phi;
#endif
#endif
#endif
        
#if defined(RT_EVOLVE_NGAMMA)
        for(kf=0;kf<N_RT_FREQ_BINS;kf++)
        {
            SphP[i].E_gamma_Pred[kf] = SphP[i].E_gamma[kf];
#if defined(RT_EVOLVE_FLUX)
            for(k=0;k<3;k++) SphP[i].Flux_Pred[kf][k] = SphP[i].Flux[kf][k];
#endif
        }
        rt_eddington_update_calculation(i);
#endif
#ifdef RT_EVOLVE_INTENSITIES
        for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k=0;k<N_RT_INTENSITY_BINS;k++) {SphP[i].Intensity_Pred[kf][k] = SphP[i].Intensity[kf][k];}}
#endif

#ifdef EOS_ELASTIC
        for(k=0;k<3;k++) {for(kf=0;kf<3;kf++) {SphP[i].Elastic_Stress_Tensor_Pred[k][kf]=SphP[i].Elastic_Stress_Tensor[k][kf];}}
#endif
        
        SphP[i].Pressure = get_pressure(i);
#ifdef EOS_ENFORCE_ADIABAT
        SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = SphP[i].Pressure / (SphP[i].Density * GAMMA_MINUS1);
#endif
    }
}



void do_sph_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
    int j; j=0;
#ifdef MAGNETIC
#ifndef MHD_ALTERNATIVE_LEAPFROG_SCHEME
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    for(j = 0; j < 3; j++) {SphP[i].B[j] += SphP[i].DtB[j] * dt_entr * BphysVolphys_to_BcodeVolCode;} // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1 / All.cf_a3inv; // for mass-based phi-fluxes (otherwise is just "1")
    /* phi units are [vcode][Bcode]=a^3 * vphys*Bphys */
    if(SphP[i].Density > 0)
    {
        /* now we're going to check for physically reasonable phi values */
        double cs_phys = All.cf_afac3 * Particle_effective_soundspeed_i(i);
        double b_phys = 0.0;
        for(j = 0; j < 3; j++) {b_phys += Get_Particle_BField(i,j)*Get_Particle_BField(i,j);}
        b_phys = sqrt(b_phys)*All.cf_a2inv;
        double vsig1 = sqrt(cs_phys*cs_phys + b_phys*b_phys/(SphP[i].Density*All.cf_a3inv));
        double vsig2 = 0.5 * All.cf_afac3 * fabs(SphP[i].MaxSignalVel);
        double vsig_max = DMAX( DMAX(vsig1,vsig2) , All.FastestWaveSpeed );
        double phi_phys_abs = fabs(Get_Particle_PhiField(i)) * All.cf_a3inv;
        double vb_phy_abs = vsig_max * b_phys;

        if((!isnan(SphP[i].DtPhi))&&(phi_phys_abs>0)&&(vb_phy_abs>0)&&(!isnan(phi_phys_abs))&&(!isnan(vb_phy_abs)))
        {
            double phi_max_tolerance = 10.0;
            if(phi_phys_abs > 1000. * phi_max_tolerance * vb_phy_abs)
            {
                /* this indicates a serious problem! issue a warning and zero phi */
#ifndef IO_REDUCED_MODE
                printf("WARNING: MAJOR GROWTH IN PHI-FIELD: phi_phys_abs=%g vb_phy_abs=%g vsig_max=%g b_phys=%g particle_id_i=%d dtphi_code=%g Pressure=%g rho=%g x/y/z=%g/%g/%g vx/vy/vz=%g/%g/%g Bx/By/Bz=%g/%g/%g h=%g u=%g m=%g phi=%g bin=%d SigVel=%g a=%g \n",
                       phi_phys_abs,vb_phy_abs,vsig_max,b_phys,i,SphP[i].DtPhi,
                       SphP[i].Pressure,SphP[i].Density,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],
                       P[i].Vel[0],P[i].Vel[1],P[i].Vel[2],SphP[i].B[0],SphP[i].B[1],SphP[i].B[2],
                       PPP[i].Hsml,SphP[i].InternalEnergy,P[i].Mass,SphP[i].Phi,P[i].TimeBin,
                       SphP[i].MaxSignalVel,All.cf_atime);
                fflush(stdout);
#endif
                SphP[i].PhiPred = SphP[i].Phi = SphP[i].DtPhi = 0;
            } else {
                if(phi_phys_abs > phi_max_tolerance * vb_phy_abs)
                {
                    /* in this limit, only allow for decay of phi: to avoid over-shooting, we apply the force as damping */
                    if(SphP[i].Phi > 0) {SphP[i].DtPhi=DMIN(SphP[i].DtPhi,0);} else {SphP[i].DtPhi=DMAX(SphP[i].DtPhi,0);}
                    double dtphi_code = dt_entr * PhiphysVolphys_to_PhicodeVolCode * SphP[i].DtPhi;
                    if(SphP[i].Phi != 0) {SphP[i].Phi *= exp( - fabs(dtphi_code) / fabs(SphP[i].Phi) );}
                } else {
                    /* ok, in this regime, we're safe to apply the 'normal' time evolution */
                    double dtphi_code = dt_entr * PhiphysVolphys_to_PhicodeVolCode * SphP[i].DtPhi;
                    SphP[i].Phi += dtphi_code;
                }
            }
        }
        if(isnan(SphP[i].DtPhi)) {SphP[i].DtPhi=0;}
        if(isnan(SphP[i].Phi)) {SphP[i].Phi=0;}
        if(isnan(SphP[i].PhiPred)) {SphP[i].PhiPred=SphP[i].Phi;}
    } else {
        SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
    }
    /* now apply the usual damping term */
    double t_damp = Get_Particle_PhiField_DampingTimeInv(i);
    if((t_damp>0) && (!isnan(t_damp)))
    {
        SphP[i].Phi *= exp( -dt_entr * t_damp );
    }
#endif
#endif
#endif
    
    
    
#ifdef RADTRANSFER
    rt_update_driftkick(i,dt_entr,0);
#endif

#ifdef EOS_ELASTIC
    elastic_body_update_driftkick(i,dt_entr,0);
#endif
}

