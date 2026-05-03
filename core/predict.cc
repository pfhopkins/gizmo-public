#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"

/*! Routines for the drift/predict step */

/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially in detail (although the highest-level algorithm
 * structure remains essentially the same)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO, and many new
 * options and subroutines added for flexibility with different
 * hydro solvers, timestepping schemes, boundary conditions, and
 * different mesh-motion options added.
 */

void reconstruct_timebins(void)
{
    int i, bin;
    long long glob_sum;
    
    for(bin = 0; bin < TIMEBINS; bin++)
    {
        TimeBinCount[bin] = 0;
        TimeBinCountGas[bin] = 0;
        FirstInTimeBin[bin] = -1;
        LastInTimeBin[bin] = -1;
#ifdef GALSF
        TimeBinSfr[bin] = 0;
#endif
#ifdef SINK_PARTICLES
        TimeBin_Sink_mass[bin] = 0;
        TimeBin_Sink_dynamicalmass[bin] = 0;
        TimeBin_Sink_Mdot[bin] = 0;
        TimeBin_Sink_Medd[bin] = 0;
#endif
    }
    
    for(i = 0; i < NumPart; i++)
    {
        bin = P[i].TimeBin;
        
        if(TimeBinCount[bin] > 0)
        {
            PrevInTimeBin[i] = LastInTimeBin[bin];
            NextInTimeBin[i] = -1;
            NextInTimeBin[LastInTimeBin[bin]] = i;
            LastInTimeBin[bin] = i;
        }
        else
        {
            FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
            PrevInTimeBin[i] = NextInTimeBin[i] = -1;
        }
        TimeBinCount[bin]++;
        if(P[i].Type == 0)
            TimeBinCountGas[bin]++;
        
#ifdef GALSF
        if(P[i].Type == 0)
            TimeBinSfr[bin] += CellP[i].Sfr;
#endif
#ifdef SINK_PARTICLES
        if(P[i].Type == 5)
        {
            TimeBin_Sink_mass[bin] += P[i].Sink_Mass;
            TimeBin_Sink_dynamicalmass[bin] += P[i].Mass;
            TimeBin_Sink_Mdot[bin] += P[i].Sink_Mdot;
            TimeBin_Sink_Medd[bin] += P[i].Sink_Mdot / P[i].Sink_Mass;
        }
#endif
    }
    
    make_list_of_active_particles();
    
    NumForceUpdate = 0;
    for (int i : ActiveParticleList)
    {
        NumForceUpdate++;
        if(i >= NumPart)
        {
            printf("Bummer i=%d\n", i);
            terminate("inconsistent list");
        }
    }
    
    sumup_large_ints(1, &NumForceUpdate, &glob_sum);
    GlobNumForceUpdate = glob_sum;
}





void drift_particle(int i, integertime time1)
{
    int j __attribute__((unused)); double dt_drift; integertime time0 = P[i].Ti_current;
    if(time1 < time0)
    {
        printf("i=%d time0=%lld time1=%lld\n", i, (long long)time0, (long long)time1);
        terminate("no prediction into past allowed");
    }
    if(time1 == time0) {return;}
    
    dt_drift = get_drift_factor(time0, time1, i, 0);
        
#if !defined(FREEZE_HYDRO)
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    if(P[i].Type==0) {advect_mesh_point(i,dt_drift);} else {P[i].Pos += P[i].Vel * dt_drift;}
#elif (SINGLE_STAR_TIMESTEPPING > 0)
    Vec3<double> fewbody_drift_dx, fewbody_kick_dv; // if super-timestepping, the updates above account for COM motion of the binary; now we account for the internal motion
#ifdef KETJU_REGULARIZATION
    if( (P[i].Type == 5) && (P[i].SuperTimestepFlag>=2) && !P[i].KetjuIntegrated )
#else
    if( (P[i].Type == 5) && (P[i].SuperTimestepFlag>=2) )
#endif
    {
        Vec3<double> COM_Vel = P[i].Vel + P[i].comp_dv * (P[i].comp_Mass/(P[i].Mass+P[i].comp_Mass)); //center of mass velocity
        P[i].Pos += COM_Vel * dt_drift; //center of mass drift
        odeint_super_timestep(i, dt_drift, fewbody_kick_dv, fewbody_drift_dx, 1); // do_fewbody_drift
        P[i].GravAccel = P[i].COM_GravAccel; //Overwrite the acceleration with center of mass value
        P[i].Pos += fewbody_drift_dx; //Keplerian evolution
        P[i].Vel += fewbody_kick_dv; //move on binary.orbit
    } else {
       P[i].Pos += P[i].Vel * dt_drift;
    }
#else
    P[i].Pos += P[i].Vel * dt_drift;
#endif
#endif // FREEZE_HYDRO clause
#if (NUMDIMS==1)
    P[i].Pos[1]=P[i].Pos[2]=0; // force zero-ing
#endif
#if (NUMDIMS==2)
    P[i].Pos[2]=0; // force zero-ing
#endif

#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    double a_fac = return_timestep_dilation_factor(i,0);
    if(a_fac > 1.) {
        double cfac = dt_drift * (1. - 1./a_fac);
        P[i].Pos += P[i].vel_of_nearest_special * cfac; /* add back in the mean drift of the surrounding stuff to dilate only the local dynamics */
    }
#endif

    double divv_fac = P[i].Particle_DivVel * dt_drift;
    /* How far the predicted kernel radius may move in one drift. The two directions are not the
       same risk: growing follows the flow and simply predicts a larger search, while shrinking
       can empty a kernel of the neighbours the next density solve needs. So growth carries the
       looser bound, which also tracks the divergence it approximates more closely than the old
       symmetric 0.3 did -- that tight cap lowers the guess's order and fails the density guess
       more often. Adaptive softening has used the loose bound on both sides for the same reason. */
    double divv_fac_max_grow   = 4.0;
    double divv_fac_max_shrink = 1.5;
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(ags_density_isactive(i) && P[i].Type>0) {divv_fac_max_grow = 4; divv_fac_max_shrink = 4;} // can [should] allow larger changes when using adapting soft for all
#endif
    if(divv_fac > +divv_fac_max_grow)   divv_fac = +divv_fac_max_grow;
    if(divv_fac < -divv_fac_max_shrink) divv_fac = -divv_fac_max_shrink;
    
#ifdef GRAIN_FLUID
    if((1 << P[i].Type) & (GRAIN_PTYPES))
    {
        P[i].KernelRadius *= exp((double)divv_fac / ((double)NUMDIMS));
        if(P[i].KernelRadius < All.MinKernelRadius) {P[i].KernelRadius = All.MinKernelRadius;}
        if(P[i].KernelRadius > All.MaxKernelRadius) {P[i].KernelRadius = All.MaxKernelRadius;}
    }
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(ags_density_isactive(i) && (dt_drift>0)) /* particle is AGS-active */
    {
        double minsoft = ags_return_minsoft(i), maxsoft = ags_return_maxsoft(i);
        P[i].AGS_KernelRadius *= exp((double)divv_fac / ((double)NUMDIMS));
        if(P[i].AGS_KernelRadius < minsoft) {P[i].AGS_KernelRadius = minsoft;}
        if(P[i].AGS_KernelRadius > maxsoft) {P[i].AGS_KernelRadius = maxsoft;}
    } else {P[i].AGS_KernelRadius = ForceSoftening_KernelRadius(i);} /* non-AGS-active particles use fixed softening */
#endif
    
#ifdef DM_FUZZY
    do_dm_fuzzy_drift_kick(i, dt_drift, 1);
#endif
    
    if((P[i].Type == 0) && (P[i].Mass > 0))
        {
            double dt_gravkick, dt_gravkick_pm, dt_hydrokick, dt_entr;
            dt_entr = dt_hydrokick = (time1 - time0) * unit_integertime_in_physical(i);
            dt_gravkick = get_gravkick_factor(time0, time1, i, 0);
            
#ifdef PMGRID
            dt_gravkick_pm = get_gravkick_factor(time0, time1, -1, 0);
            CellP[i].VelPred += P[i].GravAccel*dt_gravkick + P[i].GravPM*dt_gravkick_pm + CellP[i].HydroAccel*(dt_hydrokick*All.cf_atime); /* make sure v is in code units */
#else
            CellP[i].VelPred += P[i].GravAccel * dt_gravkick + CellP[i].HydroAccel * (dt_hydrokick*All.cf_atime); /* make sure v is in code units */
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
	        if((P[i].Type == 5) && (P[i].SuperTimestepFlag>=2)) {CellP[i].VelPred += fewbody_kick_dv;}
#endif	    
            
#if defined(TURB_DRIVING)
            CellP[i].VelPred += CellP[i].TurbAccel * dt_gravkick;
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            CellP[i].VelPred += CellP[i].Rad_Accel * (All.cf_atime * dt_hydrokick);
#endif
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            P[i].Mass = DMAX(P[i].Mass + CellP[i].DtMass * dt_entr, 0.5 * CellP[i].MassTrue); CellP[i].Mass = P[i].Mass;
#endif
            
            CellP[i].Density *= exp(-divv_fac);
            double etmp = CellP[i].InternalEnergyPred + CellP[i].DtInternalEnergy * dt_entr;
#if defined(RADTRANSFER) && defined(RT_EVOLVE_ENERGY) /* block here to deal with tricky cases where radiation energy density is -much- larger than thermal */ 
            int kfreq; double erad_tot=0,tot_e_min=0,enew=0,int_e_min=0,dErad=0,rsol_fac=C_LIGHT_CODE_REDUCED/C_LIGHT_CODE; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad_tot+=CellP[i].Rad_E_gamma_Pred[kfreq];}
            if(erad_tot > 0)
            {
                int_e_min=0.025*CellP[i].InternalEnergyPred; tot_e_min=0.025*(erad_tot/rsol_fac+CellP[i].InternalEnergyPred*P[i].Mass);
                enew=DMAX(erad_tot/rsol_fac+etmp*P[i].Mass,tot_e_min); etmp=(enew-erad_tot/rsol_fac)/P[i].Mass; if(etmp<int_e_min) {dErad=rsol_fac*(etmp-int_e_min); etmp=int_e_min;}
                if(dErad<-0.975*erad_tot) {dErad=-0.975*erad_tot;} CellP[i].InternalEnergyPred = etmp; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {CellP[i].Rad_E_gamma_Pred[kfreq] *= 1 + dErad/erad_tot;}
            } else {
                if(etmp<0.5*CellP[i].InternalEnergyPred) {CellP[i].InternalEnergyPred *= 0.5;} else {CellP[i].InternalEnergyPred=etmp;}
            }
#else
            if(etmp<0.5*CellP[i].InternalEnergyPred) {CellP[i].InternalEnergyPred *= 0.5;} else {CellP[i].InternalEnergyPred=etmp;}
#endif
            if(CellP[i].InternalEnergyPred<All.MinEgySpec) CellP[i].InternalEnergyPred=All.MinEgySpec;
            
#ifdef HYDRO_PRESSURE_SPH
            CellP[i].EgyWtDensity *= exp(-divv_fac);
#endif
            
#if (HYDRO_FIX_MESH_MOTION > 0)
            P[i].KernelRadius *= exp((double)divv_fac / ((double)NUMDIMS));
            if(P[i].KernelRadius < All.MinKernelRadius) {P[i].KernelRadius = All.MinKernelRadius;}
            if(P[i].KernelRadius > All.MaxKernelRadius) {P[i].KernelRadius = All.MaxKernelRadius;}
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            if(1 & ADAPTIVE_GRAVSOFT_FORALL) {P[i].AGS_KernelRadius = P[i].KernelRadius;} /* gas is AGS-active, so needs to be set here to match updated KernelRadius */
#endif
#endif
            drift_extra_physics(i, time0, time1, dt_entr);

            set_eos_pressure(i, P, CellP);
        }
    
    /* check for reflecting or outflow or otherwise special boundaries: if so, do the reflection/boundary! */
    apply_special_boundary_conditions(i,P[i].Mass,0);

    P[i].Ti_current = time1;
}





void move_particles(integertime time1)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(i=0; i<NumPart; i++) {drift_particle(i, time1);}
}





void drift_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
#ifdef MAGNETIC
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    CellP[i].BPred += CellP[i].DtB * (dt_entr * BphysVolphys_to_BcodeVolCode); // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1 / All.cf_a3inv; // for mass-based phi fluxes (otherwise coefficient is 1)
    double dtphi_code = (PhiphysVolphys_to_PhicodeVolCode) * CellP[i].DtPhi;
    CellP[i].PhiPred += dtphi_code  * dt_entr;
    double t_damp = Get_Gas_PhiField_DampingTimeInv(i);
    if((t_damp>0) && (!isnan(t_damp)))
    {
        CellP[i].PhiPred *= exp( -dt_entr * t_damp );
    }
#endif
#ifdef MHD_ALTERNATIVE_LEAPFROG_SCHEME
    CellP[i].B = CellP[i].BPred;
#ifdef DIVBCLEANING_DEDNER
    CellP[i].Phi=CellP[i].PhiPred;
#endif
#endif
#endif
#ifdef COSMIC_RAY_FLUID
    CosmicRay_Update_DriftKick(i, dt_entr, 1, P, CellP);
#endif
#ifdef RADTRANSFER
    rt_update_driftkick(i, dt_entr, 1, P, CellP);
#endif
#ifdef EOS_ELASTIC
    elastic_body_update_driftkick(i,dt_entr,1);
#endif
}





/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef BOX_PERIODIC
void do_box_wrapping(void)
{
    int i, j;
    double boxsize[3];
    boxsize[0] = boxSize_X;
    boxsize[1] = boxSize_Y;
    boxsize[2] = boxSize_Z;
    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(j)
#endif
    for(i = 0; i < NumPart; i++)
    {
        for(j = 0; j < 3; j++)
        {
            while(P[i].Pos[j] < 0)
            {
                P[i].Pos[j] += boxsize[j];
#ifdef BOX_SHEARING
                if(j==0)
                {
                    P[i].Vel[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
                    P[i].dp[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset * P[i].Mass;
                    if(P[i].Type==0)
                    {
                        CellP[i].VelPred[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) // if have moving cells need to wrap them, too (if cells aren't moving, should never reach this wrap) //
                        CellP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#endif
                    }
#if (BOX_SHEARING > 1)
                    /* if we're not assuming axisymmetry, we need to shift the coordinates for the shear flow at the boundary */
                    P[i].Pos[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Pos_Offset;
#endif
                }
#endif
            }
            
            while(P[i].Pos[j] >= boxsize[j])
            {
                P[i].Pos[j] -= boxsize[j];
#ifdef BOX_SHEARING
                if(j==0)
                {
                    P[i].Vel[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
                    P[i].dp[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset * P[i].Mass;
                    if(P[i].Type==0)
                    {
                        CellP[i].VelPred[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) // if have moving cells need to wrap them, too (if cells aren't moving, should never reach this wrap) //
                        CellP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#endif
                    }
#if (BOX_SHEARING > 1)
                    /* if we're not assuming axisymmetry, we need to shift the coordinates for the shear flow at the boundary */
                    P[i].Pos[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Pos_Offset;
#endif
                }
#endif
            }
        }
    }
}
#endif




/* ====================================================================== */
/* ================== Functions for physical information ================ */
/* ====================================================================== */




double INLINE_FUNC Get_Particle_Expected_Area(double h)
{
#if (NUMDIMS == 1)
    return 2;
#endif
#if (NUMDIMS == 2)
    return 2 * M_PI * h;
#endif
#if (NUMDIMS == 3)
    return 4 * M_PI * h * h;
#endif
}


/* return the estimated local column (physical units) from a local Sobolev approximation, or using the 'treecol' approximation from the gravity tree if the relevant config flag options are enabled */
double evaluate_NH_from_GradRho(MyFloat gradrho[3], double rkern, double rho, double numngb_ndim, double include_h, int target, struct particle_data *pp)
{
    double gradrho_mag=0;
    if(rho>0)
    {
#ifdef RT_USE_TREECOL_FOR_NH
        gradrho_mag = include_h * rho * rkern / numngb_ndim; if(target>=0) {gradrho_mag += pp[target].SigmaEff;}
#else             
        gradrho_mag = sqrt(gradrho[0]*gradrho[0]+gradrho[1]*gradrho[1]+gradrho[2]*gradrho[2]);
        if(gradrho_mag > 0) {gradrho_mag = rho*rho/gradrho_mag;} else {gradrho_mag=0;}
        if(include_h > 0) if(numngb_ndim > 0) gradrho_mag += include_h * rho * rkern / numngb_ndim; // quick-and-dirty approximation to the effective neighbor number needed here
#endif        
    }
    return gradrho_mag * All.cf_a2inv; // (physical units) // *(Z/Zsolar) add metallicity dependence
}







#ifdef MAGNETIC
/* this function is needed to control volume fluxes of the normal components of B and phi in the 
    -bad- situation where the meshless method 'faces' do not properly close (usually means you are 
    using boundary conditions that you should not) */
double Get_DtB_FaceArea_Limiter(int i, struct particle_data *pp, struct gas_cell_data *cell)
{
#ifdef HYDRO_SPH
    return 1;
#else
    /* define some variables */
    double dt_entr = get_particle_timestep_in_physical(i);
    /* check the magnitude of the predicted change in B-fields, vs. B-magnitude */
    Vec3<double> dB = cell[i].DtB * (dt_entr / All.cf_atime); /* converts to code units of Vol_code*B_code = Vol_phys*B_phys/a */
    double dBmag = dB.norm(), Bmag = cell[i].BPred.norm();
    /* also make sure to check the actual pressure, since if P>>B, we will need to allow larger changes in B per timestep */
    double P_BV_units = sqrt(2.*cell[i].Pressure*All.cf_a3inv)*pp[i].Mass/cell[i].Density / All.cf_a2inv;
    /* the above should be in CODE Bcode*Vol_code units! */
    double Bmag_max = DMAX(Bmag, DMIN( P_BV_units, 10.*Bmag ));
    /* now check how accurately the cell is 'closed': the face areas are ideally zero */
    double area_sum = fabs(cell[i].Face_Area[0])+fabs(cell[i].Face_Area[1])+fabs(cell[i].Face_Area[2]);
    /* but this needs to be normalized to the 'expected' area given KernelRadius */
    double area_norm = Get_Particle_Expected_Area(pp[i].KernelRadius * All.cf_atime);
    /* ok, with that in hand, define an error tolerance based on this */
    if(area_norm>0)
    {
        double area_norm_min_threshold = 0.001;
        double area_norm_weight = 200.0;
        if(area_sum/area_norm > area_norm_min_threshold)
        {
            double tol = (All.CourantFac/0.2) * DMAX( 0.01, area_norm/(area_norm_weight * area_sum) );
            tol *= Bmag_max; /* give the limiter dimensions */
            if(dBmag > tol) {return tol/dBmag;} /* now actually check if we exceed this */
        }
    }
    return 1;
#endif
}


#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Gas_PhiField(int i_particle_id)
{
    //return CellP[i_particle_id].PhiPred * CellP[i_particle_id].Density / P[i_particle_id].Mass; // volumetric phy-flux (requires extra term compared to mass-based flux)
    return CellP[i_particle_id].PhiPred / P[i_particle_id].Mass; // mass-based phi-flux
}

double INLINE_FUNC Get_Gas_PhiField_DampingTimeInv(int i_particle_id)
{
    /* this timescale should always be returned as a -physical- time */
#ifdef HYDRO_SPH
    /* PFH: add simple damping (-phi/tau) term */
    double damping_tinv = 0.5 * All.DivBcleanParabolicSigma * (CellP[i_particle_id].MaxSignalVel / (All.cf_atime*P[i_particle_id].Get_Particle_Size()));
#else
    double damping_tinv;
#ifdef SELFGRAVITY_OFF
    damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveSpeed / (All.cf_atime*P[i_particle_id].Get_Particle_Size()); // fastest wavespeed has units of [vphys]
    //double damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveDecay * All.cf_a2inv; // no improvement over fastestwavespeed; decay has units [vphys/rphys]
#else
    // only see a small performance drop from fastestwavespeed above to maxsignalvel below, despite the fact that below is purely local (so allows more flexible adapting to high dynamic range)
    damping_tinv = 0.0;
    
    if(P[i_particle_id].KernelRadius > 0)
    {
        double h_eff = P[i_particle_id].Get_Particle_Size();
        double vsig2 = 0.5 * fabs(CellP[i_particle_id].MaxSignalVel);
        double phi_B_eff = 0.0;
        if(vsig2 > 0) {phi_B_eff = Get_Gas_PhiField(i_particle_id) / (All.cf_atime * vsig2);}
        double vsig1 = 0.0;
        if(CellP[i_particle_id].Density > 0)
        {
            vsig1 = sqrt( CellP[i_particle_id].effective_soundspeed()*CellP[i_particle_id].effective_soundspeed() +
                 (1. / All.cf_atime) *
                 (CellP[i_particle_id].Bfield().norm_sq() +
                  phi_B_eff*phi_B_eff) / CellP[i_particle_id].Density );
        }
        vsig1 = DMAX(vsig1, vsig2);
        vsig2 = 0.0;
        vsig2 = CellP[i_particle_id].Gradients.Velocity.frobenius_norm();
        vsig2 = 3.0 * h_eff * DMAX( vsig2, fabs(P[i_particle_id].Particle_DivVel)) / All.cf_atime;
        double prefac_fastest = 0.1;
        double prefac_tinv = 0.5;
        double area_0 = 0.1;
#ifdef MHD_CONSTRAINED_GRADIENT
        prefac_fastest = 1.0;
        prefac_tinv = 2.0;
        area_0 = 0.05;
        vsig2 *= 5.0;
        if(CellP[i_particle_id].FlagForConstrainedGradients <= 0) prefac_tinv *= 30;
#endif
        prefac_tinv *= sqrt(1. + CellP[i_particle_id].ConditionNumber/100.);
        double area = fabs(CellP[i_particle_id].Face_Area[0]) + fabs(CellP[i_particle_id].Face_Area[1]) + fabs(CellP[i_particle_id].Face_Area[2]);
        area /= Get_Particle_Expected_Area(P[i_particle_id].KernelRadius);
        prefac_tinv *= (1. + area/area_0)*(1. + area/area_0);
        
        double vsig_max = DMAX( DMAX(vsig1,vsig2) , prefac_fastest * All.FastestWaveSpeed );
        damping_tinv = prefac_tinv * All.DivBcleanParabolicSigma * (vsig_max / (All.cf_atime * h_eff));
    }
#endif
#endif
    return damping_tinv;
}

#endif // dedner
#endif // magnetic





/* -------------------------------------------------------------------------------------------------------------------------------------
 ------------------- the following routines are not setting the velocity, but instead are useful routines for computation of
 -------------------  various quantities needed in the mesh motion for different coordinate systems or assumed mesh shapes
 ------------------------------------------------------------------------------------------------------------------------------------- */

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
/* time-step the positions of the mesh points. this is trivial except if we are evolving the mesh points in non-cartesian coordinates
    (cylindrical or spherical) based on assumed fixed initial velocities (if HYDRO_FIX_MESH_MOTION=2 or 3),
    in which case we have to convert back and forth. */
void advect_mesh_point(int i, double dt)
{
#if (HYDRO_FIX_MESH_MOTION == 2) || (HYDRO_FIX_MESH_MOTION == 3) // cylindrical or spherical coordinates
    // define the location relative to the origin (needed in these coordinate systems)
    Vec3<double> dp = P[i].Pos; Vec3<double> dp_offset = {}; // assume center is at coordinate origin
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE) // unless we use a special anchor, to define the center
    dp_offset = P[i].Pos - P[i].Min_xyz_to_Sink;
#elif defined(BOX_PERIODIC) // or if periodic, the box mid-point is instead the center
#if (NUMDIMS==1)
    dp_offset[0] = -boxHalf_X;
#elif (NUMDIMS==2)
    dp_offset[0] = -boxHalf_X; dp_offset[1] = -boxHalf_Y;
#else
    dp_offset = Vec3<double>{-boxHalf_X, -boxHalf_Y, -boxHalf_Z};
#endif
#endif
    dp += dp_offset;
#if (HYDRO_FIX_MESH_MOTION == 2) // cylindrical
    double r2=dp[0]*dp[0]+dp[1]*dp[1], r=sqrt(r2), c0=dp[0]/r, s0=dp[1]/r, z=dp[2]; // get r, sin/cos theta, z
    double vr=c0*CellP[i].ParticleVel[0] + s0*CellP[i].ParticleVel[1], vt=s0*CellP[i].ParticleVel[0] - c0*CellP[i].ParticleVel[1], vz=CellP[i].ParticleVel[2]; // velocities in these directions
    double r_n=r+vr*dt, z_n=z+vz*dt, c_n=c0-s0*(vt/r)*dt, s_n=s0+c0*(vt/r)*dt; // updated cylindrical values
    dp[0] = c_n*r_n; dp[1] = s_n*r_n; dp[2] = z_n; // back to coordinates
    CellP[i].ParticleVel[0] = c_n*vr + s_n*vt; // re-set velocities in these coordinates //
    CellP[i].ParticleVel[1] = s_n*vr - c_n*vt;
    CellP[i].ParticleVel[2] = vz;
    return;
#elif (HYDRO_FIX_MESH_MOTION == 3) // spherical
    Vec3<double> v = CellP[i].ParticleVel; double r2=dp.norm_sq(); // assume center is at coordinate origin
    double r=sqrt(r2), rxy=sqrt(dp[0]*dp[0]+dp[1]*dp[1]), vr=dot(dp,v)/r; // updated r is easy
    double ct = 1./sqrt(1.+dp[1]*dp[1]/(dp[0]*dp[0])), st = (dp[1]/dp[0])*ct; // cos and sin theta
    double cp = sqrt(1.-dp[2]*dp[2]/(r*r)), sp = dp[2]/r; // cos and sin phi
    double t_dot = (v[0]*dp[1]-v[1]*dp[0])/(rxy*rxy), p_dot = (dp[2]*(dp[0]*v[0]+dp[1]*v[1])-rxy*rxy*v[2])/(r*r*rxy); // theta, phi derivatives
    double r_n=r+vr*dt, ct_n=ct-st*t_dot, st_n=st+ct*t_dot, cp_n=cp-sp*t_dot, sp_n=sp+cp*t_dot; // updated angles and positions in spherical
    dp[0] = r_n * ct_n * cp_n; dp[1] = r_n * st_n * cp_n; dp[2] = r_n * sp_n; // back to coordinates
    rxy = sqrt(dp[0]*dp[0] + dp[1]*dp[1]); // updated rxy
    CellP[i].ParticleVel[0] = (dp[0]/r_n) * vr + dp[1] * t_dot + dp[0]*dp[2]/rxy * p_dot; // back to cartesian velocities
    CellP[i].ParticleVel[1] = (dp[1]/r_n) * vr - dp[0] * t_dot + dp[1]*dp[2]/rxy * p_dot; // back to cartesian velocities
    CellP[i].ParticleVel[2] = (dp[2]/r_n) * vr - rxy * p_dot; // back to cartesian velocities
    return;
#endif
    // ok now have the updated x/y/z positions relative to the origin, convert these back to the simulation coordinate frame
    P[i].Pos = dp - dp_offset;
#endif // ok done with cylindrical/spherical coordinates
    
    
    // ok anything else ('normal' coordinates), does down here
    P[i].Pos += CellP[i].ParticleVel * dt; // for standard grid velocities, this is trivial //
    return;
}




/* routine to calculate the overlapping face area of two cuboids in NDIMS dimensions based on their relative positions */
double calculate_face_area_for_cartesian_mesh(const Vec3<double>& dp, double rinv, double l_side, Vec3<double>& Face_Area_Vec)
{
    Face_Area_Vec = {}; double Face_Area_Norm;
#if (NUMDIMS==1)
    Face_Area_Norm = 1; Face_Area_Vec[0] = Face_Area_Norm * dp[0]/fabs(dp[0]);
#elif (NUMDIMS==2)
    if(fabs(dp[0]) > fabs(dp[1])) {Face_Area_Vec[0] = Face_Area_Norm = std::max(0.,l_side-fabs(dp[1])) * dp[0]/fabs(dp[0]) * All.cf_atime;} else {Face_Area_Vec[1] = Face_Area_Norm = std::max(0.,l_side-fabs(dp[0])) * dp[1]/fabs(dp[1]) * All.cf_atime;}
#else
    Vec3<double> dp_abs = {fabs(dp[0]), fabs(dp[1]), fabs(dp[2])};
    int kdir;
    if((dp_abs[0]>=dp_abs[1])&&(dp_abs[0]>=dp_abs[2])) {kdir=0;} else if ((dp_abs[1]>=dp_abs[0])&&(dp_abs[1]>=dp_abs[2])) {kdir=1;} else {kdir=2;}
    Face_Area_Norm=1; for(int k=0;k<3;k++) {if(k!=kdir) {Face_Area_Norm *= std::max(0.,l_side-dp_abs[k]) * All.cf_atime*All.cf_atime;}}
    Face_Area_Vec[kdir] = Face_Area_Norm * dp[kdir]/fabs(dp[kdir]);
#endif
    return fabs(Face_Area_Norm);
}

#endif




#ifdef SPECIAL_POINT_WEIGHTED_MOTION
double weight_function_for_weighted_motion_smoothing(double r, int mode)
{
    double wt = 1, amax = 1.e2, rmax_pc = 100., slope_index = 1, r_pc = r * All.cf_atime * UNIT_LENGTH_IN_PC;
    if(r_pc < rmax_pc) {wt = DMAX(pow(r_pc / rmax_pc , slope_index) , 1./amax);}
    if(mode) {return 1 - sqrt(wt);} else {return wt;}
}
#endif
