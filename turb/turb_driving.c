#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"

/* This file contains the routines for driven turbulence/stirring; use for things
 like idealized turbulence tests, large-eddy simulations, and the like */
/*
 *  This code was originally written for GADGET3 by Andreas Bauer; it has been
 *   modified slightly by Phil Hopkins for GIZMO, but is largely intact.
 */

#if defined(TURB_DRIVING)

/* block of global variables used specifically for the set of subroutines below, which need to be carried between timesteps */
double* StOUPhases; // random fluctuating component of the amplitudes
double* StAmpl; // relative amplitude for each k
double* StAka; // phases (real part)
double* StAkb; // phases (imag part)
double* StMode; // k vectors
int StNModes; // total number of modes
integertime StTPrev; // time of last update (to determine when next will be)
gsl_rng* StRng; // random number generator key


/* routine to initialize the different modes and their relative amplitudes and other global variables needed for the turbulence driving routines */
void init_turb(void)
{
    int ikx, iky, ikz; double kx,ky,kz,k, ampl;
    
    double kmin = All.TurbDriving_Global_DrivingScaleKMinVar, kmax = All.TurbDriving_Global_DrivingScaleKMaxVar;
#if !defined(TURB_DRIVING_OLDFORMAT)
    kmin = 2.*M_PI / All.TurbDriving_Global_DrivingScaleKMinVar; kmax = 2.*M_PI / All.TurbDriving_Global_DrivingScaleKMaxVar; // convert these from spatial lengths to wavenumbers in the new convention we are using
#endif
    
    int ikxmax = boxSize_X * kmax/2./M_PI, ikymax = 0, ikzmax = 0;
#if (NUMDIMS > 1)
    ikymax = boxSize_Y * kmax/2./M_PI;
#endif
#if (NUMDIMS > 2)
    ikzmax = boxSize_Z * kmax/2./M_PI;
#endif
    
    StNModes = 0;
    for(ikx = 0;ikx <= ikxmax; ikx++)
    {
        kx = 2.*M_PI*ikx/boxSize_X;
        for(iky = 0;iky <= ikymax; iky++)
        {
            ky = 2.*M_PI*iky/boxSize_Y;
            for(ikz = 0;ikz <= ikzmax; ikz++)
            {
                kz = 2.*M_PI*ikz/boxSize_Z;
                k = sqrt(kx*kx+ky*ky+kz*kz);
                if(k>=kmin && k<=kmax)
                {
#if NUMDIMS ==1
                    StNModes+=1;
#endif
#if NUMDIMS == 2
                    StNModes+=2;
#endif
#if NUMDIMS == 3
                    StNModes+=4;
#endif
                }
            }
        }
    }
    
    PRINT_STATUS("Initializing turbulent driving: max integer mode number ikx/iky/ikz = %d %d %d",ikxmax,ikymax,ikzmax);
    StMode = (double*) mymalloc_movable(&StMode,"StModes", StNModes * 3 * sizeof(double));
    StAka = (double*) mymalloc_movable(&StAka,"StAka", StNModes * 3 * sizeof(double));
    StAkb = (double*) mymalloc_movable(&StAkb,"StAkb", StNModes * 3 * sizeof(double));
    StAmpl = (double*) mymalloc_movable(&StAmpl,"StAmpl", StNModes * sizeof(double));
    StOUPhases = (double*) mymalloc_movable(StOUPhases,"StOUPhases", StNModes * 6 * sizeof(double));
    double kc = 0.5*(kmin+kmax), amin = 0., amplitude_integrated_allmodes = 0.;
    StNModes = 0;
    
    for(ikx = 0;ikx <= ikxmax; ikx++)
    {
        kx = 2.*M_PI*ikx/boxSize_X;
        
        for(iky = 0;iky <= ikymax; iky++)
        {
            ky = 2.*M_PI*iky/boxSize_Y;
            
            for(ikz = 0;ikz <= ikzmax; ikz++)
            {
                kz = 2.*M_PI*ikz/boxSize_Z;
                
                k = sqrt(kx*kx+ky*ky+kz*kz);
                if(k>=kmin && k<=kmax)
                {
                    if(All.TurbDriving_Global_DrivingSpectrumKey == 0)
                    {
                        ampl = 1.; // uniform amplitude for all
                    }
                    else if(All.TurbDriving_Global_DrivingSpectrumKey ==  1)
                    {
                        ampl = 4.0*(amin-1.0)/((kmax-kmin)*(kmax-kmin))*((k-kc)*(k-kc))+1.0; // spike at kc = k_driving
                    }
                    else if(All.TurbDriving_Global_DrivingSpectrumKey == 2)
                    {
                        //ampl = pow(k/kmin, (1.-NUMDIMS)- 5./3. ); // because this is E[vector_k] for NUMDIMS, need extra NUMDIMS-1 power term here
                        ampl = pow(k/kmin, 1./3. - 0*0.5*NUMDIMS); // this should scale as the acceleration per eddy, ~v^2/L, crudely
                    }
                    else if(All.TurbDriving_Global_DrivingSpectrumKey == 3)
                    {
                        //ampl = pow(k/kmin, (1.-NUMDIMS)- 2. ); // because this is E[vector_k] for NUMDIMS, need extra NUMDIMS-1 power term here
                        ampl = pow(k/kmin, 0. - 0*0.5*NUMDIMS); // this should scale as the acceleration per eddy, ~v^2/L, crudely
                    }
                    else
                    {
                        terminate("unknown spectral form");
                    }
                    
                    StAmpl[StNModes] = ampl;
                    StMode[3*StNModes+0] = kx;
                    StMode[3*StNModes+1] = ky;
                    StMode[3*StNModes+2] = kz;
                    PRINT_STATUS("  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,iky,ikz,StMode[3*StNModes+0],StMode[3*StNModes+1],StMode[3*StNModes+2],StAmpl[StNModes]);
                    amplitude_integrated_allmodes += ampl*ampl;
                    StNModes++;
                    
#if (NUMDIMS > 1)
                    if(ikx>0 || iky>0) // if both of these are zero, only non-degenerate modes are the +/- z modes [ensured below]
                    {
                        StAmpl[StNModes] = ampl;
                        if(iky!=0) {StMode[3*StNModes+0] = kx;} else {StMode[3*StNModes+0] = -kx;}
                        StMode[3*StNModes+1] = -ky;
                        StMode[3*StNModes+2] = kz;
                        PRINT_STATUS("  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,-iky,ikz,StMode[3*StNModes+0],StMode[3*StNModes+1],StMode[3*StNModes+2],StAmpl[StNModes]);
                        amplitude_integrated_allmodes += ampl*ampl;
                        StNModes++;
                    }

#if (NUMDIMS > 2)
                    if((iky>0 || ikz>0) && (ikx>0 || ikz>0)) // if both of these are zero, only non-degenerate modes are the +/- x or +/- y modes [already ensured above]
                    {
                        StAmpl[StNModes] = ampl;
                        if(ikz!=0) {StMode[3*StNModes+0] = kx;} else {StMode[3*StNModes+0] = -kx;}
                        StMode[3*StNModes+1] = ky;
                        StMode[3*StNModes+2] = -kz;
                        PRINT_STATUS("  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,iky,-ikz,StMode[3*StNModes+0],StMode[3*StNModes+1],StMode[3*StNModes+2],StAmpl[StNModes]);
                        amplitude_integrated_allmodes += ampl*ampl;
                        StNModes++;
                    }

                    if((ikx>0 || iky>0) && (ikx>0 || ikz>0) && (iky>0 || ikz>0)) // if both of these are zero, only non-degenerate modes are +/- z or +/- y or +/- x modes already handled above
                    {
                        StAmpl[StNModes] = ampl;
                        if(ikz==0 || iky==0) {StMode[3*StNModes+0] = -kx;} else {StMode[3*StNModes+0] = kx;}
                        StMode[3*StNModes+1] = -ky;
                        StMode[3*StNModes+2] = -kz;
                        PRINT_STATUS("  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,-iky,-ikz,StMode[3*StNModes+0],StMode[3*StNModes+1],StMode[3*StNModes+2],StAmpl[StNModes]);
                        amplitude_integrated_allmodes += ampl*ampl;
                        StNModes++;
                    }
#endif
#endif
                }
            }
        }
    }
#if !defined(TURB_DRIVING_OLDFORMAT)
    int i; for(i=0; i<StNModes; i++) {StAmpl[i] *= sqrt(1./amplitude_integrated_allmodes);} // normalize total driving amplitude across all modes here
#endif
    StTPrev = -1; // mark some arbitrarily old time as last update of turb driving fields
    StRng = gsl_rng_alloc(gsl_rng_ranlxd1); // allocate seed variables
    gsl_rng_set(StRng, All.TurbDriving_Global_DrivingRandomNumberKey); // initialize seed
    int j; for(j=0;j<100;j++) {double tmp; tmp=st_turbdrive_get_gaussian_random_variable();} // cycle past initial seed
    st_turbdrive_init_ouseq(); // initialize variable for phases
    st_turbdrive_calc_phases(); // initialize phases
    set_turb_ampl(); // set initial amplitudes and calculate initial quantities needed for dissipation measures
    StTPrev = All.Ti_Current; // mark current time as last update of turb driving fields
}


/* initialize phase variables */
void st_turbdrive_init_ouseq(void)
{
    int i; for(i = 0;i<6*StNModes;i++) {StOUPhases[i] = st_turbdrive_get_gaussian_random_variable()*st_return_rms_acceleration();}
}


/* return the rms acceleration we expect, using either the 'dissipation rate' or 'turbulent velocity' conventions for our variables */
double st_return_rms_acceleration(void)
{
#if !defined(TURB_DRIVING_OLDFORMAT)
    return All.TurbDriving_Global_AccelerationPowerVariable / st_return_mode_correlation_time(); // new convention, hoping this is more clear re: meaning of variable
#else
    return sqrt(All.TurbDriving_Global_AccelerationPowerVariable / st_return_mode_correlation_time()); // old convention, very confusing about what this variable meant: the numerator is the "D = sigma^2 / 2" in the classic OU form, but confusing b/c this is acceleration not v for Brownian motion, so this is not an energy in any sense
#endif
}


/* return the driving scale needed for scaling some other quantities below, corresponding to our global variable convention */
double st_return_driving_scale(void)
{
#if !defined(TURB_DRIVING_OLDFORMAT)
    return All.TurbDriving_Global_DrivingScaleKMinVar; // this is now spatial scale
#else
    return 2.*M_PI / All.TurbDriving_Global_DrivingScaleKMinVar; // this was K before
#endif
}


/* return the coherence time of the driving scale modes. if global variable negative, it uses the eddy turnover time of the driving-scale modes */
double st_return_mode_correlation_time(void)
{
#if !defined(TURB_DRIVING_OLDFORMAT)
    if(All.TurbDriving_Global_DecayTime > 0) {return All.TurbDriving_Global_DecayTime;} else {return st_return_driving_scale() / All.TurbDriving_Global_AccelerationPowerVariable;}
#else
    if(All.TurbDriving_Global_DecayTime > 0) {return All.TurbDriving_Global_DecayTime;} else {return pow(st_return_driving_scale(), 2./3.) / pow(All.TurbDriving_Global_AccelerationPowerVariable, 1./3.);}
#endif
}


/* return time interval between turbulent driving field updates based on global variable. if negative, default to small interval of coherence time by default. */
double st_return_dt_between_updates(void)
{
    if(All.TurbDriving_Global_DtTurbUpdates > 0) {return All.TurbDriving_Global_DtTurbUpdates;} else {return 0.01*st_return_mode_correlation_time();}
}


/* return factor needed to renormalize below based on fraction of power projected out in our solenoidal projection, to return the correct normalization for accelerations. */
double solenoidal_frac_total_weight_renormalization(void)
{
#if (NUMDIMS >= 3)
    return sqrt(3.0/3.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*All.TurbDriving_Global_SolenoidalFraction+3.0*All.TurbDriving_Global_SolenoidalFraction*All.TurbDriving_Global_SolenoidalFraction);
#endif
#if (NUMDIMS == 2)
    return sqrt(3.0/2.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*All.TurbDriving_Global_SolenoidalFraction+2.0*All.TurbDriving_Global_SolenoidalFraction*All.TurbDriving_Global_SolenoidalFraction);
#endif
#if (NUMDIMS == 1)
    return sqrt(3.0/1.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*All.TurbDriving_Global_SolenoidalFraction+1.0*All.TurbDriving_Global_SolenoidalFraction*All.TurbDriving_Global_SolenoidalFraction);
#endif
}


/* update the Markov random variable that is the dimensional multiplier for the acceleration field, which has a correlation time specified */
void st_update_ouseq(void)
{
    int i; double damping = exp( -st_return_dt_between_updates()/st_return_mode_correlation_time());
    for(i = 0;i<6*StNModes;i++) {StOUPhases[i] = StOUPhases[i] * damping + st_return_rms_acceleration() * sqrt(1.-damping*damping)*st_turbdrive_get_gaussian_random_variable();}
}


/* routine to return gaussian random number with zero mean and unity variance */
double st_turbdrive_get_gaussian_random_variable(void)
{
    double r0 = gsl_rng_uniform(StRng), r1 = gsl_rng_uniform(StRng);
    return sqrt(2. * log(1. / r0) ) * cos(2. * M_PI * r1);
}


/* routine to calculate the projected phases/acceleration field variables, using the fourier-space solenoidal/compressible projection */
void st_turbdrive_calc_phases(void)
{
    int i,j;
    for(i = 0; i < StNModes;i++)
    {
        double ka = 0., kb = 0., kk = 0.; int dim = NUMDIMS;
        for(j = 0; j<dim;j++)
        {
            kk += StMode[3*i+j]*StMode[3*i+j];
            ka += StMode[3*i+j]*StOUPhases[6*i+2*j+1];
            kb += StMode[3*i+j]*StOUPhases[6*i+2*j+0];
        }
        for(j = 0; j<dim;j++)
        {
            double diva = StMode[3*i+j]*ka/kk;
            double divb = StMode[3*i+j]*kb/kk;
            double curla = StOUPhases[6*i+2*j+0] - divb;
            double curlb = StOUPhases[6*i+2*j+1] - diva;
            
            StAka[3*i+j] = All.TurbDriving_Global_SolenoidalFraction*curla+(1.-All.TurbDriving_Global_SolenoidalFraction)*divb;
            StAkb[3*i+j] = All.TurbDriving_Global_SolenoidalFraction*curlb+(1.-All.TurbDriving_Global_SolenoidalFraction)*diva;
        }
    }
}


/* parent routine to initialize and update turbulent driving fields and to track different variables used for analyzing power spectra of dissipation, etc. */
void set_turb_ampl(void)
{
    double delta = (All.Ti_Current - StTPrev) * UNIT_INTEGERTIME_IN_PHYSICAL, Dt_Update=st_return_dt_between_updates();
    if(delta >= Dt_Update)
    {
        if(delta > 0)
        {
            int i; double e_diss_sum=0, e_drive_sum=0, glob_diss_sum=0, glob_drive_sum=0;
            PRINT_STATUS(" ..updating fields tracked for following injected energy and dissipation");
            for(i=0; i < NumPart; i++)
            {
                if(P[i].Type == 0)
                {
                    if(P[i].Mass > 0)
                    {
                        e_diss_sum += SphP[i].EgyDiss;
                        SphP[i].DuDt_diss = (SphP[i].EgyDiss / P[i].Mass) / delta;
                        SphP[i].EgyDiss = 0;
                        e_drive_sum += SphP[i].EgyDrive;
                        SphP[i].DuDt_drive = (SphP[i].EgyDrive / P[i].Mass) / delta;
                        SphP[i].EgyDrive = 0;
                    } else {
                        SphP[i].DuDt_diss = SphP[i].EgyDiss = SphP[i].DuDt_drive = SphP[i].EgyDrive = 0;
                    }
                }
            }
            MPI_Allreduce(&e_diss_sum, &glob_diss_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&e_drive_sum, &glob_drive_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            All.TurbDissipatedEnergy += glob_diss_sum;
            All.TurbInjectedEnergy += glob_drive_sum;
        }
        PRINT_STATUS(" ..updating fourier-space phase information");
        st_update_ouseq();
        PRINT_STATUS(" ..calculating coefficients and phases following desired projection");
        st_turbdrive_calc_phases();
        StTPrev = StTPrev + Dt_Update / All.Timebase_interval;
        PRINT_STATUS(" ..updated turbulent stirring field at time %f", StTPrev * All.Timebase_interval);
    }
}


/* routine to actually calculate the turbulent acceleration 'driving field' force on every resolution element */
void add_turb_accel()
{
    set_turb_ampl();
    int i, j, m; double acc[3], fac_sol = 2.*solenoidal_frac_total_weight_renormalization();
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0)
        {
            double fx = 0, fy = 0, fz = 0;
            for(m=0; m<StNModes; m++) // calc force
            {
                double kxx = StMode[3*m+0]*P[i].Pos[0], kyy = StMode[3*m+1]*P[i].Pos[1], kzz = StMode[3*m+2]*P[i].Pos[2];
                double kdotx = kxx+kyy+kzz, ampl = StAmpl[m], realt = cos(kdotx), imagt = sin(kdotx);
                
                fx += ampl*(StAka[3*m+0]*realt - StAkb[3*m+0]*imagt);
                fy += ampl*(StAka[3*m+1]*realt - StAkb[3*m+1]*imagt);
                fz += ampl*(StAka[3*m+2]*realt - StAkb[3*m+2]*imagt);
            }
            fx *= fac_sol; fy *= fac_sol; fz *= fac_sol;
            
            if(P[i].Mass > 0.)
            {
                acc[0] = fx; acc[1] = acc[2] = 0;
#if (NUMDIMS > 1)
                acc[1] = fy;
#endif
#if (NUMDIMS > 2)
                acc[2] = fz;
#endif
                for(j=0; j<3; j++) {SphP[i].TurbAccel[j] = acc[j];}
            } else {
                SphP[i].TurbAccel[0] = SphP[i].TurbAccel[1] = SphP[i].TurbAccel[2]=0;
            }
        }
    }
    PRINT_STATUS("Finished turbulence driving (acceleration) computation");
}


/* routine to integrate the turbulent driving forces, specifically the 'TurbAccel' variables that need to be drifted and kicked: note that we actually do drifting and kicking in the normal routines, this is just to integrate the dissipation rates etc used for our tracking */
void do_turb_driving_step_first_half(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    int i, j; integertime ti_step, tstart, tend; double dvel[3], dt_gravkick;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        ti_step = GET_PARTICLE_INTEGERTIME(i); tstart = P[i].Ti_begstep; tend = P[i].Ti_begstep + ti_step / 2;	/* beginning / midpoint of step */
        if(All.ComovingIntegrationOn) {dt_gravkick = get_gravkick_factor(tstart, tend);} else {dt_gravkick = (tend - tstart) * All.Timebase_interval;}
        if(P[i].Type == 0)
        {
            double vtmp[3], ekin0 = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
            for(j=0;j<3;j++) {dvel[j] = SphP[i].TurbAccel[j] * dt_gravkick; vtmp[j] = P[i].Vel[j] + dvel[j];}
            double ekin1 = 0.5 * P[i].Mass * (vtmp[0]*vtmp[0]+vtmp[1]*vtmp[1]+vtmp[2]*vtmp[2]);
            SphP[i].EgyDrive += ekin1 - ekin0;
        }
    }
    CPU_Step[CPU_DRIFT] += measure_time();
}


/* routine to integrate the turbulent driving forces, specifically the 'TurbAccel' variables that need to be drifted and kicked: note that we actually do drifting and kicking in the normal routines, this is just to integrate the dissipation rates etc used for our tracking */
void do_turb_driving_step_second_half(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    int i, j; integertime ti_step, tstart, tend; double dvel[3], dt_gravkick;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        ti_step = GET_PARTICLE_INTEGERTIME(i); tstart = P[i].Ti_begstep + ti_step / 2; tend = P[i].Ti_begstep + ti_step;	/* midpoint/end of step */
        if(All.ComovingIntegrationOn) {dt_gravkick = get_gravkick_factor(tstart, tend);} else {dt_gravkick = (tend - tstart) * All.Timebase_interval;}
        if(P[i].Type == 0)
        {
            double vtmp[3], ekin0 = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
            for(j=0;j<3;j++) {dvel[j] = SphP[i].TurbAccel[j] * dt_gravkick; vtmp[j] = P[i].Vel[j] + dvel[j];}
            double ekin1 = 0.5 * P[i].Mass * (vtmp[0]*vtmp[0]+vtmp[1]*vtmp[1]+vtmp[2]*vtmp[2]);
            SphP[i].EgyDrive += ekin1 - ekin0;
        }
    }
    CPU_Step[CPU_DRIFT] += measure_time();
}


/* routine to record and optionally write to output files various statistics of driven turbulence here (most relevant to idealized sims with a hard-coded adiabatic EOS */
void log_turb_temp(void)
{
#ifndef IO_REDUCED_MODE
    int i; double dudt_drive = 0, dudt_diss = 0, mass = 0, ekin = 0, ethermal = 0;
    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Type == 0)
        {
            dudt_drive += P[i].Mass * SphP[i].DuDt_drive;
            dudt_diss += P[i].Mass * SphP[i].DuDt_diss;
            ekin += 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
            ethermal += P[i].Mass * SphP[i].InternalEnergy;
            mass += P[i].Mass;
        }
    }
    double glob_mass, glob_dudt_drive, glob_dudt_diss, glob_ekin, glob_ethermal;
    MPI_Allreduce(&mass, &glob_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&dudt_drive, &glob_dudt_drive, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&dudt_diss, &glob_dudt_diss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ekin, &glob_ekin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ethermal, &glob_ethermal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double mach = sqrt(2.*glob_ekin / (GAMMA_DEFAULT*(GAMMA_DEFAULT-1)*glob_ethermal));
    
    if(ThisTask == 0)
    {
        fprintf(FdTurb, "%g %g %g %g %g %g %g\n", All.Time, mach, (glob_ekin + glob_ethermal) / glob_mass, glob_dudt_drive / glob_mass,
                glob_dudt_diss / glob_mass, All.TurbInjectedEnergy / glob_mass, All.TurbDissipatedEnergy / glob_mass); fflush(FdTurb);
    }
#endif
}


#endif
