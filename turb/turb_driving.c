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


double st_grn(void);
void st_init_ouseq(void);
void st_calc_phases(void);
void st_compute_injection(void);


//Ornstein-Uhlenbeck variables
double StOUVar;
double* StOUPhases;
gsl_rng* StRng;


//forcing field in fourier space
double* StAmpl;
double* StAka; //phases (real part)
double* StAkb; //phases (imag part)
double* StMode;
int StNModes;


integertime StTPrev;
double StSolWeightNorm;


void init_turb(void)
{
    int ikx, iky, ikz;
    double kx,ky,kz,k;
    double ampl;
    
    int ikxmax = boxSize_X * All.StKmax/2./M_PI;
    int ikymax = 0;
    int ikzmax = 0;
#if (NUMDIMS > 1)
    ikymax = boxSize_Y * All.StKmax/2./M_PI;
#endif
#if (NUMDIMS > 2)
    ikzmax = boxSize_Z * All.StKmax/2./M_PI;
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
                if(k>=All.StKmin && k<=All.StKmax)
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
    
    PRINT_STATUS(" ..using %d modes, %d %d %d\n",StNModes,ikxmax,ikymax,ikzmax);
    
    StMode = (double*) mymalloc_movable(&StMode,"StModes", StNModes * 3 * sizeof(double));
    StAka = (double*) mymalloc_movable(&StAka,"StAka", StNModes * 3 * sizeof(double));
    StAkb = (double*) mymalloc_movable(&StAkb,"StAkb", StNModes * 3 * sizeof(double));
    StAmpl = (double*) mymalloc_movable(&StAmpl,"StAmpl", StNModes * sizeof(double));
    StOUPhases = (double*) mymalloc_movable(StOUPhases,"StOUPhases", StNModes * 6 * sizeof(double));
    
    StOUVar = sqrt(All.StEnergy/All.StDecay);
    double kc = 0.5*(All.StKmin+All.StKmax);
    double amin = 0.;
    
#if NUMDIMS == 3
    StSolWeightNorm = sqrt(3.0/3.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*All.StSolWeight+3.0*All.StSolWeight*All.StSolWeight);
#endif
#if NUMDIMS == 2
    StSolWeightNorm = sqrt(3.0/2.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*All.StSolWeight+2.0*All.StSolWeight*All.StSolWeight);
#endif
#if NUMDIMS == 1
    StSolWeightNorm = sqrt(3.0/1.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*All.StSolWeight+1.0*All.StSolWeight*All.StSolWeight);
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
                if(k>=All.StKmin && k<=All.StKmax)
                {
                    if(All.StSpectForm == 0)
                    {
                        ampl = 1.;
                    }
                    else if(All.StSpectForm ==  1)
                    {
                        ampl = 4.0*(amin-1.0)/((All.StKmax-All.StKmin)*(All.StKmax-All.StKmin))*((k-kc)*(k-kc))+1.0;
                    }
                    else if(All.StSpectForm == 2)
                    {
                        ampl = pow(All.StKmin,5./3)/pow(k,5./3);
                    }
                    else if(All.StSpectForm == 3)
                    {
                        ampl = pow(All.StKmin,2.)/pow(k,2.);
                    }
                    else
                    {
                        terminate("unkown spectral form");
                    }
                    
                    
                    StAmpl[StNModes] = ampl;
                    StMode[3*StNModes+0] = kx;
                    StMode[3*StNModes+1] = ky;
                    StMode[3*StNModes+2] = kz;
                    PRINT_STATUS("Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,iky,ikz,kx,ky,kz,ampl);
                    StNModes++;
                    
                    
#if NUMDIMS > 1
                    StAmpl[StNModes] = ampl;
                    StMode[3*StNModes+0] = kx;
                    StMode[3*StNModes+1] = -ky;
                    StMode[3*StNModes+2] = kz;
                    PRINT_STATUS("Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,-iky,ikz,kx,-ky,kz,ampl);
                    StNModes++;
                    
#if NUMDIMS > 2
                    StAmpl[StNModes] = ampl;
                    StMode[3*StNModes+0] = kx;
                    StMode[3*StNModes+1] = ky;
                    StMode[3*StNModes+2] = -kz;
                    PRINT_STATUS("Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,iky,-ikz,kx,ky,-kz,ampl);
                    StNModes++;
                    
                    StAmpl[StNModes] = ampl;
                    StMode[3*StNModes+0] = kx;
                    StMode[3*StNModes+1] = -ky;
                    StMode[3*StNModes+2] = -kz;
                    PRINT_STATUS("Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f",StNModes,ikx,-iky,-ikz,kx,-ky,-kz,ampl);
                    StNModes++;
#endif
#endif
                }
            }
        }
    }
    
    StTPrev = -1;
    
    StRng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(StRng, All.StSeed);
    
    st_init_ouseq();
    st_calc_phases();
    
    PRINT_STATUS(" ..calling set_turb_ampl in init_turb\n");
    set_turb_ampl();
    
    StTPrev = All.Ti_Current;
}


void st_init_ouseq(void)
{
    int i;
    for(i = 0;i<6*StNModes;i++)
    {
        StOUPhases[i] = st_grn()*StOUVar;
    }
}


void st_update_ouseq(void)
{
    int i;
    double damping = exp( -All.StDtFreq/All.StDecay);
    for(i = 0;i<6*StNModes;i++)
    {
        StOUPhases[i] = StOUPhases[i] * damping + StOUVar * sqrt(1.-damping*damping)*st_grn();
    }
}


double st_grn(void)
{
    double r0 = gsl_rng_uniform(StRng);
    double r1 = gsl_rng_uniform(StRng);
    return sqrt(2. * log(1. / r0) ) * cos(2. * M_PI * r1);
}



void st_calc_phases(void)
{
    int i,j;
    for(i = 0; i < StNModes;i++)
    {
        double ka = 0.;
        double kb = 0.;
        double kk = 0.;
        
        int dim = NUMDIMS;
        
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
            
            StAka[3*i+j] = All.StSolWeight*curla+(1.-All.StSolWeight)*divb;
            StAkb[3*i+j] = All.StSolWeight*curlb+(1.-All.StSolWeight)*diva;
        }
    }
}


void set_turb_ampl(void)
{
    int i; double delta = (All.Ti_Current - StTPrev) * All.Timebase_interval / All.cf_hubble_a;
    double e_diss_sum=0, e_drive_sum=0, glob_diss_sum=0, glob_drive_sum=0;
    if(delta >= All.StDtFreq)
    {
        if(delta > 0)
        {
            PRINT_STATUS(" ..updating dudt_*\n");
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
        PRINT_STATUS(" ..st_update_ouseq() ... \n");
        st_update_ouseq();
        PRINT_STATUS(" ..st_calc_phases() ... \n");
        st_calc_phases();
        StTPrev = StTPrev + All.StDtFreq/All.Timebase_interval;
        PRINT_STATUS(" ..updated turbulent stirring field at time %f.\n", StTPrev * All.Timebase_interval);
    }
}



void add_turb_accel()
{
    int i, j, m;
    double acc[3];
    
    set_turb_ampl();
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0)
        {
            double fx = 0;
            double fy = 0;
            double fz = 0;
            
            for(m = 0;m<StNModes;m++) //calc force
            {
                double kxx = StMode[3*m+0]*P[i].Pos[0];
                double kyy = StMode[3*m+1]*P[i].Pos[1];
                double kzz = StMode[3*m+2]*P[i].Pos[2];
                double kdotx = kxx+kyy+kzz;
                double ampl = StAmpl[m];
                
                double realt = cos(kdotx);
                double imagt = sin(kdotx);
                
                fx += ampl*(StAka[3*m+0]*realt - StAkb[3*m+0]*imagt);
                fy += ampl*(StAka[3*m+1]*realt - StAkb[3*m+1]*imagt);
                fz += ampl*(StAka[3*m+2]*realt - StAkb[3*m+2]*imagt);
            }
            
            fx *= 2.*All.StAmplFac*StSolWeightNorm;
            fy *= 2.*All.StAmplFac*StSolWeightNorm;
            fz *= 2.*All.StAmplFac*StSolWeightNorm;
            
            if(P[i].Mass > 0.)
            {
                acc[0] = fx;
                acc[1] = fy;
                acc[2] = 0;
#if (NUMDIMS > 2)
                acc[2] = fz;
#endif
                for(j = 0; j < 3; j++)
                {
                    SphP[i].TurbAccel[j] = acc[j];
                }
            } else {
                SphP[i].TurbAccel[0]=SphP[i].TurbAccel[1]=SphP[i].TurbAccel[2]=0;
            }
        }
    }
    PRINT_STATUS("Finished turbulent accel computation");
}



void do_turb_driving_step_first_half(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    add_turb_accel();
    int i, j; integertime ti_step, tstart, tend; double dvel[3], dt_gravkick;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0; tstart = P[i].Ti_begstep; tend = P[i].Ti_begstep + ti_step / 2;	/* beginning / midpoint of step */
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


void do_turb_driving_step_second_half(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    int i, j; integertime ti_step, tstart, tend; double dvel[3], dt_gravkick;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0; tstart = P[i].Ti_begstep + ti_step / 2; tend = P[i].Ti_begstep + ti_step;	/* midpoint/end of step */
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
        fprintf(FdTurb, "%g %g %g %g %g %g %g\n", All.Time, mach, (glob_ekin + glob_ethermal) / glob_mass, glob_dudt_drive / glob_mass,
                glob_dudt_diss / glob_mass, All.TurbInjectedEnergy / glob_mass, All.TurbDissipatedEnergy / glob_mass);
#endif
}

#endif
