#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/* this file handles the FIRE short-range radiation-pressure and
    photo-ionization terms. written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#if defined(GALSF_FB_FIRE_RT_LOCALRP) /* first the radiation pressure coupled in the immediate vicinity of the star */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here! -- */
void radiation_pressure_winds_consolidated(void)
{
    double age_threshold_in_gyr = 0.15; // don't bother for older populations, they contribute negligibly here //
#if defined(SINGLE_STAR_SINK_DYNAMICS) || (FLAG_NOT_IN_PUBLIC_CODE > 2)
    age_threshold_in_gyr = 1.0e10; // for the single-star problems, or updated algorithm [where it adds little expense], we want to include everything, for completeness //
#endif
    if(All.RP_Local_Momentum_Renormalization<=0) return;
    Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
    PRINT_STATUS("Local Radiation-Pressure acceleration calculation");
    MyDouble *pos; int N_MAX_KERNEL,N_MIN_KERNEL,MAXITER_FB,NITER,startnode,dummy,numngb_inbox,i,j,k,n;
    double h,wt_sum,delta_v_imparted_rp=0,total_n_wind=0,total_mom_wind=0,total_prob_kick=0,avg_v_kick=0,avg_taufac=0;


    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3))))
        {
            double star_age = evaluate_stellar_age_Gyr(i);
            if( (star_age < age_threshold_in_gyr) && (P[i].Mass > 0) && (P[i].DensAroundStar > 0) )
            {
                /* calculate some basic luminosity properties of the stars */
                double lm_ssp = evaluate_light_to_mass_ratio(star_age, i); // light-to-mass ratio in solar
                double lum_cgs = (lm_ssp * SOLAR_LUM_CGS) * (P[i].Mass*UNIT_MASS_IN_SOLAR); // total L in CGS of star particle
                double f_lum_ion = particle_ionizing_luminosity_in_cgs(i) / lum_cgs; f_lum_ion=DMAX(0.,DMIN(1.,f_lum_ion)); // fraction of luminosity in H-ionizing radiation
                double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
                double dE_over_c = All.RP_Local_Momentum_Renormalization * lum_cgs * (dt*UNIT_TIME_IN_CGS) / C_LIGHT_CGS; // total photon momentum emitted in timestep, in CGS (= L*dt/c)
                dE_over_c /= (UNIT_MASS_IN_CGS * UNIT_VEL_IN_CGS); // total photon momentum now in code units
                total_prob_kick += dE_over_c; // sum contributions

                /* calculate some pre-amble properties */
#if defined(GALSF_FB_FIRE_RT_LOCALRP_OPTIMIZERS_TEST)
                double RtauMax = DMIN( 10./(UNIT_LENGTH_IN_KPC*All.cf_atime) , 4.*P[i].Hsml );
#else
                double RtauMax = P[i].Hsml * (5. + 2.0 * rt_kappa(i,RT_FREQ_BIN_FIRE_UV) * P[i].Hsml*P[i].DensAroundStar*All.cf_a2inv); // guess search radius which is a few H, plus larger factor if optically thick //
                RtauMax = DMAX( 1./(UNIT_LENGTH_IN_KPC*All.cf_atime) , DMIN( 10./(UNIT_LENGTH_IN_KPC*All.cf_atime) , RtauMax )); // restrict to 1-10 kpc here
#endif
                
#ifndef GALSF_FB_FIRE_RT_CONTINUOUSRP
                /* if kicks are stochastic, we don't want to waste time doing a neighbor search every timestep; it can be much faster to pre-estimate the kick probabilities */
                double v_wind_threshold = 15. / UNIT_VEL_IN_KMS; // unit velocity for kicks
#if defined(SINGLE_STAR_SINK_DYNAMICS) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
                v_wind_threshold = 0.2 / UNIT_VEL_IN_KMS; // for this module use lower unit mass for kicks
#endif
                double rho_phys=P[i].DensAroundStar*All.cf_a3inv, h_phys=P[i].Hsml*All.cf_atime; // density and h in -physical- units
                double v_grav_guess; v_grav_guess = DMIN( 1.82*(65.748/UNIT_VEL_IN_KMS)*pow(1.+rho_phys*UNIT_DENSITY_IN_NHCGS,-0.25) , sqrt(All.G*(P[i].Mass + NORM_COEFF*rho_phys*h_phys*h_phys*h_phys)/h_phys) ); // don't want to 'under-kick' if there are small local characteristic velocities in the region of interest
                delta_v_imparted_rp = v_wind_threshold; // always couple this 'discrete' kick, to avoid having to couple every single timestep for every single star particle
                double dv_imparted_perpart_guess = (dE_over_c/P[i].Mass); // estimate of summed dv_imparted [in code units] from single-scattering: = momentum/mass of particle
                double tau_IR_guess = rt_kappa(i,RT_FREQ_BIN_FIRE_IR) * rho_phys*h_phys; // guess of IR optical depth. everything in physical code units //
                dv_imparted_perpart_guess += (dE_over_c/P[i].Mass) * tau_IR_guess; // estimate of additional IR term [1+tau_IR]*L/c assumed here as coupling //
                double prob = dv_imparted_perpart_guess / delta_v_imparted_rp; prob *= 2000.; // need to include a buffer for errors in the estimates above
                double p_random = get_random_number(P[i].ID+ThisTask+i+2); // random number for use below
                if(p_random <= prob) // alright, its worth doing the loop!
#endif
                { // within loop
                    /* ok, now open the neighbor list for the star particle */
#if defined(GALSF_FB_FIRE_RT_LOCALRP_OPTIMIZERS_TEST)
                    N_MIN_KERNEL=1; N_MAX_KERNEL=256; MAXITER_FB=10; h=0.5*P[i].Hsml;
#else
                    N_MIN_KERNEL=10; N_MAX_KERNEL=256; MAXITER_FB=100; h=1.0*P[i].Hsml;
#endif
                    NITER=0; wt_sum=0; startnode=All.MaxPart; dummy=0; numngb_inbox=0; pos=P[i].Pos;
                    if(h<=0) {h=All.ForceSoftening[0];} else {if(h>RtauMax) {h=RtauMax;}}
                    do {
                        numngb_inbox = ngb_treefind_variable_targeted(pos, h, -1, &startnode, 0, &dummy, &dummy, 1); // search for gas (2^0=1 for bitflag), use the 'see one way' search, since weights below are all within-kernel, for now
                        if((numngb_inbox>=N_MIN_KERNEL)&&(numngb_inbox<=N_MAX_KERNEL))
                        {
                            wt_sum=0; /* note these lines and many below assume 3D sims! */
                            for(n=0; n<numngb_inbox; n++)
                            {
                                j = Ngblist[n];
#ifdef BH_WIND_SPAWN
                                if(P[j].ID == All.AGNWindID) {continue;} // dont couple to jet cells
#endif
                                if((P[j].Mass>0) && (SphP[j].Density>0))
                                {
                                    double dp[3],r2=0; for(k=0;k<3;k++) {dp[k]=P[j].Pos[k]-P[i].Pos[k];}
                                    NEAREST_XYZ(dp[0],dp[1],dp[2],1); for(k=0;k<3;k++) {r2+=dp[k]*dp[k];} /* find the closest image in the given box size */
                                    if(r2>=h*h || r2<=0) {continue;}
                                    double h_eff_j = Get_Particle_Size(j); wt_sum += h_eff_j*h_eff_j; // weight factor for neighbors
                                } /* if( (P[j].Mass>0) && (SphP[j].Density>0) ) */
                            } /* for(n=0; n<numngb_inbox; n++) */
                            if(wt_sum <= 0) {h*= 1.2123212335; startnode=All.MaxPart;} /* wt_sum <= 0; no particles found inside corners - expand */
                        } else {
                            startnode=All.MaxPart;
                            if(numngb_inbox<N_MIN_KERNEL) {if(numngb_inbox<=0) {h*=2.0;} else {if(NITER<=5) {h*=pow((float)numngb_inbox/(float)N_MIN_KERNEL,-0.3333);} else {h*=1.26;}}} /* iterate until find appropriate > N_MIN # particles */
                            if(numngb_inbox>N_MAX_KERNEL) {if(NITER<=5) {h*=pow((float)numngb_inbox/(float)N_MAX_KERNEL,-0.3333);} else {h/=1.31;}} /* iterate until find appropriate < N_MAX # particles */
                            if((NITER>MAXITER_FB/2) && (N_MIN_KERNEL>2)) {N_MIN_KERNEL/=2; N_MAX_KERNEL*=2; if(N_MIN_KERNEL<2) {N_MIN_KERNEL=2;} if(N_MAX_KERNEL>1000) {N_MAX_KERNEL=1000;}} // expand tolerance if we are doing a lot of iterations here //
                        }
                        if(h>20.*RtauMax) {h=20.*RtauMax; if(NITER<MAXITER_FB-1) {NITER=MAXITER_FB-1;}} /* if h exceeds the maximum now, set it to that value, and set NITER to maximum to force end of iteration */
                        NITER++;
                    } while( (startnode >= 0) && (NITER<=MAXITER_FB) );

                    if(wt_sum > 0)  /* found at least one massive neighbor, can proceed */
                    {
                        for(n=0; n<numngb_inbox; n++)
                        {
                            j = Ngblist[n];
#ifdef BH_WIND_SPAWN
                            if(P[j].ID == All.AGNWindID) {continue;} // dont couple to jet cells
#endif
                            if((P[j].Mass>0) && (SphP[j].Density>0))
                            {
                                double dp[3],r2=0; for(k=0;k<3;k++) {dp[k]=P[j].Pos[k]-P[i].Pos[k];}
                                NEAREST_XYZ(dp[0],dp[1],dp[2],1); for(k=0;k<3;k++) {r2+=dp[k]*dp[k];} /* find the closest image in the given box size */
                                if(r2>=h*h || r2<=0) {continue;}
                                double h_eff_i = DMIN(h, Get_Particle_Size(i)), h_eff_j = Get_Particle_Size(j);
                                r2 += MIN_REAL_NUMBER + (h_eff_i/5.)*(h_eff_i/5.); // just a small number to prevent errors on near-overlaps
                                double wk = h_eff_j*h_eff_j / wt_sum; // dimensionless weight factor

                                /* first -- share out the UV luminosity among the local neighbors, weighted by the gas kernel */
                                double dv_imparted_singlescattering = wk * (dE_over_c / P[j].Mass); // fractional initial photon momentum seen by this neighbor
                                /* velocity imparted by IR acceleration : = kappa*flux/c, flux scales as 1/r2 from source, kappa with metallicity */
                                double kappa_ir_codeunits = rt_kappa(j,RT_FREQ_BIN_FIRE_IR); // opacity in code units
                                double dv_imparted_multiplescattering = All.RP_Local_Momentum_Renormalization * (dE_over_c / P[j].Mass) * kappa_ir_codeunits * (P[j].Mass/(4.*M_PI*r2*All.cf_atime*All.cf_atime));
#if defined(GALSF_FB_FIRE_RT_CONTINUOUSRP) || (FLAG_NOT_IN_PUBLIC_CODE > 2)
                                delta_v_imparted_rp = dv_imparted_multiplescattering + dv_imparted_singlescattering;
#else
                                prob = (dv_imparted_multiplescattering+dv_imparted_singlescattering) / delta_v_imparted_rp; if(prob>1) {delta_v_imparted_rp *= prob;}
                                if(n>0) {p_random=get_random_number(P[j].ID+P[i].ID +ThisTask+ 3);} //else p_random=0;
                                if(p_random < prob)
#endif
                                { /* open subloop with wind kick */
                                    if(delta_v_imparted_rp>1.e4/UNIT_VEL_IN_KMS) {delta_v_imparted_rp=1.e4/UNIT_VEL_IN_KMS;} /* limiter */
                                    /* collect numbers to output */
                                    total_n_wind += 1.0; total_mom_wind += P[j].Mass*delta_v_imparted_rp; avg_v_kick += delta_v_imparted_rp;
                                    avg_taufac +=  (P[j].Mass*delta_v_imparted_rp) * (dv_imparted_multiplescattering / (dE_over_c / P[j].Mass));

                                    /* determine the direction of the kick */
                                    double dir[3], norm=0;
#if defined(GALSF_FB_FIRE_RT_CONTINUOUSRP) || (FLAG_NOT_IN_PUBLIC_CODE > 2)
                                    delta_v_imparted_rp = dv_imparted_multiplescattering; // ir kick: directed along opacity gradient //
                                    for(k=0;k<3;k++) {dir[k]=-P[j].GradRho[k]; norm+=dir[k]*dir[k];} // based on density gradient near star //
#else
                                    if(dv_imparted_singlescattering > dv_imparted_multiplescattering) {for(k=0;k<3;k++) {dir[k]=dp[k]; norm+=dir[k]*dir[k];}} // if kick is primarily from uv, then orient directly //
                                        else {for(k=0;k<3;k++) {dir[k]=-P[j].GradRho[k]; norm+=dir[k]*dir[k];}} // otherwise, along opacity gradient //
#endif
                                    if(norm>0) {norm=sqrt(norm); for(k=0;k<3;k++) dir[k] /= norm;} else {dir[0]=0; dir[1]=0; dir[2]=1; norm=1;}
                                    for(k=0;k<3;k++) {P[j].Vel[k] += delta_v_imparted_rp * All.cf_atime * dir[k]; SphP[j].VelPred[k] += delta_v_imparted_rp * All.cf_atime * dir[k];} /* apply the kick [put into comoving code units as oppropriate */

#if defined(GALSF_FB_FIRE_RT_CONTINUOUSRP) || (FLAG_NOT_IN_PUBLIC_CODE > 2)
                                    /* if we're not forcing the kick orientation, need to separately apply the UV kick */
                                    delta_v_imparted_rp = dv_imparted_singlescattering; // uv kick: directed from star //
                                    norm=0; for(k=0;k<3;k++) {dir[k]=dp[k]; norm+=dir[k]*dir[k];}
                                    if(norm>0) {norm=sqrt(norm); for(k=0;k<3;k++) {dir[k] /= norm;}} else {dir[0]=0; dir[1]=0; dir[2]=1; norm=1;}
                                    for(k=0; k<3; k++) {P[j].Vel[k] += delta_v_imparted_rp * All.cf_atime * dir[k]; SphP[j].VelPred[k] += delta_v_imparted_rp * All.cf_atime * dir[k];} /* apply the kick */
#endif
                                } /* closes if(get_random_number(P[i].ID + 2) < prob) */
                            } /* if( (P[j].Mass>0) && (SphP[j].Density>0) ) */
                        } /* for(n=0; n<numngb_inbox; n++) */
                    } /* if (rho>0) */
                } // // within loop
            } // star age, mass check:: (star_age < 0.1) && (P[i].Mass > 0) && (P[i].DensAroundStar > 0)
        } // particle type check::  if((P[i].Type == 4)....
    } // main particle loop for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    myfree(Ngblist);

    double totMPI_n_wind=0,totMPI_mom_wind=0,totMPI_avg_v=0,totMPI_avg_taufac=0,totMPI_prob_kick=0;
    MPI_Reduce(&total_n_wind, &totMPI_n_wind, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_mom_wind, &totMPI_mom_wind, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_v_kick, &totMPI_avg_v, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_taufac, &totMPI_avg_taufac, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_prob_kick, &totMPI_prob_kick, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
#ifdef IO_REDUCED_MODE
        if(totMPI_n_wind>0)
#endif
        if(totMPI_prob_kick>0)
        {
            totMPI_avg_v /= MIN_REAL_NUMBER + totMPI_n_wind; totMPI_avg_taufac /= MIN_REAL_NUMBER + totMPI_mom_wind;
            fprintf(FdMomWinds, "%.16g %g %g %g %g %g \n", All.Time,totMPI_n_wind,totMPI_prob_kick,totMPI_mom_wind,totMPI_avg_v,totMPI_avg_taufac); fflush(FdMomWinds);
            PRINT_STATUS(" ..Ncells_pushed=%g (L/c)dt=%g dP_coupled=%g <dv_cell>=%g <dP_multi/dP_single>=%g",totMPI_n_wind,totMPI_prob_kick,totMPI_mom_wind,totMPI_avg_v,totMPI_avg_taufac);
        }
    } // if(ThisTask==0)
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin && ThisTask == 0) {fflush(FdMomWinds);}
    PRINT_STATUS(" ..completed local Radiation-Pressure acceleration");
    CPU_Step[CPU_LOCALWIND] += measure_time(); /* collect timings and reset clock for next timing */
} // end routine :: void radiation_pressure_winds_consolidated(void)

#endif /* closes defined(GALSF_FB_FIRE_RT_LOCALRP)  */





#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
/* Routines for simple FIRE local photo-ionization heating feedback model. This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here! -- */


void HII_heating_singledomain(void)    /* this version of the HII routine only communicates with particles on the same processor */
{
#ifdef RT_CHEM_PHOTOION
    return; // the work here is done in the actual RT routines if this switch is enabled //
#endif
    if(All.HIIRegion_fLum_Coupled<=0) {return;}
    if(All.Time<=0) {return;}
    PRINT_STATUS("Local HII-Region photo-heating/ionization calculation");
    MyDouble *pos; MyFloat h_i, dt, rho; int startnode, numngb, j, n, i, NITER_HIIFB, MAX_N_ITERATIONS_HIIFB, jnearest,already_ionized,do_ionize,dummy;
    double total_N_ionizing_part=0,total_Ndot_ionizing=0,total_m_ionized=0,total_N_ionized=0,avg_RHII=0,mionizable=0,mionized=0,mion_actual=0;
    double RHII,RHIIMAX,R_search,rnearest,stellum,prob,rho_j,prandom,m_available,m_effective,RHII_initial,RHIImultiplier;
    double uion; uion = HIIRegion_Temp / (0.59 * (5./3.-1.) * U_TO_TEMP_UNITS); /* assume fully-ionized gas with gamma=5/3; this is a global variable below */
    Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
    MAX_N_ITERATIONS_HIIFB = 5; NITER_HIIFB = 0;

    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef BH_HII_HEATING
        if((P[i].Type == 5)||(((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3))))))
#else
        if((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3))))
#endif
        {
            dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
            if(P[i].Type == 5) {dt = P[i].dt_since_last_gas_search;}
#endif
            if(dt<=0) {continue;} // don't keep going with this loop

            stellum = All.HIIRegion_fLum_Coupled * particle_ionizing_luminosity_in_cgs(i); // ionizing luminosity in cgs [will be appropriately weighted for assumed spectral shape]
#ifdef CHIMES_HII_REGIONS
            stellum = chimes_ion_luminosity(evaluate_stellar_age_Gyr(i)*1000.,P[i].Mass*UNIT_MASS_IN_SOLAR) * 4.68e-11; // chimes ionizing photon flux rescaled to mean spectrum here appropriately (~29eV per photon)
#endif
            if(stellum <= 0) {continue;}
            pos = P[i].Pos; rho = P[i].DensAroundStar; h_i = PPP[i].Hsml;
            RHII = 4.78e-9*pow(stellum,0.333)*pow(rho*All.cf_a3inv*UNIT_DENSITY_IN_CGS,-0.66667); // Stromgren radius, RHII, computed using a case B recombination coefficient at 10^4 K of 2.59e-13 cm^3 s^-1, and assuming a Hydrogen mass fraction ~0.74.
            RHII /= All.cf_atime*UNIT_LENGTH_IN_CGS; // convert to code units
            RHIIMAX = 2. * 240.0*pow(stellum,0.5) / (All.cf_atime*UNIT_LENGTH_IN_CGS); // crude estimate of where flux falls below cosmic background, x2 safety factor
#ifdef GALSF_FB_FIRE_RT_HIIHEATING_OPTIMIZERS_TEST
            if(RHIIMAX < 1.2*h_i) {RHIIMAX=1.2*h_i;} // limit max search radius: can't be below 2x kernel size
            if(RHIIMAX > 4.0*h_i) {RHIIMAX=4.*h_i;} // limit search radius to 10x kernel size
#else
            if(RHIIMAX < 2.0*h_i) {RHIIMAX=2.0*h_i;} // limit max search radius: can't be below 2x kernel size
            if(RHIIMAX > 10.0*h_i) {RHIIMAX=10.*h_i;} // limit search radius to 10x kernel size
#endif
            mionizable = NORM_COEFF*rho*RHII*RHII*RHII; // estimated ionizable gas mass in code units, based on the gas density at star location [will be rescaled]
            double M_ionizing_emitted = (3.05e10 * PROTONMASS_CGS) * stellum * (dt * UNIT_TIME_IN_CGS) ; // number of ionizing photons times proton mass, gives max mass ionized [in cgs]
            mionizable = DMIN( mionizable , M_ionizing_emitted/UNIT_MASS_IN_CGS ); // in code units
            if(RHII > RHIIMAX) {RHII = RHIIMAX;} // limit initial guess to max
            if(RHII < 0.3*h_i) {RHII=0.3*h_i;} // limit initial guess to above 1/2 kernel, so can find neighbors
            RHII_initial=RHII;
            total_N_ionizing_part += 1; total_Ndot_ionizing += stellum * (3.05e10/HYDROGEN_MASSFRAC); // increment before loop //

            prandom = get_random_number(P[i].ID + 7); // pre-calc the (eventually) needed random number
            // guesstimate if this is even close to being interesting for the particle masses of interest
#ifdef GALSF_FB_FIRE_RT_HIIHEATING_OPTIMIZERS_TEST
            if(prandom < 3.0*mionizable/P[i].Mass) // prandom > this, won't be able to ionize anything interesting
#else
            if(prandom < 5.0*mionizable/P[i].Mass) // prandom > this, won't be able to ionize anything interesting
#endif
            {
                mionized=0.0; startnode = All.MaxPart; jnearest=-1; rnearest=MAX_REAL_NUMBER; dummy=0; NITER_HIIFB=0;
                do {
                    double RHII_2 = RHII*RHII;
                    jnearest=-1; rnearest=MAX_REAL_NUMBER;
                    R_search = RHII;
                    if(h_i>0.5*R_search) {R_search=0.5*h_i;}
                    numngb = ngb_treefind_variable_targeted(pos, R_search, -1, &startnode, 0, &dummy, &dummy, 1); // search for gas (2^0=1 for bitflag), use the 'see one way' search, since weights below are all within-kernel, for now
                    if(numngb>0)
                    {
                        int ngb_list_touse[numngb]; for(n=0; n<numngb; n++) {ngb_list_touse[n]=Ngblist[n];}
                        for(n = 0; n < numngb; n++)
                        {
                            if(mionized>=mionizable) {break;}
                            j = ngb_list_touse[n];
                            if(P[j].Mass <= 0 || P[j].Type !=0) {continue;}
                            if(SphP[j].DelayTimeHII > 0) {continue;}
                            if(P[j].Type == 0 && P[j].Mass > 0)
                            {
                                double dx=pos[0]-P[j].Pos[0], dy=pos[1]-P[j].Pos[1], dz=pos[2]-P[j].Pos[2];
                                NEAREST_XYZ(dx,dy,dz,1); /*  now find the closest image in the given box size */
                                double r2 = dx * dx + dy * dy + dz * dz;
                                if(r2 > RHII_2) {continue;}
                                double r=sqrt(r2), u=0;
                                /* check whether the particle is already ionized */
                                already_ionized = 0; rho_j = Get_Gas_density_for_energy_i(j);
                                if(SphP[j].InternalEnergy<SphP[j].InternalEnergyPred) {u=SphP[j].InternalEnergy;} else {u=SphP[j].InternalEnergyPred;}
                                if(SphP[j].DelayTimeHII > 0) {already_ionized=1;}
#if !defined(CHIMES_HII_REGIONS)
                                if(u>uion) {already_ionized=1;}
#endif
                                if(already_ionized) continue;
                                /* now, if inside RHII and mionized<mionizeable and not already ionized, can be ionized! */
                                do_ionize=0; prob=0;
                                if((r<=RHII)&&(already_ionized==0)&&(mionized<mionizable))
                                {
                                    m_effective = P[j].Mass*(SphP[j].Density/rho);
                                    // weight by density b/c of how the recombination rate in each particle scales
                                    m_available = mionizable-mionized;
                                    if(m_effective<=m_available) {
                                        do_ionize=1; prob = 1.001; // we can ionize the entire cell
                                    } else {
                                        prob = m_available/m_effective; // partial ionization. determine randomly if ionized
                                        if(prandom < prob) {do_ionize=1;}
                                    } // if(m_effective<=m_available) {
                                    if(do_ionize==1)
                                    {
                                        already_ionized = do_the_local_ionization(j,dt,i);
                                        total_N_ionized += 1;
                                        mion_actual += P[j].Mass;
                                        avg_RHII += P[j].Mass*r*All.cf_atime*UNIT_LENGTH_IN_KPC;
                                    }
                                    mionized += prob*m_effective;
                                } // if((r<=RHII)&&(already_ionized==0)&&(mionized<mionizable))

                                /* if nearest un-ionized particle, mark as such */
                                if((r<rnearest)&&(already_ionized==0)) {rnearest = r; jnearest = j;}
                            } // if(P[j].Type == 0 && P[j].Mass > 0)
                        } // for(n = 0; n < numngb; n++)
                    } // if(numngb>0)

                    // if still have photons and jnearest is un-ionized
                    if((mionized<mionizable)&&(jnearest>=0))
                    {
                        j=jnearest; m_effective=P[j].Mass*(SphP[j].Density/rho); m_available=mionizable-mionized; prob=m_available/m_effective; do_ionize=0;
                        if(prandom < prob) {do_ionize=1;}
                        if(do_ionize==1)
                        {
                            already_ionized = do_the_local_ionization(j,dt,i);
                            total_N_ionized += 1;
                            mion_actual += P[j].Mass;
                            double dx=pos[0]-P[j].Pos[0],dy=pos[1]-P[j].Pos[1],dz=pos[2]-P[j].Pos[2],r2; NEAREST_XYZ(dx,dy,dz,1); r2=dx*dx+dy*dy+dz*dz;
                            avg_RHII += P[j].Mass*sqrt(r2)*All.cf_atime*UNIT_LENGTH_IN_KPC;
                        }
                        mionized += prob*m_effective;
                    } // if((mionized<mionizable)&&(jnearest>=0))

                    /* now check if we have ionized sufficient material, and if not, iterate with larger regions until we do */
                    RHIImultiplier=1.10;
                    if(mionized < 0.95*mionizable)
                    {
                        /* ok, this guy did not find enough gas to ionize, it needs to expand its search */
#ifdef GALSF_FB_FIRE_RT_HIIHEATING_OPTIMIZERS_TEST
                        if((RHII >= DMIN(30.0*RHII_initial, RHIIMAX)) || (NITER_HIIFB >= MAX_N_ITERATIONS_HIIFB))
#else
                        if((RHII >= DMAX(30.0*RHII_initial, RHIIMAX)) || (NITER_HIIFB >= MAX_N_ITERATIONS_HIIFB))
#endif
                        {
                            /* we're done looping, this is just too big an HII region */
                            mionized = 1.001*mionizable;
                        } else {
                            /* in this case we're allowed to keep expanding RHII */
                            if(mionized <= 0)
                            {
                                RHIImultiplier = 2.0;
                            } else {
                                RHIImultiplier = pow(mionized/mionizable , -0.333);
                                if(RHIImultiplier>5.0) {RHIImultiplier=5.0;}
                                if(RHIImultiplier<1.26) {RHIImultiplier=1.26;}
                            } // if(mionized <= 0)
                            RHII *= RHIImultiplier; if(RHII>1.26*RHIIMAX) {RHII=1.26*RHIIMAX;}
                            startnode=All.MaxPart; // this will trigger the while loop to continue
                        } // if((RHII >= 5.0*RHII_initial)||(RHII>=RHIIMAX)||(NITER_HIIFB >= MAX_N_ITERATIONS_HIIFB))
                    } // if(mionized < 0.95*mionizable)
                    NITER_HIIFB++;
                } while(startnode >= 0 && mionized<mionizable);
                if(mion_actual>0) {total_m_ionized += mion_actual;}
            } // if(prandom < 2.0*mionizable/P[j].Mass)
        } // if((P[i].Type == 4)||(P[i].Type == 2)||(P[i].Type == 3))
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    myfree(Ngblist);

    double totMPI_N_ionizing_part=0,totMPI_Ndot_ionizing=0,totMPI_m_ionized=0,totMPI_avg_RHII=0,totMPI_N_ionized=0;
    MPI_Reduce(&total_N_ionizing_part, &totMPI_N_ionizing_part, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_Ndot_ionizing, &totMPI_Ndot_ionizing, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_m_ionized, &totMPI_m_ionized, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_N_ionized, &totMPI_N_ionized, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_RHII, &totMPI_avg_RHII, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        if(totMPI_N_ionizing_part>0)
        {
            if(totMPI_m_ionized>0) {totMPI_avg_RHII /= totMPI_m_ionized;}
            PRINT_STATUS(" ..Nsources=%g with dN/dt=%g/s ionized N=%g (M=%g sol) cells in <R_HII>=%g kpc",totMPI_N_ionizing_part,totMPI_Ndot_ionizing,totMPI_N_ionized,totMPI_m_ionized*UNIT_MASS_IN_SOLAR,totMPI_avg_RHII);
            fprintf(FdHIIHeating,"%.16g %g %g %g %g %g \n",All.Time,totMPI_N_ionizing_part,totMPI_Ndot_ionizing,totMPI_N_ionized,totMPI_m_ionized*UNIT_MASS_IN_SOLAR,totMPI_avg_RHII); fflush(FdHIIHeating);
        }
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) {fflush(FdHIIHeating);}
    } // ThisTask == 0
    CPU_Step[CPU_HIIHEATING] += measure_time();
} // void HII_heating_singledomain(void)



int do_the_local_ionization(int target, double dt, int source)
{
#if defined(CHIMES_HII_REGIONS) // set a number of chimes-specific quantities here //
    int k,age_bin=0; double stellar_age_myr=1000.*evaluate_stellar_age_Gyr(source), log_age_Myr=log10(stellar_age_myr); // determine stellar age bin
    if(log_age_Myr<CHIMES_LOCAL_UV_AGE_LOW) {age_bin=0;} else if(log_age_Myr < CHIMES_LOCAL_UV_AGE_MID) {age_bin = (int) floor(((log_age_Myr - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW) + 1);}
    else {age_bin = (int) floor((((log_age_Myr - CHIMES_LOCAL_UV_AGE_MID) / CHIMES_LOCAL_UV_DELTA_AGE_HI) + ((CHIMES_LOCAL_UV_AGE_MID - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW)) + 1); if(age_bin > CHIMES_LOCAL_UV_NBINS - 1) {age_bin = CHIMES_LOCAL_UV_NBINS - 1;}}
    // reset all of the HII-region chimes quantities to null
    for(k=0;k<CHIMES_LOCAL_UV_NBINS;k++) {SphP[target].Chimes_fluxPhotIon_HII[k]=0; SphP[target].Chimes_G0_HII[k]=0;}
    // set the quantities desired for this age bin specifically: need a softened radius, for use here //
    double dp[3],r2=0,stellar_mass=P[source].Mass*UNIT_MASS_IN_SOLAR; for(k=0;k<3;k++) {dp[k]=P[source].Pos[k]-P[target].Pos[k];}
    NEAREST_XYZ(dp[0],dp[1],dp[2],1); for(k=0;k<3;k++) {dp[k]*=All.cf_atime*UNIT_LENGTH_IN_CGS; r2+=dp[k]*dp[k];} // separation in cgs
    double eps_cgs=KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER*ForceSoftening_KernelRadius(source)*All.cf_atime*UNIT_LENGTH_IN_CGS; // plummer equivalent softening
    r2+=eps_cgs*eps_cgs; // gravitational Softening (cgs units)
    SphP[target].Chimes_fluxPhotIon_HII[age_bin] = (1.0 - All.Chimes_f_esc_ion) * chimes_ion_luminosity(stellar_age_myr, stellar_mass) / r2; // cgs flux of H-ionising photons per second seen by the star particle
    SphP[target].Chimes_G0_HII[age_bin] = (1.0 - All.Chimes_f_esc_G0) * chimes_G0_luminosity(stellar_age_myr, stellar_mass) / r2; // cgs flux in the 6-13.6 eV band

#else

    SphP[target].InternalEnergy = DMAX(SphP[target].InternalEnergy , HIIRegion_Temp / (0.59 * (5./3.-1.) * U_TO_TEMP_UNITS)); /* assume fully-ionized gas with gamma=5/3 */
    SphP[target].InternalEnergyPred = SphP[target].InternalEnergy; /* full reset of the internal energy */
    SphP[target].Ne = 1.0 + 2.0*yhelium(target); /* set the cell to fully ionized */

#endif

    SphP[target].DelayTimeHII = DMIN(dt, 10./UNIT_TIME_IN_MYR); /* tell the code to flag this in the cooling subroutine */
    return 1;
}

#endif
