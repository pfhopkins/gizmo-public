#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/* Routines for mechanical feedback/enrichment models: stellar winds, supernovae, etc
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef GALSF_FB_MECHANICAL

int N_Gas_Couplings_ThisTask; // define variable to use below to record if and how many times a coupling to a cell occurs in a timestep per MPI task

int addFB_evaluate_active_check(int i, int fb_loop_iteration);
int addFB_evaluate_active_check(int i, int fb_loop_iteration)
{
    if(P[i].Type <= 1) {return 0;} // note quantities used here must -not- change in the loop [hence not using mass here], b/c can change offsets for return from different processors, giving a negative mass and undefined behaviors
    if(PPP[i].Hsml <= 0) {return 0;}
    if(PPP[i].NumNgb <= 0) {return 0;}
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
    if(P[i].Type == 5 && !P[i].do_gas_search_this_timestep) {return 0;}
#endif
    if(P[i].SNe_ThisTimeStep>0) {if(fb_loop_iteration<0 || fb_loop_iteration==0) {return 1;}}
    return 0;
}


void determine_where_SNe_occur(void)
{
    if(All.Time<=0) return;
    int i; double dt,star_age,npossible,nhosttotal,ntotal,ptotal,dtmean,rmean;
    npossible=nhosttotal=ntotal=ptotal=dtmean=rmean=0;
    double mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean;
    mpi_npossible=mpi_nhosttotal=mpi_ntotal=mpi_ptotal=mpi_dtmean=mpi_rmean=0;
    // loop over particles //
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        P[i].SNe_ThisTimeStep=0;


#if defined(SINGLE_STAR_SINK_DYNAMICS)
        if(P[i].Type == 0) {continue;} // any non-gas type is eligible to be a 'star' here
#else
        if(All.ComovingIntegrationOn) {if(P[i].Type != 4) {continue;}} // in cosmological simulations, 'stars' have particle type=4
        if(All.ComovingIntegrationOn==0) {if((P[i].Type<2)||(P[i].Type>4)) {continue;}} // in non-cosmological sims, types 2,3,4 are valid 'stars'
#endif
        if(P[i].Mass<=0) {continue;}
        dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
        if(P[i].Type == 5) {dt = P[i].dt_since_last_gas_search;}
#endif
        if(dt<=0) {continue;} // no time, no events
        star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if(star_age<=0) {continue;} // unphysical age, no events
        // now use a calculation of mechanical event rates to determine where/when the events actually occur //
        npossible++;
        double RSNe = mechanical_fb_calculate_eventrates(i,dt);
        rmean += RSNe; ptotal += RSNe * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR);
#ifdef GALSF_SFR_IMF_SAMPLING
        if(P[i].Type<5) {if(P[i].IMF_NumMassiveStars>0) {P[i].IMF_NumMassiveStars=DMAX(0,P[i].IMF_NumMassiveStars-P[i].SNe_ThisTimeStep);}} // lose an O-star for every SNe //
#endif
        if(P[i].SNe_ThisTimeStep>0) {ntotal+=P[i].SNe_ThisTimeStep; nhosttotal++;}
        dtmean += dt;
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //

    MPI_Reduce(&dtmean, &mpi_dtmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rmean, &mpi_rmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ptotal, &mpi_ptotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nhosttotal, &mpi_nhosttotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ntotal, &mpi_ntotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&npossible, &mpi_npossible, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(ThisTask == 0)
    {
        if(mpi_ntotal > 0 && mpi_nhosttotal > 0 && mpi_dtmean > 0)
        if(mpi_npossible>0)
        {
            mpi_dtmean /= mpi_npossible; mpi_rmean /= mpi_npossible;
            fprintf(FdSneIIHeating, "%lg %g %g %g %g %g %g \n", All.Time,mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean); fflush(FdSneIIHeating);
        }
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) {fflush(FdSneIIHeating);}
    } // if(ThisTask == 0) //

} // void determine_where_SNe_occur() //



/* this is a temporary structure for quantities used ONLY in the loops below,
 for example for ensuring conservation if there are many overlapping events */
static struct temporary_mech_fb_data_tohold
{
    int N_injected; double m_injected, p_injected[3], KE_injected, TE_injected, Z_injected[NUM_METAL_SPECIES];
}
*LocalGasMechFBInfoTemp;





#define CORE_FUNCTION_NAME addFB_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_addFB    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_addFB  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(addFB_evaluate_active_check(i,loop_iteration)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

// define kernel structure (purely for convenience, will hold variables below) //
struct kernel_addFB {double dp[3], r, wk, dwk, hinv, hinv3, hinv4;};

struct OUTPUT_STRUCT_NAME
{
    MyFloat M_coupled, Area_weighted_sum[AREA_WEIGHTED_SUM_ELEMENTS];
}
*DATARESULT_NAME, *DATAOUT_NAME;


void particle2in_addFB(struct addFB_evaluate_data_in_ *in, int i, int loop_iteration)
{
    // pre-assign various values that will be used regardless of feedback physics //
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k]; in->Vel[k] = P[i].Vel[k];}
    double heff = PPP[i].Hsml / PPP[i].NumNgb; in->V_i = heff*heff*heff; in->Hsml = PPP[i].Hsml;
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.0;}
#endif
    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {in->Area_weighted_sum[k] = P[i].Area_weighted_sum[k];}
    in->Msne = 0; in->unit_mom_SNe = 0; in->SNe_v_ejecta = 0;
    if((P[i].DensAroundStar <= 0)||(P[i].Mass <= 0)) {return;} // events not possible [catch for mass->0]
    if(loop_iteration < 0) {in->Msne=P[i].Mass; in->unit_mom_SNe=1.e-4; in->SNe_v_ejecta=1.0e-4; return;} // weighting loop
    particle2in_addFB_fromstars(in,i,loop_iteration); // subroutine that actually deals with the assignment of feedback properties
    in->unit_mom_SNe = in->Msne * in->SNe_v_ejecta;
}

void out2particle_addFB(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    if(P[i].Mass > 0)
    {
        if(loop_iteration < 0)
        {
            int k=0, kmin=0, kmax=7; if(loop_iteration == -1) {kmin=kmax; kmax=AREA_WEIGHTED_SUM_ELEMENTS;}
#ifdef GALSF_USE_SNE_ONELOOP_SCHEME
            kmin=0; kmax=AREA_WEIGHTED_SUM_ELEMENTS;
#endif
            for(k=kmin;k<kmax;k++) {ASSIGN_ADD(P[i].Area_weighted_sum[k], out->Area_weighted_sum[k], mode);}
        } else {
            P[i].Mass -= out->M_coupled; if((P[i].Mass<0)||(isnan(P[i].Mass))) {P[i].Mass=0;}
        }
    }
}


/* here we have the subroutine that is the work center of this module. does the key calculations over neighbors and actually couples the relevant feedback quantities */
/*!   -- this subroutine [in all versions] writes important conservative variables to shared memory [updating the neighbor values]:
        we need to protect these writes for openmp (thread safety) below. the read-in values are themselves modified, so we need to protect -both- the read and write operations in all cases */

#ifdef GALSF_USE_SNE_ONELOOP_SCHEME

int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double u,r2,h2,kernel_zero,wk,dM,dP,E_coupled,dP_sum,dP_boost_sum;
    struct kernel_addFB kernel; struct addFB_evaluate_data_in_ local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1); wk=0;
    if(mode == 0) {particle2in_addFB(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* Load the data for the particle injecting feedback */
    if(local.Msne<=0) {return 0;} // no SNe for the origin particle! nothing to do here //
    if(local.Hsml<=0) {return 0;} // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml; kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    double unitlength_in_kpc=UNIT_LENGTH_IN_KPC * All.cf_atime, density_to_n=All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, unit_egy_SNe = 1.0e51/UNIT_ENERGY_IN_CGS; // some units (just used below, but handy to define for clarity) //


    // now define quantities that will be used below //
    double Esne51; Esne51 = 0.5*local.SNe_v_ejecta*local.SNe_v_ejecta*local.Msne / unit_egy_SNe;
    double RsneKPC, RsneKPC_0; RsneKPC=0.; RsneKPC_0=(0.0284/unitlength_in_kpc) * pow(1+Esne51,0.286); //Cioffi: weak external pressure
    double r2max_phys = 2.0/unitlength_in_kpc; r2max_phys *= r2max_phys; // no super-long-range effects allowed! (of course this is arbitrary in code units) //

    /* Now start the actual FB computation for this particle */
    if(mode == 0) {startnode = All.MaxPart;} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;} /* root node & node opening */
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}

            E_coupled = dP_sum = dP_boost_sum = 0;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Type != 0) {continue;} // require a gas particle //
                
                double Mass_j, InternalEnergy_j, rho_j, Vel_j[3]; // initialize holders for the local variables that might change below
                #pragma omp atomic read
                Mass_j = P[j].Mass; // this can get modified below, so we need to read it thread-safe now
                
                // quick block of checks to make sure it's actually worth continuing!
                if(Mass_j <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2], 1); // find the closest image in the given box size  //
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) {continue;} // same particle //
                double h2j = PPP[j].Hsml * PPP[j].Hsml;
                if((r2>h2)&&(r2>h2j)) {continue;} // outside kernel (in both 'directions') //
#ifdef FIRE1_SNE_COUPLING
                if(r2>h2) {continue;} // only search 'one way' for particles seen by the BH
#endif
                if(r2 > r2max_phys) {continue;} // outside long-range cutoff //
                kernel.r = sqrt(r2); if(kernel.r <= 0) {continue;}
                
                // calculate kernel quantities //
                #pragma omp atomic read
                rho_j = SphP[j].Density;
                u = kernel.r * kernel.hinv;
                double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j; /* note these lines and many below assume 3D sims! */
                double wk_j = 0, dwk_j = 0, u_j = kernel.r * hinv_j, hinv4_j = hinv_j*hinv3_j, V_j = Mass_j / rho_j;
                if(u<1) {kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 1);} else {kernel.dwk=kernel.wk=0;}
                if(u_j<1) {kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);} else {wk_j=dwk_j=0;}
                if(local.V_i<0 || isnan(local.V_i)) {local.V_i=0;}
                if(V_j<0 || isnan(V_j)) {V_j=0;}
                double sph_area = fabs(local.V_i*local.V_i*kernel.dwk + V_j*V_j*dwk_j); // effective face area //
                wk = 0.5 * (1 - 1/sqrt(1 + sph_area / (M_PI*kernel.r*kernel.r))); // corresponding geometric weight //
#ifdef FIRE1_SNE_COUPLING
                if(u<1) {kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);} else {kernel.wk=kernel.dwk=0;}
                wk = (Mass_j/rho_j) * kernel.wk;
#endif
                if((wk <= 0)||(isnan(wk))) continue; // no point in going further, there's no physical weight here
                double wk_vec[AREA_WEIGHTED_SUM_ELEMENTS]={0}, wk_tmp=0;
                wk_vec[0] = wk;
                wk_tmp=wk*kernel.dp[0]/kernel.r; if(kernel.dp[0]>0) {wk_vec[1]=wk_tmp;} else {wk_vec[2]=wk_tmp;}
                wk_tmp=wk*kernel.dp[1]/kernel.r; if(kernel.dp[1]>0) {wk_vec[3]=wk_tmp;} else {wk_vec[4]=wk_tmp;}
                wk_tmp=wk*kernel.dp[2]/kernel.r; if(kernel.dp[2]>0) {wk_vec[5]=wk_tmp;} else {wk_vec[6]=wk_tmp;}

                // if loop_iteration==-1, this is a pre-calc loop to get the relevant weights for coupling //
                if(loop_iteration < 0)
                {
                    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) out.Area_weighted_sum[k] += wk_vec[k];
                    continue;
                }
                // NOW do the actual feedback calculation //
                double wk_norm = 1. / (MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[0])); // normalization for scalar weight sum
                wk *= wk_norm; // this way wk matches the value summed above for the weighting //
                if((wk <= 0)||(isnan(wk))) continue;

                // ok worth initializing other variables we will use below
                #pragma omp atomic read
                InternalEnergy_j = SphP[j].InternalEnergy; // this can get modified below, so we need to read it thread-safe now
                for(k=0;k<3;k++) {
                    #pragma omp atomic read
                    Vel_j[k] = P[j].Vel[k]; // this can get modified below, so we need to read it thread-safe now
                }
                double InternalEnergy_j_0 = InternalEnergy_j, Mass_j_0 = Mass_j, rho_j_0 = rho_j, Vel_j_0[3]; for(k=0;k<3;k++) {Vel_j_0[k]=Vel_j[k];} // save initial values to use below
#ifdef METALS
                double Metallicity_j[NUM_METAL_SPECIES], Metallicity_j_0[NUM_METAL_SPECIES];
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic read
                    Metallicity_j[k] = P[j].Metallicity[k]; // this can get modified below, so we need to read it thread-safe now
                    Metallicity_j_0[k] = Metallicity_j[k]; // save initial values to  use below
                }
#endif
                
                /* define initial mass and ejecta velocity in this 'cone' */
                double v_bw[3]={0}, e_shock=0, pnorm = 0, pvec[3]={0};
                for(k=0; k<3; k++)
                {
                    double q; q = 0; int i1=2*k+1, i2=i1+1;
                    double q_i1 = fabs(local.Area_weighted_sum[i1]);
                    double q_i2 = fabs(local.Area_weighted_sum[i2]);
#ifdef FIRE1_SNE_COUPLING
                    q_i1=q_i2=1;
#endif
                    if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                    {
                        double rr = q_i2/q_i1;
                        double rr2 = rr * rr;
                        if(wk_vec[i1] != 0)
                        {
                            q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                        } else {
                            q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                        }
                    } else {
                        q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                    }
                    pvec[k] = -q;
                    pnorm += pvec[k]*pvec[k];
                }
                pnorm = sqrt(pnorm);

                wk = pnorm; // this (vector norm) is the new 'weight function' for our purposes
                dM = wk * local.Msne;

                /* now, add contribution from relative star-gas particle motion to shock energy */
                for(k=0;k<3;k++)
                {
                    v_bw[k] = local.SNe_v_ejecta*pvec[k]/pnorm + (local.Vel[k]-Vel_j[k])/All.cf_atime;
                    e_shock += v_bw[k]*v_bw[k];
                }
                double mj_preshock, dM_ejecta_in, massratio_ejecta, mu_j;
                mj_preshock = Mass_j;
                dM_ejecta_in = dM;
                massratio_ejecta = dM_ejecta_in / (dM_ejecta_in + Mass_j);
                mu_j = Mass_j / (dM + Mass_j);
                e_shock *= pnorm * 0.5*local.Msne * mu_j;

                if((wk <= 0)||(isnan(wk))) continue;

                RsneKPC = RsneKPC_0;
                double n0 = rho_j*density_to_n;
                /* this is tedious, but is a fast approximation (essentially a lookup table) for the -0.429 power above */
                if(n0 < 1.e-3) {RsneKPC *= 19.4;} else {
                    if(n0 < 1.e-2) {RsneKPC *= 1.9 + 23./(1.+333.*n0);} else {
                        if(n0 < 1.e-1) {RsneKPC *= 0.7 + 8.4/(1.+33.3*n0);} else {
                            if(n0 < 1) {RsneKPC *= 0.08 + 3.1/(1.+2.5*n0);} else {
                                if(n0 < 10) {RsneKPC *= 0.1 + 1.14/(1.+0.333*n0);} else {
                                    if(n0 < 100) {RsneKPC *= 0.035 + 0.43/(1.+0.0333*n0);} else {
                                        if(n0 < 1000) {RsneKPC *= 0.017 + 0.154/(1.+0.00333*n0);} else {
                                            if(n0 < 1.e4) {RsneKPC *= 0.006 + 0.057/(1.+0.000333*n0);} else {
                                                RsneKPC *= pow(n0, -0.429); }}}}}}}}


                /* below expression is again just as good a fit to the simulations, and much faster to evaluate */
                double z0 = Metallicity_j[0]/All.SolarAbundances[0];
                if(z0 < 0.01) {RsneKPC *= 2.0;} else {
                    if(z0 < 1) {RsneKPC *= 0.93 + 0.0615 / (0.05 + 0.8*z0);} else {RsneKPC *= 0.8 + 0.4 / (1 + z0);}}
                /* calculates cooling radius given density and metallicity in this annulus into which the ejecta propagate */

                /* if coupling radius > R_cooling, account for thermal energy loss in the post-shock medium:
                 from Thornton et al. thermal energy scales as R^(-6.5) for R>R_cool */
                double r_eff_ij = sqrt(r2) - Get_Particle_Size(j);
                if(r_eff_ij > RsneKPC) {e_shock *= RsneKPC*RsneKPC*RsneKPC/(r_eff_ij*r_eff_ij*r_eff_ij);}

                /* now we have the proper energy to couple */
                E_coupled += e_shock;

                /* inject actual mass from mass return */
                int couple_anything_but_scalar_mass_and_metals = 1; // key to indicate whether or not we actually need to do the next set of steps beyond pure scalar mass+metal couplings //
                if(P[j].Hsml<=0) {if(rho_j>0){rho_j*=(1+dM_ejecta_in/Mass_j);} else {rho_j=dM_ejecta_in*kernel.hinv3;}} else {rho_j+=kernel_zero*dM_ejecta_in*hinv3_j;}
                rho_j *= 1 + dM_ejecta_in/Mass_j; // inject mass at constant particle volume //
                Mass_j += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
#if defined(METALS) /* inject metals */
                for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {Metallicity_j[k]=(1-massratio_ejecta)*Metallicity_j[k] + massratio_ejecta*local.yields[k];}
#ifdef GALSF_FB_FIRE_AGE_TRACERS
                if(loop_iteration == 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k] += pnorm*local.yields[k]/Mass_j;}} // add age tracers in taking yields to mean MASS, so we can make it large without actually exchanging large masses
#ifndef GALSF_FB_FIRE_AGE_TRACERS_DISABLE_SURFACE_YIELDS
                if(loop_iteration != 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k]=(1-massratio_ejecta)*Metallicity_j[k] + massratio_ejecta*local.yields[k];}} // treat like any other yield when doing stellar mass exchange
#endif
#endif
#endif
                if(couple_anything_but_scalar_mass_and_metals)
                {
                    
                /* inject the post-shock energy and momentum (convert to specific units as needed first) */
                e_shock *= 1 / Mass_j;
                InternalEnergy_j += e_shock;
                /* inject momentum */
                double m_ej_input = pnorm * local.Msne;
                /* appropriate factor for the ejecta being energy-conserving inside the cooling radius (or Hsml, if thats smaller) */
                double m_cooling = 4.18879*pnorm*rho_j*RsneKPC*RsneKPC*RsneKPC;
                /* apply limiter for energy conservation */
                double mom_boost_fac = 1 + sqrt(DMIN(mj_preshock , m_cooling) / m_ej_input);
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) && (FLAG_NOT_IN_PUBLIC_CODE > 2)) || defined(SINGLE_STAR_SINK_DYNAMICS) 
                if(loop_iteration > 0) {mom_boost_fac=1;} /* no unresolved PdV component for winds+r-process */
#endif
                /* save summation values for outputs */
                dP = local.unit_mom_SNe / Mass_j * pnorm;
                dP_sum += dP;
                dP_boost_sum += dP * mom_boost_fac;

                /* actually do the injection */
                double q0 = All.cf_atime * (pnorm*local.Msne/Mass_j) * mom_boost_fac;
                for(k=0; k<3; k++)
                {
                    double q = q0 * v_bw[k];
                    Vel_j[k] += q;
                }
                    
                } // couple_anything_but_scalar_mass_and_metals
                
                /* we updated variables that need to get assigned to element 'j' -- let's do it here in a thread-safe manner */
                #pragma omp atomic
                P[j].Mass += Mass_j - Mass_j_0; // finite mass update [delta difference added here, allowing for another element to update in the meantime]. done this way to ensure conservation.
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                #pragma omp atomic
                SphP[j].MassTrue += Mass_j - Mass_j_0; // finite mass update
#endif
                if(rho_j_0 > 0) {
                    #pragma omp atomic
                    SphP[j].Density *= rho_j / rho_j_0; // inject mass at constant particle volume [no need to be exactly conservative here] //
                }
                for(k=0;k<3;k++) {
                    #pragma omp atomic
                    P[j].Vel[k] += Vel_j[k] - Vel_j_0[k]; // delta-update
                    #pragma omp atomic
                    SphP[j].VelPred[k] += Vel_j[k] - Vel_j_0[k]; // delta-update
                    #pragma omp atomic
                    P[j].dp[k] += Mass_j*Vel_j[k] - Mass_j_0*Vel_j_0[k]; // discrete momentum change
                }
                #pragma omp atomic
                SphP[j].InternalEnergy += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                #pragma omp atomic
                SphP[j].InternalEnergyPred += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic
                    P[j].Metallicity[k] += Metallicity_j[k] - Metallicity_j_0[k]; // delta-update
                }
                
                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}}}    /* open it */
    } // while(startnode >= 0)

    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addFB(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
    return 0;
} // int addFB_evaluate



#else // un-protected [updated, more fixed energy-injecting SNe scheme]


int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double u,r2,kernel_zero,wk,dM_ejecta_in,dP,E_coupled,dP_sum,dP_boost_sum;
    struct kernel_addFB kernel; struct addFB_evaluate_data_in_ local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));

    /* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_addFB(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(local.Msne<=0) {return 0;} // no SNe for the origin particle! nothing to do here //
    if(local.Hsml<=0) {return 0;} // zero-extent kernel, no particles //

    // some units (just used below, but handy to define for clarity) //
    double h2 = local.Hsml*local.Hsml; kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1); wk=0; // define the kernel zero-point value, needed to prevent some nasty behavior when no neighbors found
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4); // define kernel quantities
    double unitlength_in_kpc= UNIT_LENGTH_IN_KPC * All.cf_atime, density_to_n=All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, unit_egy_SNe = 1.0e51/UNIT_ENERGY_IN_CGS;

    // now define quantities that will be used below //
    double psi_cool=1, psi_egycon=1, v_ejecta_eff_init=local.SNe_v_ejecta, v_ejecta_eff=v_ejecta_eff_init, residual_thermal_frac=0; // separate initial [pre-shock] ejecta velocity, which defines energy, and post-shock
    double wk_norm = 1. / (MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[0])); // normalization for scalar weight sum
    double pnorm_sum = 1./(MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[10])); // re-normalization after second pass for normalized "pnorm" (should be close to ~1)
    double sedov_phase_thermal_to_kinetic_ratio = 2.54, f_sedov_kin; // define ratio of thermal-to-kinetic energy for sedov-taylor phase
    int retain_thermal_flag = 1; // flag to determine whether additional post-shock thermal energy should be considered to have radiated away already
    int feedback_type_is_SNe = 0; // variable to know if this is a sne or continuous wind or other form of mechanical feedback (discrete injection treated differently from continuous in some ways below)
    if(loop_iteration == 0) {feedback_type_is_SNe = 1;} // assume, for now, that loop 0 represents SNe, for purposes of energy-momentum switch below //
    if(feedback_type_is_SNe == 0) {sedov_phase_thermal_to_kinetic_ratio = 1.e-2;} // negligible excess thermal component for continuous sources
    f_sedov_kin=1./(1.+sedov_phase_thermal_to_kinetic_ratio); // fraction of total energy in kinetic after reverse shock, to assume for sedov-taylor phase
    if((local.Area_weighted_sum[0] > MIN_REAL_NUMBER) && (loop_iteration >= 0))
    {
        double vba_2_eff = pnorm_sum * local.Area_weighted_sum[7]; // phi term for energy: weighted mass-deposited KE for ejecta neighbors
        v_ejecta_eff_init = sqrt(local.SNe_v_ejecta*local.SNe_v_ejecta + vba_2_eff); // account for all terms to get the revised KE term here
        v_ejecta_eff = v_ejecta_eff_init * sqrt(f_sedov_kin); // effective velocity after reverse shock in sedov phase
        double beta_egycon = sqrt(pnorm_sum / local.Msne) * (1./v_ejecta_eff) * local.Area_weighted_sum[8]; // beta term for re-normalization for energy [can be positive or negative]
        double beta_cool = pnorm_sum * local.Area_weighted_sum[9]; // beta term if all particles in terminal-momentum-limit
        psi_egycon = sqrt(1. + beta_egycon*beta_egycon) - beta_egycon; // exact solution for energy equation for constant psi
        if(beta_egycon > 20.) {psi_egycon = 1./(2.*beta_egycon);} // replace with series expansion to avoid roundoff error at high beta
        //if(beta_egycon < 0) {psi_egycon=1; residual_thermal_frac=DMAX(0., sqrt(1.+beta_egycon*beta_egycon)-beta_egycon-1.) * f_sedov_kin;} // in this case (blastwave in a converging flow) we don't boost the momentum beyond the 'normal' maximum but assign the residual energy to thermal, since the timescale to convert this additional post-shock thermal energy to kinetic is actually quite long for this situation
        if(beta_egycon < 0) {psi_egycon=1; residual_thermal_frac=DMAX(0., -2.*beta_egycon) * f_sedov_kin;} // in this case (blastwave in a converging flow) we don't boost the momentum beyond the 'normal' maximum but assign the residual energy to thermal, since the timescale to convert this additional post-shock thermal energy to kinetic is actually quite long for this situation
        if(beta_cool > 0.5) {psi_cool = 1./(2.*beta_cool);} // for cooling limit, only need upper limit to psi, all else will use less energy
    }
    psi_egycon = DMIN(1 , psi_egycon); // this should be gauranteed by the above checks, but enforce it regardless


    double Energy_injected_codeunits = 0.5 * local.Msne * v_ejecta_eff_init * v_ejecta_eff_init;
    double Esne51 = Energy_injected_codeunits / unit_egy_SNe;
    double RsneKPC = 0., RsneKPC_3 = 0., m_cooling = 0., v_cooling = 210./UNIT_VEL_IN_KMS;
    double RsneKPC_0 = (0.0284/unitlength_in_kpc);
    if(feedback_type_is_SNe == 1) // check for SNe specifically
    {
        RsneKPC_0 *= pow(1+Esne51,0.286); //SNe: using scaling from Cioffi with weak external pressure
    } else {
        RsneKPC_0 *= pow(Esne51,0.286); // ensures smooth conservation for winds and tracers as mass-loading goes to vanishingly small values
    }
    double r2max_phys = 2.0/unitlength_in_kpc; // no super-long-range effects allowed! (of course this is arbitrary in code units) //
    if(local.Hsml >= r2max_phys) {psi_egycon=DMIN(psi_egycon,1); psi_cool=DMIN(psi_cool,1);}
    r2max_phys *= r2max_phys;

    /* Now start the actual FB computation for this particle */
    if(mode == 0) {startnode = All.MaxPart;} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;}    /* start at root node, open it */
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}

            E_coupled = dP_sum = dP_boost_sum = 0;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Type != 0) {continue;} // require a gas particle //
                
                double Mass_j, InternalEnergy_j, rho_j, Vel_j[3]; // initialize holders for the local variables that might change below
                #pragma omp atomic read
                Mass_j = P[j].Mass; // this can get modified below, so we need to read it thread-safe now

                // now consider a block of conditions we will use to evaluate whether its worth opening this loop at all //
                if(Mass_j <= 0) {continue;} // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) {continue;} // same particle //
                double h2j = PPP[j].Hsml * PPP[j].Hsml;
                if((r2>h2)&&(r2>h2j)) {continue;} // outside kernel (in both 'directions') //
                if(r2 > r2max_phys) {continue;} // outside long-range cutoff //
                kernel.r = sqrt(r2); if(kernel.r <= 0) {continue;}

                
                // calculate kernel quantities //
                #pragma omp atomic read
                rho_j = SphP[j].Density;
                u = kernel.r * kernel.hinv;
                double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j;
                double wk_j = 0, dwk_j = 0, u_j = kernel.r * hinv_j, hinv4_j = hinv_j*hinv3_j, V_j = Mass_j / rho_j;
                if(u<1) {kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 1);} else {kernel.wk=kernel.dwk=0;}
                if(u_j<1) {kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);} else {wk_j=dwk_j=0;}
                if(local.V_i<0 || isnan(local.V_i)) {local.V_i=0;}
                if(V_j<0 || isnan(V_j)) {V_j=0;}
                double sph_area = fabs(local.V_i*local.V_i*kernel.dwk + V_j*V_j*dwk_j); // effective face area //
                wk = 0.5 * (1 - 1/sqrt(1 + sph_area / (M_PI*kernel.r*kernel.r))); // corresponding geometric weight //

                if((wk <= 0)||(isnan(wk))) {continue;} // no point in going further, there's no physical weight here

                double wk_vec[AREA_WEIGHTED_SUM_ELEMENTS] = {0};
                wk_vec[0] = wk;
                if(kernel.dp[0]>0) {wk_vec[1]=wk*kernel.dp[0]/kernel.r; wk_vec[2]=0;} else {wk_vec[1]=0; wk_vec[2]=wk*kernel.dp[0]/kernel.r;}
                if(kernel.dp[1]>0) {wk_vec[3]=wk*kernel.dp[1]/kernel.r; wk_vec[4]=0;} else {wk_vec[3]=0; wk_vec[4]=wk*kernel.dp[1]/kernel.r;}
                if(kernel.dp[2]>0) {wk_vec[5]=wk*kernel.dp[2]/kernel.r; wk_vec[6]=0;} else {wk_vec[5]=0; wk_vec[6]=wk*kernel.dp[2]/kernel.r;}
                if(loop_iteration==-2) {
                    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {out.Area_weighted_sum[k] += wk_vec[k];} // normal summation on the first loop
                    continue;
                }
                // ok worth initializing other variables we will use below
                #pragma omp atomic read
                InternalEnergy_j = SphP[j].InternalEnergy; // this can get modified below, so we need to read it thread-safe now
                for(k=0;k<3;k++) {
                    #pragma omp atomic read
                    Vel_j[k] = P[j].Vel[k]; // this can get modified below, so we need to read it thread-safe now
                }
                double InternalEnergy_j_0,Mass_j_0,rho_j_0,Vel_j_0[3]; InternalEnergy_j_0=InternalEnergy_j; Mass_j_0=Mass_j; rho_j_0=rho_j; for(k=0;k<3;k++) {Vel_j_0[k]=Vel_j[k];} // save initial values to use below
#ifdef METALS
                double Metallicity_j[NUM_METAL_SPECIES], Metallicity_j_0[NUM_METAL_SPECIES];
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic read
                    Metallicity_j[k] = P[j].Metallicity[k]; // this can get modified below, so we need to read it thread-safe now
                    Metallicity_j_0[k] = Metallicity_j[k]; // save initial values to  use below
                }
#endif
                
                RsneKPC = RsneKPC_0;
                /* calculate cooling radius given density and metallicity in this annulus into which the ejecta propagate */
                if(loop_iteration < 2)
                {
                    double e0 = Esne51;
                    if(loop_iteration < 0) {e0=1;} else {if(feedback_type_is_SNe == 1) {e0+=1;} else {e0=0.1;}} // set to small number for non-SNe feedback
                    double n0 = rho_j*density_to_n; if(n0 < 0.001) {n0=0.001;}
                    double z0 = Metallicity_j[0]/All.SolarAbundances[0];
                    if(z0 < 0.01) {z0 = 0.01;}
                    double z0_term=1.; if(z0 < 1.) {z0_term = z0*sqrt(z0);} else {z0_term = z0;}
                    double nz_dep  = pow(n0 * z0_term , 0.14);
                    v_cooling = 210. * DMAX(nz_dep,0.5) / UNIT_VEL_IN_KMS;
                    m_cooling = 4.56e36 * e0 / (nz_dep*nz_dep * UNIT_MASS_IN_CGS);
                    if(loop_iteration >= 0 && feedback_type_is_SNe == 0) {v_cooling *= 1.e10; m_cooling *= 1.e10;} // for non-SNe, ignore finite cooling radii and directly couple; wont matter unless choose to include boost term below, with fixes we've added
                    RsneKPC = pow( 0.238732 * m_cooling/rho_j , 1./3. );
                }
                RsneKPC_3 = RsneKPC*RsneKPC*RsneKPC;
                // if loop_iteration==-1, this is a pre-calc loop to get the relevant weights for coupling //
                if(loop_iteration < 0)
                {
                    if(loop_iteration==-1) // the Area_weighted_sum quantities are computed on loop=-2; these quantities must be computed on loop=-1 (after Area_weighted_sums are computed)
                    {
                        /* calculate the corrected momentum vectors that we will actually use in the coupling proper */
                        double pnorm=0, pvec[3]={0}, vel_ba_2=0, cos_vel_ba_pcoupled=0;
                        for(k=0;k<3;k++)
                        {
                            double q = 0; int i1=2*k+1, i2=i1+1;
                            double q_i1 = fabs(local.Area_weighted_sum[i1]);
                            double q_i2 = fabs(local.Area_weighted_sum[i2]);
                            if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                            {
                                double rr = q_i2/q_i1;
                                double rr2 = rr * rr;
                                if(wk_vec[i1] != 0)
                                {
                                    q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                                } else {
                                    q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                                }
                            } else {
                                q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                            }
                            pvec[k] = -q;
                            pnorm += pvec[k]*pvec[k];
                        }
                        pnorm = sqrt(pnorm);
                        /* now calculate the additional weights that are needed for energy terms */
                        for(k=0;k<3;k++)
                        {
                            double v_ba = (Vel_j[k] - local.Vel[k]) / All.cf_atime; // relative gas-star velocity //
                            vel_ba_2 += v_ba*v_ba; // magnitude of velocity vector (for corrected post-shock energies to distribute)
                            cos_vel_ba_pcoupled += v_ba * pvec[k]/pnorm; // direction of ejecta [after correction loop]
                        }
                        wk_vec[7] = pnorm * vel_ba_2; // phi_0 term : residual KE term from mass-coupling for {small, second-order} energy correction
                        wk_vec[8] = sqrt(pnorm * Mass_j) * cos_vel_ba_pcoupled; // beta_0 term : cross-term for momentum coupling effect on energy-coupling
                        wk_vec[9] = pnorm * cos_vel_ba_pcoupled / v_cooling; // calculate the beta term as if all particles hit terminal: more accurate result in that limit
                        wk_vec[10] = pnorm; // normalization (so that we can divide by its sum to properly normalize the beta_egy and beta_cool quantities)
                        for(k=7;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {out.Area_weighted_sum[k] += wk_vec[k];} // pass these, which have not yet been calculated
                    }
                    continue;
                }
                // NOW do the actual feedback calculation //
                wk *= wk_norm; // this way wk matches the value summed above for the weighting //

                if((wk <= 0)||(isnan(wk))) continue;
            
                /* define initial mass and ejecta velocity in this 'cone' */
                double pnorm = 0, pvec[3] = {0};
                for(k=0; k<3; k++)
                {
                    double q = 0; int i1=2*k+1, i2=i1+1;
                    double q_i1 = fabs(local.Area_weighted_sum[i1]);
                    double q_i2 = fabs(local.Area_weighted_sum[i2]);
                    if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                    {
                        double rr = q_i2/q_i1;
                        double rr2 = rr * rr;
                        if(wk_vec[i1] != 0)
                        {
                            q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                        } else {
                            q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                        }
                    } else {
                        q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                    }
                    pvec[k] = -q;
                    pnorm += pvec[k]*pvec[k];
                }
                pnorm = sqrt(pnorm); // this (vector norm) is the new 'weight function' for our purposes
                pnorm *= pnorm_sum; for(k=0;k<3;k++) {pvec[k] *= pnorm_sum;} // normalize following sum [10] to ensure faces sum to unity
                dM_ejecta_in = pnorm * local.Msne;
                double mj_preshock, massratio_ejecta;
                mj_preshock = Mass_j;
                massratio_ejecta = dM_ejecta_in / (dM_ejecta_in + Mass_j);

                /* inject actual mass from mass return */
                int couple_anything_but_scalar_mass_and_metals = 1; // key to indicate whether or not we actually need to do the next set of steps beyond pure scalar mass+metal couplings //
                if(P[j].Hsml<=0) {if(rho_j>0){rho_j*=(1+dM_ejecta_in/Mass_j);} else {rho_j=dM_ejecta_in*kernel.hinv3;}} else {rho_j+=kernel_zero*dM_ejecta_in*hinv3_j;}
                rho_j *= 1 + dM_ejecta_in/Mass_j; // inject mass at constant particle volume //
                Mass_j += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
                
#ifdef METALS   /* inject metals */
                for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {Metallicity_j[k]=(1-massratio_ejecta)*Metallicity_j[k] + massratio_ejecta*local.yields[k];}
#ifdef GALSF_FB_FIRE_AGE_TRACERS
                if(loop_iteration == 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k] += pnorm*local.yields[k]/Mass_j;}} // add age tracers in taking yields to mean MASS, so we can make it large without actually exchanging large masses
#ifndef GALSF_FB_FIRE_AGE_TRACERS_DISABLE_SURFACE_YIELDS
                if(loop_iteration != 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k]=(1-massratio_ejecta)*Metallicity_j[k] + massratio_ejecta*local.yields[k];}} // treat like any other yield when doing stellar mass exchange
#endif
#endif
#endif
                
                double KE_initial=0, KE_final=0;
                if(couple_anything_but_scalar_mass_and_metals)
                {
                /* inject momentum: account for ejecta being energy-conserving inside the cooling radius (or Hsml, if thats smaller) */
                double wk_m_cooling = pnorm * m_cooling; // effective cooling mass for this particle
                double boost_terminal = sqrt(1 + wk_m_cooling / dM_ejecta_in); // terminal momentum boost-factor
                double boost_egycon = sqrt(1 + mj_preshock / dM_ejecta_in); // energy-conserving limit for coupling through neighbors
                double mom_boost_fac = 1; // set default value of momentum term and dummy variable to check which [energy or momentum-conserving] solution to apply
                if(feedback_type_is_SNe == 1) // question here is whether wind-type feedback should assume kinetic energy conserved here
                {
                    boost_terminal *= psi_cool; // appropriately re-weight boost to avoid energy conservation errors [cooling-limit]
                    boost_egycon *= psi_egycon; // appropriately re-weight boost to avoid energy conservation errors [energy-conserving-limit]
                    mom_boost_fac = boost_egycon; // default to energy-conserving solution, then modify if needed in lines below
                    double dv_dp_phys = 0; for(k=0;k<3;k++) {dv_dp_phys += (1-massratio_ejecta) * (kernel.dp[k]/kernel.r) * ((local.Vel[k] - Vel_j[k])/All.cf_atime);} // recession velocity of cell from SNe
                    double vcool_eff = v_cooling / psi_cool; // effective shell speed when the cooling radius is reached
#if 0
                    double dv_eff = vcool_eff + 2.*dv_dp_phys; // effective relative velocity for determining if you can reach that shell speed. i.e. when recession velocity equals nominal cooling mass for a real solution, assuming you cool when the post-shock temperature reaches a Tcool that corresponds to the post-shock velocity for some outward vcool, giving e.g. half the desired cooling mass is obtained when you have outward v = vcool
                    double mcool_mod = wk_m_cooling; // this should be our starting point for applying the cooling solution, unless we include the compressive term, but don't want to over-shoot there
                    if(dv_eff > 0) {mcool_mod = wk_m_cooling; if(dv_dp_phys > 0) {mcool_mod *= vcool_eff / (1.e-20*vcool_eff + dv_eff);}} // use the above recession velocity information to modify the cooling mass compared to the total cell mass to determine which solution to use
                    if(mcool_mod < mj_preshock) {mom_boost_fac = DMIN(boost_terminal,boost_egycon); residual_thermal_frac=0; retain_thermal_flag=0;} // if swept mass where reach the terminal solution is less than the cell mass, apply it, otherwise apply the conservative solution
#else
                    double dv_eff = vcool_eff + dv_dp_phys; // effective relative velocity for determining if you can reach that shell speed. i.e. when recession velocity equals nominal cooling mass for a real solution, assuming you cool when the post-shock temperature reaches a Tcool that corresponds to the post-shock velocity for some outward vcool, giving e.g. half the desired cooling mass is obtained when you have outward v = vcool
                    double m_terminal=wk_m_cooling, m_cooling=1.e10*mj_preshock; // this should be our starting point for applying the cooling solution, unless we include the compressive term, but don't want to over-shoot there
                    if(dv_eff > 1.e-10*vcool_eff) {double fac_ve=vcool_eff/(MIN_REAL_NUMBER+1.e-10*vcool_eff + dv_eff); m_cooling*=fac_ve; if(fac_ve<1.) {m_terminal*=fac_ve;}} // use the above recession velocity information to modify the cooling mass compared to the total cell mass to determine which solution to use
                    if(mj_preshock > m_cooling) {residual_thermal_frac=0; retain_thermal_flag=0;} // figure out if we've passed the cooling time/distance (which is -not- the same as reaching the terminal momentum for an arbitrary inflow structure around the explosion)
                    if(mj_preshock > m_terminal) {mom_boost_fac = DMIN(boost_terminal,boost_egycon);} // if swept mass where reach the terminal solution is less than the cell mass, apply it, otherwise apply the conservative solution
#endif
                } else {mom_boost_fac=DMIN(1,boost_egycon*psi_egycon); residual_thermal_frac=0;} // prevent energy conservation issues when coupling mass-loss
                
                /* save summation values for outputs */
                dP = local.unit_mom_SNe / Mass_j * pnorm;
                dP_sum += dP; dP_boost_sum += dP * mom_boost_fac;

                /* actually do the injection */
                double mom_prefactor =  mom_boost_fac * massratio_ejecta * (All.cf_atime*v_ejecta_eff) / pnorm; // this gives the appropriately-normalized tap-able momentum from the energy-conserving solution
                if(mom_prefactor > 0)
                {
                    for(k=0; k<3; k++)
                    {
                        double d_vel = mom_prefactor * pvec[k] + massratio_ejecta*(local.Vel[k] - Vel_j[k]); // local.Vel term from extra momentum of moving star, Vel_j term from going from momentum to velocity boost with added mass
                        KE_initial += Vel_j_0[k]*Vel_j_0[k]; Vel_j[k] += d_vel; KE_final += Vel_j[k]*Vel_j[k]; // calculate initial and final kinetic energies (here v^2, multiply by mass below)
                    }
                    KE_initial *= 0.5*mj_preshock*All.cf_a2inv; KE_final *= 0.5*Mass_j*All.cf_a2inv;
                }
                double E_sne_initial = pnorm * Energy_injected_codeunits;
                double d_Egy_internal = KE_initial + E_sne_initial - KE_final; /* now calculate the residual energy with option to add it as thermal */
#if !defined(SINGLE_STAR_FB_WINDS) /* (for single-star modules we ignore this b/c assume always trying to resolve R_cool) */
                if(feedback_type_is_SNe == 1) /* if coupling radius > R_cooling, account for thermal energy loss in the post-shock medium: from Thornton et al. thermal energy scales as R^(-6.5) for R>R_cool. only use for SNe b/c scalings [like momentum] only apply there. over-cooling if code wants to do it will easily occur next timestep. */
                {
                    d_Egy_internal = (1.-f_sedov_kin+residual_thermal_frac)*E_sne_initial; // the thermal energy component is constant (proportional to the injected area) and determined -as part of the energy-conserving solution- here, you should not artificially decrease it or renormalize it cell-by-cell for this form of the solutions
                } else {d_Egy_internal = DMAX(DMIN(d_Egy_internal , 2.*E_sne_initial),0);}
#endif
                if(retain_thermal_flag==0) {d_Egy_internal=0;} // use flag to determined if we should retain this residual thermal energy for this stage
                d_Egy_internal /= Mass_j; // convert to specific internal energy, finally //
#ifndef MECHANICAL_FB_MOMENTUM_ONLY
                if(d_Egy_internal > 0) {InternalEnergy_j += d_Egy_internal; E_coupled += d_Egy_internal;}
#endif
                    
                } // couple_anything_but_scalar_mass_and_metals
                    
                /* we updated variables that need to get assigned to element 'j' -- let's do it here in a thread-safe manner */
                #pragma omp atomic
                N_Gas_Couplings_ThisTask++; // note that a cell recieved some feedback
                #pragma omp atomic
                LocalGasMechFBInfoTemp[j].N_injected++; // flag that -this- cell received feedback
                #pragma omp atomic
                LocalGasMechFBInfoTemp[j].m_injected += Mass_j - Mass_j_0; // finite mass update [delta difference added here, allowing for another element to update in the meantime]. done this way to ensure conservation.
                #pragma omp atomic
                LocalGasMechFBInfoTemp[j].TE_injected += Mass_j*InternalEnergy_j - Mass_j_0*InternalEnergy_j_0; // delta-update of conserved quantity (total internal energy)
                #pragma omp atomic
                LocalGasMechFBInfoTemp[j].KE_injected += KE_final - KE_initial; // delta-update of conserved quantity (total kinetic energy)
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic
                    LocalGasMechFBInfoTemp[j].Z_injected[k] += Mass_j*Metallicity_j[k] - Mass_j_0*Metallicity_j_0[k]; // delta-update of conserved quantity (total metal mass)
                }
                for(k=0;k<3;k++) {
                    #pragma omp atomic
                    LocalGasMechFBInfoTemp[j].p_injected[k] += (Mass_j*Vel_j[k] - Mass_j_0*Vel_j_0[k]) / All.cf_atime; // delta-update of conserved quantity (total momentum), converted to physical units
                }
                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}}}    /* open it */
    } // while(startnode >= 0)

    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addFB(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
    return 0;
} // int addFB_evaluate

#endif // GALSF_USE_SNE_ONELOOP_SCHEME else



/* subroutine to check for total kinetic energy change and thermal energy change after integrating the effects of all SNe over all cells, and coupling this to particles,  */
void verify_and_assign_local_mechfb_integrals(void)
{
    if(N_Gas_Couplings_ThisTask <= 0) {return;} /* no cells had feedback deposited */
    int j,k,ndone=0; for(j=0;j<N_gas;j++)
    {
        if(LocalGasMechFBInfoTemp[j].N_injected <= 0) {continue;} /* all mechanisms deposit non-zero mass, so skip if this is not >0*/
        if(P[j].Type==0 && P[j].Mass>0)
        {
            double m0=P[j].Mass, dm=LocalGasMechFBInfoTemp[j].m_injected; P[j].Mass += dm; /* update mass */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            m0=SphP[j].MassTrue; SphP[j].MassTrue += dm; /* update conserved mass */
#endif
            double mf=m0+dm; /* save for below */
            for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k] = (m0/mf)*P[j].Metallicity[k] + (1./mf)*LocalGasMechFBInfoTemp[j].Z_injected[k];} /* update metallicity */
            SphP[j].Density *= mf/m0; /* update density [semi-drift] */
            double dTE=LocalGasMechFBInfoTemp[j].TE_injected;
            if(dTE > 0)
            {
                double TE_0=m0*SphP[j].InternalEnergy; dTE=DMAX(-TE_0,dTE); /* ensure against non-negative values */
                double dU = (-dm/mf)*SphP[j].InternalEnergy + (1./mf)*dTE; /* using new mass get updated internal energy */
                double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j), implied_heating_cgs=(dU*UNIT_SPECEGY_IN_CGS*PROTONMASS_CGS)/(dt*UNIT_TIME_IN_CGS), typical_cooling_cgs=1.e-23*(SphP[j].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS);
                if((implied_heating_cgs < 0.3*typical_cooling_cgs) && (dt > MIN_REAL_NUMBER) && ((dU < 4.*SphP[j].InternalEnergy) || ((dU < 1000.*SphP[j].InternalEnergy) && ((dU+SphP[j].InternalEnergy)*U_TO_TEMP_UNITS*2./3.*1.28 < 5.e5)))) {SphP[j].DtInternalEnergy += dU/dt;} else {SphP[j].InternalEnergy += dU; SphP[j].InternalEnergyPred += dU;}
                //SphP[j].InternalEnergy += dU; SphP[j].InternalEnergyPred += dU; /* update internal energy; simpler (old) way to do it - less accurate phase diagrams at high density, however */
            }
            double dKE=LocalGasMechFBInfoTemp[j].KE_injected, dp[3];
            if(dKE != 0 || LocalGasMechFBInfoTemp[j].p_injected[0] != 0 || LocalGasMechFBInfoTemp[j].p_injected[1] != 0 || LocalGasMechFBInfoTemp[j].p_injected[2] != 0 )
            {
                double KE_0=0, p0[3], p2=0, dp2=0, pdp=0; /* define variables for momentum update (the non-trivial update term) */
                for(k=0;k<3;k++) {p0[k]=m0*P[j].Vel[k]/All.cf_atime; dp[k]=LocalGasMechFBInfoTemp[j].p_injected[k]; p2+=p0[k]*p0[k]; dp2+=dp[k]*dp[k]; pdp+=p0[k]*dp[k];} /* define variables */
                KE_0=p2/(2.*m0); dKE=DMAX(-KE_0,dKE); /* limit to retain positive definite kinetic energy */
                double a0=dp2, b0=2.*pdp, c0=p2*dm/m0+2.*mf*dKE, sfac=b0*b0+4.*a0*c0, f0=0; /* assume coupled delta_p = f0*delta_p, and solve for the multiplier f0 that gives the desired total kinetic energy change */
                if(sfac<=0 || !isfinite(sfac)) {if((fabs(a0)<1.e-40) || !isfinite(b0/a0)) {f0=0;} else {f0=-b0/(2.*a0);}} // catches for floating-point error
                    else {if((b0>0) && fabs(4.*a0*c0)<1.e-3*b0*b0) {f0=c0/b0;} // catches for floating-point error
                    else {if(fabs(a0)<1.e-40 || !isfinite(a0)) {f0=0;} else {f0=(-b0+sqrt(b0*b0+4.*a0*c0))/(2.*a0);}}} // catches for floating-point error, if pass all of them, use exact solution for desired KE
                if(f0<0 || !isfinite(f0)) {f0=0;} else {if(f0>1.) {f0=1.;}} /* limit to physical values (should never be an issue but again because of float error it could be) */
                for(k=0;k<3;k++) {
                    double dv = (-dm/mf)*P[j].Vel[k] + f0*(1./mf)*dp[k]*All.cf_atime; /* calculate total momentum change and mass change and therefore final velocity (in code units) */
                    P[j].Vel[k] += dv; SphP[j].VelPred[k] += dv; P[j].dp[k] += f0*dp[k]; /* update velocities */
                }
            }
            ndone++; /* note another cell accounted for */
            if(ndone >= N_Gas_Couplings_ThisTask) {break;} /* we have done all cells (note is possible if the same cell is hit many times, N_Gas_Couplings_ThisTask can be much larger than ndone after the full loop. but if its smaller, then we -must- be done with this loop */
        }
    }
    return;
}



/* parent routine which calls all of the mechanical fb loops and verifies coupling, if needed */
void mechanical_fb_calc_toplevel(void)
{
    PRINT_STATUS("Start mechanical feedback computation...");
#ifndef GALSF_USE_SNE_ONELOOP_SCHEME
    /* allocate temporary stucture which will hold the total change, to compare when done to check for non-linear effects if too many cells act at once */
    LocalGasMechFBInfoTemp = (struct temporary_mech_fb_data_tohold *) mymalloc("LocalGasMechFBInfoTemp",N_gas * sizeof(struct temporary_mech_fb_data_tohold)); /* allocate */
    N_Gas_Couplings_ThisTask = 0; /* initialize this to zero [default to assume no coupled feedback] */
    int i; for(i=0;i<N_gas;i++) {if(P[i].Type==0) {memset(&LocalGasMechFBInfoTemp[i], 0, sizeof(struct temporary_mech_fb_data_tohold));}} /* zero it out before loops */
    mechanical_fb_calc(-2); /* compute weights for coupling [first weight-calculation pass] */
#endif
    mechanical_fb_calc(-1); /* compute weights for coupling [second weight-calculation pass] */
    mechanical_fb_calc(0); /* actually do the mechanical feedback coupling */
#ifndef GALSF_USE_SNE_ONELOOP_SCHEME
    verify_and_assign_local_mechfb_integrals();
    myfree(LocalGasMechFBInfoTemp); /* free the structure */
#endif
}


/* sub-routine which calls each of the relevant loops */
void mechanical_fb_calc(int fb_loop_iteration)
{
    PRINT_STATUS(" ..mechanical feedback loop: iteration %d",fb_loop_iteration);
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    loop_iteration = fb_loop_iteration; /* sets the appropriate feedback type for the calls below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    CPU_Step[CPU_SNIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */


#endif /* GALSF_FB_MECHANICAL */
