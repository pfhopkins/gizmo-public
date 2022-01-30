#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>

#include "./allvars.h"
#include "./proto.h"
#include "./kernel.h"
#ifdef BH_WIND_SPAWN
#define MASS_THRESHOLD_FOR_WINDPROMO (DMAX(5.*All.BAL_wind_particle_mass,0.25*All.MaxMassForParticleSplit))
#endif /* define a mass threshold for this model above which a 'hyper-element' has accreted enough to be treated as 'normal' */


/*! This file contains the operations needed for merging/splitting gas particles/cells on-the-fly in the simulations.
    If more complicated routines, etc. are to be added to determine when (and how) splitting/merging occurs, they should also be
    added here. The split routine should also be the template for spawning new gas particles (collisionless particles are spawned
    much more easily; for those, see the star formation routines). */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/*! Here we can insert any desired criteria for particle mergers: by default, this will occur
    when particles fall below some minimum mass threshold */
int does_particle_need_to_be_merged(int i)
{
    if(P[i].Mass <= 0) return 0;
#ifdef PREVENT_PARTICLE_MERGE_SPLIT
    return 0;
#else
#ifdef GRAIN_RDI_TESTPROBLEM
    return 0;
#endif
#ifdef BH_WIND_SPAWN
    if(P[i].ID == All.AGNWindID)
    {
#ifdef BH_DEBUG_SPAWN_JET_TEST
        MyFloat vr2 = (P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2]) * All.cf_a2inv; // physical
        if(vr2 <= 0.01 * All.BAL_v_outflow*All.BAL_v_outflow) {return 1;} else {return 0;} // merge only if velocity condition satisfied, even if surrounded by more massive particles //
#else
        if(P[i].Mass < (All.MaxMassForParticleSplit*target_mass_renormalization_factor_for_mergesplit(i))) {return 1;}
        if(P[i].Mass >= MASS_THRESHOLD_FOR_WINDPROMO*target_mass_renormalization_factor_for_mergesplit(i)) {return 1;}
#endif
    }
#endif
#ifdef PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT
    if(P[i].Type==0)
    {
        double lambda_J = Get_Gas_Fast_MHD_wavespeed_i(i) * sqrt(M_PI / (All.G * SphP[i].Density * All.cf_a3inv));
        if((lambda_J > 4. * PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT * Get_Particle_Size(i)*All.cf_atime) && (P[i].Mass < All.MaxMassForParticleSplit)) {return 1;} // de-refine
    }
#endif
    if((P[i].Type>0) && (P[i].Mass > 0.5*All.MinMassForParticleMerger*target_mass_renormalization_factor_for_mergesplit(i))) {return 0;}
    if(P[i].Mass <= (All.MinMassForParticleMerger*target_mass_renormalization_factor_for_mergesplit(i))) {return 1;}
    return 0;
#endif
}


/*! Here we can insert any desired criteria for particle splitting: by default, this will occur
    when particles become too massive, but it could also be done when Hsml gets very large, densities are high, etc */
int does_particle_need_to_be_split(int i)
{
    if(P[i].Type != 0) {return 0;} // default behavior: only gas particles split //
#ifdef PREVENT_PARTICLE_MERGE_SPLIT
    return 0;
#else
#ifdef BH_DEBUG_SPAWN_JET_TEST
    if(P[i].ID == All.AGNWindID) {return 0;}
#endif
    if(P[i].Mass >= (All.MaxMassForParticleSplit*target_mass_renormalization_factor_for_mergesplit(i))) {return 1;}
#ifdef PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT
    if(P[i].Type == 0)
    {
        double lambda_J = Get_Gas_Fast_MHD_wavespeed_i(i) * sqrt(M_PI / (All.G * SphP[i].Density * All.cf_a3inv));
        if((lambda_J < PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT * Get_Particle_Size(i)*All.cf_atime) && (P[i].Mass > 2*All.MinMassForParticleMerger)) {return 1;} // refine
    }
#endif
    return 0;
#endif
}

/*! A multiplicative factor that determines the target mass of a particle for the (de)refinement routines */
double target_mass_renormalization_factor_for_mergesplit(int i)
{
    double ref_factor=1.0;
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
    if(P[i].Type==0)
    {
        return 1; // need to determine appropriate desired refinement criterion, if resolution is not strictly pre-defined //
    }
#endif
/*!
 #if defined(BH_CALC_DISTANCES) && !defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE) && !defined(SINGLE_STAR_SINK_DYNAMICS)
    ref_factor = DMIN(1.,sqrt(P[i].min_dist_to_bh + 0.0001)); // this is an example of the kind of routine you could use to scale resolution with BH distance //
#endif
 */
    return ref_factor;
}


/*! This is the parent routine to actually determine if mergers/splits need to be performed, and if so, to do them
  modified by Takashi Okamoto (t.t.okamoto@gmail.com) on 20/6/2019
 */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here!!! -- */
void merge_and_split_particles(void)
{
    struct flags_merg_split {
        int flag; // 0 for nothing, -1 for clipping, 1 for merging, 2 for splitting, and 3 marked as merged
        int target_index;
    } *Ptmp;

    int target_for_merger,dummy=0,numngb_inbox,startnode,i,j,n; double threshold_val;
    int n_particles_merged,n_particles_split,n_particles_gas_split,MPI_n_particles_merged,MPI_n_particles_split,MPI_n_particles_gas_split;
    Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
    Gas_split=0; n_particles_merged=0; n_particles_split=0; n_particles_gas_split=0; MPI_n_particles_merged=0; MPI_n_particles_split=0; MPI_n_particles_gas_split=0;
    Ptmp = (struct flags_merg_split *) mymalloc("Ptmp", NumPart * sizeof(struct flags_merg_split));

    // TO: need initialization
    for (i = 0; i < NumPart; i++) {
      Ptmp[i].flag = 0;
      Ptmp[i].target_index = -1;
    }

    for (i = 0; i < NumPart; i++)
    {
        int Pi_BITFLAG = (1 << (int)P[i].Type); // bitflag for particles of type matching "i", used for restricting neighbor search
        if (P[i].Mass <= 0) continue;

#ifdef PM_HIRES_REGION_CLIPDM
        /* here we need to check whether a low-res DM particle is surrounded by all high-res particles,
            in which case we clip its mass down or split it to prevent the most problematic contamination artifacts */
        if(((P[i].Type == 2)||(P[i].Type == 3)||(P[i].Type == 5))&&(TimeBinActive[P[i].TimeBin]))
        {
#ifdef BLACK_HOLES
            if(P[i].Type == 5) continue;
#endif
            /* do a neighbor loop ON THE SAME DOMAIN to determine the neighbors */
            int n_search_min = 32;
            int n_search_max = 320;
            double h_search_max = 10. * All.ForceSoftening[P[i].Type];
            double h_search_min = 0.1 * All.ForceSoftening[P[i].Type];
            double h_guess; numngb_inbox=0; int NITER=0, NITER_MAX=30;
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
            h_guess = PPP[i].AGS_Hsml; if(h_guess > h_search_max) {h_search_max=h_guess;} if(h_guess < h_search_min) {h_search_min=h_guess;}
#else
            h_guess = 5.0 * All.ForceSoftening[P[i].Type];
#endif
            startnode=All.MaxPart;
            do {
                numngb_inbox = ngb_treefind_variable_targeted(P[i].Pos,h_guess,-1,&startnode,0,&dummy,&dummy,62); // search for all particle types -except- gas: 62=2^1+2^2+2^3+2^4+2^5
                if((numngb_inbox < n_search_min) && (h_guess < h_search_max) && (NITER < NITER_MAX))
                {
                    h_guess *= 1.27;
                    startnode=All.MaxPart; // this will trigger the while loop to continue
                }
                if((numngb_inbox > n_search_max) && (h_guess > h_search_min) && (NITER < NITER_MAX))
                {
                    h_guess /= 1.25;
                    startnode=All.MaxPart; // this will trigger the while loop to continue
                }
                NITER++;
            } while(startnode >= 0);
            int do_clipping = 0;
            if(numngb_inbox >= n_search_min-1) // if can't find enough neighbors, don't clip //
            {
                do_clipping = 1;
                for(n=0; n<numngb_inbox; n++)
                {
                    j = Ngblist[n];
                    if(j == i) {if(numngb_inbox > 1) continue;}
#ifdef BLACK_HOLES
                    if((P[j].Type == 2) || (P[j].Type == 3))
#else
                    if((P[j].Type == 2) || (P[j].Type == 3) || (P[j].Type == 5))
#endif
                    {
                        /* found a neighbor with a low-res particle type, so don't clip this particle */
                        do_clipping = 0;
                        break;
                    }
                } // for(n=0; n<numngb_inbox; n++)
            }
            //printf("Particle %d clipping %d low/hi-res DM: neighbors=%d h_search=%g soft=%g iterations=%d \n",i,do_clipping,numngb_inbox,h_guess,All.ForceSoftening[P[i].Type],NITER);
            if(do_clipping)
            {
                /* ok, the particle has neighbors but is completely surrounded by high-res particles, it should be clipped */
                printf("Particle %d clipping low/hi-res DM: neighbors=%d h_search=%g soft=%g iterations=%d \n",i,numngb_inbox,h_guess,All.ForceSoftening[P[i].Type],NITER);
                Ptmp[i].flag = -1;
            }
        }
#endif
#if defined(GALSF)
        if(((P[i].Type==0)||(P[i].Type==4))&&(TimeBinActive[P[i].TimeBin])) /* if SF active, allow star particles to merge if they get too small */
#else
        if((P[i].Type==0)&&(TimeBinActive[P[i].TimeBin])) /* default mode, only gas particles merged */
#endif
        {
            /* we have a gas particle, ask if it needs to be merged */
            if(does_particle_need_to_be_merged(i))
            {
                /* if merging: do a neighbor loop ON THE SAME DOMAIN to determine the target */
                startnode=All.MaxPart;
                numngb_inbox = ngb_treefind_variable_targeted(P[i].Pos,PPP[i].Hsml,-1,&startnode,0,&dummy,&dummy,Pi_BITFLAG); // search for particles of matching type
                if(numngb_inbox>0)
                {
                    target_for_merger = -1;
                    threshold_val = MAX_REAL_NUMBER;
                    for(n=0; n<numngb_inbox; n++) /* loop over neighbors */
                    {
                        j = Ngblist[n]; double m_eff = P[j].Mass; int do_allow_merger = 0; // boolean flag to check
                        if((P[j].Mass >= P[i].Mass) && (P[i].Mass+P[j].Mass < All.MaxMassForParticleSplit)) {do_allow_merger = 1;}
#ifdef BH_WIND_SPAWN
                        if(P[i].ID==All.AGNWindID)
                        {
                            if(P[i].Mass>=MASS_THRESHOLD_FOR_WINDPROMO)
                            {
                                if((P[j].ID!=All.AGNWindID) || (P[j].Mass>=MASS_THRESHOLD_FOR_WINDPROMO)) {do_allow_merger=1;}
                            } else if(do_allow_merger) {
                                double v2_tmp=0,vr_tmp=0; int ktmp=0; for(ktmp=0;ktmp<3;ktmp++) {v2_tmp+=(P[i].Vel[ktmp]-P[j].Vel[ktmp])*(P[i].Vel[ktmp]-P[j].Vel[ktmp]); vr_tmp+=(P[i].Vel[ktmp]-P[j].Vel[ktmp])*(P[i].Pos[ktmp]-P[j].Pos[ktmp]);}
                                if(vr_tmp > 0) {do_allow_merger=0;}
                                if(v2_tmp > 0) {v2_tmp=sqrt(v2_tmp*All.cf_a2inv);} else {v2_tmp=0;}
#if defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS)
                                if(v2_tmp >  DMIN(Get_Gas_effective_soundspeed_i(i),Get_Gas_effective_soundspeed_i(j))) {do_allow_merger = 0;}
                                if(P[j].ID == All.AGNWindID) {do_allow_merger = 0;} // wind particles can't intermerge
#else
                                if((v2_tmp > 0.25*All.BAL_v_outflow) && (v2_tmp > 0.9*Get_Gas_effective_soundspeed_i(j)*All.cf_afac3)) {do_allow_merger=0;}
#endif
                            }
                        }
                        if(P[j].ID == All.AGNWindID) {m_eff *= 1.0e10;} /* boost this enough to ensure the spawned element will never chosen if 'real' candidate exists */
#endif
                        /* make sure we're not taking the same particle (and that its available to be merged into)! and that its the least-massive available candidate for merging onto */
                        if((j<0)||(j==i)||(P[j].Type!=P[i].Type)||(P[j].Mass<=0)||(Ptmp[j].flag!=0)||(m_eff>=threshold_val)) {do_allow_merger=0;}
                        if(do_allow_merger) {threshold_val=m_eff; target_for_merger=j;} /* tell the code this can be merged! */
                    }
                    if (target_for_merger >= 0) { /* mark as merging pairs */
                        Ptmp[i].flag = 1; Ptmp[target_for_merger].flag = 3; Ptmp[i].target_index = target_for_merger;
                    }
                }

            }
            /* now ask if the particle needs to be split */
            else if(does_particle_need_to_be_split(i) && (Ptmp[i].flag == 0)) {
                /* if splitting: do a neighbor loop ON THE SAME DOMAIN to determine the nearest particle (so dont overshoot it) */
                startnode=All.MaxPart;
                numngb_inbox = ngb_treefind_variable_targeted(P[i].Pos,PPP[i].Hsml,-1,&startnode,0,&dummy,&dummy,Pi_BITFLAG); // search for particles of matching type
                if(numngb_inbox>0)
                {
                    target_for_merger = -1;
                    threshold_val = MAX_REAL_NUMBER;
                    /* loop over neighbors */
                    for(n=0; n<numngb_inbox; n++)
                    {
                        j = Ngblist[n];
                        /* make sure we're not taking the same particle */
                        if((j>=0)&&(j!=i)&&(P[j].Type==P[i].Type) && (P[j].Mass > 0) && (Ptmp[j].flag == 0)) {
                            double dp[3]; int k; double r2=0;
                            for(k=0;k<3;k++) {dp[k]=P[i].Pos[k]-P[j].Pos[k];}
                            NEAREST_XYZ(dp[0],dp[1],dp[2],1);
                            for(k=0;k<3;k++) {r2+=dp[k]*dp[k];}
                            if(r2<threshold_val) {threshold_val=r2; target_for_merger=j;} // position-based //
                        }
                    }
                    if (target_for_merger >= 0) {
                        Ptmp[i].flag = 2; // mark for splitting
                        Ptmp[i].target_index = target_for_merger;
                    }
                }
            }
        }
    }

    // actual merge-splitting loop loop
    // No tree-walk is allowed below here
    for (i = 0; i < NumPart; i++) {
#ifdef PM_HIRES_REGION_CLIPDM
        if (Ptmp[i].flag == -1) {
            // clipping
            P[i].Type = 1; // 'graduate' to high-res DM particle
            P[i].Mass = All.MassOfClippedDMParticles; // set mass to the 'safe' mass of typical high-res particles
        }
#endif
        if (Ptmp[i].flag == 1) { // merge this particle
            merge_particles_ij(i, Ptmp[i].target_index);
            n_particles_merged++;
        }
        if (Ptmp[i].flag == 2) {
            split_particle_i(i, n_particles_split, Ptmp[i].target_index);
            n_particles_split++;
            if(P[i].Type==0) {n_particles_gas_split++;}
        }
    }

#ifdef BOX_PERIODIC
    /* map the particles back onto the box (make sure they get wrapped if they go off the edges). this is redundant here,
     because we only do splits in the beginning of a domain decomposition step, where this will be called as soon as
     the particle re-order is completed. but it is still useful to keep here in case this changes (and to note what needs
     to be done for any more complicated splitting operations */
    do_box_wrapping();
#endif
    myfree(Ptmp);
    myfree(Ngblist);
    MPI_Allreduce(&n_particles_merged, &MPI_n_particles_merged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_particles_split, &MPI_n_particles_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_particles_gas_split, &MPI_n_particles_gas_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        if(MPI_n_particles_merged > 0 || MPI_n_particles_split > 0)
        {
            printf("Particle split/merge check: %d particles merged, %d particles split (%d gas) \n", MPI_n_particles_merged,MPI_n_particles_split,MPI_n_particles_gas_split);
        }
    }
    /* the reduction or increase of n_part by MPI_n_particles_merged will occur in rearrange_particle_sequence, which -must- be called immediately after this routine! */
    All.TotNumPart += (long long)MPI_n_particles_split;
    All.TotN_gas += (long long)MPI_n_particles_gas_split;
    Gas_split = n_particles_gas_split; // specific to the local processor //
    NumPart += (n_particles_split - n_particles_gas_split); // specific to the local processor; note the gas split number will be added below, this is just non-gas splits //
}




/*! This is the routine that does the particle splitting. Note this is a tricky operation if we're not using meshes to divide the volume,
    so care needs to be taken modifying this so that it's done in a way that is (1) conservative, (2) minimizes perturbations to the
    volumetric quantities of the flow, and (3) doesn't crash the tree or lead to particle 'overlap'
    Modified by Takashi Okamoto on 20/6/2019.  */
//void split_particle_i(int i, int n_particles_split, int i_nearest, double r2_nearest)
void split_particle_i(int i, int n_particles_split, int i_nearest)
{
    double mass_of_new_particle;
    if( ((P[i].Type==0) && (NumPart + n_particles_split >= All.MaxPartSph)) || ((P[i].Type!=0) && (NumPart + n_particles_split >= All.MaxPart)) )
    {
        printf ("On Task=%d with NumPart=%d we tried to split a particle, but there is no space left...(All.MaxPart=%d). Try using more nodes, or raising PartAllocFac, or changing the split conditions to avoid this.\n", ThisTask, NumPart, All.MaxPart);
        fflush(stdout);
        endrun(8888);
    }
#ifndef SPAWN_PARTICLES_VIA_SPLITTING
    if(P[i].Type != 0) {printf("SPLITTING NON-GAS-PARTICLE: i=%d ID=%llu Type=%d \n",i,(unsigned long long) P[i].ID,P[i].Type);} //fflush(stdout); endrun(8889);
#endif

    /* here is where the details of the split are coded, the rest is bookkeeping */
    mass_of_new_particle = 0.5;

    int k; double phi,cos_theta;
    k=0;
    phi = 2.0*M_PI*get_random_number(i+1+ThisTask); // random from 0 to 2pi //
    cos_theta = 2.0*(get_random_number(i+3+2*ThisTask)-0.5); // random between 1 to -1 //
    double d_r = 0.25 * KERNEL_CORE_SIZE*PPP[i].Hsml; // needs to be epsilon*Hsml where epsilon<<1, to maintain stability //
    double dp[3], r_near=0; for(k = 0; k < 3; k++) {dp[k] =P[i].Pos[k] - P[i_nearest].Pos[k];}
    NEAREST_XYZ(dp[0],dp[1],dp[2],1);
    for(k = 0; k < 3; k++) {r_near += dp[k]*dp[k];}
    r_near = 0.35 * sqrt(r_near);
    d_r = DMIN(d_r , r_near); // use a 'buffer' to limit to some multiple of the distance to the nearest particle //
    /*
    double r_near = sqrt(r2_nearest);
    double hsml = Get_Particle_Size(i);
    if(hsml < r_near) {hsml = r_near;}
    r_near *= 0.35;
    double d_r = 0.25 * hsml; // needs to be epsilon*Hsml where epsilon<<1, to maintain stability //
    d_r = DMAX( DMAX(0.1*r_near , 0.005*hsml) , DMIN(d_r , r_near) ); // use a 'buffer' to limit to some multiple of the distance to the nearest particle //
    */ // the change above appears to cause some numerical instability //
#ifndef SELFGRAVITY_OFF
    d_r = DMAX(d_r , 2.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[P[i].Type]);
#endif
#ifdef BOX_BND_PARTICLES
    if(P[i].Type != 0 && P[i].ID == 0) {d_r *= 1.e-3;}
#endif

    /* find the first non-gas particle and move it to the end of the particle list */
    long j = NumPart + n_particles_split;
    /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */
    P[j] = P[i];
    //memcpy(P[j],P[i],sizeof(struct particle_data)); // safer copy to make sure we don't just end up with a pointer re-direct

#ifdef CHIMES
    int abunIndex;
    ChimesGasVars[j] = ChimesGasVars[i];
    allocate_gas_abundances_memory(&(ChimesGasVars[j]), &ChimesGlobalVars);
    for (abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++) {ChimesGasVars[j].abundances[abunIndex] = ChimesGasVars[i].abundances[abunIndex];}
#endif

    // need to assign new particle a unique ID:
    // new method: preserve the original "ID" field, but assign a unique -child- ID: this is unique up to ~32 *GENERATIONS* of repeated splitting!
    P[j].ID_child_number = P[i].ID_child_number + (MyIDType)(1 << ((int)P[i].ID_generation)); // particle 'i' retains its child number; this ensures uniqueness
    P[i].ID_generation = P[i].ID_generation + 1;
    if(P[i].ID_generation > 30) {P[i].ID_generation=0;} // roll over at 32 generations (unlikely to ever reach this)
    P[j].ID_generation = P[i].ID_generation; // ok, all set!

    /* assign masses to both particles (so they sum correctly) */
    P[j].Mass = mass_of_new_particle * P[i].Mass;
    P[i].Mass -= P[j].Mass;
#ifdef BOX_BND_PARTICLES
    if(P[i].ID==0)  {P[j].ID=1; double m0=P[i].Mass+P[j].Mass; P[i].Mass=P[j].Mass=m0;}
#endif

    /* prepare to shift the particle locations according to the random number we drew above */
    double dx, dy, dz;
#if (NUMDIMS == 1)
    dy=dz=0; dx=d_r; // here the split direction is trivial //
#else
    /* in 2D and 3D its not so trivial how to split the directions */
    double sin_theta = sqrt(1 - cos_theta*cos_theta);
    dx = d_r * sin_theta * cos(phi);
    dy = d_r * sin_theta * sin(phi);
    dz = d_r * cos_theta;
#if (NUMDIMS == 2)
    dz=0; dx=d_r*cos(phi); dy=d_r*sin(phi);
#endif
#endif

    if(P[i].Type==0)
    {
        /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */
        SphP[j] = SphP[i];
        //memcpy(SphP[j],SphP[i],sizeof(struct sph_particle_data)); // safer copy to make sure we don't just end up with a pointer re-direct
        /* boost the condition number to be conservative, so we don't trigger madness in the kernel */
        SphP[i].ConditionNumber *= 10.0;
        SphP[j].ConditionNumber = SphP[i].ConditionNumber;
#ifdef MAGNETIC
        /* we evolve the -conserved- VB and Vphi, so this must be partitioned */
        for(k=0;k<3;k++)
        {
            SphP[j].B[k] = mass_of_new_particle * SphP[i].B[k]; SphP[i].B[k] -= SphP[j].B[k];
            SphP[j].BPred[k] = mass_of_new_particle * SphP[i].BPred[k]; SphP[i].BPred[k] -= SphP[j].BPred[k];
            SphP[j].DtB[k] = mass_of_new_particle * SphP[i].DtB[k]; SphP[i].DtB[k] -= SphP[j].DtB[k];
        }
        SphP[j].divB = mass_of_new_particle * SphP[i].divB; SphP[i].divB -= SphP[j].divB;
#ifdef DIVBCLEANING_DEDNER
        SphP[j].Phi = mass_of_new_particle * SphP[i].Phi; SphP[i].Phi -= SphP[j].Phi;
        SphP[j].DtPhi = mass_of_new_particle * SphP[i].DtPhi; SphP[i].DtPhi -= SphP[j].DtPhi;
        SphP[j].PhiPred = mass_of_new_particle * SphP[i].PhiPred; SphP[i].PhiPred -= SphP[j].PhiPred;
#endif
        /* ideally, particle-splits should be accompanied by a re-partition of the density via the density() call
         for the particles affected, after the tree-reconstruction, with quantities like B used to re-calculate after */
#endif
#ifdef RADTRANSFER
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            int k_dir; k_dir=0;
            SphP[j].Rad_E_gamma[k] = mass_of_new_particle * SphP[i].Rad_E_gamma[k]; SphP[i].Rad_E_gamma[k] -= SphP[j].Rad_E_gamma[k];
#if defined(RT_EVOLVE_ENERGY)
            SphP[j].Rad_E_gamma_Pred[k] = mass_of_new_particle * SphP[i].Rad_E_gamma_Pred[k]; SphP[i].Rad_E_gamma_Pred[k] -= SphP[j].Rad_E_gamma_Pred[k];
            SphP[j].Dt_Rad_E_gamma[k] = mass_of_new_particle * SphP[i].Dt_Rad_E_gamma[k]; SphP[i].Dt_Rad_E_gamma[k] -= SphP[j].Dt_Rad_E_gamma[k];
#endif
#if defined(RT_EVOLVE_FLUX)
            for(k_dir=0;k_dir<3;k_dir++)
            {
                SphP[j].Rad_Flux[k][k_dir] = mass_of_new_particle * SphP[i].Rad_Flux[k][k_dir]; SphP[i].Rad_Flux[k][k_dir] -= SphP[j].Rad_Flux[k][k_dir];
                SphP[j].Rad_Flux_Pred[k][k_dir] = mass_of_new_particle * SphP[i].Rad_Flux_Pred[k][k_dir]; SphP[i].Rad_Flux_Pred[k][k_dir] -= SphP[j].Rad_Flux_Pred[k][k_dir];
                SphP[j].Dt_Rad_Flux[k][k_dir] = mass_of_new_particle * SphP[i].Dt_Rad_Flux[k][k_dir]; SphP[i].Dt_Rad_Flux[k][k_dir] -= SphP[j].Dt_Rad_Flux[k][k_dir];
            }
#endif
#ifdef RT_EVOLVE_INTENSITIES
            for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++)
            {
                SphP[j].Rad_Intensity[k][k_dir] = mass_of_new_particle * SphP[i].Rad_Intensity[k][k_dir]; SphP[i].Rad_Intensity[k][k_dir] -= SphP[j].Rad_Intensity[k][k_dir];
                SphP[j].Rad_Intensity_Pred[k][k_dir] = mass_of_new_particle * SphP[i].Rad_Intensity_Pred[k][k_dir]; SphP[i].Rad_Intensity_Pred[k][k_dir] -= SphP[j].Rad_Intensity_Pred[k][k_dir];
                SphP[j].Dt_Rad_Intensity[k][k_dir] = mass_of_new_particle * SphP[i].Dt_Rad_Intensity[k][k_dir]; SphP[i].Dt_Rad_Intensity[k][k_dir] -= SphP[j].Dt_Rad_Intensity[k][k_dir];
            }
#endif
        }
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        double dmass = mass_of_new_particle * SphP[i].DtMass;
        SphP[j].DtMass = dmass;
        SphP[i].DtMass -= dmass;
        dmass = mass_of_new_particle * SphP[i].dMass;
        SphP[j].dMass = dmass;
        SphP[i].dMass -= dmass;
        for(k=0;k<3;k++)
        {
            SphP[j].GravWorkTerm[k] =0;//= mass_of_new_particle * SphP[i].GravWorkTerm[k];//appears more stable with this zero'd
            SphP[i].GravWorkTerm[k] =0;//-= SphP[j].GravWorkTerm[k];//appears more stable with this zero'd
        }
        SphP[j].MassTrue = mass_of_new_particle * SphP[i].MassTrue;
        SphP[i].MassTrue -= SphP[j].MassTrue;
#endif
#ifdef COSMIC_RAY_FLUID
        int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {
            SphP[j].CosmicRayEnergy[k_CRegy] = mass_of_new_particle * SphP[i].CosmicRayEnergy[k_CRegy]; SphP[i].CosmicRayEnergy[k_CRegy] -= SphP[j].CosmicRayEnergy[k_CRegy];
            SphP[j].CosmicRayEnergyPred[k_CRegy] = mass_of_new_particle * SphP[i].CosmicRayEnergyPred[k_CRegy]; SphP[i].CosmicRayEnergyPred[k_CRegy] -= SphP[j].CosmicRayEnergyPred[k_CRegy];
            SphP[j].DtCosmicRayEnergy[k_CRegy] = mass_of_new_particle * SphP[i].DtCosmicRayEnergy[k_CRegy]; SphP[i].DtCosmicRayEnergy[k_CRegy] -= SphP[j].DtCosmicRayEnergy[k_CRegy];
#ifdef CRFLUID_M1
            for(k=0;k<3;k++)
            {
                SphP[j].CosmicRayFlux[k_CRegy][k] = mass_of_new_particle * SphP[i].CosmicRayFlux[k_CRegy][k]; SphP[i].CosmicRayFlux[k_CRegy][k] -= SphP[j].CosmicRayFlux[k_CRegy][k];
                SphP[j].CosmicRayFluxPred[k_CRegy][k] = mass_of_new_particle * SphP[i].CosmicRayFluxPred[k_CRegy][k]; SphP[i].CosmicRayFluxPred[k_CRegy][k] -= SphP[j].CosmicRayFluxPred[k_CRegy][k];
            }
#endif
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(k=0;k<2;k++)
            {
                SphP[j].CosmicRayAlfvenEnergy[k_CRegy][k] = mass_of_new_particle * SphP[i].CosmicRayAlfvenEnergy[k_CRegy][k]; SphP[i].CosmicRayAlfvenEnergy[k_CRegy][k] -= SphP[j].CosmicRayAlfvenEnergy[k_CRegy][k];
                SphP[j].CosmicRayAlfvenEnergyPred[k_CRegy][k] = mass_of_new_particle * SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][k]; SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][k] -= SphP[j].CosmicRayAlfvenEnergyPred[k_CRegy][k];
                SphP[j].DtCosmicRayAlfvenEnergy[k_CRegy][k] = mass_of_new_particle * SphP[i].DtCosmicRayAlfvenEnergy[k_CRegy][k]; SphP[i].DtCosmicRayAlfvenEnergy[k_CRegy][k] -= SphP[j].DtCosmicRayAlfvenEnergy[k_CRegy][k];
            }
#endif
        }
#endif

        /* use a better particle shift based on the moment of inertia tensor to place new particles in the direction which is less well-sampled */
#if (NUMDIMS > 1)        
        double norm=0, dp[3]; int m; dp[0]=dp[1]=dp[2]=0;

        // get the eigenvector of NV_T that has the smallest eigenvalue (= sparsest sampling direction)
        double nvt[NUMDIMS*NUMDIMS]={0}; for(k=0;k<NUMDIMS;k++) {for(m=0;m<NUMDIMS;m++) {nvt[NUMDIMS*k + m]=SphP[i].NV_T[k][m];}} // auxiliary array to store NV_T in for feeding to GSL eigen routine
        gsl_matrix_view M = gsl_matrix_view_array(nvt,NUMDIMS,NUMDIMS); gsl_vector *eigvals = gsl_vector_alloc(NUMDIMS); gsl_matrix *eigvecs = gsl_matrix_alloc(NUMDIMS,NUMDIMS);
        gsl_eigen_symmv_workspace *v = gsl_eigen_symmv_alloc(NUMDIMS); gsl_eigen_symmv(&M.matrix, eigvals, eigvecs, v);
        int min_eigvec_index = 0; double min_eigval = MAX_REAL_NUMBER;
        for(k=0;k<NUMDIMS;k++) {if(gsl_vector_get(eigvals,k) < min_eigval){min_eigval = gsl_vector_get(eigvals,k); min_eigvec_index=k;}}
        for(k=0;k<NUMDIMS;k++) {dp[k] = gsl_matrix_get(eigvecs, k, min_eigvec_index);}
        gsl_eigen_symmv_free(v); gsl_vector_free(eigvals); gsl_matrix_free(eigvecs);
        for(k=0;k<NUMDIMS;k++) {norm += dp[k] * dp[k];}
        if(norm > 0)
        {
            norm = 1/sqrt(norm); for(k=0;k<NUMDIMS;k++) {dp[k] *= norm;}
            dx=d_r*dp[0]; dy=d_r*dp[1]; dz=d_r*dp[2];
            /* rotate to 90-degree offset from above orientation, if using the density gradient */
            // if(dp[2]==1) {dx=d_r; dy=0; dz=0;} else {dz = sqrt(dp[1]*dp[1] + dp[0]*dp[0]); dx = -d_r * dp[1]/dz; dy = d_r * dp[0]/dz; dz = 0.0;}
        }
#endif
#ifdef WAKEUP  /* TO: rather conservative. But we want to update Density and Hsml after the particle masses were changed */
        PPPZ[i].wakeup = 1; PPPZ[j].wakeup = 1; NeedToWakeupParticles_local = 1;
#endif

    } // closes special operations required only of gas particles

    /* this is allowed to push particles over the 'edges' of periodic boxes, because we will call the box-wrapping routine immediately below.
     but it is important that the periodicity of the box be accounted for in relative positions and that we correct for this before allowing
     any other operations on the particles */
    P[i].Pos[0] += dx; P[j].Pos[0] -= dx; P[i].Pos[1] += dy; P[j].Pos[1] -= dy; P[i].Pos[2] += dz; P[j].Pos[2] -= dz;

    /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
#ifdef PARTICLE_MERGE_SPLIT_EVERY_TIMESTEP    
    long bin = P[i].TimeBin;
    if(FirstInTimeBin[bin] < 0) {FirstInTimeBin[bin]=j; LastInTimeBin[bin]=j; NextInTimeBin[j]=-1; PrevInTimeBin[j]=-1;} /* only particle in this time bin on this task */
    else {NextInTimeBin[j]=FirstInTimeBin[bin]; PrevInTimeBin[j]=-1; PrevInTimeBin[FirstInTimeBin[bin]]=j; FirstInTimeBin[bin]=j;} /* there is already at least one particle; add this one "to the front" of the list */
    force_add_star_to_tree(i, j);
#endif    
    /* we solve this by only calling the merge/split algorithm when we're doing the new domain decomposition */
}



/*! Routine to merge particle 'i' into particle 'j'
    The volumetric quantities (density, pressure, gradients, kernel lengths)
    can be re-estimated after the merger operation. but we need to make sure
    all conserved quantities are appropriately dealt with. This also requires some care, to be
    done appropriately, but is a little bit less sensitive and more well-defined compared to
    particle splitting */
void merge_particles_ij(int i, int j)
{
    int k;
    if(P[i].Mass <= 0)
    {
        P[i].Mass = 0;
        return;
    }
    if(P[j].Mass <= 0)
    {
        P[j].Mass = 0;
        return;
    }
    double mtot = P[j].Mass + P[i].Mass;
    double wt_i = P[i].Mass / mtot;
    double wt_j = P[j].Mass / mtot;

    // block for merging non-gas particles (much simpler, assume collisionless)
    if((P[i].Type>0)&&(P[j].Type>0))
    {
        double pos_new_xyz[3], dp[3];
        for(k=0;k<3;k++) {dp[k]=P[j].Pos[k]-P[i].Pos[k];}
        NEAREST_XYZ(dp[0],dp[1],dp[2],-1);
        for(k=0;k<3;k++) {pos_new_xyz[k] = P[i].Pos[k] + wt_j * dp[k];}

        double p_old_i[3],p_old_j[3];
        for(k=0;k<3;k++)
        {
            p_old_i[k] = P[i].Mass * P[i].Vel[k];
            p_old_j[k] = P[j].Mass * P[j].Vel[k];
        }
        for(k=0;k<3;k++)
        {
            P[j].Pos[k] = pos_new_xyz[k]; // center-of-mass conserving //
            P[j].Vel[k] = wt_j*P[j].Vel[k] + wt_i*P[i].Vel[k]; // momentum-conserving //
            P[j].GravAccel[k] = wt_j*P[j].GravAccel[k] + wt_i*P[i].GravAccel[k]; // force-conserving //
#ifdef PMGRID
            P[j].GravPM[k] = wt_j*P[j].GravPM[k] + wt_i*P[i].GravPM[k]; // force-conserving //
#endif
        }
        PPP[j].Hsml = pow(pow(PPP[j].Hsml,NUMDIMS)+pow(PPP[i].Hsml,NUMDIMS),1.0/NUMDIMS);
#ifdef METALS
        for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k] = wt_j*P[j].Metallicity[k] + wt_i*P[i].Metallicity[k];} /* metal-mass conserving */
#endif
        /* finally zero out the particle mass so it will be deleted */
        P[i].Mass = 0;
        P[j].Mass = mtot;
        for(k=0;k<3;k++)
        {
            /* momentum shift for passing to tree (so we know how to move it) */
            P[i].dp[k] += P[i].Mass*P[i].Vel[k] - p_old_i[k];
            P[j].dp[k] += P[j].Mass*P[j].Vel[k] - p_old_j[k];
        }
        return;
    } // closes merger of non-gas particles, only gas particles will see the blocks below //


    // now we have to deal with gas particle mergers //
    if(P[i].TimeBin < P[j].TimeBin)
    {
#ifdef WAKEUP
        PPPZ[j].wakeup = 1; NeedToWakeupParticles_local = 1;
#endif
    }
    double dm_i=0,dm_j=0,de_i=0,de_j=0,dp_i[3],dp_j[3],dm_ij,de_ij,dp_ij[3];
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    dm_i += SphP[i].DtMass;
    dm_j += SphP[j].DtMass;
#endif
    de_i = P[i].Mass * SphP[i].DtInternalEnergy + dm_i*SphP[i].InternalEnergy;
    de_j = P[j].Mass * SphP[j].DtInternalEnergy + dm_j*SphP[j].InternalEnergy;
    for(k=0;k<3;k++)
    {
        dp_i[k] = P[i].Mass * SphP[i].HydroAccel[k] + dm_i * SphP[i].VelPred[k] / All.cf_atime;
        dp_j[k] = P[j].Mass * SphP[j].HydroAccel[k] + dm_j * SphP[j].VelPred[k] / All.cf_atime;
        de_i += dp_i[k] * SphP[i].VelPred[k] / All.cf_atime - 0.5 * dm_i * SphP[i].VelPred[k] * SphP[i].VelPred[k] * All.cf_a2inv;
        de_j += dp_j[k] * SphP[j].VelPred[k] / All.cf_atime - 0.5 * dm_j * SphP[j].VelPred[k] * SphP[j].VelPred[k] * All.cf_a2inv;
        dp_ij[k] = dp_i[k] + dp_j[k];
    }
    dm_ij = dm_i+dm_j;
    de_ij = de_i+de_j;

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    SphP[j].MassTrue += SphP[i].MassTrue;
    SphP[i].MassTrue = 0;
    SphP[j].DtMass=dm_ij;
    SphP[i].DtMass=0;
    SphP[j].dMass = SphP[i].dMass + SphP[j].dMass;
    SphP[i].dMass = 0;
#endif

    /* make sure to update the conserved variables correctly: mass and momentum are easy, energy is non-trivial */
    double egy_old = 0;
    egy_old += mtot * (wt_j*SphP[j].InternalEnergy + wt_i*SphP[i].InternalEnergy); // internal energy //
    double pos_new_xyz[3], dp[3];
    /* for periodic boxes, we need to (arbitrarily) pick one position as our coordinate center. we pick i. then everything defined in
        position differences relative to i. the final position will be appropriately box-wrapped after these operations are completed */
    for(k=0;k<3;k++) {dp[k]=P[j].Pos[k]-P[i].Pos[k];}
    NEAREST_XYZ(dp[0],dp[1],dp[2],-1);
    for(k=0;k<3;k++) {pos_new_xyz[k] = P[i].Pos[k] + wt_j * dp[k];}

    for(k=0;k<3;k++)
    {
        egy_old += mtot*wt_j * 0.5 * P[j].Vel[k]*P[j].Vel[k]*All.cf_a2inv; // kinetic energy (j) //
        egy_old += mtot*wt_i * 0.5 * P[i].Vel[k]*P[i].Vel[k]*All.cf_a2inv; // kinetic energy (i) //
        // gravitational energy terms need to be added (including work for moving particles 'together') //
        // Egrav = m*g*h = m * (-grav_acc) * (position relative to zero point) //
        egy_old += mtot*wt_j * (P[i].Pos[k]+dp[k] - pos_new_xyz[k])*All.cf_atime * (-P[j].GravAccel[k])*All.cf_a2inv; // work (j) //
        egy_old += mtot*wt_i * (P[i].Pos[k] - pos_new_xyz[k])*All.cf_atime * (-P[i].GravAccel[k])*All.cf_a2inv; // work (i) //
#ifdef PMGRID
        egy_old += mtot*wt_j * (P[i].Pos[k]+dp[k] - pos_new_xyz[k])*All.cf_atime * (-P[j].GravPM[k])*All.cf_a2inv; // work (j) [PMGRID] //
        egy_old += mtot*wt_i * (P[i].Pos[k] - pos_new_xyz[k])*All.cf_atime * (-P[i].GravPM[k])*All.cf_a2inv; // work (i) [PMGRID] //
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].GravWorkTerm[k] = 0; // since we're accounting for the work above and dont want to accidentally double-count //
#endif
    }


    SphP[j].InternalEnergy = wt_j*SphP[j].InternalEnergy + wt_i*SphP[i].InternalEnergy;
    SphP[j].InternalEnergyPred = wt_j*SphP[j].InternalEnergyPred + wt_i*SphP[i].InternalEnergyPred;
    double p_old_i[3],p_old_j[3];
    for(k=0;k<3;k++)
    {
        p_old_i[k] = P[i].Mass * P[i].Vel[k];
        p_old_j[k] = P[j].Mass * P[j].Vel[k];
    }
    for(k=0;k<3;k++)
    {
        P[j].Pos[k] = pos_new_xyz[k]; // center-of-mass conserving //
        P[j].Vel[k] = wt_j*P[j].Vel[k] + wt_i*P[i].Vel[k]; // momentum-conserving //
        SphP[j].VelPred[k] = wt_j*SphP[j].VelPred[k] + wt_i*SphP[i].VelPred[k]; // momentum-conserving //
        P[j].GravAccel[k] = wt_j*P[j].GravAccel[k] + wt_i*P[i].GravAccel[k]; // force-conserving //
#ifdef PMGRID
        P[j].GravPM[k] = wt_j*P[j].GravPM[k] + wt_i*P[i].GravPM[k]; // force-conserving //
#endif
    }
#ifdef MAGNETIC
    // we evolve the conservative variables VB and Vpsi, these should simply add in particle-merge operations //
    for(k=0;k<3;k++)
    {
        SphP[j].B[k] += SphP[i].B[k];
        SphP[j].BPred[k] += SphP[i].BPred[k];
        SphP[j].DtB[k] += SphP[i].DtB[k];
    }
#ifdef DIVBCLEANING_DEDNER
    SphP[j].Phi += SphP[i].Phi;
    SphP[j].PhiPred += SphP[i].PhiPred;
    SphP[j].DtPhi += SphP[i].DtPhi;
#endif
#endif

    /* correct our 'guess' for the internal energy with the residual from exact energy conservation */
    double egy_new = mtot * SphP[j].InternalEnergy;
    for(k=0;k<3;k++) {egy_new += mtot * 0.5*P[j].Vel[k]*P[j].Vel[k]*All.cf_a2inv;}
    egy_new = (egy_old - egy_new) / mtot; /* this residual needs to be put into the thermal energy */
    if(egy_new < -0.5*SphP[j].InternalEnergy) egy_new = -0.5 * SphP[j].InternalEnergy;
    //SphP[j].InternalEnergy += egy_new; SphP[j].InternalEnergyPred += egy_new;//test during splits
    if(SphP[j].InternalEnergyPred<0.5*SphP[j].InternalEnergy) SphP[j].InternalEnergyPred=0.5*SphP[j].InternalEnergy;


    // now use the conserved variables to correct the derivatives to primitive variables //
    de_ij -= dm_ij * SphP[j].InternalEnergyPred;
    for(k=0;k<3;k++)
    {
        SphP[j].HydroAccel[k] = (dp_ij[k] - dm_ij * SphP[j].VelPred[k]/All.cf_atime) / mtot;
        de_ij -= mtot * SphP[j].VelPred[k]/All.cf_atime * SphP[j].HydroAccel[k] + 0.5 * dm_ij * SphP[j].VelPred[k]*SphP[j].VelPred[k]*All.cf_a2inv;
    }
    SphP[j].DtInternalEnergy = de_ij;
    // to be conservative adopt the maximum signal velocity and kernel length //
    SphP[j].MaxSignalVel = sqrt(SphP[j].MaxSignalVel*SphP[j].MaxSignalVel + SphP[i].MaxSignalVel*SphP[i].MaxSignalVel); /* need to be conservative */
    PPP[j].Hsml = pow(pow(PPP[j].Hsml,NUMDIMS)+pow(PPP[i].Hsml,NUMDIMS),1.0/NUMDIMS); /* sum the volume of the two particles */
    SphP[j].ConditionNumber = SphP[j].ConditionNumber + SphP[i].ConditionNumber; /* sum to be conservative */
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
    SphP[j].MaxKineticEnergyNgb = DMAX(SphP[j].MaxKineticEnergyNgb,SphP[i].MaxKineticEnergyNgb); /* for the entropy/energy switch condition */
#endif

    // below, we need to take care of additional physics //
#if defined(RADTRANSFER)
    for(k=0;k<N_RT_FREQ_BINS;k++)
    {
        int k_dir;
        for(k_dir=0;k_dir<6;k_dir++) SphP[j].ET[k][k_dir] = wt_j*SphP[j].ET[k][k_dir] + wt_i*SphP[i].ET[k][k_dir];
        SphP[j].Rad_E_gamma[k] = SphP[j].Rad_E_gamma[k] + SphP[i].Rad_E_gamma[k]; /* this is a photon number, so its conserved (we simply add) */
#if defined(RT_EVOLVE_ENERGY)
        SphP[j].Rad_E_gamma_Pred[k] = SphP[j].Rad_E_gamma_Pred[k] + SphP[i].Rad_E_gamma_Pred[k];
        SphP[j].Dt_Rad_E_gamma[k] = SphP[j].Dt_Rad_E_gamma[k] + SphP[i].Dt_Rad_E_gamma[k];
#endif
#if defined(RT_EVOLVE_FLUX)
        for(k_dir=0;k_dir<3;k_dir++)
        {
            SphP[j].Rad_Flux[k][k_dir] = SphP[j].Rad_Flux[k][k_dir] + SphP[i].Rad_Flux[k][k_dir];
            SphP[j].Rad_Flux_Pred[k][k_dir] = SphP[j].Rad_Flux_Pred[k][k_dir] + SphP[i].Rad_Flux_Pred[k][k_dir];
            SphP[j].Dt_Rad_Flux[k][k_dir] = SphP[j].Dt_Rad_Flux[k][k_dir] + SphP[i].Dt_Rad_Flux[k][k_dir];
        }
#endif
#ifdef RT_EVOLVE_INTENSITIES
        for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++)
        {
            SphP[j].Rad_Intensity[k][k_dir] = SphP[j].Rad_Intensity[k][k_dir] + SphP[i].Rad_Intensity[k][k_dir];
            SphP[j].Rad_Intensity_Pred[k][k_dir] = SphP[j].Rad_Intensity_Pred[k][k_dir] + SphP[i].Rad_Intensity_Pred[k][k_dir];
            SphP[j].Dt_Rad_Intensity[k][k_dir] = SphP[j].Dt_Rad_Intensity[k][k_dir] + SphP[i].Dt_Rad_Intensity[k][k_dir];
        }
#endif
    }
#endif
#ifdef CHIMES
    double wt_h_i, wt_h_j; // Ratio of hydrogen mass fractions.
#ifdef COOL_METAL_LINES_BY_SPECIES
    wt_h_i = 1.0 - (P[i].Metallicity[0] + P[i].Metallicity[1]);
    wt_h_j = 1.0 - (P[j].Metallicity[0] + P[j].Metallicity[1]);
#else
    wt_h_i = 1.0;
    wt_h_j = 1.0;
#endif
    for (k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++) {ChimesGasVars[j].abundances[k] = (ChimesFloat) ((ChimesGasVars[j].abundances[k] * wt_j * wt_h_j) + (ChimesGasVars[i].abundances[k] * wt_i * wt_h_i));}
#endif // CHIMES
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k] = wt_j*P[j].Metallicity[k] + wt_i*P[i].Metallicity[k];} /* metal-mass conserving */
#endif
#ifdef COSMIC_RAY_FLUID
    int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        SphP[j].CosmicRayEnergy[k_CRegy] += SphP[i].CosmicRayEnergy[k_CRegy];
        SphP[j].CosmicRayEnergyPred[k_CRegy] += SphP[i].CosmicRayEnergyPred[k_CRegy];
        SphP[j].DtCosmicRayEnergy[k_CRegy] += SphP[i].DtCosmicRayEnergy[k_CRegy];
#ifdef CRFLUID_M1
        for(k=0;k<3;k++)
        {
            SphP[j].CosmicRayFlux[k_CRegy][k] += SphP[i].CosmicRayFlux[k_CRegy][k];
            SphP[j].CosmicRayFluxPred[k_CRegy][k] += SphP[i].CosmicRayFluxPred[k_CRegy][k];
        }
#endif
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
        for(k=0;k<3;k++)
        {
            SphP[j].CosmicRayAlfvenEnergy[k_CRegy][k] += SphP[i].CosmicRayAlfvenEnergy[k_CRegy][k];
            SphP[j].CosmicRayAlfvenEnergyPred[k_CRegy][k] += SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][k];
            SphP[j].DtCosmicRayAlfvenEnergy[k_CRegy][k] += SphP[i].DtCosmicRayAlfvenEnergy[k_CRegy][k];
        }
#endif
    }
#endif

    /* finally zero out the particle mass so it will be deleted */
    P[i].Mass = 0;
    P[j].Mass = mtot;
    for(k=0;k<3;k++)
    {
        /* momentum shift for passing to tree (so we know how to move it) */
        P[i].dp[k] += P[i].Mass*P[i].Vel[k] - p_old_i[k];
        P[j].dp[k] += P[j].Mass*P[j].Vel[k] - p_old_j[k];
    }
    /* call the pressure routine to re-calculate pressure (and sound speeds) as needed */
    SphP[j].Pressure = get_pressure(j);
    return;
}


/*!
  This routine swaps the location of two pointers/indices (either to a node or to a particle) in the treewalk needed for neighbor searches and gravity.
  This should be run if you are messing around with the indices of things but don't intend to do a whole treebuild after. - MYG
 */

void swap_treewalk_pointers(int i, int j){
    // walk the tree and any time we see a nextnode or sibling set to i, swap it to j and vice versa
    int no, next, pre_sibling_i=-1, pre_sibling_j=-1, previous_node_i, previous_node_j;
    no = All.MaxPart;

    while(no >= 0){ // walk the whole tree, starting from the root node (=All.MaxPart)
        if(no < All.MaxPart) { // we got a particle
            next=Nextnode[no];
            if(no != i && no != j){ // don't mess with Nextnodes if looking at i or j - handle that separately
                if(next == i) {Nextnode[no] = j; previous_node_i = no;}
                else if(next == j) { Nextnode[no] = i; previous_node_j = no;}
            }
            no = next;
        } else if(no < All.MaxPart+MaxNodes)  { // we have a node
            next = Nodes[no].u.d.nextnode;
            if(next == i) { previous_node_i = no; Nodes[no].u.d.nextnode = j;}
            else if(next == j) { previous_node_j = no; Nodes[no].u.d.nextnode = i;}
            if(Nodes[no].u.d.sibling == i) {Nodes[no].u.d.sibling = j; pre_sibling_i = no;}
            else if(Nodes[no].u.d.sibling == j) { Nodes[no].u.d.sibling = i; pre_sibling_j = no;}
            no = next;
        } else { // pseudoparticle
            next = Nextnode[no - MaxNodes];
            if(next==i) {Nextnode[no-MaxNodes] = j;}
            else if(next == j) {Nextnode[no-MaxNodes] = i;}
            no = next;
        }
    }

    if(Nextnode[i] == j){ // handle case if i->j
        Nextnode[i] = Nextnode[j];
        Nextnode[j] = i;
    } else if (Nextnode[j] == i) { // if j->i
        Nextnode[j] = Nextnode[i];
        Nextnode[i] = j;
    } else { // neither i->j nor j->i
        no = Nextnode[i];
        Nextnode[i] = Nextnode[j];
        Nextnode[j] = no;
    }
    no = Father[i];
    Father[i] = Father[j];
    Father[j] = no;
}


/*!
  This routine deletes a particle from the linked list for the treewalk, preserving the lists's integrity. This must be run if you are deleting particles but don't want to do a while treebuild after. - MYG
*/
void remove_particle_from_treewalk(int i){
    int no, next;
    no = All.MaxPart;
    while(no >= 0){ // walk the tree to find anything that might point to i and redirect it to i's nextnode
        if(no < All.MaxPart){
            next = Nextnode[no];
            if(Nextnode[no] == i) {Nextnode[no] = Nextnode[i];}
        } else if (no < All.MaxPart+MaxNodes){
            next = Nodes[no].u.d.nextnode;
            if(next == i) {Nodes[no].u.d.nextnode = Nextnode[i];}
            if(Nodes[no].u.d.sibling == i) {Nodes[no].u.d.sibling = Nextnode[i];}
        } else {
            next = Nextnode[no - MaxNodes];
            if(next == i) {Nextnode[no - MaxNodes] = Nextnode[i];}
        }
        no = next;
    }
}

/*! This is an important routine used throughout -- any time particle masses are variable OR particles can
    be created/destroyed: it goes through the particle list, makes sure they are in the appropriate order (gas
    must all come before collisionless particles, though the collisionless particles can be blocked into any order
    we like), swaps particles as needed to restore the correct ordering, eliminates particles from the main list
    if they have zero or negative mass (i.e. this does the actual 'deletion' operation), and then reconstructs the
    list of particles in each timestep (if particles had to be re-ordered) so that the code will not crash looking for
    the 'wrong' particles (or non-existent particles). In general, if you do any operations that involve particle
    creation, this needs to be called before anything is 'done' with those particles. If you do anything involving particle
    'deletion', the standard procedure should be to set the deleted particle mass to zero, and then let this routine
    (when it is called in standard sequence) do its job and 'clean up' the particle
 */
void rearrange_particle_sequence(void)
{
    int i, j, flag = 0, flag_sum;
    int count_elim, count_gaselim, count_bhelim, tot_elim, tot_gaselim, tot_bhelim;
    struct particle_data psave;
    struct sph_particle_data sphsave;
#ifdef CHIMES
    struct gasVariables gasVarsSave;
#endif

    int do_loop_check = 0;
    if(Gas_split>0)
    {
        N_gas += Gas_split;
        NumPart += Gas_split;
        Gas_split = 0;
        do_loop_check = 1;
    }
#ifdef GALSF
    if(Stars_converted)
    {
        N_gas -= Stars_converted;
        Stars_converted = 0;
        do_loop_check = 1;
    }
#endif
    if(NumPart <= N_gas) do_loop_check=0;
    if(N_gas <= 0) do_loop_check=0;

    /* if more gas than stars, need to be sure the block ordering is correct (gas first, then stars) */
    if(do_loop_check)
    {
        for(i = 0; i < N_gas; i++) /* loop over the gas block */
            if(P[i].Type != 0) /* and look for a particle converted to non-gas */
            {
                /* ok found a non-gas particle: */
                for(j = N_gas; j < NumPart; j++) /* loop from N_gas to Numpart, to find first labeled as gas */
                    if(P[j].Type == 0) break; /* break on that to record the j of interest */
                if(j >= NumPart) endrun(181170); /* if that j is too large, exit with error */

                psave = P[i]; /* otherwise, save the old pointer */
                P[i] = P[j]; /* now set the pointer equal to this new P[j] */
                P[j] = psave; /* now set the P[j] equal to the old, saved pointer */
                /* so we've swapped the two P[i] and P[j] */
                sphsave = SphP[i];
                SphP[i] = SphP[j];
                SphP[j] = sphsave;  /* have the gas particle take its sph pointer with it */
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(i,j);
#endif
#ifdef CHIMES /* swap chimes-specific 'gasvars' structure which is separate from SphP */
                gasVarsSave = ChimesGasVars[i]; ChimesGasVars[i] = ChimesGasVars[j]; ChimesGasVars[j] = gasVarsSave;
                /* Old particle (now at position j) is no longer a gas particle, so delete its abundance array. */
                free_gas_abundances_memory(&(ChimesGasVars[j]), &ChimesGlobalVars);
                ChimesGasVars[j].abundances = NULL; ChimesGasVars[j].isotropic_photon_density = NULL; ChimesGasVars[j].G0_parameter = NULL; ChimesGasVars[j].H2_dissocJ = NULL;
#endif /* CHIMES */

                /* ok we've now swapped the ordering so the gas particle is still inside the block */
                flag = 1;
            }
    }

    count_elim = 0;
    count_gaselim = 0;
    count_bhelim = 0;
    /* loop over entire block looking for things with zero mass, which need to be eliminated */
    for(i = 0; i < NumPart; i++)
        if(P[i].Mass <= 0)
        {
            P[i].Mass = 0;
            TimeBinCount[P[i].TimeBin]--;
            if(TimeBinActive[P[i].TimeBin]) {NumForceUpdate--;}

            if(P[i].Type == 0)
            {
                TimeBinCountSph[P[i].TimeBin]--;

                P[i] = P[N_gas - 1];
                SphP[i] = SphP[N_gas - 1];
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(i, N_gas-1);
#endif
                /* swap with properties of last gas particle (i-- below will force a check of this so its ok) */
#ifdef CHIMES
                free_gas_abundances_memory(&(ChimesGasVars[i]), &ChimesGlobalVars);
                ChimesGasVars[i] = ChimesGasVars[N_gas - 1];
                ChimesGasVars[N_gas - 1].abundances = NULL; ChimesGasVars[N_gas - 1].isotropic_photon_density = NULL; ChimesGasVars[N_gas - 1].G0_parameter = NULL; ChimesGasVars[N_gas - 1].H2_dissocJ = NULL;
#endif

                P[N_gas - 1] = P[NumPart - 1]; /* redirect the final gas pointer to go to the final particle (BH) */
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(N_gas - 1, NumPart-1);
                remove_particle_from_treewalk(NumPart - 1);
#endif
                N_gas--; /* shorten the total N_gas count */
                count_gaselim++; /* record that a BH was eliminated */
            }
            else
            {
                if(P[i].Type == 5) {count_bhelim++;} /* record elimination if BH */
                P[i] = P[NumPart - 1]; /* re-directs pointer for this particle to pointer at final particle -- so we
                                        swap the two; note that ordering -does not- matter among the non-SPH particles
                                        so its fine if this mixes up the list ordering of different particle types */
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(i, NumPart - 1);
                remove_particle_from_treewalk(NumPart - 1);
#endif
            }

            NumPart--;
            i--;
            count_elim++;
        }

    MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&count_bhelim, &tot_bhelim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(count_elim) {flag = 1;}

    if(ThisTask == 0) {if(tot_elim > 0) {printf("Rearrange: Eliminated %d/%d gas/star particles and merged away %d black holes.\n", tot_gaselim, tot_elim - tot_gaselim - tot_bhelim, tot_bhelim);}}

    All.TotNumPart -= tot_elim;
    All.TotN_gas -= tot_gaselim;
#ifdef BLACK_HOLES
    All.TotBHs -= tot_bhelim;
#endif

    MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(flag_sum) {reconstruct_timebins();}
}


/* function to apply -optional- cell excision for special cases where e.g. cells go far outside of the desired 'zoom-in region' or target region of a multi-scale simulation */
void apply_pm_hires_region_clipping_selection(int i)
{
#ifdef PM_HIRES_REGION_CLIPPING
    int clip_flag = 0; // flag for clipping
    if(All.Time <= All.TimeBegin) {return;} // no clips before run properly starts
    if(P[i].Type == 5) {return;} // no clips for sinks
    if(P[i].Type == 0 && density_isactive(i)) {if((SphP[i].Density <= 0) || (PPP[i].NumNgb <= 0)) {clip_flag=1;}} // undefined density behavior
    if(density_isactive(i)) {if(PPP[i].Hsml >= PM_HIRES_REGION_CLIPPING) {clip_flag=1;}} // far too big a kernel, outside valid domain, clip
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    if(ags_density_isactive(i)) {if(PPP[i].AGS_Hsml >= PM_HIRES_REGION_CLIPPING) {clip_flag=1;}} // far too big a kernel, outside valid domain, clip
#endif
#ifdef GALSF
    if((All.ComovingIntegrationOn) && (P[i].Type==0) && (P[i].Mass>0)) // clip material outside of a hires zoom-in region [unphysically well below cosmic mean density]
        if((SphP[i].Density>0) && (PPP[i].Hsml>0))
        {
            double rho_igm = COSMIC_BARYON_DENSITY_CGS * DMIN(1., 1000./All.cf_a3inv); /* density of IGM: cap scaling with z at z=10, so that we don't accidentally rule out very dense real stuff b/c IGM is also very dense */
            double rho_gas = DMAX( SphP[i].Density , All.DesNumNgb*P[i].Mass/(4.*M_PI/3.*PPP[i].Hsml*PPP[i].Hsml*PPP[i].Hsml) )* All.cf_a3inv * UNIT_DENSITY_IN_CGS;
            if(rho_gas < 1.e-6*rho_igm) {clip_flag=1;} // clip
        }
#endif
    if(clip_flag==1) {P[i].Mass=0;} // clip
#endif
    return; // done
}
