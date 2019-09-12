/*! \file blackhole_feed.c
 *  \brief This is where particles are marked for gas accretion.
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   on 1/9/15 by Paul Torrey (ptorrey@mit.edu) for clarity by parsing the existing code into
 *   smaller files and routines. Some communication and black hole structures were modified
 *   to reduce memory usage. Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
#include "blackhole_local.h"


void blackhole_feed_loop(void)
{
    int i, j, k, ndone_flag, ndone, ngrp, recvTask, place, nexport, nimport, dummy; MPI_Status status; MyFloat dt;
    /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) + sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackholedata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    /* Let's determine which particles may be swallowed by whom, and the weights for feedback */
    i = FirstActiveParticle;
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        
        /* do local particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])                      // DAA: can this be replaced by a loop over N_active_loc_BHs
            if(P[i].Type == 5)
                if(blackhole_feed_evaluate(i, 0, &nexport, Send_count) < 0)
                    break;
        
        MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            nimport += Recv_count[j];
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        BlackholeDataGet = (struct blackholedata_in *) mymalloc("BlackholeDataGet", nimport * sizeof(struct blackholedata_in));
        BlackholeDataIn = (struct blackholedata_in *) mymalloc("BlackholeDataIn", nexport * sizeof(struct blackholedata_in));
        
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            for(k = 0; k < 3; k++)
            {
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                BlackholeDataIn[j].Jgas_in_Kernel[k] = P[place].BH_Specific_AngMom[k];
#else
                BlackholeDataIn[j].Jgas_in_Kernel[k] = BlackholeTempInfo[P[place].IndexMapToTempStruc].Jgas_in_Kernel[k];
#endif
#endif
            }
#if defined(BH_GRAVCAPTURE_GAS)
            BlackholeDataIn[j].mass_to_swallow_edd = BlackholeTempInfo[P[place].IndexMapToTempStruc].mass_to_swallow_edd;
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
            BlackholeDataIn[j].SinkRadius = PPP[place].SinkRadius;
#endif	    
#endif
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
            BlackholeDataIn[j].AGS_Hsml = PPP[place].AGS_Hsml;
#endif	    
            BlackholeDataIn[j].Mass = P[place].Mass;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) 	|| defined(BH_WIND_CONTINUOUS)
            BlackholeDataIn[j].BH_disk_hr = P[place].BH_disk_hr;
#endif
#ifdef BH_ACCRETE_NEARESTFIRST
            BlackholeDataIn[j].BH_dr_to_NearestGasNeighbor = P[place].BH_dr_to_NearestGasNeighbor;
#endif
            BlackholeDataIn[j].Density = BPP(place).DensAroundStar;
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
#ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (((integertime) 1) << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
            BlackholeDataIn[j].ID = P[place].ID;
            memcpy(BlackholeDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        
        /* exchange particle data */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE, recvTask, TAG_BH_E,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_E, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);
        BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult",nimport * sizeof(struct blackholedata_out));
        BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));
        
        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++) {blackhole_feed_evaluate(j, 1, &dummy, &dummy);}
        
        if(i < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++) /* get the result */
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_F,
                                 &BlackholeDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_F, MPI_COMM_WORLD, &status);
                }
            }
        } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++) //
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
#ifdef BH_REPOSITION_ON_POTMIN
            if(BPP(place).BH_MinPot > BlackholeDataOut[j].BH_MinPot) {BPP(place).BH_MinPot = BlackholeDataOut[j].BH_MinPot; for(k = 0; k < 3; k++) {BPP(place).BH_MinPotPos[k] = BlackholeDataOut[j].BH_MinPotPos[k];}}
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
            BlackholeTempInfo[P[place].IndexMapToTempStruc].BH_angle_weighted_kernel_sum += BlackholeDataOut[j].BH_angle_weighted_kernel_sum;
#endif
        }
        myfree(BlackholeDataOut);
        myfree(BlackholeDataResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
}







/* do loop over neighbors to get quantities for accretion */
int blackhole_feed_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, listindex = 0;
    MyIDType id;
    MyFloat *pos, *velocity, h_i, dt, mdot, rho, mass, bh_mass;
    MyFloat ags_h_i = All.ForceSoftening[5];
    double h_i2, r2, r, u, hinv, hinv3, wk, dwk, vrel, vesc, dpos[3], dvel[3], sink_radius, f_accreted; sink_radius=0; f_accreted=1;
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
    double meddington, medd_max_accretable, mass_to_swallow_edd, eddington_factor;
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
    double norm, BH_disk_hr, J_dir[3];
    double BH_angle_weighted_kernel_sum=0;
#endif
#ifdef BH_REPOSITION_ON_POTMIN
    MyFloat minpotpos[3] = { 0, 0, 0 }, minpot = BHPOTVALUEINIT;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat bh_mass_alphadisk;
#endif
#ifdef BH_ACCRETE_NEARESTFIRST
    MyFloat BH_dr_to_NearestGasNeighbor;
#endif
#if defined(BH_SWALLOWGAS)
    double w,p,mass_markedswallow,bh_mass_withdisk;
    w=0; p=0; mass_markedswallow=0; bh_mass_withdisk=0;
#endif
    
    /* these are the BH properties */
    if(mode == 0)
    {
        pos = P[target].Pos;
        rho = P[target].DensAroundStar;       // DAA: DensAroundStar is not defined in BHP->BPP...
        mdot = BPP(target).BH_Mdot;
#ifndef WAKEUP
        dt = (P[target].TimeBin ? (((integertime) 1) << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[target].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
        h_i = PPP[target].Hsml;
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
	    ags_h_i = PPP[target].AGS_Hsml;
#endif	
        mass = P[target].Mass;
        bh_mass = BPP(target).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BPP(target).BH_Mass_AlphaDisk;
#endif
        velocity = P[target].Vel;
        id = P[target].ID;
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
        for(k=0;k<3;k++) {J_dir[k] = P[target].BH_Specific_AngMom[k];}
#else
        for(k=0;k<3;k++) {J_dir[k] = BlackholeTempInfo[P[target].IndexMapToTempStruc].Jgas_in_Kernel[k];}
#endif
        BH_disk_hr = P[target].BH_disk_hr;
#endif
#ifdef BH_ACCRETE_NEARESTFIRST
        BH_dr_to_NearestGasNeighbor = P[target].BH_dr_to_NearestGasNeighbor;
#endif
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
        mass_to_swallow_edd = BlackholeTempInfo[P[target].IndexMapToTempStruc].mass_to_swallow_edd;
#endif
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
        sink_radius = P[target].SinkRadius;
#endif
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        rho = BlackholeDataGet[target].Density;
        mdot = BlackholeDataGet[target].Mdot;
        dt = BlackholeDataGet[target].Dt;
        h_i = BlackholeDataGet[target].Hsml;
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
	    ags_h_i = BlackholeDataGet[target].AGS_Hsml;
#endif	
        mass = BlackholeDataGet[target].Mass;
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
        sink_radius = BlackholeDataGet[target].SinkRadius;
#endif
        bh_mass = BlackholeDataGet[target].BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BlackholeDataGet[target].BH_Mass_AlphaDisk;
#endif
        velocity = BlackholeDataGet[target].Vel;
        id = BlackholeDataGet[target].ID;
#if defined(FLAG_NOT_IN_PUBLIC_CODE)  || defined(BH_WIND_CONTINUOUS)
        for(k=0;k<3;k++) {J_dir[k] = BlackholeDataGet[target].Jgas_in_Kernel[k];}
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
#endif
#ifdef BH_ACCRETE_NEARESTFIRST
        BH_dr_to_NearestGasNeighbor = BlackholeDataGet[target].BH_dr_to_NearestGasNeighbor;
#endif
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
        mass_to_swallow_edd = BlackholeDataGet[target].mass_to_swallow_edd;
#endif
    }
    if((mass<0)||(h_i<=0)) return -1;
    
    /* initialize variables before SPH loop is started */
    h_i2 = h_i * h_i;
    hinv = 1 / h_i;
    hinv3 = hinv * hinv * hinv;
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
    meddington = bh_eddington_mdot(bh_mass);
    medd_max_accretable = All.BlackHoleEddingtonFactor * meddington * dt;
    eddington_factor = mass_to_swallow_edd / medd_max_accretable;   /* if <1 no problem, if >1, need to not set some swallowIDs */
#endif
#if defined(BH_SWALLOWGAS)
    bh_mass_withdisk = bh_mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += bh_mass_alphadisk;
#endif
#endif
#if defined(BH_SWALLOW_SMALLTIMESTEPS)
    double dt_min_to_accrete = DMAX(0.001 * sqrt(pow(All.SofteningTable[5],3.0)/(All.G * mass)), 20.0*DMAX(All.MinSizeTimestep,All.Timebase_interval) ); //0.001 = tolerance factor for dt_min, defined in part (ii) of 2.3.5 in Hubber 2013.
#endif
#if defined(BH_WIND_KICK) && !defined(BH_GRAVCAPTURE_GAS)
    /* DAA: increase the effective mass-loading of BAL winds to reach the desired momentum flux given the outflow velocity "All.BAL_v_outflow" chosen
       --> appropriate for cosmological simulations where particles are effectively kicked from ~kpc scales
           (i.e. we need lower velocity and higher mass outflow rates compared to accretion disk scales) - */
    f_accreted = All.BAL_f_accretion;
    if((All.BlackHoleFeedbackFactor > 0) && (All.BlackHoleFeedbackFactor != 1.)) {f_accreted /= All.BlackHoleFeedbackFactor;} else {if(All.BAL_v_outflow > 0) {f_accreted = 1./(1. + fabs(1.*BH_WIND_KICK)*All.BlackHoleRadiativeEfficiency*(C/All.UnitVelocity_in_cm_per_s)/(All.BAL_v_outflow));}}
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
    norm=0; for(k=0;k<3;k++) {norm+=J_dir[k]*J_dir[k];}
    if(norm>0) {norm=1/sqrt(norm); for(k=0;k<3;k++) {J_dir[k]*=norm;}} else {J_dir[0]=J_dir[1]=0; J_dir[2]=1;}
#endif
    
    /* Now start the actual SPH computation for this BH particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    //int particles_swallowed_this_bh_this_process = 0;
    //int particles_swallowed_this_bh_this_process_max = 1;
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
    BH_angle_weighted_kernel_sum = 0;
#endif
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb = ngb_treefind_pairs_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
            if(numngb < 0) return -1;
            
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                if(P[j].Mass > 0)
                {
                    for(k=0;k<3;k++) {dpos[k] = P[j].Pos[k] - pos[k];}
#ifdef BOX_PERIODIC
                    NEAREST_XYZ(dpos[0],dpos[1],dpos[2],-1);
#endif
                    dvel[0] = P[j].Vel[0]-velocity[0]; dvel[1] = P[j].Vel[1]-velocity[1]; dvel[2] = P[j].Vel[2]-velocity[2];
#ifdef BOX_SHEARING
                    if(pos[0] - P[j].Pos[0] > +boxHalf_X) {dvel[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                    if(pos[0] - P[j].Pos[0] < -boxHalf_X) {dvel[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
                    r2=0; for(k=0;k<3;k++) {r2 += dpos[k]*dpos[k];}
                    if(r2 < h_i2 || r2 < PPP[j].Hsml*PPP[j].Hsml)
                    {
                        vrel=0; for(k=0;k<3;k++) {vrel += dvel[k]*dvel[k];}
                        r=sqrt(r2); vrel=sqrt(vrel)/All.cf_atime;  /* do this once and use below */
                        vesc=bh_vesc(j,mass,r, ags_h_i);
#if defined(BH_SWALLOW_SMALLTIMESTEPS)
                        if(P[j].Type != 5) {if(vrel<vesc) {if(P[j].dt_step*All.Timebase_interval/All.cf_hubble_a<dt_min_to_accrete) {if(P[j].SwallowID<id) {P[j].SwallowID = id;}}}} /* Bound particles with very small timestep get eaten to avoid issues. */
#endif

#ifdef BH_REPOSITION_ON_POTMIN
                        /* check if we've found a new potential minimum which is not moving too fast to 'jump' to */
                        double boundedness_function, potential_function; boundedness_function = P[j].Potential + 0.5 * vrel*vrel * All.cf_atime; potential_function = P[j].Potential;
#if (BH_REPOSITION_ON_POTMIN == 2)
                        if( boundedness_function < 0 )
                        {
                            double wt_rsoft = r / (3.*All.ForceSoftening[5]); // normalization arbitrary here, just using for convenience for function below
                            boundedness_function *= 1./(1. + wt_rsoft*wt_rsoft); // this down-weights particles which are very far away, relative to the user-defined force softening scale, which should define some 'confidence radius' of resolution around the BH particle
                        }
                        potential_function = boundedness_function; // jumps based on -most bound- particle, not just deepest potential (down-weights fast-movers)
#endif
                        if(potential_function < minpot)
#if (BH_REPOSITION_ON_POTMIN == 1)
                        if( P[j].Type == 4 && vrel <= vesc )   // DAA: only if it is a star particle & bound
#endif
#if (BH_REPOSITION_ON_POTMIN == 2)
                        if( (P[j].Type != 0) && (P[j].Type != 5) )   // allow stars or dark matter but exclude gas, it's too messy! also exclude BHs, since we don't want to over-merge them
#endif
                        {
                            minpot=potential_function; for(k=0;k<3;k++) {minpotpos[k] = P[j].Pos[k];}
                        }
#endif
			
                        /* check_for_bh_merger.  Easy.  No Edd limit, just a pos and vel criteria. */
#if !defined(BH_DEBUG_DISABLE_MERGERS)
                        if((id != P[j].ID) && (P[j].Mass > 0) && (P[j].Type == 5))	/* we may have a black hole merger */
#ifdef SINGLE_STAR_SINK_DYNAMICS
                        if((r < 1.0001*P[j].min_dist_to_bh) && (P[j].Mass < mass) && (r < PPP[j].Hsml) && (P[j].Mass < 3*All.MinMassForParticleMerger) && (r < All.ForceSoftening[5])) /* only merge away stuff that is within the softening radius, and is no more massive that a few gas particles */
#endif
                        {
                            if(id != P[j].ID) /* check its not the same bh  (DAA: this is duplicated here...) */
                            {
                                if((vrel < BH_CSND_FRAC_BH_MERGE * vesc) && (bh_check_boundedness(j,vrel,vesc,r,sink_radius)==1))
                                {
#ifndef IO_REDUCED_MODE
                                    printf("MARKING_BH_MERGER: P[j.]ID=%llu to be swallowed by id=%llu \n", (unsigned long long) P[j].ID, (unsigned long long) id);
#endif
                                    if((P[j].SwallowID == 0) && (BPP(j).BH_Mass < bh_mass)) {P[j].SwallowID = id;} // most massive BH swallows the other - simplifies analysis
                                }
                                else
                                {
#ifndef IO_REDUCED_MODE
#ifdef BH_OUTPUT_MOREINFO           // DAA: BH merger info will be saved in a separate output file
                                    printf("ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n", ThisTask, All.Time, id, P[j].ID, vrel, vesc);
#else
                                    fprintf(FdBlackHolesDetails, "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n", ThisTask, All.Time, id, P[j].ID, vrel, vesc);
#endif
#endif
                                }
                            }
                        } // if(P[j].Type == 5) //
#endif                        
                        
                        
                        /* This is a similar loop to what we already did in blackhole_environment, but here we stochastically
                         reduce GRAVCAPT events in order to (statistically) obey the eddington limit */
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVCAPTURE_NONGAS)
                        if(P[j].Type != 5)
                        {
#ifdef SINGLE_STAR_SINK_DYNAMICS
			                double eps = DMAX(P[j].Hsml/2.8, DMAX(ags_h_i/2.8, r));			    
			                if(eps*eps*eps /(P[j].Mass + mass) <= P[j].SwallowTime)
#endif
#if defined(BH_ALPHADISK_ACCRETION)
                            if(bh_mass_alphadisk < BH_ALPHADISK_ACCRETION*bh_mass)
#endif
#if defined(BH_ACCRETE_NEARESTFIRST)
                            if((P[j].Type != 0) || (r<=1.0001*BH_dr_to_NearestGasNeighbor))
#endif
                            if((vrel < vesc)) // && (particles_swallowed_this_bh_this_process < particles_swallowed_this_bh_this_process_max))
                            { /* bound */
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                                double spec_mom=0; for(k=0;k<3;k++) {spec_mom += dvel[k]*dpos[k];} // delta_x.delta_v
                                spec_mom = (r2*vrel*vrel - spec_mom*spec_mom*All.cf_a2inv); // specific angular momentum^2 = r^2(delta_v)^2 - (delta_v.delta_x)^2;
				                if(spec_mom < All.G * (mass + P[j].Mass) * sink_radius)  // check Bate 1995 angular momentum criterion (in addition to bounded-ness)
#endif
                                if( bh_check_boundedness(j,vrel,vesc,r,sink_radius)==1 ) { /* apocenter within target distance */
#ifdef BH_GRAVCAPTURE_NONGAS        /* simply swallow non-gas particle if BH_GRAVCAPTURE_NONGAS enabled */
                                    if((P[j].Type != 0) && (P[j].SwallowID < id)) P[j].SwallowID = id;
#endif
#if defined(BH_GRAVCAPTURE_GAS)
                                    /* now deal with gas */
                                    if (P[j].Type == 0){
#if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION) /* if Eddington-limited and NO alpha-disk, do this stochastically */
                                        p = 1. / eddington_factor;
#if defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
                                        p /= All.BAL_f_accretion; // we need to accrete more, then remove the mass in winds
#endif
                                        w = get_random_number(P[j].ID);
                                        if(w < p)
                                        {
#ifndef IO_REDUCED_MODE
                                            printf("MARKING_BH_FOOD: P[j.]ID=%llu to be swallowed by id=%llu \n", (unsigned long long) P[j].ID, (unsigned long long) id);
#endif
                                            if(P[j].SwallowID < id) P[j].SwallowID = id;
                                        }
#else //if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                        if(P[j].SwallowID < id) {P[j].SwallowID = id;} /* in other cases, just swallow the particle */  //particles_swallowed_this_bh_this_process++;
#endif //else defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                    } //if (P[j].Type == 0)
#endif //ifdef BH_GRAVCAPTURE_GAS
                                } // if( apocenter in tolerance range )
                            } // if(vrel < vesc)
                        } //if(P[j].Type != 5)
#endif // defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVCAPTURE_NONGAS)
                        
                        
                        
                        /* now is the more standard accretion only of gas, according to the mdot calculated before */
                        if(P[j].Type == 0) /* here we have a gas particle */
                        {
                            u=r*hinv; kernel_main(u,hinv3,hinv*hinv3,&wk,&dwk,-1);
#if defined(BH_SWALLOWGAS) && !defined(BH_GRAVCAPTURE_GAS) /* compute accretion probability, this below is only meaningful if !defined(BH_GRAVCAPTURE_GAS)... */
                            double dm_toacc = bh_mass_withdisk - (mass + mass_markedswallow); if(dm_toacc>0) {p=dm_toacc*wk/rho;} else {p=0;}
#ifdef BH_WIND_KICK /* DAA: for stochastic winds (BH_WIND_KICK) we remove a fraction of mass from gas particles prior to kicking --> need to increase the probability here to balance black hole growth */
                            if(f_accreted>0) {p /= f_accreted; if((bh_mass_withdisk - mass) < 0) {p = ( (1-f_accreted)/f_accreted ) * mdot * dt * wk / rho;}} /* DAA: compute outflow probability when "bh_mass_withdisk < mass" - we don't need to enforce mass conservation in this case, relevant only in low-res sims where the BH seed mass is much lower than the gas particle mass */
#endif
#ifdef BH_ACCRETE_NEARESTFIRST /* put all the weight on the single nearest gas particle, instead of spreading it in a kernel-weighted fashion */
                            p=0; if(dm_toacc>0 && P[j].Mass>0 && r<1.0001*BH_dr_to_NearestGasNeighbor) {p=dm_toacc/P[j].Mass;}
#endif
                            w = get_random_number(P[j].ID);
                            if(w < p)
                            {
#ifndef IO_REDUCED_MODE
                                printf("MARKING_BH_FOOD: j %d w %g p %g TO_BE_SWALLOWED \n",j,w,p);
#endif
                                if(P[j].SwallowID < id) {P[j].SwallowID = id; mass_markedswallow += P[j].Mass*f_accreted;}
                            } // if(w < p)
#endif // BH_SWALLOWGAS

#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS) /* calculate the angle-weighting for the photon momentum */
                            if((mdot>0)&&(dt>0)&&(r>0)&&(P[j].SwallowID==0)&&(P[j].Mass>0)&&(P[j].Type==0))
                            { /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                                norm=0; for(k=0;k<3;k++) {norm+=(dpos[k]/r)*J_dir[k];}
                                BH_angle_weighted_kernel_sum += bh_angleweight_localcoupling(j,BH_disk_hr,norm,r,h_i);
                            }
#endif

#ifdef BH_THERMALFEEDBACK
                            double energy = bh_lum_bol(mdot, bh_mass, -1) * dt;
                            if(rho > 0) {SphP[j].Injected_BH_Energy += (wk/rho) * energy * P[j].Mass;}
#endif                            
                        } // if(P[j].Type == 0)

                    } // if(r2 < h_i2)
                } // if(P[j].Mass > 0)
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}	/* open it */
            }
        } // mode==1
    } // while(startnode >= 0) (outer of the double-loop)
    
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
        BlackholeTempInfo[P[target].IndexMapToTempStruc].BH_angle_weighted_kernel_sum += BH_angle_weighted_kernel_sum;  /* need to correct target index */
#endif
#ifdef BH_REPOSITION_ON_POTMIN
        BPP(target).BH_MinPot = minpot; for(k = 0; k < 3; k++) {BPP(target).BH_MinPotPos[k] = minpotpos[k];}
#endif
    }
    else
    {
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
        BlackholeDataResult[target].BH_angle_weighted_kernel_sum = BH_angle_weighted_kernel_sum;
#endif
#ifdef BH_REPOSITION_ON_POTMIN
        BlackholeDataResult[target].BH_MinPot = minpot; for(k = 0; k < 3; k++) {BlackholeDataResult[target].BH_MinPotPos[k] = minpotpos[k];}
#endif
    }
    return 0;
} /* closes bh_evaluate routine */
