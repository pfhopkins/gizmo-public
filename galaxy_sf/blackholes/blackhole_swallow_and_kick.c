/*! \file blackhole_swallow_and_kick.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
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

#ifdef BH_ALPHADISK_ACCRETION
#define accreted_BH_mass_alphaornot accreted_BH_mass_alphadisk
#else
#define accreted_BH_mass_alphaornot accreted_BH_mass
#endif


static int N_gas_swallowed, N_star_swallowed, N_dm_swallowed, N_BH_swallowed;

void blackhole_swallow_and_kick_loop(void)
{
    int i, j, k;
    int ndone_flag, ndone;
    int ngrp, recvTask, place, nexport, nimport, dummy;
    MPI_Status status;
    
    int Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed;
    
    /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) + sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackholedata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    N_gas_swallowed = N_star_swallowed = N_dm_swallowed = N_BH_swallowed = 0;
    Ntot_gas_swallowed = Ntot_star_swallowed = Ntot_dm_swallowed = Ntot_BH_swallowed = 0;
    
    i = FirstActiveParticle;	/* first particle for this task */
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 5)
                if(P[i].SwallowID == 0)     /* this particle not being swallowed */
                    if(blackhole_swallow_and_kick_evaluate(i, 0, &nexport, Send_count) < 0)
                        break;
        
        qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
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
        
        
        /* populate the struct to be exported */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            for(k = 0; k < 3; k++)
            {
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                BlackholeDataIn[j].Jgas_in_Kernel[k] = P[place].BH_Specific_AngMom[k];
#else
                BlackholeDataIn[j].Jgas_in_Kernel[k] = BlackholeTempInfo[P[place].IndexMapToTempStruc].Jgas_in_Kernel[k];
#endif
#endif
            }
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].ID = P[place].ID;
            BlackholeDataIn[j].Mass = P[place].Mass;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
            BlackholeDataIn[j].BH_disk_hr = P[place].BH_disk_hr;
            BlackholeDataIn[j].BH_angle_weighted_kernel_sum = BlackholeTempInfo[P[place].IndexMapToTempStruc].BH_angle_weighted_kernel_sum;
#endif
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
#ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (((integertime) 1) << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
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
                                 Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_G,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_G, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);
        BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult", nimport * sizeof(struct blackholedata_out));
        BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));
        
        /* do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_swallow_and_kick_evaluate(j, 1, &dummy, &dummy);  /* set BlackholeDataResult based on BlackholeDataGet */
        
        if(i < 0)
            ndone_flag = 1;
        else
            ndone_flag = 0;
        
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /* get the result */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_H,
                                 &BlackholeDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_H, MPI_COMM_WORLD, &status);
                }
            }
        }
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_Mass += BlackholeDataOut[j].accreted_Mass;
            BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_BH_Mass += BlackholeDataOut[j].accreted_BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_BH_mass_alphadisk += BlackholeDataOut[j].accreted_BH_mass_alphadisk;
#endif
            for(k = 0; k < 3; k++)
            {
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_momentum[k] += BlackholeDataOut[j].accreted_momentum[k];
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_centerofmass[k] += BlackholeDataOut[j].accreted_centerofmass[k];
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_J[k] += BlackholeDataOut[j].accreted_J[k];
#endif		
            }
#ifdef BH_COUNTPROGS
            BPP(place).BH_CountProgs += BlackholeDataOut[j].BH_CountProgs;
#endif
#ifdef GALSF
            if(P[place].StellarAge > BlackholeDataOut[j].Accreted_Age)
                P[place].StellarAge = BlackholeDataOut[j].Accreted_Age;
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
    
    
    MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_star_swallowed, &Ntot_star_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_dm_swallowed, &Ntot_dm_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if((ThisTask == 0)&&(Ntot_gas_swallowed+Ntot_star_swallowed+Ntot_dm_swallowed+Ntot_BH_swallowed>0))
    {
        printf("Accretion done: swallowed %d gas, %d star, %d dm, and %d BH particles\n",
               Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed);
    }
    
}




int blackhole_swallow_and_kick_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, bin, listindex = 0;
    MyIDType id;
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
    MyLongDouble accreted_momentum[3]={0};
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
    MyLongDouble accreted_centerofmass[3]={0};
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
    MyLongDouble accreted_J[3]={0};
#endif
    MyFloat *pos, *vel, h_i, bh_mass;
#if (defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK))
    MyFloat hinv, hinv3;
#endif
    MyFloat f_accreted=0;
#if defined(BH_WIND_KICK) || defined(BH_OUTPUT_GASSWALLOW)
    MyFloat mass;
#ifdef BH_WIND_KICK
    MyFloat v_kick=0;
    MyFloat bh_mass_withdisk;
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat bh_mass_alphadisk;     // DAA: we need bh_mass_alphadisk for BH_WIND_KICK winds below
#endif
#endif
#endif
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS))
    MyFloat mdot,dt;
#endif
    
    MyFloat dir[3], norm, mom;
    mom=0; norm=0; dir[0]=0;
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK) 
    MyFloat J_dir[3];
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
    double BH_angle_weighted_kernel_sum, mom_wt;
    MyFloat BH_disk_hr,kernel_zero,dwk;
    kernel_main(0.0,1.0,1.0,&kernel_zero,&dwk,-1);
#endif
#ifdef GALSF
    double accreted_age = 1;
#endif
    
    int mod_index = 0;
    
    if(mode == 0)
    {
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = PPP[target].Hsml;
        id = P[target].ID;
#if defined(BH_WIND_KICK) || defined(BH_OUTPUT_GASSWALLOW)
        mass = P[target].Mass;    
#endif
#if defined(BH_ALPHADISK_ACCRETION) && defined(BH_WIND_KICK)
        bh_mass_alphadisk = BPP(target).BH_Mass_AlphaDisk;
#endif
        bh_mass = BPP(target).BH_Mass;
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS))
        mdot = BPP(target).BH_Mdot;
#ifndef WAKEUP
        dt = (P[target].TimeBin ? (((integertime) 1) << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[target].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
        for(k=0;k<3;k++) {J_dir[k] = P[target].BH_Specific_AngMom[k];}
#else
        for(k=0;k<3;k++) {J_dir[k] = BlackholeTempInfo[P[target].IndexMapToTempStruc].Jgas_in_Kernel[k];}
#endif
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
        BH_disk_hr = P[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeTempInfo[P[target].IndexMapToTempStruc].BH_angle_weighted_kernel_sum;
#endif
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        vel = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
#if defined(BH_WIND_KICK) || defined(BH_OUTPUT_GASSWALLOW)
        mass = BlackholeDataGet[target].Mass;
#if defined(BH_ALPHADISK_ACCRETION) && defined(BH_WIND_KICK)
        bh_mass_alphadisk = BlackholeDataGet[target].BH_Mass_AlphaDisk;      
#endif
#endif
        bh_mass = BlackholeDataGet[target].BH_Mass;
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
        mdot = BlackholeDataGet[target].Mdot;
        dt = BlackholeDataGet[target].Dt;
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
        for(k=0;k<3;k++) {J_dir[k] = BlackholeDataGet[target].Jgas_in_Kernel[k];}
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeDataGet[target].BH_angle_weighted_kernel_sum;
#endif
    }

#ifdef BH_WIND_KICK
    bh_mass_withdisk = bh_mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += bh_mass_alphadisk;
#endif
#endif
    MyLongDouble accreted_mass = 0, accreted_BH_mass = 0, accreted_BH_mass_alphadisk = 0;
#ifdef BH_COUNTPROGS
    int accreted_BH_progs = 0;
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
    mom = bh_lum_bol(mdot, bh_mass, -1) * dt / (C / All.UnitVelocity_in_cm_per_s); mom_wt = 0;
#endif
#if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)
    hinv=h_i; hinv3=hinv*hinv*hinv;
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
    norm=0; for(k=0;k<3;k++) {norm+=J_dir[k]*J_dir[k];}
    if(norm>0) {norm=1/sqrt(norm); for(k=0;k<3;k++) {J_dir[k]*=norm;}} else {J_dir[0]=J_dir[1]=0; J_dir[2]=1;}
#endif
    
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb = ngb_treefind_pairs_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
            if(numngb < 0) return -1;
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n]; MyIDType OriginallyMarkedSwallowID = P[j].SwallowID; // record this to help prevent double-counting below
                double dP[3]={0},dv[3]={0}; for(k=0;k<3;k++) {dP[k]=P[j].Pos[k]-pos[k];}
#ifdef BOX_PERIODIC
                NEAREST_XYZ(dP[0],dP[1],dP[2],-1); /*  find the closest image in the given box size  */
#endif
                dv[0] = P[j].Vel[0]-vel[0]; dv[1] = P[j].Vel[1]-vel[1]; dv[2] = P[j].Vel[2]-vel[2];
#ifdef BOX_SHEARING
                if(pos[0] - P[j].Pos[0] > +boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                if(pos[0] - P[j].Pos[0] < -boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif

                
                
                /* we've found a particle to be swallowed.  This could be a BH merger, DM particle, or baryon w/ feedback */
                if(P[j].SwallowID == id && P[j].Mass > 0)
                {   /* accreted quantities to be added [regardless of particle type] */
                    f_accreted = 1; /* default to accreting entire particle */
#ifdef BH_WIND_KICK
                    if(P[j].Type == 0)
                    {
                        f_accreted = All.BAL_f_accretion; /* if particle is gas, only a fraction gets accreted in these particular modules */
#ifndef BH_GRAVCAPTURE_GAS
                        if((All.BlackHoleFeedbackFactor > 0) && (All.BlackHoleFeedbackFactor != 1.)) {f_accreted /= All.BlackHoleFeedbackFactor;} else {if(All.BAL_v_outflow > 0) f_accreted = 1./(1. + fabs(1.*BH_WIND_KICK)*All.BlackHoleRadiativeEfficiency*(C/All.UnitVelocity_in_cm_per_s)/(All.BAL_v_outflow));}
                        if((bh_mass_withdisk - mass) <= 0) {f_accreted=0;} // DAA: no need to accrete gas particle to enforce mass conservation (we will simply kick),  note that here the particle mass P.Mass is larger than the physical BH mass P.BH_Mass
#endif
                    }
#endif

                    /* handle accretion/conservation of certain conserved quantities, depending on whether we are intending our sub-grid model to follow them */
                    double mcount_for_conserve = f_accreted * P[j].Mass;
#if (BH_FOLLOW_ACCRETED_ANGMOM == 1) /* in this case we are only counting this if its coming from BH particles */
                    if(P[j].Type!=5) {mcount_for_conserve=0;} else {mcount_for_conserve=BPP(j).BH_Mass;}
#ifdef BH_ALPHADISK_ACCRETION
                    if(P[j].Type==5) {mcount_for_conserve += BPP(j).BH_Mass_AlphaDisk;}
#endif
#endif
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
                    for(k=0;k<3;k++) {accreted_momentum[k] += FLT( mcount_for_conserve * dv[k]);}
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
                    for(k=0;k<3;k++) {accreted_centerofmass[k] += FLT(mcount_for_conserve * dP[k]);}
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                    accreted_J[0] += FLT(mcount_for_conserve * ( dP[1]*dv[2] - dP[2]*dv[1] ));
                    accreted_J[1] += FLT(mcount_for_conserve * ( dP[2]*dv[0] - dP[0]*dv[2] ));
                    accreted_J[2] += FLT(mcount_for_conserve * ( dP[0]*dv[1] - dP[1]*dv[0] ));
                    if(P[j].Type==5) {for(k=0;k<3;k++) {accreted_J[k] += FLT(mcount_for_conserve * BPP(j).BH_Specific_AngMom[k]);}}
#endif

                    

                    if(P[j].Type == 5)  /* this is a BH-BH merger */
                    {
#ifdef BH_OUTPUT_MOREINFO
                        fprintf(FdBhMergerDetails,"%g  %u %g %2.7f %2.7f %2.7f  %u %g %2.7f %2.7f %2.7f\n", All.Time,  id,bh_mass,pos[0],pos[1],pos[2],  P[j].ID,BPP(j).BH_Mass,P[j].Pos[0],P[j].Pos[1],P[j].Pos[2]);
#else
#ifndef IO_REDUCED_MODE
                        fprintf(FdBlackHolesDetails,"ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n", ThisTask, All.Time, id, P[j].ID, bh_mass, BPP(j).BH_Mass);
#endif
#endif
#ifdef BH_INCREASE_DYNAMIC_MASS
                        /* the true dynamical mass of the merging BH is P[j].Mass/BH_INCREASE_DYNAMIC_MASS unless exceeded by physical growth
                         - in the limit BPP(j).BH_Mass > BH_INCREASE_DYNAMIC_MASS x m_b, then bh_mass=P[j].Mass on average and we are good as well  */
                        accreted_mass    += FLT( DMAX(BPP(j).BH_Mass, P[j].Mass/BH_INCREASE_DYNAMIC_MASS) );
#else
                        accreted_mass    += FLT(P[j].Mass);
#endif
                        accreted_BH_mass += FLT(BPP(j).BH_Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        accreted_BH_mass_alphadisk += FLT(BPP(j).BH_Mass_AlphaDisk);
#endif
#ifdef BH_WIND_SPAWN
                        accreted_BH_mass_alphaornot += FLT(BPP(j).unspawned_wind_mass);
#endif
#ifdef BH_COUNTPROGS
                        accreted_BH_progs += BPP(j).BH_CountProgs;
#endif
                        bin = P[j].TimeBin; TimeBin_BH_mass[bin] -= BPP(j).BH_Mass; TimeBin_BH_dynamicalmass[bin] -= P[j].Mass; TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
                        if(BPP(j).BH_Mass > 0) {TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;}
                        P[j].Mass = 0; BPP(j).BH_Mass = 0; BPP(j).BH_Mdot = 0;
#ifdef GALSF
                        accreted_age = P[j].StellarAge;
#endif
                        N_BH_swallowed++;
                    } // if(P[j].Type == 5) -- BH + BH merger


#ifdef BH_GRAVCAPTURE_NONGAS /* DM and star particles can only be accreted ifdef BH_GRAVCAPTURE_NONGAS */
                    if((P[j].Type == 1) || (All.ComovingIntegrationOn && (P[j].Type==2||P[j].Type==3)) )
                    {   /* this is a DM particle: In this case, no kick, so just zero out the mass and 'get rid of' the particle (preferably by putting it somewhere irrelevant) */
#ifndef IO_REDUCED_MODE
                        printf("BH_swallow_DM: j %d Type(j) %d  M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g \n", j,P[j].Type,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2]);
#endif
                        accreted_mass += FLT(P[j].Mass); accreted_BH_mass += FLT(P[j].Mass); P[j].Mass = 0;
                        N_dm_swallowed++;
                    }
                    if((P[j].Type==4) || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn) ))
                    {   /* this is a star particle: If there is an alpha-disk, we let them go to the disk. If there is no alpha-disk, stars go to the BH directly and won't affect feedback. (Can be simply modified if we need something different.) */
                        accreted_mass += FLT(P[j].Mass); accreted_BH_mass_alphaornot += FLT(P[j].Mass); P[j].Mass = 0;
                        N_star_swallowed++;
                    }
#endif // #ifdef BH_GRAVCAPTURE_NONGAS -- BH + DM or Star merger


                    /* this is a gas particle: DAA: we need to see if the gas particle has to be accreted in full or not, depending on BH_WIND_KICK
                     the only difference with BH_ALPHADISK_ACCRETION should be that the mass goes first to the alphadisk */
                    if(P[j].Type == 0)                    
                    {
                        accreted_mass += FLT(f_accreted*P[j].Mass);
#ifdef BH_GRAVCAPTURE_GAS
                        accreted_BH_mass_alphaornot += FLT(f_accreted*P[j].Mass);
#endif
                        P[j].Mass *= (1-f_accreted);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue *= (1-f_accreted);
#endif

#ifdef BH_WIND_KICK     /* BAL kicking operations. NOTE: we have two separate BAL wind models, particle kicking and smooth wind model. This is where we do the particle kicking BAL model. This should also work when there is alpha-disk. */
                        v_kick=All.BAL_v_outflow; for(k=0;k<3;k++) {dir[k]=dP[k];} // DAA: default direction is radially outwards
#if defined(BH_COSMIC_RAYS) /* inject cosmic rays alongside wind injection */
                        double dEcr = All.BH_CosmicRay_Injection_Efficiency * P[j].Mass * (All.BAL_f_accretion/(1.-All.BAL_f_accretion)) * (C / All.UnitVelocity_in_cm_per_s)*(C / All.UnitVelocity_in_cm_per_s); SphP[j].CosmicRayEnergy+=dEcr; SphP[j].CosmicRayEnergyPred+=dEcr;
#ifdef COSMIC_RAYS_M1
                        dEcr*=COSMIC_RAYS_M1; for(k=0;k<3;k++) {SphP[j].CosmicRayFlux[k]+=dEcr*dir[k]; SphP[j].CosmicRayFluxPred[k]+=dEcr*dir[k];}
#endif
#endif
#if (BH_WIND_KICK < 0)  /* DAA: along polar axis defined by angular momentum within Kernel (we could add finite opening angle) work out the geometry w/r to the plane of the disk */
                        if((dir[0]*J_dir[0] + dir[1]*J_dir[1] + dir[2]*J_dir[2]) > 0){for(k=0;k<3;k++) {dir[k]=J_dir[k];}} else {for(k=0;k<3;k++) {dir[k]=-J_dir[k];}}
#endif
                        for(k=0,norm=0;k<3;k++) {norm+=dir[k]*dir[k];} if(norm<=0) {dir[0]=0;dir[1]=0;dir[2]=1;norm=1;} else {norm=sqrt(norm); dir[0]/=norm;dir[1]/=norm;dir[2]/=norm;}
                        for(k=0;k<3;k++) {P[j].Vel[k]+=v_kick*All.cf_atime*dir[k]; SphP[j].VelPred[k]+=v_kick*All.cf_atime*dir[k];}
#ifdef GALSF_SUBGRID_WINDS // if sub-grid galactic winds are decoupled from the hydro, we decouple the BH kick winds as well
                        SphP[j].DelayTime = All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
#endif
#ifndef IO_REDUCED_MODE
                        printf("BAL kick: P[j].ID %llu ID %llu Type(j) %d f_acc %g M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g v_out %g \n",(unsigned long long) P[j].ID, (unsigned long long) P[j].SwallowID,P[j].Type, All.BAL_f_accretion,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2],v_kick);
#endif
#ifdef BH_OUTPUT_MOREINFO
                        fprintf(FdBhWindDetails,"%g  %u %g  %2.7f %2.7f %2.7f  %2.7f %2.7f %2.7f  %g %g %g  %u  %2.7f %2.7f %2.7f\n",All.Time, P[j].ID, P[j].Mass,  P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],  P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],dir[0]/norm,dir[1]/norm,dir[2]/norm, id, pos[0],pos[1],pos[2]);
#endif
#endif // #ifdef BH_WIND_KICK
                        N_gas_swallowed++;
#ifdef BH_OUTPUT_GASSWALLOW
                        MyDouble tempB[3]={0,0,0};
#ifdef MAGNETIC
                        tempB[0]=SphP[j].B[0];tempB[1]=SphP[j].B[1];tempB[2]=SphP[j].B[2]; //use particle magnetic field
#endif
                        fprintf(FdBhSwallowDetails,"%g  %u %g %2.7f %2.7f %2.7f  %u %g %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f\n", All.Time,  id,mass,pos[0],pos[1],pos[2],  P[j].ID, P[j].Mass, (P[j].Pos[0]-pos[0]),(P[j].Pos[1]-pos[1]),(P[j].Pos[2]-pos[2]), (P[j].Vel[0]-vel[0]),(P[j].Vel[1]-vel[1]),(P[j].Vel[2]-vel[2]), SphP[j].InternalEnergy, tempB[0], tempB[1], tempB[2]);
#endif
                    }  // if(P[j].Type == 0)

                    P[j].SwallowID = 0; /* DAA: make sure it is not accreted (or ejected) by the same BH again if inactive in the next timestep */
                } // if(P[j].SwallowID == id)  -- particles being entirely or partially swallowed!!!

#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)                
                /* now, do any other feedback "kick" operations (which used the previous loops to calculate weights) */
                if(mom>0 && mdot>0 && dt>0 && OriginallyMarkedSwallowID==0 && P[j].SwallowID==0 && P[j].Mass>0 && P[j].Type==0) // particles NOT being swallowed!
                {
                    double r=0; for(k=0;k<3;k++) {dir[k]=dP[k]; r+=dir[k]*dir[k];} // should be away from BH
                    if(r>0)
                    {
                        r=sqrt(r); for(k=0;k<3;k++) {dir[k]/=r;} /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                        for(norm=0,k=0;k<3;k++) {norm+=dir[k]*J_dir[k];}
                        mom_wt = bh_angleweight_localcoupling(j,BH_disk_hr,norm,r,h_i) / BH_angle_weighted_kernel_sum;
                        if(BH_angle_weighted_kernel_sum<=0) mom_wt=0;
                                
#if defined(BH_COSMIC_RAYS) && defined(BH_WIND_CONTINUOUS) /* inject cosmic rays alongside continuous wind injection */
                        double dEcr = All.BH_CosmicRay_Injection_Efficiency * mom_wt * (C / All.UnitVelocity_in_cm_per_s)*(C / All.UnitVelocity_in_cm_per_s) * mdot*dt;
                        SphP[j].CosmicRayEnergy+=dEcr; SphP[j].CosmicRayEnergyPred+=dEcr;
#ifdef COSMIC_RAYS_M1
                        dEcr*=COSMIC_RAYS_M1; for(k=0;k<3;k++) {SphP[j].CosmicRayFlux[k]+=dEcr*dir[k]; SphP[j].CosmicRayFluxPred[k]+=dEcr*dir[k];}
#endif
#endif
#if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK) /* inject BAL winds, this is the more standard smooth feedback model */
                        double m_wind = mom_wt * (1-All.BAL_f_accretion)/(All.BAL_f_accretion) * mdot*dt; /* mass to couple */
                        if(BH_angle_weighted_kernel_sum<=0) m_wind=0;
                        //1. check if (Vw-V0)*rhat <= 0   [ equivalently, check if   |Vw| <= V0*rhat ]
                        //2. if (1) is False, the wind will catch the particle, couple mass, momentum, energy, according to the equations above
                        //3. if (1) is True, the wind will not catch the particle, or will only asymptotically catch it. For the sake of mass conservation in the disk, I think it is easiest to treat this like the 'marginal' case where the wind barely catches the particle. In this case, add the mass normally, but no momentum, and no energy, giving:
                        //dm = m_wind, dV = 0, du = -mu*u0   [decrease the thermal energy slightly to account for adding more 'cold' material to it]
                        double dvr_gas_to_bh, dr_gas_to_bh;
                        for(dvr_gas_to_bh=dr_gas_to_bh=0, k=0;k<3;k++) {dvr_gas_to_bh += dv[k]*dP[k]; dr_gas_to_bh  += dP[k]*dP[k];}
                        dvr_gas_to_bh /= dr_gas_to_bh ;
                        
                        /* add wind mass to particle, correcting density as needed */
                        if(P[j].Hsml<=0)
                        {
                            if(SphP[j].Density>0){SphP[j].Density*=(1+m_wind/P[j].Mass);} else {SphP[j].Density=m_wind*hinv3;}
                        } else {
                            SphP[j].Density += kernel_zero * m_wind/(P[j].Hsml*P[j].Hsml*P[j].Hsml);
                        }
                        P[j].Mass += m_wind;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue += m_wind;
#endif
                        /* now add wind momentum to particle */
                        if(dvr_gas_to_bh < All.BAL_v_outflow)   // gas moving away from BH at v < BAL speed
                        {
                            double e_wind = 0;
                            for(k=0;k<3;k++)
                            {
                                norm = All.cf_atime*All.BAL_v_outflow*dir[k] - dv[k]; // relative wind-particle velocity (in code units) including BH-particle motion;
                                P[j].Vel[k] += All.BlackHoleFeedbackFactor * norm * m_wind/P[j].Mass; // momentum conservation gives updated velocity
                                SphP[j].VelPred[k] += All.BlackHoleFeedbackFactor * norm * m_wind/P[j].Mass;
                                e_wind += (norm/All.cf_atime)*(norm/All.cf_atime); // -specific- shocked wind energy
                            }
                            e_wind *= 0.5*m_wind/P[j].Mass; // make total wind energy, add to particle as specific energy of -particle-
                            SphP[j].InternalEnergy += e_wind; SphP[j].InternalEnergyPred += e_wind;
                        } else {    // gas moving away from BH at wind speed (or faster) already.
                            if(SphP[j].InternalEnergy * ( P[j].Mass - m_wind ) / P[j].Mass > 0) {SphP[j].InternalEnergy = SphP[j].InternalEnergy * ( P[j].Mass - m_wind ) / P[j].Mass;}
                        }
#endif // if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)
                    } // r > 0
                } // (check if valid gas neighbor of interest)
#endif // defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)                
                

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
        }
    } // while(startnode >= 0)
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
        BlackholeTempInfo[mod_index].accreted_Mass = accreted_mass;
        BlackholeTempInfo[mod_index].accreted_BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BlackholeTempInfo[mod_index].accreted_BH_mass_alphadisk = accreted_BH_mass_alphadisk;
#endif
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
        for(k=0;k<3;k++) {BlackholeTempInfo[mod_index].accreted_momentum[k] = accreted_momentum[k];}
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
        for(k=0;k<3;k++) {BlackholeTempInfo[mod_index].accreted_centerofmass[k] = accreted_centerofmass[k];}
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
        for(k=0;k<3;k++) {BlackholeTempInfo[mod_index].accreted_J[k] = accreted_J[k];}
#endif
#ifdef BH_COUNTPROGS
        BPP(target).BH_CountProgs += accreted_BH_progs;
#endif
#ifdef GALSF
        if(P[target].StellarAge > accreted_age) {P[target].StellarAge = accreted_age;}
#endif
    }
    else
    {
        BlackholeDataResult[target].accreted_Mass = accreted_mass;
        BlackholeDataResult[target].accreted_BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BlackholeDataResult[target].accreted_BH_mass_alphadisk = accreted_BH_mass_alphadisk;
#endif
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
        for(k=0;k<3;k++) {BlackholeDataResult[target].accreted_momentum[k] = accreted_momentum[k];}
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
        for(k=0;k<3;k++) {BlackholeDataResult[target].accreted_centerofmass[k] = accreted_centerofmass[k];}
#endif	    
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
        for(k=0;k<3;k++) {BlackholeDataResult[target].accreted_J[k] = accreted_J[k];}
#endif
#ifdef BH_COUNTPROGS
        BlackholeDataResult[target].BH_CountProgs = accreted_BH_progs;
#endif
#ifdef GALSF
        BlackholeDataResult[target].Accreted_Age = accreted_age;
#endif
    }
    
    return 0;
} /* closes bh_evaluate_swallow */



#ifdef BH_WIND_SPAWN
void spawn_bh_wind_feedback(void)
{
    int i, n_particles_split = 0, MPI_n_particles_split, dummy_gas_tag=0;
    for(i = 0; i < NumPart; i++)
        if(P[i].Type==0)
        {
            dummy_gas_tag=i;
            break;
        }
    
    /* don't loop or go forward if there are no gas particles in the domain, or the code will crash */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        //long nmax = (int)(0.9*All.MaxPart); if(All.MaxPart-1000 < nmax) nmax=All.MaxPart-1000; /* stricter criterion for allowing spawns, more relaxed below */
        //if((NumPart+n_particles_split+(int)(2.*(BH_WIND_SPAWN+0.1)) < nmax) && (n_particles_split<1) && (P[i].Type==5))

        long nmax = (int)(0.99*All.MaxPart); if(All.MaxPart-20 < nmax) nmax=All.MaxPart-20;
        if((NumPart+n_particles_split+(int)(2.*(BH_WIND_SPAWN+0.1)) < nmax) && (P[i].Type==5))
        {
            if(BPP(i).unspawned_wind_mass >= (BH_WIND_SPAWN)*All.BAL_wind_particle_mass)
            {
                int j; dummy_gas_tag=-1; double r2=MAX_REAL_NUMBER;
                for(j=0; j<N_gas; j++) /* find the closest gas particle on the domain to act as the dummy */
                {
                    if(P[j].Type==0)
                    {
                        double dx2=(P[j].Pos[0]-P[i].Pos[0])*(P[j].Pos[0]-P[i].Pos[0]) + (P[j].Pos[1]-P[i].Pos[1])*(P[j].Pos[1]-P[i].Pos[1]) + (P[j].Pos[2]-P[i].Pos[2])*(P[j].Pos[2]-P[i].Pos[2]);
                        if(dx2 < r2) {r2=dx2; dummy_gas_tag=j;}
                    }
                }
                if(dummy_gas_tag >= 0)
                {
                    n_particles_split += blackhole_spawn_particle_wind_shell( i , dummy_gas_tag, n_particles_split);
                }
            }
        }
    }
    MPI_Allreduce(&n_particles_split, &MPI_n_particles_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(MPI_n_particles_split>0){TreeReconstructFlag = 1;}
    if(ThisTask == 0) {printf("Particle BH spawn check: %d particles spawned \n", MPI_n_particles_split);}

    /* rearrange_particle_sequence -must- be called immediately after this routine! */
    All.TotNumPart += (long long)MPI_n_particles_split;
    All.TotN_gas   += (long long)MPI_n_particles_split;
    Gas_split       = n_particles_split;                    // specific to the local processor //
    
    //rearrange_particle_sequence();
}




/*! this code copies what was used in merge_split.c for the gas particle split case */
int blackhole_spawn_particle_wind_shell( int i, int dummy_sph_i_to_clone, int num_already_spawned )
{
    double total_mass_in_winds = BPP(i).unspawned_wind_mass;
    int n_particles_split   = floor( total_mass_in_winds / All.BAL_wind_particle_mass );
    if( (n_particles_split == 0) || (n_particles_split < 1) ) {return 0;}
    int n0max = DMAX(20 , (int)(3.*(BH_WIND_SPAWN)+0.1));
    if(n_particles_split > n0max) {n_particles_split = n0max;}

    /* here is where the details of the split are coded, the rest is bookkeeping */
    //double mass_of_new_particle = total_mass_in_winds / n_particles_split; /* don't do this, as can produce particles with extremely large masses; instead wait to spawn */
    double mass_of_new_particle = All.BAL_wind_particle_mass;
#ifndef IO_REDUCED_MODE
    printf("Task %d wants to create %g mass in wind with %d new particles each of mass %g \n", ThisTask,total_mass_in_winds, n_particles_split, mass_of_new_particle);
    printf(" splitting BH %d using SphP particle %d\n", i, dummy_sph_i_to_clone);
#endif
    int k=0; long j;
    if(NumPart + num_already_spawned + n_particles_split >= All.MaxPart)
    {
        printf ("On Task=%d with NumPart=%d (+N_spawned=%d) we tried to split a particle, but there is no space left...(All.MaxPart=%d). Try using more nodes, or raising PartAllocFac, or changing the split conditions to avoid this.\n", ThisTask, NumPart, num_already_spawned, All.MaxPart);
        fflush(stdout); endrun(8888);
    }
    
#ifndef WAKEUP
    double dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
    double dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
    double d_r = 0.25 * KERNEL_CORE_SIZE*PPP[i].Hsml; // needs to be epsilon*Hsml where epsilon<<1, to maintain stability //
    double r2=0; for(k=0;k<3;k++) {r2+=(P[dummy_sph_i_to_clone].Pos[k]-P[i].Pos[k])*(P[dummy_sph_i_to_clone].Pos[k]-P[i].Pos[k]);}
    d_r = DMIN(d_r, 0.5*sqrt(r2));
#ifndef SELFGRAVITY_OFF
    d_r = DMAX(d_r , 2.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0]);
#endif
#ifdef BH_DEBUG_SPAWN_JET_TEST
    d_r = DMIN(d_r , 0.01); /* PFH: need to write this in a way that does not make assumptions about units/problem structure */
#endif
    double jz[3]={0,0,1},jy[3]={0,1,0},jx[3]={1,0,0};  /* set up a coordinate system [xyz if we don't have any other information */
#ifdef BH_FOLLOW_ACCRETED_ANGMOM  /* use local angular momentum to estimate preferred directions/coordinates for spawning */
    double Jtot=0; for(k=0;k<3;k++) {Jtot+=P[i].BH_Specific_AngMom[k]*P[i].BH_Specific_AngMom[k];}
    if(Jtot>0) {Jtot=1/sqrt(Jtot); for(k=0;k<3;k++) {jz[k]=P[i].BH_Specific_AngMom[k]*Jtot;}}
#ifdef JET_DIRECTION_FROM_KERNEL_AND_SINK //direction from the mass weighted average of the sink and the gas kernel angular momentum
    Jtot=0; for(k=0;k<3;k++) {Jtot+=P[i].Jgas_in_Kernel[k]*P[i].Jgas_in_Kernel[k];}
    if(Jtot>0) {Jtot=1/sqrt(Jtot); for(k=0;k<3;k++) {jz[k]=jz[k]*P[i].Mass + P[i].Jgas_in_Kernel[k]*Jtot*P[i].Mgas_in_Kernel;}}
#endif
    Jtot=jz[1]*jz[1]+jz[2]*jz[2]; if(Jtot>0) {Jtot=1/sqrt(Jtot); jy[1]=jz[2]*Jtot; jy[2]=-jz[1]*Jtot; for(k=0;k<3;k++) {jz[k]*=Jtot;}}
    jx[0]=jz[1]*jy[2]-jz[2]*jy[1]; jx[1]=jz[2]*jy[0]-jz[0]*jy[2]; jx[2]=jz[0]*jy[1]-jz[1]*jy[0];
#endif
    long bin, bin_0; for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinCount[bin] > 0) break;} /* gives minimum active timebin of any particle */
    bin_0 = bin; int i0 = i; /* save minimum timebin, also save ID of BH particle for use below */    
    bin = P[i0].TimeBin; /* make this particle active on the BH/star timestep */
#ifdef BH_DEBUG_SPAWN_JET_TEST
    bin = bin_0; i0 = dummy_sph_i_to_clone; /* make this particle active on the minimum timestep, and order with respect to the cloned particle */
#endif

    
    /* create the  new particles to be added to the end of the particle list :
        i is the BH particle tag, j is the new "spawed" particle's location, dummy_sph_i_to_clone is a dummy SPH particle's tag to be used to init the wind particle */
    for(j = NumPart + num_already_spawned; j < NumPart + num_already_spawned + n_particles_split; j++)
    {   /* first, clone the 'dummy' particle so various fields are set appropriately */
        P[j] = P[dummy_sph_i_to_clone]; SphP[j] = SphP[dummy_sph_i_to_clone]; /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */

        /* now we need to make sure everything is correctly placed in timebins for the tree */
        P[j].TimeBin = bin; P[j].dt_step = bin ? (((integertime) 1) << bin) : 0; // put this particle into the appropriate timebin
        NextActiveParticle[j] = FirstActiveParticle; FirstActiveParticle = j; NumForceUpdate++;
        TimeBinCount[bin]++; TimeBinCountSph[bin]++; PrevInTimeBin[j] = i0; /* likewise add it to the counters that register how many particles are in each timebin */
#ifndef BH_DEBUG_SPAWN_JET_TEST
        NextInTimeBin[j] = NextInTimeBin[i0]; if(NextInTimeBin[i0] >= 0) {PrevInTimeBin[NextInTimeBin[i0]] = j;}
        NextInTimeBin[i0] = j; if(LastInTimeBin[bin] == i0) {LastInTimeBin[bin] = j;}
#else
        if(FirstInTimeBin[bin] < 0) {FirstInTimeBin[bin]=j; LastInTimeBin[bin]=j; NextInTimeBin[j]=-1; PrevInTimeBin[j]=-1;} /* only particle in this time bin on this task */
            else {NextInTimeBin[j]=FirstInTimeBin[bin]; PrevInTimeBin[j]=-1; PrevInTimeBin[FirstInTimeBin[bin]]=j; FirstInTimeBin[bin]=j;} /* there is already at least one particle; add this one "to the front" of the list */
#endif
        P[j].Ti_begstep = All.Ti_Current; P[j].Ti_current = All.Ti_Current;
#ifdef WAKEUP
        PPPZ[j].wakeup = 1;
#endif
        /* this is a giant pile of variables to zero out. dont need everything here because we cloned a valid particle, but handy anyways */
        P[j].Particle_DivVel = 0; SphP[j].DtInternalEnergy = 0; for(k=0;k<3;k++) {SphP[j].HydroAccel[k] = 0; P[j].GravAccel[k] = 0;}
        P[j].NumNgb=All.DesNumNgb;
#ifdef PMGRID
        for(k=0;k<3;k++) {P[j].GravPM[k] = 0;}
#endif
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        SphP[j].MaxKineticEnergyNgb = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].dMass = 0; SphP[j].DtMass = 0; SphP[j].MassTrue = P[j].Mass; for(k=0;k<3;k++) {SphP[j].GravWorkTerm[k] = 0;}
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        PPPZ[j].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        PPP[j].AGS_Hsml = PPP[j].Hsml;
#endif
#endif
#ifdef CONDUCTION
        SphP[j].Kappa_Conduction = 0;
#endif
#ifdef MHD_NON_IDEAL
        SphP[j].Eta_MHD_OhmicResistivity_Coeff = 0; SphP[j].Eta_MHD_HallEffect_Coeff = 0; SphP[j].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
#ifdef VISCOSITY
        SphP[j].Eta_ShearViscosity = 0; SphP[j].Zeta_BulkViscosity = 0;
#endif
#ifdef TURB_DIFFUSION
        SphP[j].TD_DiffCoeff = 0;
#endif
#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==1)
        SphP[j].HostHaloMass = 0;
#endif
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        SphP[j].DelayTimeCoolingSNe = 0;
#endif
#ifdef GALSF
        SphP[j].Sfr = 0;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        SphP[j].alpha = 0.0;
#endif
#if defined(BH_THERMALFEEDBACK)
        SphP[j].Injected_BH_Energy = 0;
#endif
#ifdef RADTRANSFER
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            SphP[j].E_gamma[k] = 0;
#if defined(RT_EVOLVE_NGAMMA)
            SphP[j].E_gamma_Pred[k] = 0; SphP[j].Dt_E_gamma[k] = 0;
#endif
        }
#endif
        /* note, if you want to use this routine to inject magnetic flux or cosmic rays, do this below */
#ifdef MAGNETIC
        SphP[j].divB = 0; for(k=0;k<3;k++) {SphP[j].B[k]*=1.e-10; SphP[j].BPred[k]*=1.e-10; SphP[j].DtB[k]=0;} /* add magnetic flux here if desired */
#ifdef DIVBCLEANING_DEDNER
        SphP[j].DtPhi = SphP[j].PhiPred = SphP[j].Phi = 0;
#endif
#endif
        
        /* now set the real hydro variables. */
        /* set the particle ID */ // unsigned int bits; int SPLIT_GENERATIONS = 4; for(bits = 0; SPLIT_GENERATIONS > (1 << bits); bits++); /* the particle needs an ID: we give it a bit-flip from the original particle to signify the split */
        P[j].ID = All.AGNWindID; /* update:  We are using a fixed wind ID, to allow for trivial wind particle identification */
        P[j].ID_child_number = P[i].ID_child_number; P[i].ID_child_number +=1; P[j].ID_generation = P[i].ID; // this allows us to track spawned particles by giving them unique sub-IDs
        P[j].Mass = mass_of_new_particle; /* assign masses to both particles (so they sum correctly) */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].MassTrue = P[j].Mass;
#endif
#ifndef BH_DEBUG_FIX_MASS
        P[i].Mass -= P[j].Mass; /* make sure the operation is mass conserving! */
#endif
        BPP(i).unspawned_wind_mass -= P[j].Mass; /* remove the mass successfully spawned, to update the remaining unspawned mass */
        /* positions: uniformly sample unit sphere, and rotate into preferred coordinate system for use below */
        double phi=2.*M_PI*get_random_number(j+1+ThisTask), cos_theta=2.*(get_random_number(j+3+2*ThisTask)-0.5), sin_theta=sqrt(1-cos_theta*cos_theta), dx[3], sin_phi=sin(phi), cos_phi=cos(phi);
        for(k=0;k<3;k++) {P[j].Pos[k]=P[i].Pos[k] + (sin_theta*cos_phi*jx[k] + sin_theta*sin_phi*jy[k] + cos_theta*jz[k])*d_r;} // actually lay down position (in code coordinates)

        /* velocities (determined by wind velocity) */
        double veldir[3]; veldir[0]=sin_theta*cos_phi; veldir[1]=sin_theta*sin_phi; veldir[2]=cos_theta; // default to velocity pointed radially away from BH
#if defined(BH_DEBUG_SPAWN_JET_TEST) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(JET_DIRECTION_FROM_KERNEL_AND_SINK) || defined(BH_FB_COLLIMATED)
        double theta0=0.01, thetamax=80.*(M_PI/180.); // "flattening parameter" and max opening angle of jet velocity distribution from Matzner & McKee 1999, sets the collimation of the jets
        double jet_theta=atan(theta0*tan(get_random_number(j+7+5*ThisTask)*atan(sqrt(1+theta0*theta0)*tan(thetamax)/theta0))/sqrt(1+theta0*theta0)); // biased sampling to get collimation
        if(cos_theta<0) {jet_theta=M_PI-jet_theta;} // determines 'up' or 'down' based on which hemisphere particle is in
        veldir[0]=sin(jet_theta)*cos_phi; veldir[1]=sin(jet_theta)*sin_phi; veldir[2]=cos(jet_theta);//relative direction of velocity compared to BH_Specific_AngMom
#endif
        double v_magnitude = All.BAL_v_outflow * All.cf_atime; // velocity of the jet
        for(k=0;k<3;k++) {P[j].Vel[k]=P[i].Vel[k] + (veldir[0]*jx[k]+veldir[1]*jy[k]+veldir[2]*jz[k])*v_magnitude; SphP[j].VelPred[k]=P[j].Vel[k];}
        
        /* condition number, smoothing length, and density */
        SphP[j].ConditionNumber *= 100.0; /* boost the condition number to be conservative, so we don't trigger madness in the kernel */
        //SphP[j].Density *= 1e-10; SphP[j].Pressure *= 1e-10; PPP[j].Hsml = All.SofteningTable[0];  /* set dummy values: will be re-generated anyways [actually better to use nearest-neighbor values to start] */
#ifdef BH_DEBUG_SPAWN_JET_TEST
        PPP[j].Hsml=5.*d_r; SphP[j].Density=mass_of_new_particle/pow(KERNEL_CORE_SIZE*PPP[j].Hsml,NUMDIMS); /* PFH: need to write this in a way that does not make assumptions about units/problem structure */
#endif
        /* internal energy, determined by desired wind temperature */
        SphP[j].InternalEnergy = All.BAL_internal_temperature / (  PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g  ); SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;

#if defined(BH_COSMIC_RAYS) /* inject cosmic rays alongside wind injection */
        double dEcr = All.BH_CosmicRay_Injection_Efficiency * P[j].Mass * (All.BAL_f_accretion/(1.-All.BAL_f_accretion)) * (C / All.UnitVelocity_in_cm_per_s)*(C / All.UnitVelocity_in_cm_per_s);
        SphP[j].CosmicRayEnergy=dEcr; SphP[j].CosmicRayEnergyPred=dEcr;
#ifdef COSMIC_RAYS_M1
        dEcr*=COSMIC_RAYS_M1; for(k=0;k<3;k++) {SphP[j].CosmicRayFlux[k]=dEcr*dx[k]; SphP[j].CosmicRayFluxPred[k]=SphP[j].CosmicRayFlux[k];}
#endif
#endif
        /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
        force_add_star_to_tree(i0, j);// (buggy) /* we solve this by only calling the merge/split algorithm when we're doing the new domain decomposition */
    }    
    if(BPP(i).unspawned_wind_mass < 0) {BPP(i).unspawned_wind_mass=0;}
    return n_particles_split;
}
#endif
