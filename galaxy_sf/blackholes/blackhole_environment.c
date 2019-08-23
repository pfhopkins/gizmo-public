#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"

/*! \file blackhole_environment.c
 *  \brief routines for evaluating black hole environment
 */
/*
 * This file was largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   by Paul Torrey (ptorrey@mit.edu) on 1/9/15 for clairity.  The main functional difference is that BlackholeTempInfo
 *   is now allocated only for N_active_loc_BHs, rather than NumPart (as was done before).  Some
 *   extra index gymnastics are required to follow this change through in the MPI comm routines.
 *   Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */



/* quantities that pass IN to the 'blackhole_environment_evaluate' routines */
static struct blackholedata_in
{
#if defined(BH_GRAVCAPTURE_GAS) || (BH_GRAVACCRETION == 8)
    MyDouble Mass;
#endif
#if defined(BH_GRAVCAPTURE_FIXEDSINKRADIUS)
    MyFloat SinkRadius;
#endif  
    MyDouble Pos[3];
    MyFloat Vel[3];
#ifdef BH_GRAVACCRETION
    MyFloat Jgas[3];
    MyFloat Jstar[3];
#endif
    MyFloat Hsml;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    MyFloat AGS_Hsml;
#endif    
    MyIDType ID;
    int NodeList[NODELISTLENGTH];
#ifdef BH_WAKEUP_GAS
    MyFloat TimeBin;
#endif
}
*BlackholeDataIn, *BlackholeDataGet;



void blackhole_environment_loop(void)
{
    int i, j, k, nexport, nimport, place, ngrp, recvTask, dummy;
    int ndone_flag, ndone;
    MPI_Status status;
    
    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackhole_temp_particle_data) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackhole_temp_particle_data))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    /* Scan gas particles for angular momentum, boundedness, etc */
    i = FirstActiveParticle;	/* first particle for this task */
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local active BH particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 5)
                if(blackhole_environment_evaluate(i, 0, &nexport, Send_count) < 0)
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
            }
#if defined(BH_GRAVCAPTURE_GAS) || (BH_GRAVACCRETION == 8)
            BlackholeDataIn[j].Mass = P[place].Mass;
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
            BlackholeDataIn[j].SinkRadius = P[place].SinkRadius;
#endif	    
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            BlackholeDataIn[j].AGS_Hsml = PPP[place].AGS_Hsml;
#endif	    
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;	    
            BlackholeDataIn[j].ID = P[place].ID;
#ifdef BH_WAKEUP_GAS
	        BlackholeDataIn[j].TimeBin = P[place].TimeBin;
#endif
            memcpy(BlackholeDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
                                 recvTask, TAG_BH_A,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_A, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);
   
        /*  DAA: PasserResult cycles over nimport!  PasserOut cycles over nexport! */
        BlackholeDataPasserResult = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserResult", nimport * sizeof(struct blackhole_temp_particle_data));
        BlackholeDataPasserOut = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserOut", nexport * sizeof(struct blackhole_temp_particle_data));
        
        memset( &BlackholeDataPasserResult[0], 0, nimport * sizeof(struct blackhole_temp_particle_data)  );
        memset( &BlackholeDataPasserOut[0],    0, nexport * sizeof(struct blackhole_temp_particle_data)  );

        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_environment_evaluate(j, 1, &dummy, &dummy);
        
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
                    MPI_Sendrecv(&BlackholeDataPasserResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                 MPI_BYTE, recvTask, TAG_BH_B,
                                 &BlackholeDataPasserOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                 MPI_BYTE, recvTask, TAG_BH_B, MPI_COMM_WORLD, &status);
                }
            }
        } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
#if defined(BH_ACCRETE_NEARESTFIRST) && defined(BH_GRAVCAPTURE_GAS)
            if(P[place].BH_dr_to_NearestGasNeighbor > BlackholeDataPasserOut[j].BH_dr_to_NearestGasNeighbor) {P[place].BH_dr_to_NearestGasNeighbor  = BlackholeDataPasserOut[j].BH_dr_to_NearestGasNeighbor;}
#endif
            out2particle_blackhole(&BlackholeDataPasserOut[j], P[place].IndexMapToTempStruc, 1);
        } // for(j = 0; j < nexport; j++)
        myfree(BlackholeDataPasserOut);
        myfree(BlackholeDataPasserResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);

    /* DAA: normalize/finalized some of the environment variables */
    for(i=0; i<N_active_loc_BHs; i++) {normalize_temp_info_struct(i);}        
}




/* routine to return the values we need of the properties of the gas, stars, etc in the vicinity of the BH -- these all factor into the BHAR */
int blackhole_environment_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    /* initialize variables before SPH loop is started */
    int startnode, numngb, j, k, n, listindex=0, mod_index;
    MyFloat *pos, h_i, *vel, hinv, mass, hinv3; mass=0; MyFloat ags_h_i = All.ForceSoftening[5]; MyIDType id;
    MyDouble dP[3],dv[3],wt;
    struct blackhole_temp_particle_data out;
    memset(&out, 0, sizeof(struct blackhole_temp_particle_data));
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS) || (BH_GRAVACCRETION == 8)
    MyFloat wk, dwk, u; u=wk=dwk=0;
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    MyFloat vrel, vbound, r2, sink_radius; sink_radius=0;
#ifdef BH_ACCRETE_NEARESTFIRST
    out.BH_dr_to_NearestGasNeighbor = MAX_REAL_NUMBER; // initialize large value
#endif
#endif
    
    /* these are the BH properties */
    if(mode == 0)
    {
#if defined(BH_GRAVCAPTURE_GAS) || (BH_GRAVACCRETION == 8)
        mass = P[target].Mass;
#if defined(BH_GRAVCAPTURE_FIXEDSINKRADIUS)
        sink_radius = P[target].SinkRadius;
#endif	
#endif
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = PPP[target].Hsml;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
	    ags_h_i = PPP[target].AGS_Hsml;
#endif	
        id = P[target].ID;
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
#ifdef BH_WAKEUP_GAS
	    bh_timebin = P[target].TimeBin;
#endif
    }
    else
    {
#if defined(BH_GRAVCAPTURE_GAS) || (BH_GRAVACCRETION == 8)
        mass = BlackholeDataGet[target].Mass;
#if defined(BH_GRAVCAPTURE_FIXEDSINKRADIUS)
        sink_radius = BlackholeDataGet[target].SinkRadius;
#endif		
#endif
        pos = BlackholeDataGet[target].Pos;
        vel = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
	    ags_h_i = BlackholeDataGet[target].AGS_Hsml;
#endif	
        id = BlackholeDataGet[target].ID;
        mod_index = 0;                              /* this is not used for mode==1, but this avoids compiler error */
#ifdef BH_WAKEUP_GAS
	    bh_timebin = BlackholeDataGet[target].TimeBin;
#endif 
    }
    
    if(h_i < 0) return -1;
    hinv = 1./h_i; hinv3 = hinv*hinv*hinv;
    if(mode == 0)
    {
        startnode = All.MaxPart;  /* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;        /* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb = ngb_treefind_pairs_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG);	   
            if(numngb < 0) return -1;
            
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];

#ifdef BH_WAKEUP_GAS
		        if (bh_timebin < P[j].LowestBHTimeBin) {P[j].LowestBHTimeBin = bh_timebin;}
#endif
                if( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
                {
                    wt = P[j].Mass;
                    dP[0] = P[j].Pos[0]-pos[0]; dP[1] = P[j].Pos[1]-pos[1]; dP[2] = P[j].Pos[2]-pos[2];
#ifdef BOX_PERIODIC
                    NEAREST_XYZ(dP[0],dP[1],dP[2],-1); /*  find the closest image in the given box size  */
#endif
                    dv[0] = P[j].Vel[0]-vel[0]; dv[1] = P[j].Vel[1]-vel[1]; dv[2] = P[j].Vel[2]-vel[2];
#ifdef BOX_SHEARING
                    if(pos[0] - P[j].Pos[0] > +boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                    if(pos[0] - P[j].Pos[0] < -boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif

#ifdef BH_DYNFRICTION
#if (BH_DYNFRICTION == 1)    // DAA: dark matter + stars
                    if( !(P[j].Type==0) )
#if (BH_REPOSITION_ON_POTMIN == 2)
                    if( (P[j].Type != 5) )
#endif
#elif (BH_DYNFRICTION == 2)  // DAA: stars only
                    if( P[j].Type==4 || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn)) )
#endif
                    {
                        double wtfac = wt;
#if (BH_REPOSITION_ON_POTMIN == 2)
                        double rfac = (dP[0]*dP[0] + dP[1]*dP[1] + dP[2]*dP[2]) * (10./(h_i*h_i) + 0.1/(All.ForceSoftening[5]*All.ForceSoftening[5]));
                        wtfac = wt / (1. + rfac); // simple function scaling ~ 1/r^2 for large r, to weight elements closer to the BH, so doesnt get 'pulled' by far-away elements //
#endif
                        if(P[j].Mass>out.DF_mmax_particles) out.DF_mmax_particles=P[j].Mass;
                        for (k=0;k<3;k++)
                        {
                            out.DF_mean_vel[k] += wt*dv[k];
#if (BH_REPOSITION_ON_POTMIN == 2)
                            out.DF_rms_vel += wt;
#else
                            out.DF_rms_vel += wt*dv[k]*dv[k];
#endif
                        }
                    }
#endif
                    
                    /* DAA: compute mass/angular momentum for GAS/STAR/DM components within BH kernel
                            this is done always now (regardless of the specific BH options used) */
                    if(P[j].Type==0)
                    {
                        /* we found gas in BH's kernel */
                        out.Mgas_in_Kernel += wt;
                        out.Sfr_in_Kernel += SphP[j].Sfr;
                        out.BH_InternalEnergy += wt*SphP[j].InternalEnergy;
                        out.Jgas_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jgas_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jgas_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
                        u=0; for(k=0;k<3;k++) {u+=dP[k]*dP[k];}
                        u=sqrt(u)/h_i; kernel_main(u,hinv3,hinv3*hinv,&wk,&dwk,1);
                        dwk /= u*h_i; for(k=0;k<3;k++) out.GradRho_in_Kernel[k] += wt * dwk * fabs(dP[k]);
#endif
#if defined(BH_BONDI) || defined(BH_DRAG) || (BH_GRAVACCRETION >= 5)
                        for(k=0;k<3;k++) {out.BH_SurroundingGasVel[k] += wt*dv[k];}
#endif
#if (BH_GRAVACCRETION == 8)
                        u=0; for(k=0;k<3;k++) {u+=dP[k]*dP[k];}
                        u=sqrt(u)/h_i; kernel_main(u,hinv3,hinv3*hinv,&wk,&dwk,-1);
                        double rj=u*h_i*All.cf_atime; double csj=Particle_effective_soundspeed_i(j);
                        double vdotrj=0; for(k=0;k<3;k++) {vdotrj+=-dP[k]*dv[k];}
                        double vr_mdot = 4*M_PI * wt*(wk*All.cf_a3inv) * rj*vdotrj;
			
                        if(rj < All.ForceSoftening[5]*All.cf_atime)
                        {
                            double bondi_mdot = 4*M_PI*All.G*All.G * mass*mass / pow(csj*csj + (dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2])*All.cf_a2inv, 1.5) * wt * (wk*All.cf_a3inv);
                            vr_mdot = DMAX(vr_mdot , bondi_mdot);
                            out.hubber_mdot_bondi_limiter += bondi_mdot;
                        }
                        out.hubber_mdot_vr_estimator += vr_mdot; /* physical */
                        out.hubber_mdot_disk_estimator += wt*wk * sqrt(rj) / (SphP[j].Density * csj*csj); /* physical */
#endif
                    }
                    else if( P[j].Type==4 || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn)) ) 
                    {
                        /* stars */
                        out.Mstar_in_Kernel += wt;
                        out.Jstar_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jstar_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jstar_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                    }
                    else 
                    { 
                        /* dark matter */ // DAA: Jalt_in_Kernel and Malt_in_Kernel are updated in normalize_temp_info_struct() to be TOTAL angular momentum and mass 
                        out.Malt_in_Kernel += wt;
                        out.Jalt_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jalt_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jalt_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                    }
                    

#if defined(BH_GRAVCAPTURE_GAS)
                    /* XM: I formally distinguish BH_GRAVCAPTURE_GAS and BH_GRAVCAPTURE_NONGAS. The former applies to
                     gas ONLY, as an accretion model. The later can be combined with any accretion model.
                     Currently, I only allow gas accretion to contribute to BH_Mdot (consistent with the energy radiating away).
                     For star particles, if there is an alpha-disk, they are captured to the disk. If not, they directly go
                     to the hole, without any contribution to BH_Mdot and feedback. This can be modified in the swallow loop
                     for other purposes. The goal of the following part is to estimate BH_Mdot, which will be used to evaluate feedback strength.
                     Therefore, we only need it when we enable BH_GRAVCAPTURE_GAS as gas accretion model. */
                    if( (P[j].Mass > 0) && (P[j].Type == 0))
                    {
                        vrel=r2=0; for(k=0;k<3;k++) {vrel+=dv[k]*dv[k]; r2+=dP[k]*dP[k];}
                        double dr_code = sqrt(r2); vrel = sqrt(vrel) / All.cf_atime;
                        vbound = bh_vesc(j, mass, dr_code, ags_h_i);
                        if(vrel < vbound) { /* bound */
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                            double spec_mom=0; for(k=0;k<3;k++) {spec_mom += dv[k]*dP[k];} // delta_x.delta_v
                            spec_mom = (r2*vrel*vrel - spec_mom*spec_mom*All.cf_a2inv);  // specific angular momentum^2 = r^2(delta_v)^2 - (delta_v.delta_x)^2;
                            if(spec_mom < All.G * (mass + P[j].Mass) * sink_radius) // check Bate 1995 angular momentum criterion (in addition to bounded-ness)
#endif
                            if( bh_check_boundedness(j,vrel,vbound,dr_code,sink_radius)==1 ) 
                            { /* apocenter within 2.8*epsilon (softening length) */
#ifdef SINGLE_STAR_SINK_DYNAMICS
                                double eps = DMAX(P[j].Hsml/2.8, DMAX(dr_code, ags_h_i/2.8));
			                    double tff = eps*eps*eps / (mass + P[j].Mass);
			                    if(tff < P[j].SwallowTime) {P[j].SwallowTime = tff;}
#endif
#if defined(BH_ACCRETE_NEARESTFIRST)
                                if((out.BH_dr_to_NearestGasNeighbor > dr_code) && (P[j].SwallowID < id)) {out.BH_dr_to_NearestGasNeighbor = dr_code; out.mass_to_swallow_edd = P[j].Mass;}
#else
                                if(P[j].SwallowID < id) {out.mass_to_swallow_edd += P[j].Mass;} /* mark as 'will be swallowed' on next loop, to correct accretion rate */
#endif
                            } /* if( apocenter in tolerance range ) */
                        } /* if(vrel < vbound) */
                    } /* type check */
#endif // BH_GRAVCAPTURE_GAS
                } // ( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
            } // for(n = 0; n < numngb; n++)

            if(mode == 0) /* local -> send directly to local temp struct */
            {
#if defined(BH_ACCRETE_NEARESTFIRST) && defined(BH_GRAVCAPTURE_GAS)
                P[target].BH_dr_to_NearestGasNeighbor = out.BH_dr_to_NearestGasNeighbor;
#endif
                out2particle_blackhole(&out, mod_index, 0);     /* target -> mod_index for reduced size struc */
            }
            else
            {
                BlackholeDataPasserResult[target] = out;        /* this is okay because target cycles over nimport only */
            }
            
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];   /* non-local target o.k. */
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            } // if(listindex < NODELISTLENGTH)
        } // if(mode == 1)
    } // while(startnode >= 0)
    return 0;
}





/* -----------------------------------------------------------------------------------------------------
 * DAA: modified versions of blackhole_environment_loop and blackhole_environment_evaluate for a second
 * environment loop. Here we do a Bulge-Disk kinematic decomposition for gravitational torque accretion
 * ----------------------------------------------------------------------------------------------------- 
 */


#if defined(BH_GRAVACCRETION) && (BH_GRAVACCRETION == 0)

void blackhole_environment_second_loop(void)
{
    int i, j, k, nexport, nimport, place, ngrp, recvTask, dummy;
    int ndone_flag, ndone;
    int mod_index;
    MPI_Status status;

    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackhole_temp_particle_data) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackhole_temp_particle_data))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    /* Scan gas particles for angular momentum, boundedness, etc */
    i = FirstActiveParticle;    /* first particle for this task */
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local active BH particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 5)
                if(blackhole_environment_second_evaluate(i, 0, &nexport, Send_count) < 0)
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
            mod_index = P[place].IndexMapToTempStruc;  /* DAA: the index of the BlackholeTempInfo structure */
            for(k = 0; k < 3; k++)
            {
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
                BlackholeDataIn[j].Jgas[k] = BlackholeTempInfo[mod_index].Jgas_in_Kernel[k];    // DAA: Jgas is available after first environment loop
                BlackholeDataIn[j].Jstar[k] = BlackholeTempInfo[mod_index].Jstar_in_Kernel[k];
            }
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            BlackholeDataIn[j].AGS_Hsml = PPP[place].AGS_Hsml;
#endif	    
            BlackholeDataIn[j].ID = P[place].ID;
            memcpy(BlackholeDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
                                 recvTask, TAG_BH_C,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_C, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);

        //DAA: PasserResult cycles over nimport, PasserOut cycles over nexport:
        BlackholeDataPasserResult = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserResult", nimport * sizeof(struct blackhole_temp_particle_data));
        BlackholeDataPasserOut = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserOut", nexport * sizeof(struct blackhole_temp_particle_data));

 /* -------------- NOTE
 *         DAA: note that we actually don't need the full BlackholeDataPasserResult/BlackholeDataPasserOut structures here... 
 *         --> we could only pass Mbulge_in_Kernel !
 *  -------------- NOTE
 */

        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_environment_second_evaluate(j, 1, &dummy, &dummy);

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
                                MPI_Sendrecv(&BlackholeDataPasserResult[Recv_offset[recvTask]],
                                             Recv_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                             MPI_BYTE, recvTask, TAG_BH_D,
                                             &BlackholeDataPasserOut[Send_offset[recvTask]],
                                             Send_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                             MPI_BYTE, recvTask, TAG_BH_D, MPI_COMM_WORLD, &status);
                            }
                        }
                    } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++)

        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            mod_index = P[place].IndexMapToTempStruc;
            /* DAA: we only need to update Mbulge_in_Kernel */
            BlackholeTempInfo[mod_index].MgasBulge_in_Kernel += BlackholeDataPasserOut[j].MgasBulge_in_Kernel;
            BlackholeTempInfo[mod_index].MstarBulge_in_Kernel += BlackholeDataPasserOut[j].MstarBulge_in_Kernel;
        } // for(j = 0; j < nexport; j++)
        myfree(BlackholeDataPasserOut);
        myfree(BlackholeDataPasserResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);

}




int blackhole_environment_second_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    /* initialize variables before SPH loop is started */
    int startnode, numngb, j, n, listindex=0, mod_index;
    MyFloat *pos, h_i, *vel, *Jgas, *Jstar;
    MyIDType id;

    double dP[3],dv[3],J_tmp[3],wt;
    double MgasBulge_tmp, MstarBulge_tmp;
    MgasBulge_tmp=0; MstarBulge_tmp=0;

    /* these are the BH properties */
    if(mode == 0)
    {
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = PPP[target].Hsml;	
        id = P[target].ID;
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
        Jgas = BlackholeTempInfo[mod_index].Jgas_in_Kernel;      // DAA: Jgas/Jstar available after first environment loop
        Jstar = BlackholeTempInfo[mod_index].Jstar_in_Kernel;
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        vel = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
        mod_index = 0;                              /* this is not used for mode==1, but this avoids compiler error */
        Jgas = BlackholeDataGet[target].Jgas;
        Jstar = BlackholeDataGet[target].Jstar;
    }

    if(h_i < 0) return -1;

    if(mode == 0)
    {
        startnode = All.MaxPart;  /* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;        /* open it */
    }

    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb = ngb_treefind_pairs_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
            if(numngb < 0) return -1;

            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];


                if( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
                {
                    wt = P[j].Mass;
                    dP[0] = P[j].Pos[0]-pos[0];
                    dP[1] = P[j].Pos[1]-pos[1];
                    dP[2] = P[j].Pos[2]-pos[2];
#ifdef BOX_PERIODIC           
                    NEAREST_XYZ(dP[0],dP[1],dP[2],-1); /*  find the closest image in the given box size  */
#endif
                    dv[0] = P[j].Vel[0]-vel[0];
                    dv[1] = P[j].Vel[1]-vel[1];
                    dv[2] = P[j].Vel[2]-vel[2];
#ifdef BOX_SHEARING
                    if(pos[0] - P[j].Pos[0] > +boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                    if(pos[0] - P[j].Pos[0] < -boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif

                    // DAA: jx,jy,jz, are independent of 'a' because ~ m*r*v, vphys=v/a, rphys=r*a //
                    J_tmp[0] = wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                    J_tmp[1] = wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                    J_tmp[2] = wt*(dP[0]*dv[1] - dP[1]*dv[0]);

                    if(P[j].Type==0)
                    {
                        /* GAS */
                        if( (J_tmp[0]*Jgas[0] + J_tmp[1]*Jgas[1] + J_tmp[2]*Jgas[2]) < 0 )  
                        {
                            /* DAA: assume the bulge component contains as many particles with positive azimuthal velocities  
                               as with negative azimuthal velocities relative to the angular momentum vector */
                            MgasBulge_tmp += 2.0 * wt;
                        }
                    }
                    else if( P[j].Type==4 || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn)) )
                    {
                        /* STARS */
                        if( (J_tmp[0]*Jstar[0] + J_tmp[1]*Jstar[1] + J_tmp[2]*Jstar[2]) < 0 )   
                        {
                            MstarBulge_tmp += 2.0 * wt;
                        }
                    }


                } // ( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
            } // for(n = 0; n < numngb; n++)

            if(mode == 0) /* local -> send directly to local temp struct */
            {
                BlackholeTempInfo[mod_index].MgasBulge_in_Kernel = MgasBulge_tmp;
                BlackholeTempInfo[mod_index].MstarBulge_in_Kernel = MstarBulge_tmp;
            }
            else
            {
                BlackholeDataPasserResult[target].MgasBulge_in_Kernel = MgasBulge_tmp;
                BlackholeDataPasserResult[target].MstarBulge_in_Kernel = MstarBulge_tmp;
            }

        } // while(startnode >= 0)

        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];   /* non-local target o.k. */
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;  /* open it */
            } // if(listindex < NODELISTLENGTH)
        } // if(mode == 1)
    } // while(startnode >= 0)
    return 0;
}


#endif   //BH_GRAVACCRETION == 0
