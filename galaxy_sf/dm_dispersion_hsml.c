#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file dm_dispersion_hsml
 *  \brief smoothing length and velocity dispersion calculation for dark matter particles around gas particles
 *
 *  This file contains a loop modeled on the gas density computation which
 *    determines kernel lengths for dark matter particles around a given set of gas particles;
 *    this is used by the flag GALSF_SUBGRID_WIND_SCALING==2 to estimate the local dark
 *    matter velocity dispersion around a given gas particle, which (in turn) is used to set the
 *    sub-grid wind velocity and mass loading. The loop here needs to be called for these models (note this
 *    in general will require a different smoothing length from, say, the dm-dm force softening, or the
 *    gas kernel length for hydro, hence it requires a whole additional loop, even though the loop is
 *    functionally identical - modulo which particles are used - to the loop for adaptive gravitational softening
 *
 * This file was written by Qirong Zhu, for GIZMO, based on Phil Hopkins's adaptive gravitational softening
 *    routine. It has been modified by Phil on re-merger into the main branch of GIZMO with various optimizations
 *
 */

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)

#define CORE_FUNCTION_NAME disp_density_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME disp_particle2in_density    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME disp_out2particle_density  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(disp_density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/* this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3];
    MyFloat HsmlDM;
    int NodeList[NODELISTLENGTH];
}
*DATAIN_NAME, *DATAGET_NAME;

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
void disp_particle2in_density(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    in->HsmlDM = SphP[i].HsmlDM;
}


/* this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble Ngb, DM_Vel_Disp, DM_Vx, DM_Vy, DM_Vz;
}
*DATARESULT_NAME, *DATAOUT_NAME;

/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
void disp_out2particle_density(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    ASSIGN_ADD(SphP[i].DM_Vx, out->DM_Vx, mode);
    ASSIGN_ADD(SphP[i].DM_Vy, out->DM_Vy, mode);
    ASSIGN_ADD(SphP[i].DM_Vz, out->DM_Vz, mode);
    ASSIGN_ADD(SphP[i].DM_VelDisp, out->DM_Vel_Disp, mode);
    ASSIGN_ADD(SphP[i].NumNgbDM, out->Ngb, mode);
}


/* routine to determine if we need to use disp_density to calculate Hsml */
int disp_density_isactive(int i);
int disp_density_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0;
    if(P[i].Type > 0) return 0; // only gas particles //
    if(P[i].Mass <= 0) return 0;
    return 1;
}


/*! This function represents the core of the density computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine contains no writes to shared memory -- */
int disp_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, n; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* define variables and zero memory and import data for local target*/
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* imports the data to the correct place and names */
    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.HsmlDM, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, 2); // search for high-res DM particles only: 2^1 = 2
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Mass <= 0) continue;
                out.DM_Vx += P[j].Vel[0]; out.DM_Vy += P[j].Vel[1]; out.DM_Vz += P[j].Vel[2];
                out.DM_Vel_Disp += (P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);
                out.Ngb++;
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


void disp_density(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    MyFloat *Left, *Right; double desnumngb=64, desnumngbdev=48; long long ntot; 
    int i, npleft, iter=0, redo_particle;
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
    /* initialize anything we need to about the active particles before their loop */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {if(disp_density_isactive(i)) {SphP[i].NumNgbDM = 0; Left[i] = Right[i] = 0;}}
    
    /* allocate buffers to arrange communication */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */

        /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        double tstart = my_second(), tend;
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(disp_density_isactive(i))
            {
                redo_particle = 0; /* now check whether we have enough neighbours, and are below the maximum search radius */
                double maxsoft = DMIN(All.MaxHsml, 10.0*PPP[i].Hsml);
                if(((SphP[i].NumNgbDM < desnumngb - desnumngbdev) || (SphP[i].NumNgbDM > (desnumngb + desnumngbdev)))
                   && (Right[i]-Left[i] > 0.001*Left[i] || Left[i]==0 || Right[i]==0))
                {
                    redo_particle = 1;
                }
                if(SphP[i].HsmlDM >= maxsoft)
                {
                    SphP[i].HsmlDM = maxsoft;
                    redo_particle = 0;
                }

                if(redo_particle)
                {
                    /* need to redo this particle */
                    npleft++;
                    
                    if(SphP[i].NumNgbDM < desnumngb-desnumngbdev) {Left[i]=DMAX(SphP[i].HsmlDM, Left[i]);}
                    if(SphP[i].NumNgbDM > desnumngb+desnumngbdev) {if(Right[i]>0) {Right[i]=DMIN(Right[i],SphP[i].HsmlDM);} else {Right[i]=SphP[i].HsmlDM;}}
                    
                    if(iter >= MAXITER - 10)
                    {
                        PRINT_WARNING("DM disp: i=%d task=%d ID=%llu Type=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                         i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, SphP[i].HsmlDM, Left[i], Right[i], (float) SphP[i].NumNgbDM, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                    }
                    
                    // right/left define upper/lower bounds from previous iterations
                    if(Right[i] > 0 && Left[i] > 0)
                    {
                        // geometric interpolation between right/left //
                        if(SphP[i].NumNgbDM > 1)
                        {
                            SphP[i].HsmlDM *= pow( desnumngb / SphP[i].NumNgbDM , 1./NUMDIMS );
                        } else {
                            SphP[i].HsmlDM *= 2.0;
                        }
                        if((SphP[i].HsmlDM<Right[i])&&(SphP[i].HsmlDM>Left[i]))
                        {
                            SphP[i].HsmlDM = pow(SphP[i].HsmlDM*SphP[i].HsmlDM*SphP[i].HsmlDM*SphP[i].HsmlDM * Left[i]*Right[i] , 1./6.);
                        } else {
                            if(SphP[i].HsmlDM>Right[i]) SphP[i].HsmlDM=Right[i];
                            if(SphP[i].HsmlDM<Left[i]) SphP[i].HsmlDM=Left[i];
                            SphP[i].HsmlDM = pow(SphP[i].HsmlDM * Left[i] * Right[i] , 1.0/3.0);
                        }
                    }
                    else
                    {
                        if(Right[i] == 0 && Left[i] == 0)
                        {
                            char buf[1000]; sprintf(buf, "DM disp: Right[i] == 0 && Left[i] == 0 && SphP[i].HsmlDM=%g\n", SphP[i].HsmlDM); terminate(buf);
                        }
                        double fac;
                        if(Right[i] == 0 && Left[i] > 0)
                        {
                            if(SphP[i].NumNgbDM > 1) {fac = log( desnumngb / SphP[i].NumNgbDM ) / NUMDIMS;} else {fac=1.4;}
                            if((SphP[i].NumNgbDM < 2*desnumngb)&&(SphP[i].NumNgbDM > 0.1*desnumngb)) {SphP[i].HsmlDM *= exp(fac);} else {SphP[i].HsmlDM *= 1.26;}
                        }
                        
                        if(Right[i] > 0 && Left[i] == 0)
                        {
                            if(SphP[i].NumNgbDM > 1) {fac = log( desnumngb / SphP[i].NumNgbDM ) / NUMDIMS;} else {fac=1.4;}
                            fac = DMAX(fac,-1.535);
                            if((SphP[i].NumNgbDM < 2*desnumngb)&&(SphP[i].NumNgbDM > 0.1*desnumngb)) {SphP[i].HsmlDM *= exp(fac);} else {SphP[i].HsmlDM /= 1.26;}
                        }
                    }
                }
                else {P[i].TimeBin = -P[i].TimeBin - 1;}	/* Mark as inactive */
            } //  if(disp_density_isactive(i))
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])

        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 0 && ThisTask == 0) {if(iter > 10) printf("DM disp: ngb iteration %d: need to repeat for %d%09d particles.\n", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("DM disp: failed to converge in neighbour iteration in disp_density()\n"); fflush(stdout); endrun(1155);}
        }
    }
    while(ntot > 0);
    
    /* iteration is done - de-malloc everything now */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    myfree(Right); myfree(Left);

    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0) {P[i].TimeBin = -P[i].TimeBin - 1;}
    }
    
    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(disp_density_isactive(i))
        {
            if(SphP[i].NumNgbDM > 0)
            {
                SphP[i].DM_Vx /= SphP[i].NumNgbDM;
                SphP[i].DM_Vy /= SphP[i].NumNgbDM;
                SphP[i].DM_Vz /= SphP[i].NumNgbDM;
                SphP[i].DM_VelDisp /= SphP[i].NumNgbDM;
                SphP[i].DM_VelDisp = (1./All.cf_atime) * sqrt(SphP[i].DM_VelDisp - SphP[i].DM_Vx * SphP[i].DM_Vx - SphP[i].DM_Vy * SphP[i].DM_Vy - SphP[i].DM_Vz * SphP[i].DM_Vz) / 1.732;//	   1d velocity dispersion
            } else {
                if((SphP[i].DM_VelDisp <= 0) || isnan(SphP[i].DM_VelDisp)) {SphP[i].DM_VelDisp = sqrt(P[i].Vel[0]*P[i].Vel[0]+P[i].Vel[1]*P[i].Vel[1]+P[i].Vel[2]*P[i].Vel[2])/All.cf_atime;}
            }
        }
    }
    
    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm; CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */


#endif
#endif




