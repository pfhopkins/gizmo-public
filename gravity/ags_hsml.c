#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file ags_hsml.c
 *  \brief kernel length determination for non-gas particles
 *
 *  This file contains a loop modeled on the gas density computation which 
 *    determines softening lengths (and appropriate correction terms) 
 *    for all particle types, to make softenings fully adaptive
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#define AGS_DSOFT_TOL (0.5)    // amount by which softening lengths are allowed to vary in single timesteps //

/*! this routine is called by the adaptive gravitational softening neighbor search and forcetree (for application 
    of the appropriate correction terms), to determine which particle types "talk to" which other particle types 
    (i.e. which particle types you search for to determine the softening radii for gravity). For effectively volume-filling
    fluids like gas or dark matter, it makes sense for this to be 'matched' to particles of the same type. For other 
    particle types like stars or black holes, it's more ambiguous, and requires some judgement on the part of the user. 
    The routine specifically returns a bitflag which defines all valid particles to which a particle of type 'primary' 
    can 'see': i.e. SUM(2^n), where n are all the particle types desired for neighbor finding,
    so e.g. if you want particle types 0 and 4, set the bitmask = 17 = 1 + 16 = 2^0 + 2^4
 */
int ags_gravity_kernel_shared_BITFLAG(short int particle_type_primary)
{
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(!((1 << particle_type_primary) & (ADAPTIVE_GRAVSOFT_FORALL))) {return 0;} /* particle is NOT one of the designated 'adaptive' types */
#endif
    
    if(particle_type_primary == 0) {return 1;} /* gas particles see gas particles */

#if (ADAPTIVE_GRAVSOFT_FORALL & 32) && defined(BLACK_HOLES)
    if(particle_type_primary == 5) {return 1;} /* black hole particles are AGS-active, but using BH physics, they see only gas */
#endif
    
#if defined(GALSF) && ( (ADAPTIVE_GRAVSOFT_FORALL & 16) || (ADAPTIVE_GRAVSOFT_FORALL & 8) || (ADAPTIVE_GRAVSOFT_FORALL & 4) )
    if(All.ComovingIntegrationOn) /* stars [4 for cosmo runs, 2+3+4 for non-cosmo runs] are AGS-active and see baryons (any type) */
    {
        if(particle_type_primary == 4) {return 17;} // 2^0+2^4
    } else {
        if((particle_type_primary == 4)||(particle_type_primary == 2)||(particle_type_primary == 3)) {return 29;} // 2^0+2^2+2^3+2^4
    }
#endif
    
#ifdef DM_SIDM
    if((1 << particle_type_primary) & (DM_SIDM)) {return DM_SIDM;} /* SIDM particles see other SIDM particles, regardless of type/mass */
#endif
    
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
    return (1 << particle_type_primary); /* if we haven't been caught by one of the above checks, we simply return whether or not we see 'ourselves' */
#endif
    
    return 0;
}



#ifdef AGS_HSML_CALCULATION_IS_ACTIVE

#define CORE_FUNCTION_NAME ags_density_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME ags_particle2in_density    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME ags_out2particle_density  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(ags_density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

struct kernel_density 
{
    double dp[3],dv[3],r, wk, dwk, hinv, hinv3, hinv4; /*! Structure for communication during the density computation. Holds data that is sent to other processors */
};

static struct INPUT_STRUCT_NAME
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat AGS_Hsml;
  int NodeList[NODELISTLENGTH];
  int Type;
}
 *DATAIN_NAME, *DATAGET_NAME;

void ags_particle2in_density(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration);
void ags_particle2in_density(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k]; in->Vel[k]=P[i].Vel[k];}
    in->AGS_Hsml = PPP[i].AGS_Hsml;
    in->Type = P[i].Type;
}


static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble Ngb;
    MyLongDouble DhsmlNgb;
    MyLongDouble AGS_zeta;
    MyLongDouble AGS_vsig;
    MyLongDouble Particle_DivVel;
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    MyLongDouble NV_T[3][3];
#endif
}
 *DATARESULT_NAME, *DATAOUT_NAME;

void ags_out2particle_density(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration);
void ags_out2particle_density(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    ASSIGN_ADD(PPP[i].NumNgb, out->Ngb, mode);
    ASSIGN_ADD(PPPZ[i].AGS_zeta, out->AGS_zeta,   mode);
    if(out->AGS_vsig > PPP[i].AGS_vsig) {PPP[i].AGS_vsig = out->AGS_vsig;}
    ASSIGN_ADD(P[i].Particle_DivVel, out->Particle_DivVel,   mode);
    ASSIGN_ADD(PPP[i].DhsmlNgbFactor, out->DhsmlNgb, mode);
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    {int j,k; for(k = 0; k < 3; k++) {for(j = 0; j < 3; j++) {ASSIGN_ADD(P[i].NV_T[k][j], out->NV_T[k][j], mode);}}}
#endif
}


/*! This function represents the core of the density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
/*!   -- this subroutine writes to shared memory [updating the neighbor values, primarily for wakeup-type updates]: need to protect these writes for openmp below */
int ags_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, n;
    int startnode, numngb_inbox, listindex = 0;
    double r2, h2, u;
    struct kernel_density kernel;
    struct INPUT_STRUCT_NAME local;
    struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    
    if(mode == 0)
        ags_particle2in_density(&local, target, loop_iteration);
    else
        local = DATAGET_NAME[target];
    
    h2 = local.AGS_Hsml * local.AGS_Hsml;
    kernel_hinv(local.AGS_Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(local.Type); // determine allowed particle types for search for adaptive gravitational softening terms
    
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = DATAGET_NAME[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }
    
    
    double fac_mu = -3 / (All.cf_afac3 * All.cf_atime);
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.AGS_Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, AGS_kernel_shared_BITFLAG);
            if(numngb_inbox < 0) {return -2;}
            
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Mass <= 0) continue;
                
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                if(r2 < h2)
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);

                    out.Ngb += kernel.wk;
                    out.DhsmlNgb += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
                    out.AGS_zeta += P[j].Mass * kernel_gravity(u, kernel.hinv, kernel.hinv3, 0); // needs to be here, should include self-contribution

                    if(kernel.r > 0)
                    {
                        if(P[j].Type==0)
                        {
                            kernel.dv[0] = local.Vel[0] - SphP[j].VelPred[0];
                            kernel.dv[1] = local.Vel[1] - SphP[j].VelPred[1];
                            kernel.dv[2] = local.Vel[2] - SphP[j].VelPred[2];
                        } else {
                            kernel.dv[0] = local.Vel[0] - P[j].Vel[0];
                            kernel.dv[1] = local.Vel[1] - P[j].Vel[1];
                            kernel.dv[2] = local.Vel[2] - P[j].Vel[2];
                        }
                        NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,kernel.dv,1); /* wrap velocities for shearing boxes if needed */
                        double v_dot_r = kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2];
                        if(v_dot_r > 0) {v_dot_r *= 0.333333;} // receding elements don't signal strong change in forces in the same manner as approaching/converging particles
                        double vsig = 0.5 * fabs( fac_mu * v_dot_r / kernel.r );
                        short int TimeBin_j = P[j].TimeBin; if(TimeBin_j < 0) {TimeBin_j = -TimeBin_j - 1;} // need to make sure we correct for the fact that TimeBin is used as a 'switch' here to determine if a particle is active for iteration, otherwise this gives nonsense!
                        if(vsig > out.AGS_vsig) {out.AGS_vsig = vsig;}
#if defined(WAKEUP) && (defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(DM_FUZZY) || defined(FLAG_NOT_IN_PUBLIC_CODE))
                        int wakeup_condition = 0; // determine if wakeup is allowed
                        if(!(TimeBinActive[TimeBin_j]) && (All.Time > All.TimeBegin) && (vsig > WAKEUP*P[j].AGS_vsig)) {wakeup_condition = 1;}
#if defined(GALSF)
                        if((P[j].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[j].Type == 2)||(P[j].Type==3)))) {wakeup_condition = 0;} // don't wakeup star particles, or risk 2x-counting feedback events! //
#endif
                        if(wakeup_condition) // do the wakeup
                        {
                                #pragma omp atomic write
                                P[j].wakeup = 1;
                                #pragma omp atomic write
                                NeedToWakeupParticles_local = 1;
                        }
#endif
                        out.Particle_DivVel -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how Hsml will evolve */
                        
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
                        out.NV_T[0][0] +=  kernel.wk * kernel.dp[0] * kernel.dp[0];
                        out.NV_T[0][1] +=  kernel.wk * kernel.dp[0] * kernel.dp[1];
                        out.NV_T[0][2] +=  kernel.wk * kernel.dp[0] * kernel.dp[2];
                        out.NV_T[1][1] +=  kernel.wk * kernel.dp[1] * kernel.dp[1];
                        out.NV_T[1][2] +=  kernel.wk * kernel.dp[1] * kernel.dp[2];
                        out.NV_T[2][2] +=  kernel.wk * kernel.dp[2] * kernel.dp[2];
#endif
                    }
                }
            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DATAGET_NAME[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}    /* open it */
            }
        }
    }
    if(mode == 0) {ags_out2particle_density(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
    return 0;
}



void ags_density(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second(); MyFloat *Left, *Right, *AGS_Prev; double fac, fac_lim, desnumngb, desnumngbdev; long long ntot;
    int i, npleft, iter=0, redo_particle, particle_set_to_minhsml_flag = 0, particle_set_to_maxhsml_flag = 0;
    AGS_Prev = (MyFloat *) mymalloc("AGS_Prev", NumPart * sizeof(MyFloat));
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
    /* initialize anything we need to about the active particles before their loop */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if(ags_density_isactive(i)) {
            Left[i] = Right[i] = 0; AGS_Prev[i] = PPP[i].AGS_Hsml; PPP[i].AGS_vsig = 0;
#ifdef WAKEUP
            P[i].wakeup = 0;
#endif
      }}

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
            if(ags_density_isactive(i))
            {
#ifdef DM_FUZZY
                P[i].AGS_Density = P[i].Mass * PPP[i].NumNgb;
#endif
                if(PPP[i].NumNgb > 0)
                {
                    PPP[i].DhsmlNgbFactor *= PPP[i].AGS_Hsml / (NUMDIMS * PPP[i].NumNgb);
                    P[i].Particle_DivVel /= PPP[i].NumNgb;
                    /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                    PPP[i].NumNgb *= NORM_COEFF * pow(PPP[i].AGS_Hsml,NUMDIMS);
                } else {
                    PPP[i].NumNgb = PPP[i].DhsmlNgbFactor = P[i].Particle_DivVel = 0;
                }
                
                // inverse of SPH volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(PPP[i].DhsmlNgbFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
                    PPP[i].DhsmlNgbFactor = 1 / (1 + PPP[i].DhsmlNgbFactor);
                else
                    PPP[i].DhsmlNgbFactor = 1;
                P[i].Particle_DivVel *= PPP[i].DhsmlNgbFactor;
                
                /* now check whether we have enough neighbours */
                redo_particle = 0;
                
                double minsoft = ags_return_minsoft(i);
                double maxsoft = ags_return_maxsoft(i);
                if(All.Time > All.TimeBegin)
                {
                    minsoft = DMAX(minsoft , AGS_Prev[i]*AGS_DSOFT_TOL);
                    maxsoft = DMIN(maxsoft , AGS_Prev[i]/AGS_DSOFT_TOL);
                }
                desnumngb = All.AGS_DesNumNgb;
                desnumngbdev = All.AGS_MaxNumNgbDeviation;
                /* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
                double ConditionNumber = do_cbe_nvt_inversion_for_faces(i); // right now we don't do anything with this, but could use to force expansion of search, as in hydro
                if(ConditionNumber > MAX_REAL_NUMBER) {PRINT_WARNING("CNUM for CBE: ThisTask=%d i=%d ConditionNumber=%g desnumngb=%g NumNgb=%g iter=%d NVT=%g/%g/%g/%g/%g/%g AGS_Hsml=%g \n",ThisTask,i,ConditionNumber,desnumngb,PPP[i].NumNgb,iter,P[i].NV_T[0][0],P[i].NV_T[1][1],P[i].NV_T[2][2],P[i].NV_T[0][1],P[i].NV_T[0][2],P[i].NV_T[1][2],PPP[i].AGS_Hsml);}
                if(iter > 10) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*((double)iter - 9.)) );}
#else
                if(iter > 4) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*((double)iter - 3.)) );}
#endif
                if(All.Time<=All.TimeBegin) {if(desnumngbdev > 0.0005) desnumngbdev=0.0005; if(iter > 50) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*((double)iter - 49.)) );}}


                /* check if we are in the 'normal' range between the max/min allowed values */
                if((PPP[i].NumNgb < (desnumngb - desnumngbdev) && PPP[i].AGS_Hsml < 0.999*maxsoft) ||
                   (PPP[i].NumNgb > (desnumngb + desnumngbdev) && PPP[i].AGS_Hsml > 1.001*minsoft))
                    redo_particle = 1;
                
                /* check maximum kernel size allowed */
                particle_set_to_maxhsml_flag = 0;
                if((PPP[i].AGS_Hsml >= 0.999*maxsoft) && (PPP[i].NumNgb < (desnumngb - desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].AGS_Hsml == maxsoft)
                    {
                        /* iteration at the maximum value is already complete */
                        particle_set_to_maxhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        PPP[i].AGS_Hsml = maxsoft;
                        particle_set_to_maxhsml_flag = 1;
                    }
                }
                
                /* check minimum kernel size allowed */
                particle_set_to_minhsml_flag = 0;
                if((PPP[i].AGS_Hsml <= 1.001*minsoft) && (PPP[i].NumNgb > (desnumngb + desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].AGS_Hsml == minsoft)
                    {
                        /* this means we've already done an iteration with the MinHsml value, so the
                         neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
                        particle_set_to_minhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        PPP[i].AGS_Hsml = minsoft;
                        particle_set_to_minhsml_flag = 1;
                    }
                }
                
                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
                        PRINT_WARNING("AGS: i=%d task=%d ID=%llu Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)\n",
                               i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, PPP[i].AGS_Hsml, PPP[i].DhsmlNgbFactor, Left[i], Right[i],
                               (float) PPP[i].NumNgb, Right[i] - Left[i], particle_set_to_maxhsml_flag, particle_set_to_minhsml_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                    }
                    
                    /* need to redo this particle */
                    npleft++;
                    
                    if(Left[i] > 0 && Right[i] > 0)
                        if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                            continue;
                        }
                    
                    if((particle_set_to_maxhsml_flag==0)&&(particle_set_to_minhsml_flag==0))
                    {
                        if(PPP[i].NumNgb < (desnumngb - desnumngbdev))
                        {
                            Left[i] = DMAX(PPP[i].AGS_Hsml, Left[i]);
                        }
                        else
                        {
                            if(Right[i] != 0)
                            {
                                if(PPP[i].AGS_Hsml < Right[i])
                                    Right[i] = PPP[i].AGS_Hsml;
                            }
                            else
                                Right[i] = PPP[i].AGS_Hsml;
                        }
                        
                        // right/left define upper/lower bounds from previous iterations
                        if(Right[i] > 0 && Left[i] > 0)
                        {
                            // geometric interpolation between right/left //
                            double maxjump=0;
                            if(iter>1) {maxjump = 0.2*log(Right[i]/Left[i]);}
                            if(PPP[i].NumNgb > 1)
                            {
                                double jumpvar = PPP[i].DhsmlNgbFactor * log( desnumngb / PPP[i].NumNgb ) / NUMDIMS;
                                if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
                                PPP[i].AGS_Hsml *= exp(jumpvar);
                            } else {
                                PPP[i].AGS_Hsml *= 2.0;
                            }
                            if((PPP[i].AGS_Hsml<Right[i])&&(PPP[i].AGS_Hsml>Left[i]))
                            {
                                if(iter > 1)
                                {
                                    double hfac = exp(maxjump);
                                    if(PPP[i].AGS_Hsml > Right[i] / hfac) {PPP[i].AGS_Hsml = Right[i] / hfac;}
                                    if(PPP[i].AGS_Hsml < Left[i] * hfac) {PPP[i].AGS_Hsml = Left[i] * hfac;}
                                }
                            } else {
                                if(PPP[i].AGS_Hsml>Right[i]) PPP[i].AGS_Hsml=Right[i];
                                if(PPP[i].AGS_Hsml<Left[i]) PPP[i].AGS_Hsml=Left[i];
                                PPP[i].AGS_Hsml = pow(PPP[i].AGS_Hsml * Left[i] * Right[i] , 1.0/3.0);
                            }
                        }
                        else
                        {
                            if(Right[i] == 0 && Left[i] == 0)
                            {
                                char buf[1000]; sprintf(buf, "AGS: Right[i] == 0 && Left[i] == 0 && PPP[i].AGS_Hsml=%g\n", PPP[i].AGS_Hsml); terminate(buf);
                            }
                            
                            if(Right[i] == 0 && Left[i] > 0)
                            {
                                if (PPP[i].NumNgb > 1)
                                    fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (+0.231=2x desnumngb)
                                else
                                    fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                                
                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4)
                                        if(PPP[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                    
                                    if(fac < fac_lim+0.231)
                                    {
                                        PPP[i].AGS_Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                    {
                                        PPP[i].AGS_Hsml *= exp(fac_lim+0.231);
                                        // fac~0.26 leads to expected doubling of number if density is constant,
                                        //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                    }
                                }
                                else
                                    PPP[i].AGS_Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                            
                            if(Right[i] > 0 && Left[i] == 0)
                            {
                                if (PPP[i].NumNgb > 1)
                                    fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                                else
                                    fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                                
                                if (fac_lim < -1.535) fac_lim = -1.535; // decreasing N_ngb by factor ~100
                                
                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=10)
                                        if(PPP[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                    
                                    if(fac > fac_lim-0.231)
                                    {
                                        PPP[i].AGS_Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                        PPP[i].AGS_Hsml *= exp(fac_lim-0.231); // limiter to prevent --too-- far a jump in a single iteration
                                }
                                else
                                    PPP[i].AGS_Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                        } // closes if(Right[i] > 0 && Left[i] > 0) else clause
                        
                    } // closes if[particle_set_to_max/minhsml_flag]
                    /* resets for max/min values */
                    if(PPP[i].AGS_Hsml < minsoft) PPP[i].AGS_Hsml = minsoft;
                    if(particle_set_to_minhsml_flag==1) PPP[i].AGS_Hsml = minsoft;
                    if(PPP[i].AGS_Hsml > maxsoft) PPP[i].AGS_Hsml = maxsoft;
                    if(particle_set_to_maxhsml_flag==1) PPP[i].AGS_Hsml = maxsoft;
                } // closes redo_particle
                else
                    P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
            } //  if(ags_density_isactive(i))
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        
        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 10 && ThisTask == 0) {printf("AGS-ngb iteration %d: need to repeat for %d%09d particles.\n", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("ags-failed to converge in neighbour iteration in density()\n"); fflush(stdout); endrun(1155);}
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
        if(ags_density_isactive(i))
        {
            if((P[i].Mass>0)&&(PPP[i].AGS_Hsml>0)&&(PPP[i].NumNgb>0))
            {
                double minsoft = ags_return_minsoft(i);
                double maxsoft = ags_return_maxsoft(i);
                minsoft = DMAX(minsoft , AGS_Prev[i]*AGS_DSOFT_TOL);
                maxsoft = DMIN(maxsoft , AGS_Prev[i]/AGS_DSOFT_TOL);
                if(PPP[i].AGS_Hsml >= maxsoft) {PPPZ[i].AGS_zeta = 0;} /* check that we're within the 'valid' range for adaptive softening terms, otherwise zeta=0 */

                double z0 = 0.5 * PPPZ[i].AGS_zeta * PPP[i].AGS_Hsml / (NUMDIMS * P[i].Mass * PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].AGS_Hsml,NUMDIMS) )); // zeta before various prefactors
                double h_eff = 2. * (KERNEL_CORE_SIZE*All.ForceSoftening[P[i].Type]); // force softening defines where Jeans pressure needs to kick in; prefactor = NJeans [=2 here]
                double Prho = 0 * h_eff*h_eff/2.; if(P[i].Particle_DivVel>0) {Prho=-Prho;} // truelove criterion. NJeans[above] , gamma=2 for effective EOS when this dominates, rho=ma*na; h_eff here can be Hsml [P/rho~H^-1] or gravsoft_min to really enforce that, as MIN, with P/rho~H^-3; if-check makes it so this term always adds KE to the system, pumping it up
                PPPZ[i].AGS_zeta = P[i].Mass*P[i].Mass * PPP[i].DhsmlNgbFactor * ( z0 + Prho ); // force correction, including corrections for adaptive softenings and EOS terms
                PPP[i].NumNgb = pow(PPP[i].NumNgb , 1./NUMDIMS); /* convert NGB to the more useful format, NumNgb^(1/NDIMS), which we can use to obtain the corrected particle sizes */
            } else {
                PPPZ[i].AGS_zeta = 0; PPP[i].NumNgb = 0; PPP[i].AGS_Hsml = All.ForceSoftening[P[i].Type];
            }
            apply_pm_hires_region_clipping_selection(i);
        }
    }
    myfree(AGS_Prev);
    
    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm; CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */









/* routine to determine if we need to use ags_density to calculate Hsml */
int ags_density_isactive(int i)
{
    int default_to_return = 0; // default to not being active - needs to be pro-actively 'activated' by some physics
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    default_to_return = 1;
    if(!((1 << P[i].Type) & (ADAPTIVE_GRAVSOFT_FORALL))) /* particle is NOT one of the designated 'adaptive' types */
    {
        PPP[i].AGS_Hsml = All.ForceSoftening[P[i].Type];
        PPPZ[i].AGS_zeta = 0;
        default_to_return = 0;
    } else {default_to_return = 1;} /* particle is AGS-active */
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
    if(P[i].Type==0)
    {
        PPP[i].AGS_Hsml = PPP[i].Hsml; // gas sees gas, these are identical
        default_to_return = 0; // don't actually need to do the loop //
    }
#endif
#ifdef DM_SIDM
    if((1 << P[i].Type) & (DM_SIDM)) {default_to_return = 1;}
#endif
#if defined(DM_FUZZY) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    if(P[i].Type == 1) {default_to_return = 1;}
#endif
    if(P[i].TimeBin < 0) {default_to_return = 0;} /* check our 'marker' for particles which have finished iterating to an Hsml solution (if they have, dont do them again) */
    return default_to_return;
}
    

/* routine to return the maximum allowed softening */
double ags_return_maxsoft(int i)
{
    double maxsoft = All.MaxHsml; // user-specified maximum: nothing is allowed to exceed this
#ifdef PMGRID /* Maximum allowed gravitational softening when using the TreePM method. The quantity is given in units of the scale used for the force split (PM_ASMTH) */
    maxsoft = DMIN(maxsoft, 1e3 * 0.5 * All.Asmth[0]); /* no more than 1/2 the size of the largest PM cell, times a 'safety factor' which can be pretty big */
#endif
#if (ADAPTIVE_GRAVSOFT_FORALL & 32) && defined(BLACK_HOLES) && !defined(SINGLE_STAR_SINK_DYNAMICS)
    if(P[i].Type == 5) {maxsoft = All.BlackHoleMaxAccretionRadius  / All.cf_atime;}   // MaxAccretionRadius is now defined in params.txt in PHYSICAL units
#endif
    return maxsoft;
}

    
/* routine to return the minimum allowed softening */
double ags_return_minsoft(int i)
{
    double minsoft = All.ForceSoftening[P[i].Type]; // this is the user-specified minimum
#if !defined(ADAPTIVE_GRAVSOFT_FORALL)
    minsoft = DMIN(All.MinHsml, minsoft);
#endif
    return minsoft;
}


/* routine to return effective particle sizes (inter-particle separation) based on AGS_Hsml saved values */
double INLINE_FUNC Get_Particle_Size_AGS(int i)
{
    /* in previous versions of the code, we took NumNgb^(1/NDIMS) here; however, now we
     take that when NumNgb is computed (at the end of the density routine), so we
     don't have to re-compute it each time. That makes this function fast enough to
     call -inside- of loops (e.g. hydro computations) */
#if (NUMDIMS == 1)
    return 2.00000 * PPP[i].AGS_Hsml / PPP[i].NumNgb; // (2)^(1/1)
#endif
#if (NUMDIMS == 2)
    return 1.77245 * PPP[i].AGS_Hsml / PPP[i].NumNgb; // (pi)^(1/2)
#endif
#if (NUMDIMS == 3)
    return 1.61199 * PPP[i].AGS_Hsml / PPP[i].NumNgb; // (4pi/3)^(1/3)
#endif
}


/* --------------------------------------------------------------------------
 very quick sub-routine to get the particle densities from their volumes
 -------------------------------------------------------------------------- */
double get_particle_volume_ags(int j)
{
    double L_j = Get_Particle_Size_AGS(j);
#if (NUMDIMS==1)
    return L_j;
#elif (NUMDIMS==2)
    return L_j*L_j;
#else
    return L_j*L_j*L_j;
#endif
}


#ifdef AGS_FACE_CALCULATION_IS_ACTIVE

/* --------------------------------------------------------------------------
 Subroutine here exists to calculate the MFM-like effective faces for purposes of face-interaction evaluation
 -------------------------------------------------------------------------- */

/* routine to invert the NV_T matrix after neighbor pass */
double do_cbe_nvt_inversion_for_faces(int i)
{
    /* initialize the matrix to be inverted */
    MyLongDouble NV_T[3][3], Tinv[3][3]; int j,k; for(j=0;j<3;j++) {for(k=0;k<3;k++) {NV_T[j][k]=P[i].NV_T[j][k];}}
    /* fill in the missing elements of NV_T (it's symmetric, so we saved time not computing these directly) */
    NV_T[1][0]=NV_T[0][1]; NV_T[2][0]=NV_T[0][2]; NV_T[2][1]=NV_T[1][2];
    /* want to work in dimensionless units for defining certain quantities robustly, so normalize out the units */
    double dimensional_NV_T_normalizer = pow( PPP[i].Hsml , 2-NUMDIMS ); /* this has the same dimensions as NV_T here */
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {NV_T[j][k] /= dimensional_NV_T_normalizer;}} /* now NV_T should be dimensionless */
    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
        this will tell us how robust our procedure is (and let us know if we need to improve the conditioning) */
    double ConditionNumber=0, ConditionNumber_threshold = 10. * CONDITION_NUMBER_DANGER; /* set a threshold condition number - above this we will 'pre-condition' the matrix for better behavior */
    double trace_initial = NV_T[0][0] + NV_T[1][1] + NV_T[2][2]; /* initial trace of this symmetric, positive-definite matrix; used below as a characteristic value for adding the identity */
    double conditioning_term_to_add = 1.05 * (trace_initial / NUMDIMS) / ConditionNumber_threshold; /* this will be added as a test value if the code does not reach the desired condition number */
    /* now enter an iterative loop to arrive at a -well-conditioned- inversion to use */
    while(1)
    {
        /* initialize the matrix this will go into */
        ConditionNumber = matrix_invert_ndims(NV_T, Tinv); // compute the matrix inverse, and return the condition number
        if(ConditionNumber < ConditionNumber_threshold) {break;} // end loop if we have reached target conditioning for the matrix
        for(j=0;j<NUMDIMS;j++) {NV_T[j][j] += conditioning_term_to_add;} /* add the conditioning term which should make the matrix better-conditioned for subsequent use: this is a normalization times the identity matrix in the relevant number of dimensions */
        conditioning_term_to_add *= 1.2; /* multiply the conditioning term so it will grow and eventually satisfy our criteria */
    } // end of loop broken when condition number is sufficiently small
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {P[i].NV_T[j][k] = Tinv[j][k] / dimensional_NV_T_normalizer;}} // now P[i].NV_T holds the inverted matrix elements //
#ifdef CBE_DEBUG
    if((ThisTask==0)&&(ConditionNumber>1.0e10)) {printf("Condition number == %g (Task=%d i=%d)\n",ConditionNumber,ThisTask,i);}
#endif
    return ConditionNumber;
}

#endif





/* ------------------------------------------------------------------------------------------------------
 Everything below here is a giant block to define the sub-routines needed to calculate additional force
  terms for particle types that do not fall into the 'hydro' category.
 -------------------------------------------------------------------------------------------------------- */
#define CORE_FUNCTION_NAME AGSForce_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(AGSForce_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

struct kernel_AGSForce
{
    double dp[3], dv[3], r, wk_i, wk_j, dwk_i, dwk_j, h_i, hinv_i, hinv3_i, hinv4_i, h_j, hinv_j, hinv3_j, hinv4_j;
};

/* structure for variables needed in evaluation sub-routines which must be passed from particles (sent to other processors) */
struct INPUT_STRUCT_NAME
{
    double Mass;
    double AGS_Hsml;
    double Pos[3];
    double Vel[3];
    int NodeList[NODELISTLENGTH];
    int Type;
    double dtime;
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    MyLongDouble NV_T[3][3];
    double V_i;
#endif
#if defined(DM_FUZZY)
    double AGS_Gradients_Density[3], AGS_Gradients2_Density[3][3], AGS_Numerical_QuantumPotential;
#if (DM_FUZZY > 0)
    double AGS_Psi_Re, AGS_Gradients_Psi_Re[3], AGS_Gradients2_Psi_Re[3][3];
    double AGS_Psi_Im, AGS_Gradients_Psi_Im[3], AGS_Gradients2_Psi_Im[3][3];
#endif
#endif
#if defined(DM_SIDM)
    double dtime_sidm;
    MyIDType ID;
#ifdef GRAIN_COLLISIONS
    double Grain_CrossSection_PerUnitMass;
#endif
#endif
}
*DATAIN_NAME, *DATAGET_NAME;

/* routine to pass particle information to the actual evaluation sub-routines */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    in->Mass = PPP[i].Mass;
    in->AGS_Hsml = PPP[i].AGS_Hsml;
    in->Type = P[i].Type;
    in->dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
    int k,k2; k=0; k2=0;
    for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    for(k=0;k<3;k++) {in->Vel[k] = P[i].Vel[k];}
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    in->V_i = get_particle_volume_ags(i);
    for(k=0;k<3;k++) {for(k2=0;k2<3;k2++) {in->NV_T[k][k2] = P[i].NV_T[k][k2];}}
#endif
#if defined(DM_FUZZY)
    for(k=0;k<3;k++) {in->AGS_Gradients_Density[k] = P[i].AGS_Gradients_Density[k];}
    for(k=0;k<3;k++) {for(k2=0;k2<3;k2++) {in->AGS_Gradients2_Density[k][k2] = P[i].AGS_Gradients2_Density[k][k2];}}
    in->AGS_Numerical_QuantumPotential = P[i].AGS_Numerical_QuantumPotential;
#if (DM_FUZZY > 0)
    in->AGS_Psi_Re = P[i].AGS_Psi_Re_Pred * P[i].AGS_Density / P[i].Mass;
    for(k=0;k<3;k++) {in->AGS_Gradients_Psi_Re[k] = P[i].AGS_Gradients_Psi_Re[k];}
    for(k=0;k<3;k++) {for(k2=0;k2<3;k2++) {in->AGS_Gradients2_Psi_Re[k][k2] = P[i].AGS_Gradients2_Psi_Re[k][k2];}}
    in->AGS_Psi_Im = P[i].AGS_Psi_Im_Pred * P[i].AGS_Density / P[i].Mass;
    for(k=0;k<3;k++) {in->AGS_Gradients_Psi_Im[k] = P[i].AGS_Gradients_Psi_Im[k];}
    for(k=0;k<3;k++) {for(k2=0;k2<3;k2++) {in->AGS_Gradients2_Psi_Im[k][k2] = P[i].AGS_Gradients2_Psi_Im[k][k2];}}
#endif
#endif
#if defined(DM_SIDM)
    in->dtime_sidm = P[i].dtime_sidm;
    in->ID = P[i].ID;
#ifdef GRAIN_COLLISIONS
    in->Grain_CrossSection_PerUnitMass = return_grain_cross_section_per_unit_mass(i);
#endif
#endif
}


/* structure for variables which must be returned -from- the evaluation sub-routines */
struct OUTPUT_STRUCT_NAME
{
#if defined(DM_SIDM)
    double sidm_kick[3], dtime_sidm; int si_count;
#endif
#ifdef DM_FUZZY
    double acc[3], AGS_Dt_Numerical_QuantumPotential;
#if (DM_FUZZY > 0)
    double AGS_Dt_Psi_Re, AGS_Dt_Psi_Im, AGS_Dt_Psi_Mass;
#endif
#endif
}
*DATARESULT_NAME, *DATAOUT_NAME;

#define ASSIGN_ADD_PRESET(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))

static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    int k,k2; k=0; k2=0;
#if defined(DM_SIDM)
    for(k=0;k<3;k++) {P[i].Vel[k] += out->sidm_kick[k];}
    MIN_ADD(P[i].dtime_sidm, out->dtime_sidm, mode);
    P[i].NInteractions += out->si_count;
#endif
#ifdef DM_FUZZY
    for(k=0;k<3;k++) {P[i].GravAccel[k] += out->acc[k];}
    ASSIGN_ADD_PRESET(P[i].AGS_Dt_Numerical_QuantumPotential,out->AGS_Dt_Numerical_QuantumPotential,mode);
#if (DM_FUZZY > 0)
    ASSIGN_ADD_PRESET(P[i].AGS_Dt_Psi_Re,out->AGS_Dt_Psi_Re,mode);
    ASSIGN_ADD_PRESET(P[i].AGS_Dt_Psi_Im,out->AGS_Dt_Psi_Im,mode);
    ASSIGN_ADD_PRESET(P[i].AGS_Dt_Psi_Mass,out->AGS_Dt_Psi_Mass,mode);
#endif
#endif
}


/* routine to determine if we need to apply the additional AGS-Force calculation[s] */
int AGSForce_isactive(int i);
int AGSForce_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0; /* check our 'marker' for particles which have finished iterating to an Hsml solution (if they have, dont do them again) */
#ifdef DM_SIDM
    if((1 << P[i].Type) & (DM_SIDM)) return 1;
#endif
#if defined(DM_FUZZY) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    if(P[i].Type == 1) return 1;
#endif
    return 0; // default to no-action, need to affirm calculation above //
}


/*!   -- this subroutine writes to shared memory [updating the neighbor values]: need to protect these writes for openmp below. none of the modified values are read, so only the write block is protected. note the writes can occur in the called code-blocks, so need to make sure they are followed so everything can be carefully constructed */
int AGSForce_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    /* zero memory and import data for local target */
    int startnode, numngb_inbox, listindex = 0, j, k, n; double r2, u_i, u_j;
    struct kernel_AGSForce kernel; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); memset(&kernel, 0, sizeof(struct kernel_AGSForce));
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(local.Mass <= 0 || local.AGS_Hsml <= 0) return 0;
    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.AGS_Hsml; kernel_hinv(kernel.h_i, &kernel.hinv_i, &kernel.hinv3_i, &kernel.hinv4_i);
    int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(local.Type); // determine allowed particle types for search for adaptive gravitational softening terms
#if defined(DM_SIDM)
    out.dtime_sidm = local.dtime_sidm;
#endif
    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            double search_len = local.AGS_Hsml;
#if defined(DM_SIDM)
            search_len *= 3.0; // need a 'buffer' because we will consider interactions with any kernel -overlap, not just inside one or the other kernel radius
#endif
            numngb_inbox = ngb_treefind_pairs_threads_targeted(local.Pos, search_len, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, AGS_kernel_shared_BITFLAG);
            if(numngb_inbox < 0) {return -2;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if((P[j].Mass <= 0)||(PPP[j].AGS_Hsml <= 0)) continue; /* make sure neighbor is valid */
                /* calculate position relative to target */
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0]; kernel.dp[1] = local.Pos[1] - P[j].Pos[1]; kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); /*  now find the closest image in the given box size  */
                r2 = kernel.dp[0]*kernel.dp[0] + kernel.dp[1]*kernel.dp[1] + kernel.dp[2]*kernel.dp[2];
                if(r2 <= 0) continue;
                kernel.r = sqrt(r2);
                kernel.h_j = PPP[j].AGS_Hsml;
#if defined(DM_SIDM)
                if(kernel.r > kernel.h_i+kernel.h_j) continue;
#else
                if(kernel.r > kernel.h_i && kernel.r > kernel.h_j) continue;
#endif
                /* calculate kernel quantities needed below */
                kernel_hinv(kernel.h_j, &kernel.hinv_j, &kernel.hinv3_j, &kernel.hinv4_j);
                u_i = kernel.r * kernel.hinv_i; u_j = kernel.r * kernel.hinv_j;
                if(u_i < 1) {kernel_main(u_i, kernel.hinv3_i, kernel.hinv4_i, &kernel.wk_i, &kernel.dwk_i, 0);} else {kernel.wk_i=kernel.dwk_i=0;}
                if(u_j < 1) {kernel_main(u_j, kernel.hinv3_j, kernel.hinv4_j, &kernel.wk_j, &kernel.dwk_j, 0);} else {kernel.wk_j=kernel.dwk_j=0;}
                for(k=0;k<3;k++)
                {
                    kernel.dv[k] = local.Vel[k] - P[j].Vel[k];
                    if(All.ComovingIntegrationOn) {kernel.dv[k] += All.cf_hubble_a * kernel.dp[k]/All.cf_a2inv;}
                }
                
#ifdef DM_FUZZY
#include "../sidm/dm_fuzzy_flux_computation.h"
#endif
#if defined(DM_SIDM)
#include "../sidm/sidm_core_flux_computation.h"
#endif

            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}



void AGSForce_calc(void)
{
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    PRINT_STATUS(" ..entering AGS-Force calculation [as hydro loop for non-gas elements]\n");
    /* before doing any operations, need to zero the appropriate memory so we can correctly do pair-wise operations */
#if defined(DM_SIDM)
    {int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {P[i].dtime_sidm = 10.*GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);}}
#endif
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    /* do final operations on results: these are operations that can be done after the complete set of iterations */
    /* collect timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm; CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */


#endif // AGS_HSML_CALCULATION_IS_ACTIVE
