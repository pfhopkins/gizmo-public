/* This is a generic code block designed for simple neighbor loops, so that they don't have to 
    be copy-pasted and can be generically optimized in a single place. specifically this is for
    the initial loop of particles on the local processor (and determination of which need passing) 

   Two blocks need to be defined or this will crash: 
   CONDITION_FOR_EVALUATION inserts the clause that actually determines
        whether or not to pass a particle to the main evaluation routine
   EVALUATION_CALL is the actual call, and needs to be written appropriately
 */
#if !defined(CONDITION_FOR_EVALUATION) || !defined(EVALUATION_CALL)
printf("Cannot compile the primary sub-loop without both CONDITION_FOR_EVALUATION and EVALUATION_CALL defined. Exiting. \n"); fflush(stdout); exit(995533);
#endif
/* variable assignment */
int i, j, *exportflag, *exportnodecount, *exportindex, *ngblist, thread_id = *(int *) p;
/* define the pointers needed for each thread to speak back regarding what needs processing */
ngblist = Ngblist + thread_id * NumPart;
exportflag = Exportflag + thread_id * NTask;
exportnodecount = Exportnodecount + thread_id * NTask;
exportindex = Exportindex + thread_id * NTask;
/* Note: exportflag is local to each thread */
for(j = 0; j < NTask; j++) {exportflag[j] = -1;}
/* now begin the actual loop */
while(1)
{
    int exitFlag = 0;
    LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
    {
        if(BufferFullFlag != 0 || NextParticle < 0)
        {
            exitFlag = 1;
        }
        else
        {
            i = NextParticle;
            ProcessedFlag[i] = 0;
            NextParticle = NextActiveParticle[NextParticle];
        }
    }
    UNLOCK_NEXPORT;
    if(exitFlag) {break;}
    CONDITION_FOR_EVALUATION
    {
        if(EVALUATION_CALL < 0) {break;} // export buffer has filled up //
    }
    ProcessedFlag[i] = 1; /* particle successfully finished */
}
/* loop completed successfully */
return NULL;
