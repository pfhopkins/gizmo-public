/* This is a generic code block designed for simple neighbor loops, so that they don't have to 
    be copy-pasted and can be generically optimized in a single place. specifically this is for
    the secondary loop of particles on a remote processor (after the primary has been passed)

    EVALUATION_CALL is the actual call, and needs to be defined appropriately, or this will crash
 */
#if !defined(EVALUATION_CALL)
printf("Cannot compile the secondary sub-loop without EVALUATION_CALL defined. Exiting. \n"); fflush(stdout); exit(995534);
#endif
int j, dummy, *ngblist, thread_id = *(int *) p;
ngblist = Ngblist + thread_id * NumPart;
while(1)
{
    LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
    {
        j = NextJ;
        NextJ++;
    }
    UNLOCK_NEXPORT;
    if(j >= Nimport) {break;}
    EVALUATION_CALL
}
/* loop completed successfully */
return NULL;
