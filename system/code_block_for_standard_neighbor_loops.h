/* This is a generic code block designed for simple neighbor loops, so that they don't have to
 be copy-pasted and can be generically optimized in a single place */

#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

void *PRIMARY_SUBFUN_NAME(void *p, int loop_iteration);
void *SECONDARY_SUBFUN_NAME(void *p, int loop_iteration);
int EVALUATION_WORKHORSE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration);

void *PRIMARY_SUBFUN_NAME(void *p, int loop_iteration)
{
#define CONDITION_FOR_EVALUATION if(CONDITIONFUNCTION_FOR_EVALUATION(i))
#define EVALUATION_CALL EVALUATION_WORKHORSE_FUNCTION_NAME(i,0,exportflag,exportnodecount,exportindex,ngblist,loop_iteration)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}

void *SECONDARY_SUBFUN_NAME(void *p, int loop_iteration)
{
#define EVALUATION_CALL EVALUATION_WORKHORSE_FUNCTION_NAME(j, 1, &dummy, &dummy, &dummy, ngblist, loop_iteration);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}

void MASTER_FUNCTION_NAME(void)
{
    /* define the number of loop iterations needed for the physics of interest, then loop over those iterations */
    int loop_iteration, number_of_loop_iterations = 1;
    for(loop_iteration=0; loop_iteration<number_of_loop_iterations; loop_iteration++)
    {
        int i, j, k, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
        double timeall=0, timecomp1=0, timecomp2=0, timecommsumm1=0, timecommsumm2=0, timewait1=0, timewait2=0;
        double timecomp, timecomm, timewait, tstart, tend, t0, t1; long long n_exported = 0, NTaskTimesNumPart;
        /* allocate buffers to arrange communication */
        NTaskTimesNumPart = maxThreads * NumPart; size_t MyBufferSize = All.BufferSize;
        All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                               sizeof(struct INPUT_STRUCT_NAME) + sizeof(struct OUTPUT_STRUCT_NAME) + sizemax(sizeof(struct INPUT_STRUCT_NAME),sizeof(struct OUTPUT_STRUCT_NAME))));
        CPU_Step[CPU_MISC] += measure_time(); t0 = my_second();
        Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
        DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
        DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
        
        NextParticle = FirstActiveParticle;    /* begin the main loop; start with this index */
        do /* do local particles and prepare export list */
        {
            BufferFullFlag = 0; Nexport = 0; save_NextParticle = NextParticle;
            for(j = 0; j < NTask; j++) {Send_count[j] = 0; Exportflag[j] = -1;}
            tstart = my_second();
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
#ifdef _OPENMP
                int mainthreadid = omp_get_thread_num();
#else
                int mainthreadid = 0;
#endif
                PRIMARY_SUBFUN_NAME(&mainthreadid, loop_iteration);    /* do local particles and prepare export list */
            }
            tend = my_second(); timecomp1 += timediff(tstart, tend);
            if(BufferFullFlag)
            {
                int last_nextparticle = NextParticle; NextParticle = save_NextParticle;
                while(NextParticle >= 0)
                {
                    if(NextParticle == last_nextparticle) {break;}
                    if(ProcessedFlag[NextParticle] != 1) {break;}
                    ProcessedFlag[NextParticle] = 2; NextParticle = NextActiveParticle[NextParticle];
                }
                if(NextParticle == save_NextParticle) {endrun(123708);} /* in this case, the buffer is too small to process even a single particle */
                int new_export = 0;
                for(j = 0, k = 0; j < Nexport; j++)
                    if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                    {
                        if(k < j + 1) {k = j + 1;}
                        for(; k < Nexport; k++)
                            if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                            {
                                int old_index = DataIndexTable[j].Index; DataIndexTable[j] = DataIndexTable[k];
                                DataNodeList[j] = DataNodeList[k]; DataIndexTable[j].IndexGet = j; new_export++;
                                DataIndexTable[k].Index = old_index; k++; break;
                            }
                    }
                    else {new_export++;}
                Nexport = new_export;
            }
            n_exported += Nexport;
            for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
            for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
            MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
            tstart = my_second();
            MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
            tend = my_second();
            timewait1 += timediff(tstart, tend);
            for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
            {
                Nimport += Recv_count[j];
                if(j > 0)
                {
                    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                    Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }
            /* prepare particle data for export */
            DATAGET_NAME = (struct INPUT_STRUCT_NAME *) mymalloc("DATAGET_NAME", Nimport * sizeof(struct INPUT_STRUCT_NAME));
            DATAIN_NAME = (struct INPUT_STRUCT_NAME *) mymalloc("DATAIN_NAME", Nexport * sizeof(struct INPUT_STRUCT_NAME));
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index; INPUTFUNCTION_NAME(&DATAIN_NAME[j], place);
                memcpy(DATAIN_NAME[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
            }
            /* exchange particle data */
            tstart = my_second();
            for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                        /* get the particles */
                        MPI_Sendrecv(&DATAIN_NAME[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct INPUT_STRUCT_NAME), MPI_BYTE, recvTask, TAG_GRADLOOP_A,
                                     &DATAGET_NAME[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct INPUT_STRUCT_NAME), MPI_BYTE,
                                     recvTask, TAG_GRADLOOP_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second(); timecommsumm1 += timediff(tstart, tend); myfree(DATAIN_NAME);
            DATARESULT_NAME = (struct OUTPUT_STRUCT_NAME *) mymalloc("DATARESULT_NAME", Nimport * sizeof(struct OUTPUT_STRUCT_NAME));
            DATAOUT_NAME = (struct OUTPUT_STRUCT_NAME *) mymalloc("DATAOUT_NAME", Nexport * sizeof(struct OUTPUT_STRUCT_NAME));
            /* now do the particles that were sent to us */
            tstart = my_second(); NextJ = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
#ifdef _OPENMP
                int mainthreadid = omp_get_thread_num();
#else
                int mainthreadid = 0;
#endif
                SECONDARY_SUBFUN_NAME(&mainthreadid, loop_iteration);
            }
            tend = my_second(); timecomp2 += timediff(tstart, tend);
            if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
            tstart = my_second();
            MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            tend = my_second(); timewait2 += timediff(tstart, tend);
            /* get the result */
            tstart = my_second();
            for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                        /* send the results */
                        MPI_Sendrecv(&DATARESULT_NAME[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME), MPI_BYTE, recvTask, TAG_GRADLOOP_B,
                                     &DATAOUT_NAME[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME),
                                     MPI_BYTE, recvTask, TAG_GRADLOOP_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second(); timecommsumm2 += timediff(tstart, tend);
            /* add the result to the local particles */
            tstart = my_second();
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index;
                OUTPUTFUNCTION_NAME(&DATAOUT_NAME[j], place, 1, loop_iteration);
            }
            tend = my_second(); timecomp1 += timediff(tstart, tend);
            myfree(DATAOUT_NAME); myfree(DATARESULT_NAME); myfree(DATAGET_NAME);
        }
        while(ndone < NTask);
        myfree(DataNodeList); myfree(DataIndexTable); myfree(Ngblist);
        
        /* do final operations on results: these are operations that can be done after the complete set of iterations */
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        {
            FINAL_OPERATIONS_FUNCTION_NAME(i, loop_iteration);
        }
        MPI_Barrier(MPI_COMM_WORLD); // force barrier so we know the first derivatives are fully-computed //
        
        /* collect some timing information */
        t1 = WallclockTime = my_second(); timeall += timediff(t0, t1); timecomp = timecomp1 + timecomp2; timewait = timewait1 + timewait2; timecomm = timecommsumm1 + timecommsumm2;
        CPU_Step[CPU_COST_CODE_NAME] += timecomp; CPU_Step[CPU_COST_CODE_NAME] += timewait; CPU_Step[CPU_COST_CODE_NAME] += timecomm;
        CPU_Step[CPU_COST_CODE_NAME] += timeall - (timecomp + timewait + timecomm);
    } // close loop over loop_iterations (master loop)
}


int EVALUATION_WORKHORSE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, n; /* define variables */
    struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* zero memory and import data for local target */
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target);} else {local = DATAGET_NAME[target];}
    if(local.Hsml <= 0) {return 0;} /* check if we should bother doing a neighbor loop */
    int kernel_shared_BITFLAG = KERNEL_BITFLAG_DEFINITION_LOCAL;
    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.Hsml, target, &startnode, mode, exportflag,
                                                                  exportnodecount, exportindex, ngblist, kernel_shared_BITFLAG);
            if(numngb_inbox < 0) {return -1;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; if((P[j].Mass <= 0)||(P[j].Hsml <= 0)) {continue;} /* make sure neighbor is valid */
                NEIGHBOROPS_FUNCTION_NAME(&local, &out, j, loop_iteration);
            } // numngb_inbox loop
        } // while(startnode)
        /* continue to open leaves if needed */
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DATAGET_NAME[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
            }
        }
    }
    /* Collect the result at the right place */
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
    return 0;
}

#undef CPU_COST_CODE_NAME
#undef KERNEL_BITFLAG_DEFINITION_LOCAL
#undef FINAL_OPERATIONS_FUNCTION_NAME
#undef NEIGHBOROPS_FUNCTION_NAME
#undef CONDITIONFUNCTION_FOR_EVALUATION
#undef SECONDARY_SUBFUN_NAME
#undef PRIMARY_SUBFUN_NAME
#undef EVALUATION_WORKHORSE_FUNCTION_NAME
#undef OUTPUTFUNCTION_NAME
#undef DATARESULT_NAME
#undef OUTPUT_STRUCT_NAME
#undef DATAOUT_NAME
#undef INPUTFUNCTION_NAME
#undef DATAGET_NAME
#undef DATAIN_NAME
#undef INPUT_STRUCT_NAME
#undef MASTER_FUNCTION_NAME
