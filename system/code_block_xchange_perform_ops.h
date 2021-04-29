/* This is a generic code block designed for simple neighbor loops, so that they don't have to
be copy-pasted and can be generically optimized in a single place */
{
    int j, k, ndone=0, ndone_flag=0, recvTask, place, save_NextParticle; long long n_exported = 0; double tstart, tend, tstart_loop; /* define some variables used only below */
    NextParticle = FirstActiveParticle;    /* begin the main loop; start with this index */
    tstart_loop = my_second();
    do /* primary point-element loop */
    {
        BufferFullFlag = 0; Nexport = 0; save_NextParticle = NextParticle; tstart = my_second();
        for(j = 0; j < NTask; j++) {Send_count[j] = 0; Exportflag[j] = -1;} /* do local particles and prepare export list */
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
        tend = my_second(); timecomp += timediff(tstart, tend);
        if(BufferFullFlag) /* we've filled the buffer or reached the end of the list, prepare for communications */
        {
            int last_nextparticle = NextParticle; NextParticle = save_NextParticle; /* figure out where we are */
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle) {break;}
                if(ProcessedFlag[NextParticle] != 1) {break;}
                ProcessedFlag[NextParticle] = 2; NextParticle = NextActiveParticle[NextParticle];
            }
            if(NextParticle == save_NextParticle)
            {
                PRINT_WARNING("NextParticle == save_NextParticle condition (the buffer appears too small to hold a single particle): NextParticle=%d save_NextParticle=%d last_nextparticle=%d ProcessedFlag[NextParticle]=%d NextActiveParticle[NextParticle]=%d NumPart=%d N_gas=%d NTaskTimesNumPart=%llu maxThreads=%d All.BunchSize=%ld All.BufferSize=%llu Nexport=%ld ndone=%d ndone_flag=%d NTask=%d",NextParticle,save_NextParticle,last_nextparticle,ProcessedFlag[NextParticle],NextActiveParticle[NextParticle],NumPart,N_gas,(unsigned long long)NTaskTimesNumPart,maxThreads,All.BunchSize,(unsigned long long)All.BufferSize,Nexport,ndone,ndone_flag,NTask);
                if(NextParticle >= 0) {PRINT_WARNING("This is a live particle: NextParticle=%d ID=%llu Mass=%g Type=%d",NextParticle,(unsigned long long)P[NextParticle].ID,P[NextParticle].Mass,P[NextParticle].Type);}
                printf("Extended Debug: Printing Processed Flag for Entire Active Particle Chain on This Task: \n"); int nj=0; for(j=FirstActiveParticle;j>=0;j=NextActiveParticle[j]) {printf("nj=%d j=%d ProcFlag[j]=%d \n",nj,j,ProcessedFlag[j]); nj++; fflush(stdout);} fflush(stdout);
                endrun(113312);
            } /* in this case, the buffer is too small to process even a single particle */
            
            int new_export = 0; /* actually calculate exports [so we can tell other tasks] */
            for(j = 0, k = 0; j < Nexport; j++)
            {
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1) {k = j + 1;}
                    for(; k < Nexport; k++)
                        if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                        {
                            int old_index = DataIndexTable[j].Index;
                            DataIndexTable[j] = DataIndexTable[k]; DataNodeList[j] = DataNodeList[k]; DataIndexTable[j].IndexGet = j; new_export++;
                            DataIndexTable[k].Index = old_index; k++;
                            break;
                        }
                }
                else {new_export++;}
            }
            Nexport = new_export; /* counting exports... */
        }
        n_exported += Nexport;
        for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
        for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare); /* construct export count tables */
        tstart = my_second();
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD); /* broadcast import/export counts */
        tend = my_second(); timewait += timediff(tstart, tend);

        for(j = 0, Send_offset[0] = 0; j < NTask; j++) {if(j > 0) {Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];}} /* calculate export table offsets */
        DATAIN_NAME = (struct INPUT_STRUCT_NAME *) mymalloc("DATAIN_NAME", Nexport * sizeof(struct INPUT_STRUCT_NAME));
        DATAOUT_NAME = (struct OUTPUT_STRUCT_NAME *) mymalloc("DATAOUT_NAME", Nexport * sizeof(struct OUTPUT_STRUCT_NAME));
        for(j = 0; j < Nexport; j++) /* prepare particle data for export [fill in the structures to be passed] */
        {
            place = DataIndexTable[j].Index;
            INPUTFUNCTION_NAME(&DATAIN_NAME[j], place, loop_iteration);
            memcpy(DATAIN_NAME[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

        /* ok now we have to figure out if there is enough memory to handle all the tasks sending us their data, and if not, break it into sub-chunks */
        int N_chunks_for_import, ngrp_initial, ngrp;
        for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) /* sub-chunking loop opener */
        {
            int flagall;
            N_chunks_for_import = (1 << PTask) - ngrp_initial;
            do {
                int flag = 0; Nimport = 0;
                for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++)
                {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask) {if(Recv_count[recvTask] > 0) {Nimport += Recv_count[recvTask];}}
                }
                size_t space_needed = Nimport * sizeof(struct INPUT_STRUCT_NAME) + Nimport * sizeof(struct OUTPUT_STRUCT_NAME) + 16384; /* extra bitflag is a padding, to avoid overflows */
                if(space_needed > FreeBytes) {flag = 1;}
                
                MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                if(flagall) {N_chunks_for_import /= 2;} else {break;}
            } while(N_chunks_for_import > 0);
            if(N_chunks_for_import == 0) {printf("Memory is insufficient for even one import-chunk: N_chunks_for_import=%d  ngrp_initial=%d  Nimport=%ld  FreeBytes=%lld , but we need to allocate=%lld \n",N_chunks_for_import, ngrp_initial, Nimport, (long long)FreeBytes,(long long)(Nimport * sizeof(struct INPUT_STRUCT_NAME) + Nimport * sizeof(struct OUTPUT_STRUCT_NAME) + 16384)); endrun(9977);}
            if(flagall) {if(ThisTask==0) PRINT_WARNING("Splitting import operation into sub-chunks as we are hitting memory limits (check this isn't imposing large communication cost)");}

            /* now allocated the import and results buffers */
            DATAGET_NAME = (struct INPUT_STRUCT_NAME *) mymalloc("DATAGET_NAME", Nimport * sizeof(struct INPUT_STRUCT_NAME));
            DATARESULT_NAME = (struct OUTPUT_STRUCT_NAME *) mymalloc("DATARESULT_NAME", Nimport * sizeof(struct OUTPUT_STRUCT_NAME));

            tstart = my_second(); Nimport = 0; /* reset because this will be cycled below to calculate the recieve offsets (Recv_offset) */
            for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) /* exchange particle data */
            {
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) /* get the particles */
                    {
                        MPI_Sendrecv(&DATAIN_NAME[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct INPUT_STRUCT_NAME), MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_A,
                                     &DATAGET_NAME[Nimport], Recv_count[recvTask] * sizeof(struct INPUT_STRUCT_NAME), MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_A,
                                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        Nimport += Recv_count[recvTask];
                    }
                }
            }
            tend = my_second(); timecomm += timediff(tstart, tend);
            
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
            tend = my_second(); timecomp += timediff(tstart, tend); tstart = my_second();
            MPI_Barrier(MPI_COMM_WORLD); /* insert MPI Barrier here - will be forced by comms below anyways but this allows for clean timing measurements */
            tend = my_second(); timewait += timediff(tstart, tend);
            
            tstart = my_second(); Nimport = 0;
            for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) /* send the results for imported elements back to their host tasks */
            {
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                        MPI_Sendrecv(&DATARESULT_NAME[Nimport], Recv_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME), MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_B,
                                     &DATAOUT_NAME[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME), MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_B,
                                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        Nimport += Recv_count[recvTask];
                    }
                }
            }
            tend = my_second(); timecomm += timediff(tstart, tend);
            myfree(DATARESULT_NAME); myfree(DATAGET_NAME); /* free the structures used to send data back to tasks, its sent */
            
        } /* close the sub-chunking loop: for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) */

        /* we have all our results back from the elements we exported: add the result to the local elements */
        tstart = my_second();
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            OUTPUTFUNCTION_NAME(&DATAOUT_NAME[j], place, 1, loop_iteration);
        }
        tend = my_second(); timecomp += timediff(tstart, tend);
        myfree(DATAOUT_NAME); myfree(DATAIN_NAME); /* free the structures used to prepare our initial export data, we're done here! */
        
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;} /* figure out if we are done with the particular active set here */
        tstart = my_second();
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /* call an allreduce to figure out if all tasks are also done here, otherwise we need to iterate */
        tend = my_second(); timewait += timediff(tstart, tend);
    }
    while(ndone < NTask);
    timeall += timediff(tstart_loop, my_second());
    
} /* closes clause, so variables don't 'leak' */

