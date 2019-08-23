/** \file
    MPI utility functions.
*/
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part (cleaned up, some routines re-organized and consolidated and a 
 * couple others added) by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <string.h>
#include "../allvars.h"
#include "../proto.h"


/** Calculates the recv_count, send_offset, and recv_offset arrays
    based on the send_count. Returns nimportbytes, the total number of
    bytes to be received. If an identical set of copies are to be
    sent to all tasks, set send_identical=1 and the send_offset will
    be zero for all tasks.

    All arrays should be allocated with NTask size. */
int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical)
{
  // Exchange the send/receive counts
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);
  
  int nimportbytes = 0;
  recv_offset[0] = 0;
  send_offset[0] = 0;
  int j;
  for(j = 0; j < NTask; j++)
    {
      nimportbytes += recv_count[j];

      if(j > 0)
	    {
	      send_offset[j] = send_offset[j - 1] + (send_identical ? 0 : send_count[j - 1]);
	      recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
	    }
    }
  return nimportbytes;
}


/** Compare function used to sort an array of int pointers into order
    of the pointer targets. */
int intpointer_compare(const void *a, const void *b)
{
  if((**(int **) a) < (**(int **) b))
    return -1;

  if((**(int **) a) > (**(int **) b))
    return +1;

  return 0;
}


/** Sort an opaque array into increasing order of an int field, given
    by the specified offset. (This would typically be field indicating
    the task.) Returns a sorted copy of the data array, that needs to
    be myfreed.

    We do this by sorting an array of pointers to the task field, and
    then using this array to deduce the reordering of the data
    array. Unfortunately this means making a copy of the data, but
    this just replaces the copy after the mpi_exchange_buffers
    anyway.  */
void sort_based_on_field(void *data, int field_offset, int n_items, int item_size, void **data2ptr)
{
  int i;
  int **perm;
  *data2ptr = (char *)mymalloc_movable(data2ptr, "data2", n_items * item_size);
  perm = (int **)mymalloc("perm", n_items * sizeof(*perm));

  for(i = 0; i < n_items; ++i)
    perm[i] = (int *) ((char *) data + i * item_size + field_offset);

  qsort(perm, n_items, sizeof(*perm), intpointer_compare);
  // reorder data into data2
  for(i = 0; i < n_items; ++i)
    {
      size_t orig_pos = ((char *) perm[i] - ((char *) data + field_offset)) / item_size;
      if (!(((char *) perm[i] - ((char *) data + field_offset)) % item_size == 0)) terminate("something wrong here!");
      memcpy((char*)*data2ptr + item_size * i, (char *) data + item_size * orig_pos, item_size);
    }

  myfree(perm);
}


/** This function distributes the members in an opaque structure to
    the tasks based on a task field given by a specified offset into
    the opaque struct. The task field must have int type. n_items is
    updated to the new size of data. max_n is the allocated size of
    the data array, and is updated if a realloc is necessary.  */
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size)
{
  int i;

  for(i = 0; i < NTask; i++) Send_count[i] = 0;

  for(i = 0; i < *n_items; i++)   /* find Send_count for each task */
    {
      int task = *((int *) ((char *) data + i * item_size + task_offset));
      if (!(task >= 0 && task < NTask)) terminate("wrong task number!");
      Send_count[task] += item_size;
    }

  void *data2;
  sort_based_on_field(data, task_offset, *n_items, item_size, &data2);   /* sort data based on target task */

  int nimportbytes = mpi_calculate_offsets(Send_count, Send_offset, Recv_count, Recv_offset, 0);   /* calculate offsets of send and receive buffers */
  if (nimportbytes % item_size !=0) terminate("should not happen!");
  int nimport = nimportbytes/item_size;
  
  if(*max_n < nimport)   /* realloc data if too small */
    {
      data = (char *)myrealloc_movable(data, nimportbytes);
      *max_n = nimport;
    }

  MPI_Alltoallv(data2, Send_count, Send_offset, MPI_BYTE, data, Recv_count, Recv_offset, MPI_BYTE, MPI_COMM_WORLD);   /* exchange data */
  
  myfree(data2);

  *n_items = nimport;
}



#ifdef MPISENDRECV_CHECKSUM

#undef MPI_Sendrecv

int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status)
{
    int checksumtag = 1000, errtag = 2000;
    int i, iter = 0, err_flag, err_flag_imported, size_sendtype, size_recvtype, Local_ThisTask, Local_NTask;
    long long sendCheckSum, recvCheckSum, importedCheckSum;
    unsigned char *p, *buf, *recvbuf;
    
    if(dest != source)
        endrun(3);
    
    MPI_Comm_rank(comm, &Local_ThisTask);
    MPI_Comm_size(comm, &Local_NTask);
    
    MPI_Type_size(sendtype, &size_sendtype);
    MPI_Type_size(recvtype, &size_recvtype);
    
    if(dest == Local_ThisTask)
    {
        memcpy(recvbufreal, sendbuf, recvcount * size_recvtype);
        return 0;
    }
    
    
    if(!(buf = mymalloc("buf", recvcount * size_recvtype + 1024)))
        endrun(6);
    
    for(i = 0, p = buf; i < recvcount * size_recvtype + 1024; i++)
        *p++ = 255;
    
    recvbuf = buf + 512;
    
    MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
                 recvbuf, recvcount, recvtype, source, recvtag, comm, status);
    
    for(i = 0, p = buf; i < 512; i++, p++)
    {
        if(*p != 255)
        {
            printf
            ("MPI-ERROR: Task=%d/%s: Recv occured before recv buffer. message-size=%d from %d, i=%d c=%d\n",
             Local_ThisTask, getenv("HOST"), recvcount, dest, i, *p);
            fflush(stdout);
            endrun(6);
        }
    }
    
    for(i = 0, p = recvbuf + recvcount * size_recvtype; i < 512; i++, p++)
    {
        if(*p != 255)
        {
            printf
            ("MPI-ERROR: Task=%d/%s: Recv occured after recv buffer. message-size=%d from %d, i=%d c=%d\n",
             Local_ThisTask, getenv("HOST"), recvcount, dest, i, *p);
            fflush(stdout);
            endrun(6);
        }
    }
    
    
    for(i = 0, p = sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
        sendCheckSum += *p;
    
    importedCheckSum = 0;
    
    if(dest > Local_ThisTask)
    {
        if(sendcount > 0)
            MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
        if(recvcount > 0)
            MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, comm, status);
    }
    else
    {
        if(recvcount > 0)
            MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, comm, status);
        if(sendcount > 0)
            MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
    }
    
    checksumtag++;
    
    for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount * size_recvtype; i++, p++)
        recvCheckSum += *p;
    
    
    err_flag = err_flag_imported = 0;
    
    if(recvCheckSum != importedCheckSum)
    {
        printf
        ("MPI-ERROR: Receive error on task=%d/%s from task=%d, message size=%d, sendcount=%d checksums= %d %d  %d %d. Try to fix it...\n",
         Local_ThisTask, getenv("HOST"), source, recvcount, sendcount, (int) (recvCheckSum >> 32),
         (int) recvCheckSum, (int) (importedCheckSum >> 32), (int) importedCheckSum);
        fflush(stdout);
        
        err_flag = 1;
    }
    
    if(dest > Local_ThisTask)
    {
        MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
        MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
    }
    else
    {
        MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
        MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
    }
    errtag++;
    
    if(err_flag > 0 || err_flag_imported > 0)
    {
        printf("Task=%d is on %s, wants to send %d and has checksum=%d %d of send data\n",
               Local_ThisTask, getenv("HOST"), sendcount, (int) (sendCheckSum >> 32), (int) sendCheckSum);
        
        do
        {
            sendtag++;
            recvtag++;
            
            for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
                *p = 0;
            
            if((iter & 1) == 0)
            {
                if(dest > Local_ThisTask)
                {
                    if(sendcount > 0)
                        MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, comm);
                    if(recvcount > 0)
                        MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, comm, status);
                }
                else
                {
                    if(recvcount > 0)
                        MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, comm, status);
                    if(sendcount > 0)
                        MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, comm);
                }
            }
            else
            {
                if(iter > 5)
                {
                    printf("we're trying to send each byte now on task=%d (iter=%d)\n", Local_ThisTask, iter);
                    if(dest > Local_ThisTask)
                    {
                        for(i = 0, p = sendbuf; i < sendcount * size_sendtype; i++, p++)
                            MPI_Ssend(p, 1, MPI_BYTE, dest, i, comm);
                        for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
                            MPI_Recv(p, 1, MPI_BYTE, dest, i, comm, status);
                    }
                    else
                    {
                        for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
                            MPI_Recv(p, 1, MPI_BYTE, dest, i, comm, status);
                        for(i = 0, p = sendbuf; i < sendcount * size_sendtype; i++, p++)
                            MPI_Ssend(p, 1, MPI_BYTE, dest, i, comm);
                    }
                }
                else
                {
                    MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
                                 recvbuf, recvcount, recvtype, source, recvtag, comm, status);
                }
            }
            
            importedCheckSum = 0;
            
            for(i = 0, p = sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
                sendCheckSum += *p;
            
            printf("Task=%d gas send_checksum=%d %d\n", Local_ThisTask, (int) (sendCheckSum >> 32), (int) sendCheckSum);
#ifndef IO_REDUCED_MODE
            fflush(stdout);
#endif
            if(dest > Local_ThisTask)
            {
                if(sendcount > 0)
                    MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
                if(recvcount > 0)
                    MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag,
                             comm, status);
            }
            else
            {
                if(recvcount > 0)
                    MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag,
                             comm, status);
                if(sendcount > 0)
                    MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
            }
            
            for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount; i++, p++)
                recvCheckSum += *p;
            
            err_flag = err_flag_imported = 0;
            
            if(recvCheckSum != importedCheckSum)
            {
                printf
                ("MPI-ERROR: Again (iter=%d) a receive error on task=%d/%s from task=%d, message size=%d, checksums= %d %d  %d %d. Try to fix it...\n",
                 iter, Local_ThisTask, getenv("HOST"), source, recvcount, (int) (recvCheckSum >> 32),
                 (int) recvCheckSum, (int) (importedCheckSum >> 32), (int) importedCheckSum);
                fflush(stdout);
                err_flag = 1;
            }
            
            if(dest > Local_ThisTask)
            {
                MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
                MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
            }
            else
            {
                MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
                MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
            }
            
            if(err_flag == 0 && err_flag_imported == 0)
                break;
            
            errtag++;
            checksumtag++;
            iter++;
        }
        while(iter < 10);
        
#ifndef IO_REDUCED_MODE
        if(iter >= 10)
        {
            char buf[1000];
            int length;
            FILE *fd;
            
            sprintf(buf, "send_data_%d.dat", Local_ThisTask);
            fd = fopen(buf, "w");
            length = sendcount * size_sendtype;
            fwrite(&length, 1, sizeof(int), fd);
            fwrite(sendbuf, sendcount, size_sendtype, fd);
            fclose(fd);
            
            sprintf(buf, "recv_data_%d.dat", Local_ThisTask);
            fd = fopen(buf, "w");
            length = recvcount * size_recvtype;
            fwrite(&length, 1, sizeof(int), fd);
            fwrite(recvbuf, recvcount, size_recvtype, fd);
            fclose(fd);
            
            printf("MPI-ERROR: Even 10 trials proved to be insufficient on task=%d/%s. Stopping\n", Local_ThisTask, getenv("HOST"));
            fflush(stdout);
            endrun(10);
        }
#endif
    }
    
    memcpy(recvbufreal, recvbuf, recvcount * size_recvtype);
    
    myfree(buf);
    
    return 0;
}

#endif


#ifdef MPISENDRECV_SIZELIMIT

#undef MPI_Sendrecv


int MPI_Sizelimited_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             int dest, int sendtag, void *recvbuf, int recvcount,
                             MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                             MPI_Status * status)
{
    int iter = 0, size_sendtype, size_recvtype, send_now, recv_now;
    int count_limit;
    
    
    if(dest != source)
        endrun(3);
    
    MPI_Type_size(sendtype, &size_sendtype);
    MPI_Type_size(recvtype, &size_recvtype);
    
    if(dest == ThisTask)
    {
        memcpy(recvbuf, sendbuf, recvcount * size_recvtype);
        return 0;
    }
    
    count_limit = (int) ((((long long) MPISENDRECV_SIZELIMIT) * 1024 * 1024) / size_sendtype);
    
    while(sendcount > 0 || recvcount > 0)
    {
        if(sendcount > count_limit)
        {
            send_now = count_limit;
            if(iter == 0)
            {
                printf("imposing size limit on MPI_Sendrecv() on task=%d (send of size=%d)\n", ThisTask, sendcount * size_sendtype);
                //fflush(stdout);
            }
            iter++;
        }
        else
            send_now = sendcount;
        
        if(recvcount > count_limit)
            recv_now = count_limit;
        else
            recv_now = recvcount;
        
        MPI_Sendrecv(sendbuf, send_now, sendtype, dest, sendtag,
                     recvbuf, recv_now, recvtype, source, recvtag, comm, status);
        
        sendcount -= send_now;
        recvcount -= recv_now;
        
        sendbuf += send_now * size_sendtype;
        recvbuf += recv_now * size_recvtype;
    }
    
    return 0;
}

#endif
