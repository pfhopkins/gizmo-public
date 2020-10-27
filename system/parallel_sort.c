#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>

#include "../allvars.h"
#include "../proto.h"

/*
 * This file was originally part of the GADGET3 code by Volker Springel.
 */

#define TRANSFER_SIZE_LIMIT  1000000000

#define TAG_TRANSFER  100

static void serial_sort(char *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
static void msort_serial_with_tmp(char *base, size_t n, size_t s, int (*compar) (const void *, const void *),
				  char *t);
static void get_local_rank(char *element,
			   size_t tie_braking_rank,
			   char *base,
			   size_t nmemb, size_t size, size_t noffs_thistask,
			   long long left, long long right,
			   size_t * loc, int (*compar) (const void *, const void *));



void parallel_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *))
{
  parallel_sort_comm(base, nmemb, size, compar, MPI_COMM_WORLD);
}

void parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *), MPI_Comm comm)
{
  int i, j, max_task, ranks_not_found, Local_ThisTask, Local_NTask, Local_PTask, Color;
  MPI_Comm MPI_CommLocal;

  /* we create a communicator that contains just those tasks with nmemb > 0. This makes 
   *  it easier to deal with CPUs that do not hold any data.
   */

  if(nmemb)
    {
      Color = 1;
      serial_sort((char*)base, nmemb, size, compar);
    }
  else
    Color = 0;

  MPI_Comm_split(comm, Color, 0, &MPI_CommLocal);
  MPI_Comm_rank(MPI_CommLocal, &Local_ThisTask);
  MPI_Comm_size(MPI_CommLocal, &Local_NTask);

  if(Local_NTask > 1 && Color == 1)
    {
      for(Local_PTask = 0; Local_NTask > (1 << Local_PTask); Local_PTask++);

      size_t *nlist = (size_t *) mymalloc("nlist", Local_NTask * sizeof(size_t));
      size_t *noffs = (size_t *) mymalloc("noffs", Local_NTask * sizeof(size_t));

      MPI_Allgather(&nmemb, sizeof(size_t), MPI_BYTE, nlist, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

      for(i = 1, noffs[0] = 0; i < Local_NTask; i++) noffs[i] = noffs[i - 1] + nlist[i - 1];

      char *element_guess = (char *) mymalloc("element_guess", (Local_NTask - 1) * size);
      size_t *element_tie_braking_rank = (size_t *) mymalloc("element_tie_braking_rank", (Local_NTask - 1) * sizeof(size_t));
      size_t *desired_glob_rank = (size_t *) mymalloc("desired_glob_rank", (Local_NTask - 1) * sizeof(size_t));
      size_t *current_glob_rank = (size_t *) mymalloc("current_glob_rank", (Local_NTask - 1) * sizeof(size_t));
      size_t *current_loc_rank = (size_t *) mymalloc("current_loc_rank", (Local_NTask - 1) * sizeof(size_t));
      long long *range_left = (long long *) mymalloc("range_left", (Local_NTask - 1) * sizeof(long long));
      long long *range_right = (long long *) mymalloc("range_right", (Local_NTask - 1) * sizeof(long long));
      int *max_loc = (int *) mymalloc("max_loc", (Local_NTask - 1) * sizeof(int));
      size_t max_nmemb;

      /* find the largest nmemb value and initialize the guesses with the median element found there */
      for(i = 0, max_task = max_nmemb = 0; i < Local_NTask; i++)
	if(max_nmemb < nlist[i])
	  {
	    max_nmemb = nlist[i];
	    max_task = i;
	  }

      if(Local_ThisTask == max_task) memcpy(element_guess, (char *) base + size * (max_nmemb / 2), size);

      MPI_Bcast(element_guess, size, MPI_BYTE, max_task, MPI_CommLocal);

      for(i = 1; i < Local_NTask - 1; i++) memcpy(element_guess + i * size, element_guess, size);

      for(i = 0; i < Local_NTask - 1; i++)
	{
	  desired_glob_rank[i] = noffs[i + 1];
	  element_tie_braking_rank[i] = (max_nmemb / 2) + noffs[max_task];
	  max_loc[i] = max_task;

	  range_left[i] = 0;	/* first element that it can be */
	  range_right[i] = nmemb;	/* first element that it can not be */

	  current_glob_rank[i] = 0;
	}

      int iter = 0;

      do
	{
	  for(i = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])
		{
		  get_local_rank(element_guess + i * size, element_tie_braking_rank[i],
				 (char *) base, nmemb, size, noffs[Local_ThisTask],
				 range_left[i], range_right[i], &current_loc_rank[i], compar);
		}
	    }

	  /* now compute the global ranks by summing the local ranks */
	  size_t *list = (size_t *) mymalloc("list", Local_NTask * sizeof(size_t));
	  for(i = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])
		{
		  MPI_Allgather(&current_loc_rank[i], sizeof(size_t), MPI_BYTE, list, sizeof(size_t),
				MPI_BYTE, MPI_CommLocal);

		  for(j = 0, current_glob_rank[i] = 0; j < Local_NTask; j++)
		    current_glob_rank[i] += list[j];
		}
	    }
	  myfree(list);

	  for(i = 0, ranks_not_found = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])	/* here we're not yet done */
		{
		  ranks_not_found++;

		  if(current_glob_rank[i] < desired_glob_rank[i])
		    {
		      range_left[i] = current_loc_rank[i];

		      if(Local_ThisTask == max_loc[i])	// && range_right[i] == range_left[i]+1)
			{
			  range_left[i]++;
			}
		    }

		  if(current_glob_rank[i] > desired_glob_rank[i])
		    range_right[i] = current_loc_rank[i];
		}
	    }

	  /* now we need to determine new element guesses */

	  long long *range_len_list =
	    (long long *) mymalloc("range_len_list", Local_NTask * sizeof(long long));

	  for(i = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])	/* here we're not yet done */
		{
		  long long max_range_len, range_len = range_right[i] - range_left[i] + 1;
		  MPI_Allgather(&range_len, sizeof(long long), MPI_BYTE, range_len_list, sizeof(long long),
				MPI_BYTE, MPI_CommLocal);

		  /* find the largest one, and the cpu that holds it */
		  for(j = 0, max_loc[i] = max_range_len = 0; j < Local_NTask; j++)
		    if(max_range_len < range_len_list[j])
		      {
			max_range_len = range_len_list[j];
			max_loc[i] = j;
		      }

		  if(Local_ThisTask == max_loc[i])
		    {
		      long long mid = (range_left[i] + range_right[i]) / 2;

		      memcpy(element_guess + i * size, (char *) base + mid * size, size);
		      element_tie_braking_rank[i] = mid + noffs[Local_ThisTask];
		    }

		  MPI_Bcast(element_guess + i * size, size, MPI_BYTE, max_loc[i], MPI_CommLocal);
		  MPI_Bcast(&element_tie_braking_rank[i], sizeof(size_t), MPI_BYTE, max_loc[i],
			    MPI_CommLocal);
		}
	    }

	  myfree(range_len_list);

	  iter++;

	  if(iter > 800 + 3*NTask && Local_ThisTask == 0)
	    {
	      printf("iter=%d: ranks_not_found=%d  Local_NTask=%d\n", iter, ranks_not_found, Local_NTask); fflush(stdout);
	      if(iter > 900 + 3*NTask) terminate("can't find the split points. That's odd");
	    }
	}
      while(ranks_not_found);


      /* At this point we have found all the elements corresponding to the desired split points */
      /* we can now go ahead and determine how many elements of the local CPU have to go to each other CPU */

      size_t *Send_count = (size_t *) mymalloc("Send_count", Local_NTask * sizeof(size_t));
      size_t *Recv_count = (size_t *) mymalloc("Recv_count", Local_NTask * sizeof(size_t));
      size_t *Send_offset = (size_t *) mymalloc("Send_offset", Local_NTask * sizeof(size_t));
      size_t *Recv_offset = (size_t *) mymalloc("Recv_offset", Local_NTask * sizeof(size_t));

      for(i = 0; i < Local_NTask; i++)
	Send_count[i] = 0;

      int target = 0;

      for(i = 0; i < (int)nmemb; i++)
	{
	  while(target < Local_NTask - 1)
	    {
	      int cmp = compar((char *) base + i * size, element_guess + target * size);
	      if(cmp == 0)
		{
		  if(i + noffs[Local_ThisTask] < element_tie_braking_rank[target])
		    cmp = -1;
		  else if(i + noffs[Local_ThisTask] > element_tie_braking_rank[target])
		    cmp = +1;
		}
	      if(cmp >= 0)
		target++;
	      else
		break;
	    }
	  Send_count[target]++;
	}

      MPI_Alltoall(Send_count, sizeof(size_t), MPI_BYTE, Recv_count, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

      size_t nimport;

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < Local_NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      if(nimport != nmemb)
	terminate("nimport != nmemb");

      char *basetmp = (char *) mymalloc("basetmp", nmemb * size);

      int ngrp, recvTask;

      /* exchange the data */
      for(ngrp = 0; ngrp < (1 << Local_PTask); ngrp++)
	{
	  recvTask = Local_ThisTask ^ ngrp;

	  if(recvTask < Local_NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
#ifndef MPISENDRECV_SIZELIMIT
		  if(Send_count[recvTask] > TRANSFER_SIZE_LIMIT || Recv_count[recvTask] > TRANSFER_SIZE_LIMIT)
		    terminate("we are above TRANSFER_SIZE_LIMIT");
#endif
		  MPI_Sendrecv((char *) base + Send_offset[recvTask] * size,
			       Send_count[recvTask] * size, MPI_BYTE,
			       recvTask, TAG_TRANSFER,
			       basetmp + Recv_offset[recvTask] * size,
			       Recv_count[recvTask] * size, MPI_BYTE,
			       recvTask, TAG_TRANSFER, MPI_CommLocal, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* instead of doing another serial sort at the end, we could also exploit that the different incoming pieces
         are already individually sorted. This can be used in a variant of merge sort to do the final step a bit quicker */

      memcpy(base, basetmp, nmemb * size);
      myfree(basetmp);
      serial_sort((char*)base, nmemb, size, compar);


      myfree(Recv_offset);
      myfree(Send_offset);
      myfree(Recv_count);
      myfree(Send_count);

      myfree(max_loc);
      myfree(range_right);
      myfree(range_left);
      myfree(current_loc_rank);
      myfree(current_glob_rank);
      myfree(desired_glob_rank);
      myfree(element_tie_braking_rank);
      myfree(element_guess);
      myfree(noffs);
      myfree(nlist);
    }

  MPI_Comm_free(&MPI_CommLocal);
}



static void get_local_rank(char *element,	/* element of which we want the rank */
			   size_t tie_braking_rank,	/* the inital global rank of this element (needed for braking ties) */
			   char *base,	/* base address of local data */
			   size_t nmemb, size_t size,	/* number and size of local data */
			   size_t noffs_thistask,	/* cumulative length of data on lower tasks */
			   long long left, long long right,	/* range of elements on local task that may hold the element */
			   size_t * loc,	/* output: local rank of the element */
			   int (*compar) (const void *, const void *))	/* user-specified  comparison function */
{
  if(right < left)
    terminate("right < left");

  if(left == 0 && right == (int)(nmemb + 1))
    {
      if(compar(base + (nmemb - 1) * size, element) < 0)
	{
	  *loc = nmemb;
	  return;
	}
      else if(compar(base, element) > 0)
	{
	  *loc = 0;
	  return;
	}
    }


  if(right == left)		/* looks like we already converged to the proper rank */
    {
      *loc = left;
    }
  else
    {
      if(compar(base + (right - 1) * size, element) < 0)	/* the last element is smaller, hence all elements are on the left */
	*loc = (right - 1) + 1;
      else if(compar(base + left * size, element) > 0)	/* the first element is already larger, hence no element is on the left */
	*loc = left;
      else
	{
	  while(right > left)
	    {
	      long long mid = ((right - 1) + left) / 2;

	      int cmp = compar(base + mid * size, element);
	      if(cmp == 0)
		{
		  if(mid + noffs_thistask < tie_braking_rank)
		    cmp = -1;
		  else if(mid + noffs_thistask > tie_braking_rank)
		    cmp = +1;
		}

	      if(cmp == 0)	/* element has exactly been found */
		{
		  *loc = mid;
		  break;
		}

	      if((right - 1) == left)	/* elements is not on this CPU */
		{
		  if(cmp < 0)
		    *loc = mid + 1;
		  else
		    *loc = mid;
		  break;
		}

	      if(cmp < 0)
		{
		  left = mid + 1;
		}
	      else
		{
		  if((right - 1) == left + 1)
		    {
		      if(mid != left)
			{
			  printf("-->left=%lld  right=%lld\n", left, right);
			  terminate("can't be");
			}
		      *loc = left;
		      break;
		    }

		  right = mid;
		}
	    }
	}
    }
}

static void serial_sort(char *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *))
{
  const size_t storage = nmemb * size;

  char *tmp = (char *) mymalloc("char *tmp", storage);

  msort_serial_with_tmp(base, nmemb, size, compar, tmp);

  myfree(tmp);
}


static void msort_serial_with_tmp(char *base, size_t n, size_t s, int (*compar) (const void *, const void *), char *t)
{
  char *tmp;
  char *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = base;
  b2 = base + n1 * s;

  msort_serial_with_tmp(b1, n1, s, compar, t);
  msort_serial_with_tmp(b2, n2, s, compar, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(compar(b1, b2) < 0)
	{
	  --n1;
	  memcpy(tmp, b1, s);
	  tmp += s;
	  b1 += s;
	}
      else
	{
	  --n2;
	  memcpy(tmp, b2, s);
	  tmp += s;
	  b2 += s;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * s);

  memcpy(base, t, (n - n2) * s);
}





int compare_densities_for_sort(const void *a, const void *b)
{
    double x = SphP[ *(int *) a].Density;
    double y = SphP[ *(int *) b].Density;
    if (x < y) {return -1;} else if(x > y) {return 1;}
    return 0;
}





/* simple recursive function to efficiently find index in arbitrarily spaced bin-list */
int binarySearch(const double * arr, const double x, const int l, const int r, const int total)
{
  if(r<l){ endrun(7777);}
  const int w     = r-l;
  const int mid   = l + w/2;
  if(w <= 1) {
    if(mid < total-1) {if(x >= arr[mid] && x < arr[mid+1]) {return mid;}}
    if(mid >0) {if(x >= arr[mid-1] && x < arr[mid]) {return mid-1;}}
  } // endif w < 1
  if(arr[mid] > x) {
     return binarySearch(arr,x,l,mid-1,total);
  } else if(arr[mid]<x) {
    return binarySearch(arr,x,mid+1,r,total);
  } else {
    return mid;
  }
  return -1;
}
