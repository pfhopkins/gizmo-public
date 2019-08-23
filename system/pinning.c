#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <unistd.h>

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org).
 */

#ifdef IMPOSE_PINNING

#define __USE_GNU

#include <sched.h>

#include "../allvars.h"
#include "../proto.h"

#ifndef SOCKETS
#define SOCKETS   4           /* Setting for our 4 x AMD Magny-Cours nodes */
#define MAX_CORES 48          /* Setting for our 4 x AMD Magny-Cours nodes */
#endif

 /* SOCKETS   = 12 for 2 x Intel 6-core nodes with hypethreading  */
 /* MAX_CORES = 24 for 2 x Intel 6-core nodes with hypethreading  */

static cpu_set_t cpuset;


void get_core_set(void)
{
  CPU_ZERO(&cpuset);
  sched_getaffinity(getpid(), sizeof(cpuset), &cpuset);
}


void pin_to_core_set(void)
{
  int core, task, num_threads;
  static cpu_set_t cpuset_new;

  char *p = getenv("OMP_NUM_THREADS");
  if(p)
    num_threads = atoi(p);
  else
    num_threads = 1;

#ifdef PTHREADS_NUM_THREADS
  if(num_threads != PTHREADS_NUM_THREADS)
    {
      terminate("You have activated PTHREADS_NUM_THREADS, but the value of the environment variable OMP_NUM_THREADS is not equal to PTHREADS_NUM_THREADS!"); 
    }
#endif  


  CPU_ZERO(&cpuset_new);

  int corestart = 0;

  for(task=0, core= corestart; task < NTask * num_threads; task++)
    {
      while(!CPU_ISSET(core, &cpuset))
	{
	  core += SOCKETS;
	  if(core >= MAX_CORES)
	    {
	      corestart++;
	      if(corestart >= SOCKETS)
		corestart = 0;

	      core = corestart;
	    }
	}

      if((task / num_threads) == ThisTask)
	CPU_SET(core, &cpuset_new);

      core += SOCKETS;
      if(core >= MAX_CORES)
	{
	  corestart++;
	  if(corestart >= SOCKETS)
	    corestart = 0;

	  core = corestart;
	}
    }

  sched_setaffinity (getpid(), sizeof(cpuset_new), &cpuset_new);
}


void report_pinning(void)
{
  cpu_set_t cpuset;
  int i;
  char buf[MAX_CORES+1];

  CPU_ZERO (&cpuset);
  sched_getaffinity (getpid(), sizeof(cpuset), &cpuset);

  for(i=0;i<MAX_CORES;i++)
    if(CPU_ISSET(i, &cpuset))
      buf[i]='1';
    else
      buf[i]='-';
  buf[MAX_CORES]=0;

#ifndef IO_REDUCED_MODE
  for(i=0; i<NTask; i++)
    {
      if(ThisTask == i)
	printf("Task=%02d: %s\n", ThisTask, buf);
      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
}

#endif
