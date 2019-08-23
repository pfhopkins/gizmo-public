#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"



/*! \file main.c
 *  \brief start of the program
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

/*!
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0. Then begrun() is called, which sets up
 *  the simulation either from IC's or from restart files.  Finally,
 *  run() is started, the main simulation loop, which iterates over
 *  the timesteps.
 */
int main(int argc, char **argv)
{
  int i;

#ifdef IMPOSE_PINNING
  get_core_set();
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

#ifdef IMPOSE_PINNING
  pin_to_core_set();
#endif

  mpi_report_comittable_memory(0);
  MPI_Barrier(MPI_COMM_WORLD);

  /* initialize OpenMP thread pool and bind (implicitly though OpenMP runtime) */
  if(ThisTask == 0)
    {
      char *username = getenv("USER");
      char hostname[201]; hostname[200] = '\0';
      int have_hn = gethostname(hostname,200);
      time_t rawtime;
      struct tm * timeinfo;
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );

      printf("\nSystem time: %s", asctime(timeinfo) );
      printf("This is GIZMO, version %s, running on %s as %s.\n",
              GIZMO_VERSION,
              have_hn == 0 ? hostname : "?",
              username ? username : "?"
      );
#ifdef BUILDINFO
      printf(BUILDINFO", " __DATE__ " " __TIME__ "\n");
#endif
      printf("\nCode was compiled with settings:\n\n");
      output_compile_time_options();
   }

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp master
    {
      maxThreads = omp_get_num_threads();
      if(ThisTask == 0)
	printf("Using %d OpenMP threads\n", maxThreads);
    }
  }
#elif defined(PTHREADS_NUM_THREADS)
  if(ThisTask == 0)
    printf("Using %d POSIX threads\n", maxThreads);
#endif

  for(PTask = 0; NTask > (1 << PTask); PTask++);

  if(argc < 2)
    {
      if(ThisTask == 0)
	{
	  printf("Parameters are missing.\n");
	  printf("Call with <ParameterFile> [<RestartFlag>] [<RestartSnapNum>]\n");
	  printf("\n");
	  printf("   RestartFlag    Action\n");
	  printf("       0          Read initial conditions and start simulation\n");
	  printf("       1          Read restart files and resume simulation\n");
	  printf("       2          Restart from specified snapshot dump and continue simulation\n");
	  printf("       3          Run FOF and optionally SUBFIND if enabled\n");
	  printf("       4          Convert snapshot file to different format\n");
	  printf("       5          Calculate power spectrum and two-point function\n");
	  printf("       6          Calculate velocity power spectrum for the gas particles\n");
	  printf("\n");
	}
      endrun(0);
    }

  strcpy(ParameterFile, argv[1]);

  if(argc >= 3)
    RestartFlag = atoi(argv[2]);
  else
    RestartFlag = 0;

  if(argc >= 4)
    RestartSnapNum = atoi(argv[3]);
  else
    RestartSnapNum = -1;

  /* initialize CPU-time/Wallclock-time measurement */
  for(i = 0; i < CPU_PARTS; i++)
    All.CPU_Sum[i] = CPU_Step[i] = 0;

  CPUThisRun = 0;
  WallclockTime = my_second();

  begrun();			/* set-up run  */

  run();			/* main simulation loop */

  MPI_Finalize();		/* clean up & finalize MPI */

  return 0;
}
