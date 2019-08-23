#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>


#include "../allvars.h"
#include "../proto.h"

/* various routines are collected here, needed to communicate with the 
 *  actual system, stdin/out, abort runs, etc
 */

/*
 * This file (and the relevant routines) were originally part of the GADGET3 
 * code developed by Volker Springel (volker.springel@h-its.org). The routines 
 * have been collected and re-arranged a bit, but not substantially modified, by 
 * Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/*  This function aborts the simulations. If a single processors
 *  wants an immediate termination,  the function needs to be
 *  called with ierr>0. A bunch of MPI-error messages will also
 *  appear in this case.
 *  For ierr=0, MPI is gracefully cleaned up, but this requires
 *  that all processors call endrun().
 */
void endrun(int ierr)
{
    if(ierr)
    {
        printf("task %d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, ierr);
        exit(0);
    }
    
    MPI_Finalize();
    exit(0);
}


#ifdef DEBUG
#include <fenv.h>
void enable_core_dumps_and_fpu_exceptions(void)
{
  struct rlimit rlim;
  extern int feenableexcept(int __excepts);

  /* enable floating point exceptions */

  /*
     feenableexcept(FE_DIVBYZERO | FE_INVALID);
   */

  /* Note: FPU exceptions appear not to work properly
   * when the Intel C-Compiler for Linux is used
   */

  /* set core-dump size to infinity */
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  /* MPICH catches the signales SIGSEGV, SIGBUS, and SIGFPE....
   * The following statements reset things to the default handlers,
   * which will generate a core file.
   */
  /*
     signal(SIGSEGV, catch_fatal);
     signal(SIGBUS, catch_fatal);
     signal(SIGFPE, catch_fatal);
     signal(SIGINT, catch_fatal);
   */

  signal(SIGSEGV, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGFPE, SIG_DFL);
  signal(SIGINT, SIG_DFL);

  /* Establish a handler for SIGABRT signals. */
  signal(SIGABRT, catch_abort);
}


void catch_abort(int sig)
{
  MPI_Finalize();
  exit(0);
}

void catch_fatal(int sig)
{
  terminate_processes();
  MPI_Finalize();

  signal(sig, SIG_DFL);
  raise(sig);
}


void terminate_processes(void)
{
  pid_t my_pid;
  char buf[500], hostname[500], *cp;
  char commandbuf[500];
  FILE *fd;
  int i, pid;

  sprintf(buf, "%s%s", All.OutputDir, "PIDs.txt");

  my_pid = getpid();

  if((fd = fopen(buf, "r")))
    {
      for(i = 0; i < NTask; i++)
	{
	  int ret;

	  ret = fscanf(fd, "%s %d", hostname, &pid);

	  cp = hostname;
	  while(*cp)
	    {
	      if(*cp == '.')
		*cp = 0;
	      else
		cp++;
	    }

	  if(my_pid != pid)
	    {
	      sprintf(commandbuf, "ssh %s kill -ABRT %d", hostname, pid);
	      printf("--> %s\n", commandbuf);
	      fflush(stdout);
#ifndef NOCALLSOFSYSTEM
	      ret = system(commandbuf);
#endif
	    }
	}

      fclose(fd);
    }
}

void write_pid_file(void)
{
  pid_t my_pid;
  char mode[8], buf[500];
  FILE *fd;
  int i;

  my_pid = getpid();

  sprintf(buf, "%s%s", All.OutputDir, "PIDs.txt");

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  for(i = 0; i < NTask; i++)
    {
      if(ThisTask == i)
	{
	  if(ThisTask == 0)
	    sprintf(mode, "w");
	  else
	    sprintf(mode, "a");

	  if((fd = fopen(buf, mode)))
	    {
	      fprintf(fd, "%s %d\n", getenv("HOST"), (int) my_pid);
	      fclose(fd);
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }
}
#endif

#ifdef PAUSE_RUN_TO_ATTACH_DEBUGGER
void pause_run_to_attach_debugger()
{
  int continue_run = 0, my_pid;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  my_pid = getpid();
  printf("PID %d on %s ready for debugger to attach\n", my_pid, hostname);
  fflush(stdout);
  
  int *all_pids;
  if (ThisTask == 0) all_pids = (int *) malloc(sizeof(int)*NTask);
  
  MPI_Gather(&my_pid, 1, MPI_INT, all_pids, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (ThisTask == 0)
  {
    FILE *fd = fopen("pid_list_for_debugger.txt","w");
    int j;
    for (j=0; j<NTask; ++j) fprintf(fd, "%i \n", all_pids[j]);
    fclose(fd);   
    free(all_pids);
    printf("PID file written\n");
    fflush(stdout);
  }
  
  while (continue_run == 0) /* continue_run needsto be set to 1 ("set var continue_run = 1") by debugger to continue */
    sleep(1);
}
#endif


double get_random_number(MyIDType id)
{
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
    return RndTable[(int) (id % RNDTABLE)];
#else
    return gsl_rng_uniform(random_generator);
#endif
}

void set_random_numbers(void)
{
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
  int i;
  for(i = 0; i < RNDTABLE; i++)
    RndTable[i] = gsl_rng_uniform(random_generator);
#endif
}


/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double my_second(void)
{
#ifdef WALLCLOCK
  return MPI_Wtime();
#else
  return ((double) clock()) / CLOCKS_PER_SEC;
#endif

  /* note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}

double measure_time(void)	/* strategy: call this at end of functions to account for time in this function, and before another (nontrivial) function is called */
{
  double t, dt;

  t = my_second();
  dt = t - WallclockTime;
  WallclockTime = t;

  return dt;
}

double report_time(void)       /* strategy: call this to measure sub-times of functions*/
{
  double t, dt;

  t = my_second();
  dt = t - WallclockTime;

  return dt;
}


/* returns the time difference between two measurements
 * obtained with my_second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)			/* overflow has occured (for systems with 32bit tick counter) */
    {
#ifdef WALLCLOCK
      dt = 0;
#else
      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}




#ifdef X86FIX

#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#define _FPU_EXTENDED 0x0300
#define _FPU_DOUBLE   0x0200

void x86_fix(void)
{
  unsigned short dummy, new_cw;
  unsigned short *old_cw;

  old_cw = &dummy;

  _FPU_GETCW(*old_cw);
  new_cw = (*old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
}

#endif


void minimum_large_ints(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
                MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = src[j];

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      if(res[j] > numlist[i * n + j])
        res[j] = numlist[i * n + j];

  myfree(numlist);
}

void sumup_large_ints(int n, int *src, long long *res)
{
  int i, j, *numlist;

  numlist = (int *) mymalloc("numlist", NTask * n * sizeof(int));
  MPI_Allgather(src, n, MPI_INT, numlist, n, MPI_INT, MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}

void sumup_longs(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
		MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}

size_t sizemax(size_t a, size_t b)
{
  if(a < b)
    return b;
  else
    return a;
}



long long report_comittable_memory(long long *MemTotal,
				   long long *Committed_AS, long long *SwapTotal, long long *SwapFree)
{
  FILE *fd;
  char buf[1024];

  if((fd = fopen("/proc/meminfo", "r")))
    {
      while(1)
	{
	  if(fgets(buf, 500, fd) != buf)
	    break;

	  if(bcmp(buf, "MemTotal", 8) == 0)
	    {
	      *MemTotal = atoll(buf + 10);
	    }
	  if(strncmp(buf, "Committed_AS", 12) == 0)
	    {
	      *Committed_AS = atoll(buf + 14);
	    }
	  if(strncmp(buf, "SwapTotal", 9) == 0)
	    {
	      *SwapTotal = atoll(buf + 11);
	    }
	  if(strncmp(buf, "SwapFree", 8) == 0)
	    {
	      *SwapFree = atoll(buf + 10);
	    }
	}
      fclose(fd);
    }

  return (*MemTotal - *Committed_AS);
}

void mpi_report_comittable_memory(long long BaseMem)
{
  long long *sizelist, maxsize[6], minsize[6];
  double avgsize[6];
  int i, imem, mintask[6], maxtask[6];
  long long Mem[6];
  char label[512];

  Mem[0] = report_comittable_memory(&Mem[1], &Mem[2], &Mem[3], &Mem[4]);
  Mem[5] = Mem[1] - Mem[0];

  for(imem = 0; imem < 6; imem++)
    {
      sizelist = (long long *) malloc(NTask * sizeof(long long));
      MPI_Allgather(&Mem[imem], sizeof(long long), MPI_BYTE, sizelist, sizeof(long long), MPI_BYTE,
		    MPI_COMM_WORLD);

      for(i = 1, mintask[imem] = 0, maxtask[imem] = 0, maxsize[imem] = minsize[imem] =
	  sizelist[0], avgsize[imem] = sizelist[0]; i < NTask; i++)
	{
	  if(sizelist[i] > maxsize[imem])
	    {
	      maxsize[imem] = sizelist[i];
	      maxtask[imem] = i;
	    }
	  if(sizelist[i] < minsize[imem])
	    {
	      minsize[imem] = sizelist[i];
	      mintask[imem] = i;
	    }
	  avgsize[imem] += sizelist[i];
	}

      free(sizelist);
    }

  if(ThisTask == 0)
    {
      printf("-------------------------------------------------------------------------------------------\n");
      for(imem = 0; imem < 6; imem++)
	{
	  switch (imem)
	    {
	    case 0:
	      sprintf(label, "AvailMem");
	      break;
	    case 1:
	      sprintf(label, "Total Mem");
	      break;
	    case 2:
	      sprintf(label, "Committed_AS");
	      break;
	    case 3:
	      sprintf(label, "SwapTotal");
	      break;
	    case 4:
	      sprintf(label, "SwapFree");
	      break;
	    case 5:
	      sprintf(label, "AllocMem");
	      break;
	    }
	  printf
	    ("%s:\t Largest = %10.2f Mb (on task=%d), Smallest = %10.2f Mb (on task=%d), Average = %10.2f Mb\n",
	     label, maxsize[imem] / (1024.), maxtask[imem], minsize[imem] / (1024.), mintask[imem],
	     avgsize[imem] / (1024. * NTask));
	}
      printf("-------------------------------------------------------------------------------------------\n");
    }
  if(ThisTask == maxtask[2])
    {
      printf("Task with the maximum commited memory");
      system("echo $HOST");
    }


  fflush(stdout);
}
