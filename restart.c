#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/file.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"

static FILE *fd;

#ifdef CHIMES 
static ChimesFloat *sphAbundancesBuf;
#endif 

static void in(int *x, int modus);
static void byten(void *x, size_t n, int modus);

int old_MaxPart = 0, new_MaxPart;


/* This function reads or writes the restart files.
 * Each processor writes its own restart file, with the
 * I/O being done in parallel. To avoid congestion of the disks
 * you can tell the program to restrict the number of files
 * that are simultaneously written to NumFilesWrittenInParallel.
 *
 * If modus>0  the restart()-routine reads, 
 * if modus==0 it writes a restart file. 
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part (adding/removing read/write items, allowing for different variables 
 * to be changed or re-initialized on restarts, and changing variable units as necessary)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void restart(int modus)
{
    char buf[200], buf_bak[200], buf_mv[500];
    double save_PartAllocFactor;
    int nprocgroup, primaryTask, groupTask;
    struct global_data_all_processes all_task0;
    int nmulti = MULTIPLEDOMAINS, regular_restarts_are_valid = 1, backup_restarts_are_valid = 1;
    

#ifdef CHIMES 
    int partIndex, abunIndex; 
#endif 
    
    if(ThisTask == 0 && modus == 0) // writing re-start files: move old files to .bak
    {
        sprintf(buf, "%s/restartfiles", All.OutputDir);
        mkdir(buf, 02755);
#ifndef NOCALLSOFSYSTEM
        int i_Task_iter;
        for(i_Task_iter=0; i_Task_iter<NTask; i_Task_iter++)
        {
            sprintf(buf, "%s/restartfiles/%s.%d", All.OutputDir, All.RestartFile, i_Task_iter);
            sprintf(buf_bak, "%s/restartfiles/%s.%d.bak", All.OutputDir, All.RestartFile, i_Task_iter);
#ifdef IO_REDUNDANT_BACKUP_RESTARTFILE_FREQUENCY
            if( (((int)(CPUThisRun/All.CpuTimeBetRestartFile)) % ((int)IO_REDUNDANT_BACKUP_RESTARTFILE_FREQUENCY)) == 0)
            {
                char buf_bak2[200]; sprintf(buf_bak2, "%s/restartfiles/%s.%d.bak2", All.OutputDir, All.RestartFile, i_Task_iter);
                rename(buf_bak,buf_bak2); // move old backup restart files to .bak2 files //
            }
#endif
            rename(buf,buf_bak); // move old restart files to .bak files //
        }
#endif
    }
    if(modus == 1) // reading re-start files. make sure to check all the files to read exist!
    {
#ifndef NOCALLSOFSYSTEM
        int i_Task_iter;
        for(i_Task_iter=0; i_Task_iter<NTask; i_Task_iter++)
        {
            sprintf(buf, "%s/restartfiles/%s.%d", All.OutputDir, All.RestartFile, i_Task_iter);
            sprintf(buf_bak, "%s/restartfiles/%s.%d.bak", All.OutputDir, All.RestartFile, i_Task_iter);
            if(!(fd = fopen(buf, "r"))) {regular_restarts_are_valid=0;} else {fclose(fd);} // check if regular restart exists
            if(!(fd = fopen(buf_bak, "r"))) {backup_restarts_are_valid=0;} else {fclose(fd);} // check if backup restart exists
        }
#endif
        if(ThisTask == 0)
        {
            if((regular_restarts_are_valid == 0) && (backup_restarts_are_valid == 0))
            {
                printf("Fatal error. Full set of restart files ('%s' or '%s') not found - check the restarts are uncorrupted and your MPI process number has not changed.\n", buf, buf_bak);
                endrun(7871);
            }
            if((regular_restarts_are_valid == 0) && (backup_restarts_are_valid == 1))
            {
                printf("Default restartfiles ('%s') not found - they are incomplete or corrupted [number of files matching MPI process number not found. But apparently valid set of backup restartfiles ('%s') found. Attempting to use those.\n", buf, buf_bak);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    sprintf(buf, "%s/restartfiles/%s.%d", All.OutputDir, All.RestartFile, ThisTask);
    if((modus == 1) && (regular_restarts_are_valid == 0) && (backup_restarts_are_valid == 1)) {sprintf(buf, "%s/restartfiles/%s.%d.bak", All.OutputDir, All.RestartFile, ThisTask);}
    sprintf(buf_bak, "%s/restartfiles/%s.%d.bak", All.OutputDir, All.RestartFile, ThisTask);
    sprintf(buf_mv, "mv %s %s", buf, buf_bak);
    
    if((NTask < All.NumFilesWrittenInParallel))
    {
        printf("Fatal error.\nNumber of processors must be greater than or equal to `NumFilesWrittenInParallel'.\n");
        endrun(2131);
    }
    
    nprocgroup = NTask / All.NumFilesWrittenInParallel;
    
    if((NTask % All.NumFilesWrittenInParallel))
    {
        nprocgroup++;
    }

  primaryTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (primaryTask + groupTask))	/* ok, it's this processor's turn */
	{
	  if(modus)
	    {
	      if(!(fd = fopen(buf, "r")))
		{
		  if(!(fd = fopen(buf_bak, "r")))
		    {
		      printf("Restart file '%s' nor '%s' found.\n", buf, buf_bak);
		      endrun(7870);
		    }
		}
	    }
	  else
	    {
	      if(!(fd = fopen(buf, "w")))
		{
		  printf("Restart file '%s' cannot be opened.\n", buf);
		  endrun(7878);
		}
	    }


	  save_PartAllocFactor = All.PartAllocFactor;

	  /* common data  */
	  byten(&All, sizeof(struct global_data_all_processes), modus);

	  if(ThisTask == 0 && modus > 0)
	    all_task0 = All;

	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }


	  if(modus)		/* read */
	    {
	      if(All.PartAllocFactor != save_PartAllocFactor)
		{
		  old_MaxPart = All.MaxPart;	/* old MaxPart */

		  if(ThisTask == 0)
		    printf("PartAllocFactor changed: %f/%f , adapting bounds ...\n",
			   All.PartAllocFactor, save_PartAllocFactor);

		  All.PartAllocFactor = save_PartAllocFactor;
		  All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));
		  All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));
#ifdef ALLOW_IMBALANCED_GASPARTICLELOAD
          All.MaxPartSph = All.MaxPart; // PFH: increasing All.MaxPartSph according to this line can allow better load-balancing in some cases. however it leads to more memory problems
#endif
          new_MaxPart = All.MaxPart;

		  save_PartAllocFactor = -1;
		}

	      if(all_task0.Time != All.Time)
		{
		  printf("The restart file on task=%d is not consistent with the one on task=0\n", ThisTask);
		  fflush(stdout);
		  endrun(16);
		}

	      allocate_memory();
	    }

	  in(&NumPart, modus);
	  if(NumPart > All.MaxPart)
	    {
	      printf
		("it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		 NumPart / (((double) All.TotNumPart) / NTask));
	      printf("fatal error\n");
	      endrun(22);
	    }

	  if(modus)		/* read */
	    {
	      if(old_MaxPart)
		All.MaxPart = old_MaxPart;	/* such that tree is still valid */
	    }


	  /* Particle data  */
	  byten(&P[0], NumPart * sizeof(struct particle_data), modus);

	  in(&N_gas, modus);
	  if(N_gas > 0)
	    {
	      if(N_gas > All.MaxPartSph)
		{
		  printf
		    ("SPH: it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		     N_gas / (((double) All.TotN_gas) / NTask));
		  printf("fatal error\n");
		  endrun(222);
		}
	      /* Sph-Particle data  */
	      byten(&SphP[0], N_gas * sizeof(struct sph_particle_data), modus);

#ifdef CHIMES 
	      sphAbundancesBuf = (ChimesFloat *) malloc(N_gas * ChimesGlobalVars.totalNumberOfSpecies * sizeof(ChimesFloat));

	      if (!modus) /* write */
		{
		  /* Read abundance arrays into buffer */
		  for (partIndex = 0; partIndex < N_gas; partIndex++)
		    {
		      for (abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
			sphAbundancesBuf[(partIndex * ChimesGlobalVars.totalNumberOfSpecies) + abunIndex] = ChimesGasVars[partIndex].abundances[abunIndex];
		    }
		}

	      /* Abundance buffer */
	      byten(&sphAbundancesBuf[0], N_gas * ChimesGlobalVars.totalNumberOfSpecies * sizeof(ChimesFloat), modus);
	      /* GasVars */
	      byten(&ChimesGasVars[0], N_gas * sizeof(struct gasVariables), modus);
			  
	      if (modus) /* read */
		{
		  for (partIndex = 0; partIndex < N_gas; partIndex++)
		    {
		      /* Allocate memory for abundance arrays */
		      allocate_gas_abundances_memory(&(ChimesGasVars[partIndex]), &ChimesGlobalVars);
		      
		      /* Read abundances from buffer */
		      for (abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
			ChimesGasVars[partIndex].abundances[abunIndex] = sphAbundancesBuf[(partIndex * ChimesGlobalVars.totalNumberOfSpecies) + abunIndex];

#ifdef CHIMES_TURB_DIFF_IONS 
		      chimes_update_turbulent_abundances(partIndex, 1); 
#endif 
		    }
		}
	      free(sphAbundancesBuf); 
#endif
	    }

	  /* write state of random number generator */
	  byten(gsl_rng_state(random_generator), gsl_rng_size(random_generator), modus);
	  byten(&SelRnd, sizeof(SelRnd), modus);

#ifdef TURB_DRIVING
      byten(gsl_rng_state(StRng), gsl_rng_size(StRng), modus);
	  byten(&StNModes, sizeof(StNModes), modus);
	  byten(StOUPhases, StNModes*6*sizeof(double),modus);
	  byten(StAmpl, StNModes*3*sizeof(double),modus);
	  byten(StAka, StNModes*3*sizeof(double),modus);
	  byten(StAkb, StNModes*3*sizeof(double),modus);
	  byten(StMode, StNModes*3*sizeof(double),modus);
	  byten(&StTPrev, sizeof(StTPrev),modus);
#endif

	  /* write flags for active timebins */
	  byten(TimeBinActive, TIMEBINS * sizeof(int), modus);

	  /* now store relevant data for tree */
        in(&Gas_split, modus);
#ifdef GALSF
        in(&Stars_converted, modus);
#endif


	  /* now store relevant data for tree */

	  in(&nmulti, modus);
        in(&NTopleaves, modus);
        in(&NTopnodes, modus);
	  if(modus != 0 && nmulti != MULTIPLEDOMAINS)
	    {
	      if(ThisTask == 0)
		printf
		  ("Looks like you changed MULTIPLEDOMAINS from %d to %d.\nWe will need to discard tree stored in restart files and construct a new one.\n",
		   nmulti, (int) MULTIPLEDOMAINS);

	      /* In this case we must do a new domain decomposition! */
	    }
	  else
	    {

	      if(modus)		/* read */
		{
		  domain_allocate();
		  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
		}

	      in(&Numnodestree, modus);

	      if(Numnodestree > MaxNodes)
		{
		  printf
		    ("Tree storage: it seems you have reduced(!) 'PartAllocFactor' below the value needed to load the restart file (task=%d). "
		     "Numnodestree=%d  MaxNodes=%d\n", ThisTask, Numnodestree, MaxNodes);
		  endrun(221);
		}

	      byten(Nodes_base, Numnodestree * sizeof(struct NODE), modus);
	      byten(Extnodes_base, Numnodestree * sizeof(struct extNODE), modus);

	      byten(Father, NumPart * sizeof(int), modus);

	      byten(Nextnode, NumPart * sizeof(int), modus);
	      byten(Nextnode + All.MaxPart, NTopnodes * sizeof(int), modus);

	      byten(DomainStartList, NTask * MULTIPLEDOMAINS * sizeof(int), modus);
	      byten(DomainEndList, NTask * MULTIPLEDOMAINS * sizeof(int), modus);
	      byten(TopNodes, NTopnodes * sizeof(struct topnode_data), modus);
	      byten(DomainTask, NTopnodes * sizeof(int), modus);
	      byten(DomainNodeIndex, NTopleaves * sizeof(int), modus);

	      byten(DomainCorner, 3 * sizeof(double), modus);
	      byten(DomainCenter, 3 * sizeof(double), modus);
	      byten(&DomainLen, sizeof(double), modus);
	      byten(&DomainFac, sizeof(double), modus);
	    }

	  fclose(fd);
	}
      else			/* wait inside the group */
	{
	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }


  if(modus != 0 && nmulti != MULTIPLEDOMAINS)	/* in this case we must force a domain decomposition */
    {
        if(ThisTask == 0) {printf("Doing extra domain decomposition because you changed MULTIPLEDOMAINS\n"); fflush(stdout);}

      domain_Decomposition(0, 0, 0);
    }
}



/* reads/writes n bytes 
 */
void byten(void *x, size_t n, int modus)
{
  if(modus)
    my_fread(x, n, 1, fd);
  else
    my_fwrite(x, n, 1, fd);
}


/* reads/writes one int 
 */
void in(int *x, int modus)
{
  if(modus)
    my_fread(x, 1, sizeof(int), fd);
  else
    my_fwrite(x, 1, sizeof(int), fd);
}
