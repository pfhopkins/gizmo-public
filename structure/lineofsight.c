#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! compute line-of-sight integrated quantities (for e.g. Lyman-alpha forest studies) */

/*
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * It has been updated by PFH for basic compatibility with GIZMO.
 */

#define  LYMAN_ALPHA_WAVELENGTH_CM_HI    (1215.6e-8)    /* 1215.6 Angstrom */
#define  LYMAN_ALPHA_WAVELENGTH_CM_HeII  (303.8e-8)    /* 303.8 Angstrom */
#define  OSCILLATOR_STRENGTH_HI    (0.41615)
#define  OSCILLATOR_STRENGTH_HeII  (0.41615)

#ifdef OUTPUT_LINEOFSIGHT

#ifdef OUTPUT_LINEOFSIGHT_SPECTRUM
#define  PIXELS 512		/* number of bins along line of sight */
#else
#define  PIXELS 1
#endif

#define  N_LOS  10		/* number of lines of sight selected  */

static double H_a, Wmax;

struct line_of_sight
{
  int xaxis, yaxis, zaxis;
  double Xpos, Ypos;		/* relative position of line-of-sight on face of box */
  double BoxSize, Wmax, Time;

  /* total gas density */
  double Rho[PIXELS];
  double Vpec[PIXELS];
  double Temp[PIXELS];
  double Metallicity[PIXELS];

  /* neutral hydrogen */
  double RhoHI[PIXELS];
  double NHI[PIXELS];
  double VpecHI[PIXELS];
  double TempHI[PIXELS];
  double TauHI[PIXELS];

  /* HeII quantities */
  double RhoHeII[PIXELS];
  double NHeII[PIXELS];
  double VpecHeII[PIXELS];
  double TempHeII[PIXELS];
  double TauHeII[PIXELS];
}
 *Los, *LosGlobal;


struct line_of_sight_particles
{
  MyFloat Pos[3];
  MyFloat Hsml;
  MyFloat Vz;
  MyFloat Utherm;
  MyFloat Mass;
  MyFloat Metallicity;
}
 *particles;



void lineofsight_output(void)
{
  char buf[500];
  int n, s, next;
  double ti;

  next = find_next_lineofsighttime(All.Ti_nextlineofsight);

  ti = All.TimeBegin * exp(next * All.Timebase_interval);
  PRINT_STATUS("Line of sight output! ThisTask=%d Time=%g  NextTime=%g", ThisTask, All.Time, ti);
  H_a = hubble_function(All.Time);
  Wmax = All.Time * H_a * All.BoxSize;


  if(ThisTask == 0)
    {
      sprintf(buf, "%s/los", All.OutputDir);
      mkdir(buf, 02755);
    }

  Los = mymalloc("Los", sizeof(struct line_of_sight));
  LosGlobal = mymalloc("LosGlobal", sizeof(struct line_of_sight));

  for(n = 0, s = 0; n < N_LOS; n++)
    {
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
      if(s + 3 >= RNDTABLE)
	{
	  set_random_numbers();
	  s = 0;
	}
#endif

      Los->zaxis = (int) (3.0 * get_random_number(s++));
      switch (Los->zaxis)
	{
	case 2:
	  Los->xaxis = 0;
	  Los->yaxis = 1;
	  break;
	case 0:
	  Los->xaxis = 1;
	  Los->yaxis = 2;
	  break;
	case 1:
	  Los->xaxis = 2;
	  Los->yaxis = 0;
	  break;
	}

      Los->Xpos = All.BoxSize * get_random_number(s++);
      Los->Ypos = All.BoxSize * get_random_number(s++);

#ifdef OUTPUT_LINEOFSIGHT_SPECTRUM
      add_along_lines_of_sight();
      sum_over_processors_and_normalize();
      absorb_along_lines_of_sight();
      output_lines_of_sight(n);
#endif

#ifdef OUTPUT_LINEOFSIGHT_PARTICLES
      find_particles_and_save_them(n);
#endif
    }

  myfree(LosGlobal);
  myfree(Los);
}


void find_particles_and_save_them(int num)
{
  int n, k, count_local, *countlist, counttot, rep;
  MyFloat dx, dy, dz, r2;
  char fname[1000];
  MPI_Status status;
  FILE *fd = 0;

  countlist = mymalloc("countlist", sizeof(int) * NTask);
  particles = mymalloc("particles", sizeof(struct line_of_sight_particles) * N_gas);

  for(n = 0, count_local = 0; n < N_gas; n++)
    {
      if(P[n].Type == 0)
	{
        dx = P[n].Pos[Los->xaxis] - Los->Xpos;
        dy = P[n].Pos[Los->yaxis] - Los->Ypos;
        dz = 0;
        NEAREST_XYZ(dx,dy,dz,-1);

	  r2 = dx * dx + dy * dy;

	  if(r2 < PPP[n].Hsml * PPP[n].Hsml)
	    {
	      for(k = 0; k < 3; k++)
		particles[count_local].Pos[k] = P[n].Pos[k];

	      particles[count_local].Hsml = PPP[n].Hsml;
	      particles[count_local].Vz = P[n].Vel[Los->zaxis];
	      particles[count_local].Utherm = SphP[n].InternalEnergyPred;
	      particles[count_local].Mass = P[n].Mass;
	      particles[count_local].Metallicity = P[n].Metallicity[0];

	      count_local++;
	    }
	}
    }

  MPI_Gather(&count_local, 1, MPI_INT, countlist, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      sprintf(fname, "%s/los/part_los_z%05.3f_%03d.dat", All.OutputDir, 1 / All.Time - 1, num);

      if(!(fd = fopen(fname, "w")))
	{
	  printf("can't open file `%s`\n", fname);
	  endrun(112);
	}

      fwrite(&count_local, sizeof(int), 1, fd);	/* will be overwritten later */
      fwrite(&LosGlobal->xaxis, sizeof(int), 1, fd);
      fwrite(&LosGlobal->yaxis, sizeof(int), 1, fd);
      fwrite(&LosGlobal->zaxis, sizeof(int), 1, fd);
      fwrite(&LosGlobal->Xpos, sizeof(double), 1, fd);
      fwrite(&LosGlobal->Ypos, sizeof(double), 1, fd);
      fwrite(&LosGlobal->BoxSize, sizeof(double), 1, fd);
      fwrite(&LosGlobal->Wmax, sizeof(double), 1, fd);
      fwrite(&LosGlobal->Time, sizeof(double), 1, fd);
    }


  for(rep = 0, counttot = 0; rep < NTask; rep++)
    {
      if(ThisTask != 0 && rep == ThisTask && count_local > 0)
	MPI_Ssend(particles, sizeof(struct line_of_sight_particles) * count_local, MPI_BYTE, 0,
		  TAG_PDATA, MPI_COMM_WORLD);

      if(ThisTask == 0)
	{
	  if(rep > 0 && countlist[rep] > 0)
	    MPI_Recv(particles, sizeof(struct line_of_sight_particles) * countlist[rep],
		     MPI_BYTE, rep, TAG_PDATA, MPI_COMM_WORLD, &status);

	  fwrite(particles, sizeof(struct line_of_sight_particles), countlist[rep], fd);

	  counttot += countlist[rep];
	}
    }

  if(ThisTask == 0)
    {
      fclose(fd);
      if(!(fd = fopen(fname, "r+")))
	{
	  printf("can't open file `%s'\n", fname);
	  endrun(113);
	}

      fseek(fd, 0, SEEK_CUR);
      fwrite(&counttot, sizeof(int), 1, fd);
      fclose(fd);
    }

  myfree(particles);
  myfree(countlist);
}






void add_along_lines_of_sight(void)
{
  int n, bin, i, iz0, iz1, iz;
  double dx, dy, dz, r, r2, ne=1, nh0=0, nHeII=0, utherm=1, temp=1e4;
  double u, wk, dwk, weight, h3inv;
  double z0, z1;

  for(i = 0; i < PIXELS; i++)
    {
      Los->Rho[i] = 0;
      Los->Vpec[i] = 0;
      Los->Temp[i] = 0;
      Los->Metallicity[i] = 0;

      Los->RhoHI[i] = 0;
      Los->NHI[i] = 0;
      Los->VpecHI[i] = 0;
      Los->TempHI[i] = 0;
      Los->TauHI[i] = 0;

      Los->RhoHeII[i] = 0;
      Los->NHeII[i] = 0;
      Los->VpecHeII[i] = 0;
      Los->TempHeII[i] = 0;
      Los->TauHeII[i] = 0;
    }


  for(n = 0; n < N_gas; n++)
    {
      if(P[n].Type == 0)
	{
        dx = P[n].Pos[Los->xaxis] - Los->Xpos;
        dy = P[n].Pos[Los->yaxis] - Los->Ypos;
        dz = 0;
        NEAREST_XYZ(dx,dy,dz,-1);

	  r2 = dx * dx + dy * dy;

	  if(r2 < PPP[n].Hsml * PPP[n].Hsml)
	    {
	      z0 = (P[n].Pos[Los->zaxis] - PPP[n].Hsml) / All.BoxSize * PIXELS;
	      z1 = (P[n].Pos[Los->zaxis] + PPP[n].Hsml) / All.BoxSize * PIXELS;
	      iz0 = (int) z0;
	      iz1 = (int) z1;
	      if(z0 < 0)
		iz0 -= 1;

	      for(iz = iz0; iz <= iz1; iz++)
		{
            dx = P[n].Pos[Los->xaxis] - Los->Xpos;
            dy = P[n].Pos[Los->yaxis] - Los->Ypos;
            dz = (iz + 0.5) / PIXELS * All.BoxSize - P[n].Pos[Los->zaxis];
            NEAREST_XYZ(dx,dy,dz,-1);
		  r = sqrt(r2 + dz * dz);

		  if(PPP[n].Hsml > All.BoxSize)
		    {
		      printf("Here:%d  n=%d %g\n", ThisTask, n, PPP[n].Hsml);
		      endrun(89);
		    }

		  if(r < PPP[n].Hsml)
		    {
		      u = r / PPP[n].Hsml;
		      h3inv = 1.0 / (PPP[n].Hsml * PPP[n].Hsml * PPP[n].Hsml);
              kernel_main(u, h3inv, 1, &wk, &dwk, -1);
                
		      bin = iz;
		      while(bin >= PIXELS)
                  bin -= PIXELS;
		      while(bin < 0)
                  bin += PIXELS;

              utherm = DMAX(All.MinEgySpec, SphP[i].InternalEnergyPred);
              double mu_in = 1, nHe0, nHepp, nhp;
              temp = ThermalProperties(utherm, SphP[n].Density * All.cf_a3inv, n, &mu_in, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);

		      /* do total gas */
		      weight = P[n].Mass * wk;
		      Los->Rho[bin] += weight;
		      Los->Metallicity[bin] += P[n].Metallicity[0] * weight;
		      Los->Temp[bin] += temp * weight;
		      Los->Vpec[bin] += P[n].Vel[Los->zaxis] * weight;

		      /* do neutral hydrogen */
		      weight = nh0 * HYDROGEN_MASSFRAC * P[n].Mass * wk;
		      Los->RhoHI[bin] += weight;
		      Los->TempHI[bin] += temp * weight;
		      Los->VpecHI[bin] += P[n].Vel[Los->zaxis] * weight;

		      /* do HeII */
		      weight = 4 * nHeII * HYDROGEN_MASSFRAC * P[n].Mass * wk;
		      Los->RhoHeII[bin] += weight;
		      Los->TempHeII[bin] += temp * weight;
		      Los->VpecHeII[bin] += P[n].Vel[Los->zaxis] * weight;
		    }
		}
	    }
	}
    }
}


void sum_over_processors_and_normalize(void)
{
  int bin;

  MPI_Reduce(Los->Rho, LosGlobal->Rho, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->Metallicity, LosGlobal->Metallicity, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->Temp, LosGlobal->Temp, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->Vpec, LosGlobal->Vpec, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(Los->RhoHI, LosGlobal->RhoHI, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->TempHI, LosGlobal->TempHI, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->VpecHI, LosGlobal->VpecHI, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(Los->RhoHeII, LosGlobal->RhoHeII, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->TempHeII, LosGlobal->TempHeII, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Los->VpecHeII, LosGlobal->VpecHeII, PIXELS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      /* normalize results by the weights */
      for(bin = 0; bin < PIXELS; bin++)
	{
	  /* total gas density */
	  LosGlobal->Metallicity[bin] /= LosGlobal->Rho[bin];
	  LosGlobal->Temp[bin] /= LosGlobal->Rho[bin];
	  LosGlobal->Vpec[bin] /= (All.Time * LosGlobal->Rho[bin]);

	  /* neutral hydrogen quantities */
	  LosGlobal->VpecHI[bin] /= LosGlobal->RhoHI[bin];
	  LosGlobal->TempHI[bin] /= LosGlobal->RhoHI[bin];
	  LosGlobal->NHI[bin] = LosGlobal->RhoHI[bin] * (UNIT_MASS_IN_CGS / PROTONMASS_CGS);

	  /* HeII quantities */
	  LosGlobal->VpecHeII[bin] /= (All.Time * LosGlobal->RhoHeII[bin]);
	  LosGlobal->TempHeII[bin] /= LosGlobal->RhoHeII[bin];
	  LosGlobal->NHeII[bin] = LosGlobal->RhoHeII[bin] * (UNIT_MASS_IN_CGS / (4 * PROTONMASS_CGS));
	}
    }
}



void absorb_along_lines_of_sight(void)
{
  double dz, dv, b, fac, fac_HeII;
  int bin, k;


  if(ThisTask == 0)
    {
      dz = All.BoxSize / PIXELS;

      for(bin = 0; bin < PIXELS; bin++)
	{
	  LosGlobal->TauHI[bin] = 0;
	  LosGlobal->TauHeII[bin] = 0;

	  for(k = 0; k < PIXELS; k++)
	    {
	      dv = (k - bin);

	      while(dv < -PIXELS / 2)
		dv += PIXELS;
	      while(dv > PIXELS / 2)
		dv -= PIXELS;

	      dv = (dv * Wmax / PIXELS + LosGlobal->VpecHI[k]) * UNIT_VEL_IN_CGS;

	      b = sqrt(2 * BOLTZMANN_CGS * LosGlobal->TempHI[k] / PROTONMASS_CGS);

	      LosGlobal->TauHI[bin] += LosGlobal->NHI[k] * exp(-dv * dv / (b * b)) / b * dz;


	      /* now HeII */
	      dv = (k - bin);

	      while(dv < -PIXELS / 2)
		dv += PIXELS;
	      while(dv > PIXELS / 2)
		dv -= PIXELS;

	      dv = (dv * Wmax / PIXELS + LosGlobal->VpecHeII[k]) * UNIT_VEL_IN_CGS;

	      b = sqrt(2 * BOLTZMANN_CGS * LosGlobal->TempHeII[k] / (4 * PROTONMASS_CGS));

	      LosGlobal->TauHeII[bin] += LosGlobal->NHeII[k] * exp(-dv * dv / (b * b)) / b * dz;
	    }
	}


      /* multiply with correct prefactors */

      /*  to get things into cgs units */
      fac = 1 / (UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS);
      fac *= OSCILLATOR_STRENGTH_HI * M_PI * LYMAN_ALPHA_WAVELENGTH_CM_HI * sqrt(3 * (THOMPSON_CX_CGS) / (8 * M_PI)) * C_LIGHT_CGS / (All.cf_atime * All.cf_atime) / sqrt(M_PI);	/* Ly-alpha cross section, Thompson cross-section in cgs */

      /* Note: For HeII, the oscillator strength is equal to that of HI,
         and the Lyman-alpha wavelength is 4 times shorter */

      fac_HeII = fac * (OSCILLATOR_STRENGTH_HeII / OSCILLATOR_STRENGTH_HI) * (LYMAN_ALPHA_WAVELENGTH_CM_HeII / LYMAN_ALPHA_WAVELENGTH_CM_HI);

      for(bin = 0; bin < PIXELS; bin++)
	{
	  LosGlobal->TauHI[bin] *= fac;
	  LosGlobal->TauHeII[bin] *= fac_HeII;
	}

      LosGlobal->BoxSize = All.BoxSize;
      LosGlobal->Wmax = Wmax;
      LosGlobal->Time = All.Time;
    }

}



void output_lines_of_sight(int num)
{
  FILE *fd;
  int dummy;
  char fname[400];

  if(ThisTask != 0)
    return;

  sprintf(fname, "%s/los/spec_los_z%05.3f_%03d.dat", All.OutputDir, 1 / All.Time - 1, num);

  if(!(fd = fopen(fname, "w")))
    {
      printf("can't open file `%s`\n", fname);
      endrun(1);
    }

  dummy = PIXELS;
  fwrite(&dummy, sizeof(int), 1, fd);
  fwrite(&Los->BoxSize, sizeof(double), 1, fd);
  fwrite(&Los->Wmax, sizeof(double), 1, fd);
  fwrite(&Los->Time, sizeof(double), 1, fd);
  fwrite(&Los->Xpos, sizeof(double), 1, fd);
  fwrite(&Los->Ypos, sizeof(double), 1, fd);
  fwrite(&Los->xaxis, sizeof(int), 1, fd);
  fwrite(&Los->yaxis, sizeof(int), 1, fd);
  fwrite(&Los->zaxis, sizeof(int), 1, fd);

  fwrite(LosGlobal->TauHI, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->TempHI, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->VpecHI, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->NHI, sizeof(double), PIXELS, fd);

  fwrite(LosGlobal->TauHeII, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->TempHeII, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->VpecHeII, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->NHeII, sizeof(double), PIXELS, fd);

  fwrite(LosGlobal->Rho, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->Vpec, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->Temp, sizeof(double), PIXELS, fd);
  fwrite(LosGlobal->Metallicity, sizeof(double), PIXELS, fd);

  fclose(fd);
}


integertime find_next_lineofsighttime(integertime time0)
{
  double a1, a2, df1, df2, u1, u2;
  double logTimeBegin, logTimeMax;
  int i1, i2, im;
  integertime time1;

  logTimeBegin = log(All.TimeBegin);
  logTimeMax = log(All.TimeMax);

  a1 = logTimeBegin + time0 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH) {i1 = DRIFT_TABLE_LENGTH - 1;}

  if(i1 <= 1) {df1 = u1 * GravKickTable[0];} else {df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);}

  df2 = df1 + All.BoxSize / C_LIGHT_CODE;

  i2 = DRIFT_TABLE_LENGTH - 1;

  while(i2 - i1 > 0)
    {
      im = (i1 - 1 + i2) / 2;

      if(GravKickTable[im] > df2)
	i2 = im;
      else
	i1 = im + 1;
    }

  u2 = (df2 - GravKickTable[i2 - 1]) / (GravKickTable[i2] - GravKickTable[i2 - 1]) + i2;

  a2 = u2 * (logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH + logTimeBegin;

  time1 = (integertime) ((a2 - logTimeBegin) / All.Timebase_interval + 0.5);

  return time1;
}


#endif
