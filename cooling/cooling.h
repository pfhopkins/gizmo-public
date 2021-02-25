#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/*
 * This file contains the definitions for the cooling.c routines
 *
 * This file was originally part of the GADGET3 code developed by
 *   Volker Springel. The code has been modified by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

double ThermalProperties(double u, double rho, int target, double *mu_guess, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess);
double return_uvb_shieldfac(int target, double gamma_12, double nHcgs, double logT);
double return_local_gammamultiplier(int target);
double evaluate_Compton_heating_cooling_rate(int target, double T, double nHcgs, double n_elec, double shielding_factor_for_exgalbg);
void   InitCool(void);
#ifndef CHIMES 
void   InitCoolMemory(void);
void   IonizeParams(void);
void   IonizeParamsFunction(void);
void   IonizeParamsTable(void);
//double INLINE_FUNC LogTemp(double u, double ne);
void   MakeCoolingTable(void);
void   ReadIonizeParams(char *fname);
void   SetZeroIonization(void);
#endif 
void   TestCool(void);

#ifndef CHIMES 
double find_abundances_and_rates(double logT, double rho, int target, double shieldfac, int return_cooling_mode,
                                 double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess);
double convert_u_to_temp(double u, double rho, int target, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess);
double CoolingRate(double logT, double rho, double nelec, int target);
double CoolingRateFromU(double u, double rho, double ne_guess, int target);
#endif 
double DoCooling(double u_old, double rho, double dt, double ne_guess, int target);
#ifndef CHIMES 
double GetCoolingTime(double u_old, double rho,  double ne_guess, int target);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double ne_guess, int target);
#endif

#ifdef COOL_GRACKLE
void InitGrackle(void);
double CallGrackle(double u_old, double rho, double dt, double ne_guess, int target, int mode);
#endif

