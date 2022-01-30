#ifndef EOS_INTERFACE_H
#define EOS_INTERFACE_H

#include "../GIZMO_config.h"


#if (defined(EOS_TILLOTSON) || defined(EOS_ELASTIC) || defined(EOS_HELMHOLTZ) || defined(COSMIC_RAY_FLUID) || defined(RT_RADPRESSURE_IN_HYDRO) || defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE) || defined(EOS_GMC_BAROTROPIC) || defined(FLAG_NOT_IN_PUBLIC_CODE)) && !defined(EOS_GENERAL)
#define EOS_GENERAL
#endif

#ifdef EOS_HELMHOLTZ
#define EOS_TABULATED
#define EOS_USES_CGS
#define EOS_CARRIES_YE
#define EOS_CARRIES_ABAR
#define EOS_CARRIES_TEMPERATURE
#define EOS_PROVIDES_ENTROPY
#define EOS_PROVIDES_CV
#endif

struct eos_input
{
  double rho;         /* Density */
  double eps;         /* Specific internal energy */
#ifdef EOS_CARRIES_YE
  double Ye;          /* Electron fraction */
#endif
#ifdef EOS_CARRIES_ABAR
  double Abar;        /* Mean atomic weight (in atomic mass units) */
#endif
#ifdef EOS_CARRIES_TEMPERATURE
  double temp;        /* Temperature initial guess */
#endif
};

struct eos_output
{
  double press;       /* Pressure */
  double csound;      /* Sound speed */
#ifdef EOS_CARRIES_TEMPERATURE
  double temp;        /* Temperature (in Kelvin) */
#endif
#ifdef EOS_PROVIDES_ENTROPY
  double entropy;     /* Entropy (in CGS) */
#endif
#ifdef EOS_PROVIDES_CV
  double cv;          /* Specific heat at constant volume (in CGS) */
#endif
};

#ifdef EOS_TABULATED
int eos_init(char const * eos_table_fname);
int eos_cleanup();
#endif

int eos_compute(struct eos_input const * in, struct eos_output * out);

#endif
