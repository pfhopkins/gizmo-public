#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <sys/types.h>

#ifdef NUCLEAR_NETWORK // top-level flag [needs to be here to prevent compiler breaking when this is not active] //

/* options:
    NETWORK_SUPERLU
    NETWORK_PARDISO
    NETWORK_SCREENING
    NETWORK_VARIABLE=1,2 // variable temperature, or temperature+density //
 */

#if defined(NETWORK_SUPERLU) || defined(NETWORK_PARDISO)
#define NETWORK_SPARSE 1
#else
#define NETWORK_SPARSE 0
#endif

#ifdef NETWORK_SCREENING
#define NETWORK_SEP_YZ
#endif

#define NETWORK_VAR_TEMP 1
#define NETWORK_VAR_RHO_TEMP 2
#if defined(NETWORK_VARIABLE) && !(NETWORK_VARIABLE == NETWORK_VAR_TEMP || NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP)
#error "unknown value for NETWORK_VARIABLE"
#endif

#define safe_fgets( str, num, stream ) util_fgets( str, num, stream, __FILE__, __LINE__ )
#define safe_fread( ptr, size, count, stream ) util_fread( ptr, size, count, stream, __FILE__, __LINE__ )

void myprintf( const char* format, ... );
char * util_fgets( char *str, int num, FILE *stream, char *file, int line );
size_t util_fread( void *ptr, size_t size, size_t count, FILE * stream, char *file, int line );

__inline static int imin( int a, int b ) { return a < b ? a : b; }
__inline static int imax( int a, int b ) { return a > b ? a : b; }
__inline static double dmin( double a, double b ) { return a < b ? a : b; }
__inline static double dmax( double a, double b ) { return a > b ? a : b; }

#define VECT_NORM(N, V) ({double VN; int vcc; for(vcc = 0, VN = 0; vcc < (N); vcc++) VN+=((V)[vcc]*(V)[vcc]); VN;})

double SwapDouble( double Val );
float SwapFloat( float Val );
int SwapInt( int Val );
int CheckSwap( char* fname, int *swap );

typedef struct {
  double v; /* the values itself */
  double drho, dT, dYe, dYz; /* the derivatives */
  double dabar, dzbar; /* additional derivatives */
} network_var;

typedef struct {
#if NETWORK_SPARSE
  double *values;
  int *columns, *rowstart;
  int n_elements;
#else
  double *matrix;
#endif
} jacob_t;

struct network_weakrate {
  int input, output; /* input and output species */
  double q1, q2;
  double lambda1[143], lambda2[143]; /* 13 * 11 */
  int isReverse;
};

struct network_rate {
  int type;
  int index; /*!< index in the global rates array */
  int ninput, noutput; /*! number of input and output species (maximum: 4) */
  int input[4], output[4]; /*! indices of the input and output species */
  int isWeak, isReverse, isElectronCapture;
  double q; /*!< q value of the reaction */
  double data[7]; /*!< temperature parametrization of the rate as in REACLIB */
#ifdef NETWORK_SCREENING
  /* some quantities for screening have to be computed only once: z_tilde1 and z_tilde2 are from (A9) in
   * Wallace, Woosley & Weaver 1982, ApJ 258; C_const is - 0.5551 * 5./3. * log(z1*z2/(z1+z2)) - 2.996
   * tau_fac is tau_12 * T_9**(1./3.) from (A2). H_weak_fac is from (A14). Gamma_fac is the prefactor in (A13) */
  double z_tilde1, z_tilde2, C_const, tau_fac, H_w_fac, Gamma_fac;
#endif
};

struct network_solver_trajectory {
    integertime ntimesteps, timestep;
    double *timesteps;
    double *rho;
    double *energy;
    double *x;
    double time, maxtime;
    double mintemp, maxtemp;
};

/* special structures holding sub-structures. need to be defined below the initial structure defs above */
struct network_solver_data {
    int nsteps;
    int maxstep;
    int *steps;
    int maxiter;
    int matrixsize;
    int matrixsize2;
    int nelements;
    double *aion;
    double tolerance;
    double *dy;
    double *ynew;
    double *yscale;
    network_var *rhs_deriv;
    double *first_rhs, *rhs;
    jacob_t jacob;
    jacob_t mod_jacob;
    double *x;
    double *err;
    double *qcol;
    double *a;
    double *alf;
#if NETWORK_VAR_RHO_T
    int iTemp;
#endif
#if NETWORK_VAR_RHO_T == NETWORK_VAR_RHO
    int iRho;
#endif
};

struct network_nucdata {
  int na, nz, nn; /* atomic number, proton number, neutron number */
  char name[6]; /* full name, e.g. he4 */
  double part[24]; /* tabulated partition function */
  double exm;  /* mass excess */
  double q; /* binding energy */
  double spin; /* spin */
  double m; /* atomic mass */
  size_t nrates, nweakrates; /* number of rates that influence the species */
  double *w, *wweak; /* weight factor for the rate */
  /* indices of the rates in the global rates array */
  size_t *irates;
  size_t *iweakrates;
  /* pointers to the rates to allow "central" computation */
  struct network_rate **prates;
  struct network_weakrate **pweakrates;
  /* copies of the rates to allow "local" computation */
  struct network_rate *rates;
  struct network_weakrate *weakrates;
};

struct network_data {
  int initialized;
  size_t nuc_count;
  struct network_nucdata *nucdata;
  size_t rate_count;
  struct network_rate *rates;
  size_t weakrate_count;
  struct network_weakrate *weakrates;
  size_t n_matrix;
  /* indices of the temperature and density equations */
#if NETWORK_VARIABLE == NETWORK_VAR_TEMP || NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  size_t iTemp;
#endif
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  size_t iRho;
#endif
#ifdef NETWORK_SEP_YZ
  size_t iYz;
#endif
  size_t electronCapture_count;
  size_t *electronCapture_rates;
  /* some lists to save if clauses in the loops */
  size_t screened_count;
  struct network_rate **screened_rates;
  struct network_rate *tripleAlpha_rate;
  float weakTemp[13], weakRhoYe[11]; /* temperature and rho*ye grid for weak reactions */
  /* additional rates needed for screening */
  struct network_nucdata nuc_he4, nuc_be8;
  struct network_rate alpha_alpha_rate, Be8_alpha_rate;
};

struct network_workspace {
  int initialized;
  size_t nuc_count;
  network_var *gg; /*!< for the temperature dependent partition functions */
  double *y;
  double *matrix; /*!< used for constructing sparse matrices */
  double old_temp, old_rho, old_ye; /* indicates whether the rates need to be recomputed */
  double *prefac; /* prefactor of x_i = prefac * exp( .. ), used in nse solver */
  network_var *baserate; /*!< the rate without the density factor */
  network_var *unscreened_rate; /*!< baserate before applying screening */
  network_var *rates, *weakrates; /*!< the final rate */
  network_var *yrates, *yweakrates; /*! rates multiplied with the abundances */
  network_var ***prates, ***pweakrates; /*! pointers to those rates (size == nuc_count) */
  network_var temp_deriv_fac, rho_deriv_fac, *dTdYi, *drhodYi; /*!< caches for use in the Jacobian */
  struct network_solver_data *nsd;
};


/* function predefs dependent on definitions above */
 
struct network_solver_data *network_solver_init( double tolerance, int matrixsize, int nelements, const struct network_nucdata *nucdata );
void network_solver_deinit( struct network_solver_data *nsd );
void network_solver_interpolate_trajectory( struct network_solver_trajectory *traj, double time, double *rho, double *energy );
void network_solver_integrate_traj( const struct network_data *nd, struct network_workspace *nw, struct network_solver_trajectory *traj );
void network_solver_integrate( double temp, double rho, double *y, double dt, const struct network_data *nd, struct network_workspace *nw );

int network_init( char *speciesfile, char *ratesfile, char *partfile, char *massesfile, char *weakratesfile, struct network_data *nd );
int network_init_onlyweak( char *speciesfile, char *ratesfile, char *partfile, char *massesfile, char *weakratesfile, struct network_data *nd );
int network_workspace_init( const struct network_data *nd, struct network_workspace *nw );
void network_deinit( struct network_data *nd );
void network_workspace_deinit( struct network_workspace *nw );

int network_getrhs(double rho, double temp, const double y[], int compute_derivs, const struct network_data *nd, struct network_workspace *nw, double *rhs, network_var *deriv);
int network_getjacob(const double y[], const network_var *rhs, const struct network_data *nd, struct network_workspace *nw, jacob_t *jacob);
int network_part(double temp, const struct network_data *nd, struct network_workspace *nw);
int getrates(double rho, double temp, double ye, double yz, int compute_derivs, const struct network_data *nd, struct network_workspace *nw);

void network_normalize(double *x, double *e, const struct network_data *nd, struct network_workspace *nw);
int network_integrate( double temp, double rho, const double *x, double *dx, double dt, double *dedt, double *drhodt, const struct network_data *nd, struct network_workspace *nw );


#endif // top-level flag
