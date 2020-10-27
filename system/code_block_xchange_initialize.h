/* this is a 'pre-amble' block which must be included before
    the code in 'code_block_xchange...(any of the other routines).h' is called.
    It sets a number of variable and function names based on the user-specified
    parent routine name (which should of course itself be unique!). It also sets
    up some definitions for the multi-threading, and pre-defines some of the
    subroutines which will be referenced. */

/* initialize some thread information */
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

/* initialize macro and variable names: these define structures/variables with names following the value of CORE_FUNCTION_NAME with the '_data_in' and other terms appended -- this should be unique within the file defined! */
#define INPUT_STRUCT_NAME   MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _data_in_)
#define DATAIN_NAME         MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _DataIn_)
#define DATAGET_NAME        MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _DataGet_)
#define OUTPUT_STRUCT_NAME  MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _data_out_)
#define DATAOUT_NAME        MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _DataOut_)
#define DATARESULT_NAME     MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _DataResult_)
#define PRIMARY_SUBFUN_NAME MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _subfun_primary_)
#define SECONDARY_SUBFUN_NAME MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _subfun_secondary_)
#ifndef INPUTFUNCTION_NAME /* assign a default name if one is not user-specified */
#define INPUTFUNCTION_NAME MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _particle2in_)
#endif
#ifndef OUTPUTFUNCTION_NAME /* assign a default name if one is not user-specified */
#define OUTPUTFUNCTION_NAME MACRO_NAME_CONCATENATE(CORE_FUNCTION_NAME, _out2particle_)
#endif


/* define generic forms of functions used below */
int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration);
//static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration);
//static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration);
static inline void *PRIMARY_SUBFUN_NAME(void *p, int loop_iteration);
static inline void *SECONDARY_SUBFUN_NAME(void *p, int loop_iteration);

