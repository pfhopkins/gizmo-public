/* This is a code block which must be inserted at the END of the code_block_xchange...(any).h
     subroutines. This de-allocates everything to prevent memory leaks and errors, and includes
     a couple of the key subroutines (which have been defined already but can be detailed here) */


static inline void *PRIMARY_SUBFUN_NAME(void *p, int loop_iteration)
{
#define CONDITION_FOR_EVALUATION CONDITIONFUNCTION_FOR_EVALUATION
#define EVALUATION_CALL CORE_FUNCTION_NAME(i, 0, exportflag, exportnodecount, exportindex, ngblist, loop_iteration)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}

static inline void *SECONDARY_SUBFUN_NAME(void *p, int loop_iteration)
{
#define EVALUATION_CALL CORE_FUNCTION_NAME(j, 1, &dummy, &dummy, &dummy, ngblist, loop_iteration);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}

#undef CONDITIONFUNCTION_FOR_EVALUATION
#undef SECONDARY_SUBFUN_NAME
#undef PRIMARY_SUBFUN_NAME
#undef OUTPUTFUNCTION_NAME
#undef DATARESULT_NAME
#undef DATAOUT_NAME
#undef OUTPUT_STRUCT_NAME
#undef INPUTFUNCTION_NAME
#undef DATAGET_NAME
#undef DATAIN_NAME
#undef INPUT_STRUCT_NAME
#undef CORE_FUNCTION_NAME
