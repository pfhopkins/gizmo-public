/* this is just a 'pre-amble' block which must be included before 
    the code in 'code_block_for_standard_neighbor_loops.h' is called. 
    all it does is define a ton of dummy variable and function names which must be 
    unique, based on the user-specified master routine name (which should of 
    course itself be unique! */

#define INPUT_STRUCT_NAME MASTER_FUNCTION_NAME##_data_in /* dummy name - must be unique */
#define DATAIN_NAME MASTER_FUNCTION_NAME##_dataIn /* dummy name - must be unique */
#define DATAGET_NAME MASTER_FUNCTION_NAME##_DataGet /* dummy name - must be unique */
#define INPUTFUNCTION_NAME MASTER_FUNCTION_NAME##_particle2in /* dummy name - must be unique */
#define OUTPUT_STRUCT_NAME MASTER_FUNCTION_NAME##_data_out /* dummy name - must be unique */
#define DATAOUT_NAME MASTER_FUNCTION_NAME##_DataOut /* dummy name - must be unique */
#define DATARESULT_NAME MASTER_FUNCTION_NAME##_DataResult /* dummy name - must be unique */
#define OUTPUTFUNCTION_NAME MASTER_FUNCTION_NAME##_out2particle /* dummy name - must be unique */
#define EVALUATION_WORKHORSE_FUNCTION_NAME MASTER_FUNCTION_NAME##_evaluate /* dummy name - must be unique */
#define PRIMARY_SUBFUN_NAME MASTER_FUNCTION_NAME##_evaluate_primary /* dummy name - must be unique */
#define SECONDARY_SUBFUN_NAME MASTER_FUNCTION_NAME##_evaluate_secondary /* dummy name - must be unique */
#define CONDITIONFUNCTION_FOR_EVALUATION MASTER_FUNCTION_NAME##_is_active /* dummy name - must be unique */
#define NEIGHBOROPS_FUNCTION_NAME MASTER_FUNCTION_NAME##_neighbor_operations /* dummy name - must be unique */
#define FINAL_OPERATIONS_FUNCTION_NAME MASTER_FUNCTION_NAME##_final_operations /* dummy name - must be unique */
