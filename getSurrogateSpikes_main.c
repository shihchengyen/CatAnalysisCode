/*
 * MATLAB Compiler: 3.0
 * Date: Wed Jan 21 11:28:01 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "getSurrogateSpikes.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "libmatlb.h"
#include "getSurrogateSpikes.h"
#include "getOptArgs.h"
#include "libmmfile.h"

extern _mex_information _main_info;

static mexFunctionTableEntry function_table[2]
  = { { "getSurrogateSpikes", mlxGetSurrogateSpikes, -3,
        1, &_local_function_table_getSurrogateSpikes },
      { "getOptArgs", mlxGetOptArgs, -3, 1,
        &_local_function_table_getOptArgs } };

static _mexInitTermTableEntry init_term_table[3]
  = { { libmmfileInitialize, libmmfileTerminate },
      { InitializeModule_getSurrogateSpikes,
        TerminateModule_getSurrogateSpikes },
      { InitializeModule_getOptArgs, TerminateModule_getOptArgs } };

_mex_information _main_info
  = { 1, 2, function_table, 0, NULL, 0, NULL, 3, init_term_table };

/*
 * The function "main" is a Compiler-generated main wrapper, suitable for
 * building a stand-alone application.  It calls a library function to perform
 * initialization, call the main function, and perform library termination.
 */
int main(int argc, const char * * argv) {
    return mclMain(argc, argv, mlxGetSurrogateSpikes, 1, &_main_info);
}
