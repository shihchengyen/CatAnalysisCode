/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jun 26 12:13:21 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "shufflesyncsurr" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "libmatlb.h"
#include "_shufflesync_private_shufflesyncsurr.h"
#include "concatenate.h"
#include "histcie.h"
#include "vecc.h"
#include "vecr.h"
#include "getOptArgs.h"
#include "removeargs.h"
#include "libmmfile.h"

extern _mex_information _main_info;

static mexFunctionTableEntry function_table[7]
  = { { "@shufflesync/private/shufflesyncsurr",
        mlx_shufflesync_private_shufflesyncsurr, -2, 0,
        &_local_function_table__shufflesync_private_shufflesyncsurr },
      { "concatenate", mlxConcatenate, -3, 1,
        &_local_function_table_concatenate },
      { "histcie", mlxHistcie, -3, 2, &_local_function_table_histcie },
      { "vecc", mlxVecc, 1, 1, &_local_function_table_vecc },
      { "vecr", mlxVecr, 1, 1, &_local_function_table_vecr },
      { "getOptArgs", mlxGetOptArgs, -3, 2,
        &_local_function_table_getOptArgs },
      { "removeargs", mlxRemoveargs, 3, 2,
        &_local_function_table_removeargs } };

static _mexInitTermTableEntry init_term_table[8]
  = { { libmmfileInitialize, libmmfileTerminate },
      { InitializeModule__shufflesync_private_shufflesyncsurr,
        TerminateModule__shufflesync_private_shufflesyncsurr },
      { InitializeModule_concatenate, TerminateModule_concatenate },
      { InitializeModule_histcie, TerminateModule_histcie },
      { InitializeModule_vecc, TerminateModule_vecc },
      { InitializeModule_vecr, TerminateModule_vecr },
      { InitializeModule_getOptArgs, TerminateModule_getOptArgs },
      { InitializeModule_removeargs, TerminateModule_removeargs } };

_mex_information _main_info
  = { 1, 7, function_table, 0, NULL, 0, NULL, 8, init_term_table };

/*
 * The function "main" is a Compiler-generated main wrapper, suitable for
 * building a stand-alone application.  It calls a library function to perform
 * initialization, call the main function, and perform library termination.
 */
int main(int argc, const char * * argv) {
    return mclMain(
             argc,
             argv,
             mlx_shufflesync_private_shufflesyncsurr,
             0,
             &_main_info);
}
