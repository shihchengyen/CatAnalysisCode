/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "libmatlb.h"
#include "_shufflesync_private_computeSurrData.h"
#include "concat.h"
#include "concatenate.h"
#include "convmtx.h"
#include "histcie.h"
#include "jbtest.h"
#include "kstest2.h"
#include "lillietest.h"
#include "nptDir.h"
#include "getOptArgs.h"
#include "vecc.h"
#include "chi2cdf.h"
#include "chi2inv.h"
#include "cdfcalc.h"
#include "normcdf.h"
#include "nptFileParts.h"
#include "removeargs.h"
#include "distchck.h"
#include "gamcdf.h"
#include "gaminv.h"
#include "norminv.h"
#include "gampdf.h"
#include "libmmfile.h"

extern _mex_information _main_info;

static mexFunctionTableEntry function_table[22]
  = { { "@shufflesync/private/computeSurrData",
        mlx_shufflesync_private_computeSurrData, 1, 5,
        &_local_function_table__shufflesync_private_computeSurrData },
      { "concat", mlxConcat, -2, 1, &_local_function_table_concat },
      { "concatenate", mlxConcatenate, -3, 1,
        &_local_function_table_concatenate },
      { "convmtx", mlxConvmtx, 2, 1, &_local_function_table_convmtx },
      { "histcie", mlxHistcie, -3, 2, &_local_function_table_histcie },
      { "jbtest", mlxJbtest, 2, 4, &_local_function_table_jbtest },
      { "kstest2", mlxKstest2, 4, 3, &_local_function_table_kstest2 },
      { "lillietest", mlxLillietest, 2, 4, &_local_function_table_lillietest },
      { "nptDir", mlxNptDir, -1, 1, &_local_function_table_nptDir },
      { "getOptArgs", mlxGetOptArgs, -3, 2,
        &_local_function_table_getOptArgs },
      { "vecc", mlxVecc, 1, 1, &_local_function_table_vecc },
      { "chi2cdf", mlxChi2cdf, 2, 1, &_local_function_table_chi2cdf },
      { "chi2inv", mlxChi2inv, 2, 1, &_local_function_table_chi2inv },
      { "cdfcalc", mlxCdfcalc, 2, 4, &_local_function_table_cdfcalc },
      { "normcdf", mlxNormcdf, 5, 3, &_local_function_table_normcdf },
      { "nptFileParts", mlxNptFileParts, 1, 4,
        &_local_function_table_nptFileParts },
      { "removeargs", mlxRemoveargs, 3, 2, &_local_function_table_removeargs },
      { "distchck", mlxDistchck, 5, 5, &_local_function_table_distchck },
      { "gamcdf", mlxGamcdf, 3, 1, &_local_function_table_gamcdf },
      { "gaminv", mlxGaminv, 3, 1, &_local_function_table_gaminv },
      { "norminv", mlxNorminv, 5, 3, &_local_function_table_norminv },
      { "gampdf", mlxGampdf, 3, 1, &_local_function_table_gampdf } };

static _mexInitTermTableEntry init_term_table[23]
  = { { libmmfileInitialize, libmmfileTerminate },
      { InitializeModule__shufflesync_private_computeSurrData,
        TerminateModule__shufflesync_private_computeSurrData },
      { InitializeModule_concat, TerminateModule_concat },
      { InitializeModule_concatenate, TerminateModule_concatenate },
      { InitializeModule_convmtx, TerminateModule_convmtx },
      { InitializeModule_histcie, TerminateModule_histcie },
      { InitializeModule_jbtest, TerminateModule_jbtest },
      { InitializeModule_kstest2, TerminateModule_kstest2 },
      { InitializeModule_lillietest, TerminateModule_lillietest },
      { InitializeModule_nptDir, TerminateModule_nptDir },
      { InitializeModule_getOptArgs, TerminateModule_getOptArgs },
      { InitializeModule_vecc, TerminateModule_vecc },
      { InitializeModule_chi2cdf, TerminateModule_chi2cdf },
      { InitializeModule_chi2inv, TerminateModule_chi2inv },
      { InitializeModule_cdfcalc, TerminateModule_cdfcalc },
      { InitializeModule_normcdf, TerminateModule_normcdf },
      { InitializeModule_nptFileParts, TerminateModule_nptFileParts },
      { InitializeModule_removeargs, TerminateModule_removeargs },
      { InitializeModule_distchck, TerminateModule_distchck },
      { InitializeModule_gamcdf, TerminateModule_gamcdf },
      { InitializeModule_gaminv, TerminateModule_gaminv },
      { InitializeModule_norminv, TerminateModule_norminv },
      { InitializeModule_gampdf, TerminateModule_gampdf } };

_mex_information _main_info
  = { 1, 22, function_table, 0, NULL, 0, NULL, 23, init_term_table };

/*
 * The function "main" is a Compiler-generated main wrapper, suitable for
 * building a stand-alone application.  It calls a library function to perform
 * initialization, call the main function, and perform library termination.
 */
int main(int argc, const char * * argv) {
    return mclMain(
             argc,
             argv,
             mlx_shufflesync_private_computeSurrData,
             1,
             &_main_info);
}
