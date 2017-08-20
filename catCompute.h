/*
 * MATLAB Compiler: 3.0
 * Date: Thu Apr 22 15:58:01 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "catCompute.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __catCompute_h
#define __catCompute_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_catCompute(void);
extern void TerminateModule_catCompute(void);
extern _mexLocalFunctionTable _local_function_table_catCompute;

extern void mlfCatCompute(mxArray * sets, mxArray * filename);
extern void mlxCatCompute(int nlhs,
                          mxArray * plhs[],
                          int nrhs,
                          mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
