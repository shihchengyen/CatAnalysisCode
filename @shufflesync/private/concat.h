/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __concat_h
#define __concat_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_concat(void);
extern void TerminateModule_concat(void);
extern _mexLocalFunctionTable _local_function_table_concat;

extern mxArray * mlfConcat(mxArray * a, ...);
extern void mlxConcat(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
