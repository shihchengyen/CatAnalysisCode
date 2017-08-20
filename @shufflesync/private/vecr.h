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

#ifndef __vecr_h
#define __vecr_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_vecr(void);
extern void TerminateModule_vecr(void);
extern _mexLocalFunctionTable _local_function_table_vecr;

extern mxArray * mlfVecr(mxArray * a);
extern void mlxVecr(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
