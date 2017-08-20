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

#ifndef __getSurrogateSpikes_h
#define __getSurrogateSpikes_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_getSurrogateSpikes(void);
extern void TerminateModule_getSurrogateSpikes(void);
extern _mexLocalFunctionTable _local_function_table_getSurrogateSpikes;

extern mxArray * mlfGetSurrogateSpikes(mxArray * rfd, mxArray * n, ...);
extern void mlxGetSurrogateSpikes(int nlhs,
                                  mxArray * plhs[],
                                  int nrhs,
                                  mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
