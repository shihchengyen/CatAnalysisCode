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

#ifndef __chi2cdf_h
#define __chi2cdf_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_chi2cdf(void);
extern void TerminateModule_chi2cdf(void);
extern _mexLocalFunctionTable _local_function_table_chi2cdf;

extern mxArray * mlfChi2cdf(mxArray * x, mxArray * v);
extern void mlxChi2cdf(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
