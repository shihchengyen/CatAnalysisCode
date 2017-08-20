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

#ifndef __normcdf_h
#define __normcdf_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_normcdf(void);
extern void TerminateModule_normcdf(void);
extern _mexLocalFunctionTable _local_function_table_normcdf;

extern mxArray * mlfNNormcdf(int nargout,
                             mxArray * * plo,
                             mxArray * * pup,
                             mxArray * x,
                             mxArray * mu,
                             mxArray * sigma,
                             mxArray * pcov,
                             mxArray * alpha);
extern mxArray * mlfNormcdf(mxArray * * plo,
                            mxArray * * pup,
                            mxArray * x,
                            mxArray * mu,
                            mxArray * sigma,
                            mxArray * pcov,
                            mxArray * alpha);
extern void mlfVNormcdf(mxArray * x,
                        mxArray * mu,
                        mxArray * sigma,
                        mxArray * pcov,
                        mxArray * alpha);
extern void mlxNormcdf(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
