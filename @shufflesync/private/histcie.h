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

#ifndef __histcie_h
#define __histcie_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_histcie(void);
extern void TerminateModule_histcie(void);
extern _mexLocalFunctionTable _local_function_table_histcie;

extern mxArray * mlfNHistcie(int nargout,
                             mxArray * * bin,
                             mxArray * x,
                             mxArray * edges,
                             ...);
extern mxArray * mlfHistcie(mxArray * * bin, mxArray * x, mxArray * edges, ...);
extern void mlfVHistcie(mxArray * x, mxArray * edges, ...);
extern void mlxHistcie(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
