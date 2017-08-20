/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jun 26 12:13:21 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "shufflesyncsurr" 
 */
#include "vecr.h"
#include "libmatlbm.h"
static mxArray * _mxarray0_;

void InitializeModule_vecr(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
}

void TerminateModule_vecr(void) {
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mvecr(int nargout_, mxArray * a);

_mexLocalFunctionTable _local_function_table_vecr
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfVecr" contains the normal interface for the "vecr"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/vecr.m"
 * (lines 1-15). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfVecr(mxArray * a) {
    int nargout = 1;
    mxArray * vr = NULL;
    mlfEnterNewContext(0, 1, a);
    vr = Mvecr(nargout, a);
    mlfRestorePreviousContext(0, 1, a);
    return mlfReturnValue(vr);
}

/*
 * The function "mlxVecr" contains the feval interface for the "vecr"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/vecr.m"
 * (lines 1-15). The feval function calls the implementation version of vecr
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxVecr(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vecr Line: 1 Column: 1 The function \"vecr\""
            " was called with more than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vecr Line: 1 Column: 1 The function \"vecr"
            "\" was called with more than the declared number of inputs (1)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = Mvecr(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mvecr" is the implementation version of the "vecr" M-function
 * from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/vecr.m"
 * (lines 1-15). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function vr = vecr(a)
 */
static mxArray * Mvecr(int nargout_, mxArray * a) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_vecr);
    mxArray * vr = NULL;
    mxArray * as = NULL;
    mclCopyArray(&a);
    /*
     * %VECR Convert to row vector
     * %   V = VECR(A) checks A to see if it is a N by 1 vector and converts
     * %   it to a 1 by N vector. Otherwise, the function returns A.
     * 
     * % get size of a
     * as = size(a);
     */
    mlfAssign(&as, mlfSize(mclValueVarargout(), mclVa(a, "a"), NULL));
    /*
     * % take transpose only if number of columns is 1 and number of rows is
     * % greater than 1
     * if( (as(1)>1) && (as(2)==1) )
     */
    if (mclScalarToBool(mclGt(mclIntArrayRef1(mclVv(as, "as"), 1), _mxarray0_))
        && mclScalarToBool(
             mclEq(mclIntArrayRef1(mclVv(as, "as"), 2), _mxarray0_))) {
        /*
         * vr = a';
         */
        mlfAssign(&vr, mlfCtranspose(mclVa(a, "a")));
    /*
     * else
     */
    } else {
        /*
         * vr = a;
         */
        mlfAssign(&vr, mclVa(a, "a"));
    /*
     * end
     */
    }
    mclValidateOutput(vr, 1, nargout_, "vr", "vecr");
    mxDestroyArray(as);
    mxDestroyArray(a);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return vr;
}
