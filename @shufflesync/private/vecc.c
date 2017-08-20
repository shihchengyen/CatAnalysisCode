/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "vecc.h"
#include "libmatlbm.h"
static mxArray * _mxarray0_;

void InitializeModule_vecc(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
}

void TerminateModule_vecc(void) {
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mvecc(int nargout_, mxArray * a);

_mexLocalFunctionTable _local_function_table_vecc
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfVecc" contains the normal interface for the "vecc"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/vecc.m"
 * (lines 1-15). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfVecc(mxArray * a) {
    int nargout = 1;
    mxArray * vc = NULL;
    mlfEnterNewContext(0, 1, a);
    vc = Mvecc(nargout, a);
    mlfRestorePreviousContext(0, 1, a);
    return mlfReturnValue(vc);
}

/*
 * The function "mlxVecc" contains the feval interface for the "vecc"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/vecc.m"
 * (lines 1-15). The feval function calls the implementation version of vecc
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxVecc(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vecc Line: 1 Column: 1 The function \"vecc\""
            " was called with more than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vecc Line: 1 Column: 1 The function \"vecc"
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
    mplhs[0] = Mvecc(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mvecc" is the implementation version of the "vecc" M-function
 * from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/vecc.m"
 * (lines 1-15). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function vc = vecc(a)
 */
static mxArray * Mvecc(int nargout_, mxArray * a) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_vecc);
    mxArray * vc = NULL;
    mxArray * as = NULL;
    mclCopyArray(&a);
    /*
     * %VECC Convert to column vector
     * %   V = VECC(A) checks A to see if it is a 1 by N vector and converts
     * %   it to a N by 1 vector. Otherwise, the function returns A.
     * 
     * % get size of a
     * as = size(a);
     */
    mlfAssign(&as, mlfSize(mclValueVarargout(), mclVa(a, "a"), NULL));
    /*
     * % take transpose only if number of row is 1 and number of columns is
     * % greater than 1
     * if( (as(2)>1) && (as(1)==1) )
     */
    if (mclScalarToBool(mclGt(mclIntArrayRef1(mclVv(as, "as"), 2), _mxarray0_))
        && mclScalarToBool(
             mclEq(mclIntArrayRef1(mclVv(as, "as"), 1), _mxarray0_))) {
        /*
         * vc = a';
         */
        mlfAssign(&vc, mlfCtranspose(mclVa(a, "a")));
    /*
     * else
     */
    } else {
        /*
         * vc = a;
         */
        mlfAssign(&vc, mclVa(a, "a"));
    /*
     * end
     */
    }
    mclValidateOutput(vc, 1, nargout_, "vc", "vecc");
    mxDestroyArray(as);
    mxDestroyArray(a);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return vc;
}
