/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "chi2cdf.h"
#include "distchck.h"
#include "gamcdf.h"
#include "libmatlbm.h"

static mxChar _array1_[29] = { 'R', 'e', 'q', 'u', 'i', 'r', 'e', 's', ' ', 't',
                               'w', 'o', ' ', 'i', 'n', 'p', 'u', 't', ' ', 'a',
                               'r', 'g', 'u', 'm', 'e', 'n', 't', 's', '.' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;

static mxChar _array5_[47] = { 'R', 'e', 'q', 'u', 'i', 'r', 'e', 's', ' ', 'n',
                               'o', 'n', '-', 's', 'c', 'a', 'l', 'a', 'r', ' ',
                               'a', 'r', 'g', 'u', 'm', 'e', 'n', 't', 's', ' ',
                               't', 'o', ' ', 'm', 'a', 't', 'c', 'h', ' ', 'i',
                               'n', ' ', 's', 'i', 'z', 'e', '.' };
static mxArray * _mxarray4_;

void InitializeModule_chi2cdf(void) {
    _mxarray0_ = mclInitializeString(29, _array1_);
    _mxarray2_ = mclInitializeDouble(2.0);
    _mxarray3_ = mclInitializeDouble(0.0);
    _mxarray4_ = mclInitializeString(47, _array5_);
}

void TerminateModule_chi2cdf(void) {
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mchi2cdf(int nargout_, mxArray * x, mxArray * v);

_mexLocalFunctionTable _local_function_table_chi2cdf
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfChi2cdf" contains the normal interface for the "chi2cdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/chi2cdf.m"
 * (lines 1-34). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfChi2cdf(mxArray * x, mxArray * v) {
    int nargout = 1;
    mxArray * p = NULL;
    mlfEnterNewContext(0, 2, x, v);
    p = Mchi2cdf(nargout, x, v);
    mlfRestorePreviousContext(0, 2, x, v);
    return mlfReturnValue(p);
}

/*
 * The function "mlxChi2cdf" contains the feval interface for the "chi2cdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/chi2cdf.m"
 * (lines 1-34). The feval function calls the implementation version of chi2cdf
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxChi2cdf(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: chi2cdf Line: 1 Column: "
            "1 The function \"chi2cdf\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: chi2cdf Line: 1 Column:"
            " 1 The function \"chi2cdf\" was called with m"
            "ore than the declared number of inputs (2)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = Mchi2cdf(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mchi2cdf" is the implementation version of the "chi2cdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/chi2cdf.m"
 * (lines 1-34). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function p = chi2cdf(x,v)
 */
static mxArray * Mchi2cdf(int nargout_, mxArray * x, mxArray * v) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_chi2cdf);
    int nargin_ = mclNargin(2, x, v, NULL);
    mxArray * p = NULL;
    mxArray * errorcode = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&x);
    mclCopyArray(&v);
    /*
     * %CHI2CDF Chi-square cumulative distribution function.
     * %   P = CHI2CDF(X,V) returns the chi-square cumulative distribution
     * %   function with V degrees of freedom at the values in X.
     * %   The chi-square density function with V degrees of freedom,
     * %   is the same as a gamma density function with parameters V/2 and 2.
     * %
     * %   The size of P is the common size of X and V. A scalar input   
     * %   functions as a constant matrix of the same size as the other input.    
     * 
     * %   References:
     * %      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
     * %      Functions", Government Printing Office, 1964, 26.4.
     * %
     * %   Notice that we do not check if the degree of freedom parameter is integer
     * %   or not. In most cases, it should be an integer. Numerically, non-integer 
     * %   values still gives a numerical answer, thus, we keep them.
     * 
     * %   Copyright 1993-2002 The MathWorks, Inc. 
     * %   $Revision: 2.10 $  $Date: 2002/01/17 21:30:14 $
     * 
     * if   nargin < 2, 
     */
    if (nargin_ < 2) {
        /*
         * error('Requires two input arguments.');
         */
        mlfError(_mxarray0_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * [errorcode x v] = distchck(2,x,v);
     */
    mlfAssign(
      &errorcode,
      mlfDistchck(
        &x,
        &v,
        NULL,
        NULL,
        _mxarray2_,
        mclVa(x, "x"),
        mclVa(v, "v"),
        NULL,
        NULL));
    /*
     * 
     * if errorcode > 0
     */
    if (mclGtBool(mclVv(errorcode, "errorcode"), _mxarray3_)) {
        /*
         * error('Requires non-scalar arguments to match in size.');
         */
        mlfError(_mxarray4_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * % Call the gamma distribution function. 
     * p = gamcdf(x,v/2,2);
     */
    mlfAssign(
      &p,
      mlfGamcdf(
        mclVa(x, "x"), mclMrdivide(mclVa(v, "v"), _mxarray2_), _mxarray2_));
    mclValidateOutput(p, 1, nargout_, "p", "chi2cdf");
    mxDestroyArray(ans);
    mxDestroyArray(errorcode);
    mxDestroyArray(v);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return p;
}
