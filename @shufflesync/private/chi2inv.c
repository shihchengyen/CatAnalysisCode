/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "chi2inv.h"
#include "distchck.h"
#include "gaminv.h"
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
static double _ieee_nan_;
static mxArray * _mxarray6_;

void InitializeModule_chi2inv(void) {
    _mxarray0_ = mclInitializeString(29, _array1_);
    _mxarray2_ = mclInitializeDouble(2.0);
    _mxarray3_ = mclInitializeDouble(0.0);
    _mxarray4_ = mclInitializeString(47, _array5_);
    _ieee_nan_ = mclGetNaN();
    _mxarray6_ = mclInitializeDouble(_ieee_nan_);
}

void TerminateModule_chi2inv(void) {
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mchi2inv(int nargout_, mxArray * p, mxArray * v);

_mexLocalFunctionTable _local_function_table_chi2inv
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfChi2inv" contains the normal interface for the "chi2inv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/chi2inv.m"
 * (lines 1-37). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfChi2inv(mxArray * p, mxArray * v) {
    int nargout = 1;
    mxArray * x = NULL;
    mlfEnterNewContext(0, 2, p, v);
    x = Mchi2inv(nargout, p, v);
    mlfRestorePreviousContext(0, 2, p, v);
    return mlfReturnValue(x);
}

/*
 * The function "mlxChi2inv" contains the feval interface for the "chi2inv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/chi2inv.m"
 * (lines 1-37). The feval function calls the implementation version of chi2inv
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxChi2inv(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: chi2inv Line: 1 Column: "
            "1 The function \"chi2inv\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: chi2inv Line: 1 Column:"
            " 1 The function \"chi2inv\" was called with m"
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
    mplhs[0] = Mchi2inv(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mchi2inv" is the implementation version of the "chi2inv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/chi2inv.m"
 * (lines 1-37). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function x = chi2inv(p,v);
 */
static mxArray * Mchi2inv(int nargout_, mxArray * p, mxArray * v) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_chi2inv);
    int nargin_ = mclNargin(2, p, v, NULL);
    mxArray * x = NULL;
    mxArray * k = NULL;
    mxArray * errorcode = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&p);
    mclCopyArray(&v);
    /*
     * %CHI2INV Inverse of the chi-square cumulative distribution function (cdf).
     * %   X = CHI2INV(P,V)  returns the inverse of the chi-square cdf with V  
     * %   degrees of freedom at the values in P. The chi-square cdf with V 
     * %   degrees of freedom, is the gamma cdf with parameters V/2 and 2.   
     * %
     * %   The size of X is the common size of P and V. A scalar input
     * %   functions as a constant matrix of the same size as the other input.   
     * 
     * %   References:
     * %      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
     * %      Functions", Government Printing Office, 1964, 26.4.
     * %      [2] E. Kreyszig, "Introductory Mathematical Statistics",
     * %      John Wiley, 1970, section 10.2 (page 144)
     * 
     * %   Copyright 1993-2002 The MathWorks, Inc. 
     * %   $Revision: 2.10 $  $Date: 2002/01/17 21:30:15 $
     * 
     * if nargin < 2, 
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
     * [errorcode p v] = distchck(2,p,v);
     */
    mlfAssign(
      &errorcode,
      mlfDistchck(
        &p,
        &v,
        NULL,
        NULL,
        _mxarray2_,
        mclVa(p, "p"),
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
     * % Call the gamma inverse function. 
     * x = gaminv(p,v/2,2);
     */
    mlfAssign(
      &x,
      mlfGaminv(
        mclVa(p, "p"), mclMrdivide(mclVa(v, "v"), _mxarray2_), _mxarray2_));
    /*
     * 
     * % Return NaN if the degrees of freedom is not positive.
     * k = (v <= 0);
     */
    mlfAssign(&k, mclLe(mclVa(v, "v"), _mxarray3_));
    /*
     * if any(k(:))
     */
    if (mlfTobool(
          mlfAny(mclArrayRef1(mclVv(k, "k"), mlfCreateColonIndex()), NULL))) {
        /*
         * x(k) = NaN;
         */
        mclArrayAssign1(&x, _mxarray6_, mclVv(k, "k"));
    /*
     * end
     */
    }
    mclValidateOutput(x, 1, nargout_, "x", "chi2inv");
    mxDestroyArray(ans);
    mxDestroyArray(errorcode);
    mxDestroyArray(k);
    mxDestroyArray(v);
    mxDestroyArray(p);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return x;
}
