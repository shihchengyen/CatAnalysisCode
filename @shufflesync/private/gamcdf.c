/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "gamcdf.h"
#include "distchck.h"
#include "libmatlbm.h"
#include "libmmfile.h"
static mxArray * _mxarray0_;

static mxChar _array2_[38] = { 'R', 'e', 'q', 'u', 'i', 'r', 'e', 's', ' ', 'a',
                               't', ' ', 'l', 'e', 'a', 's', 't', ' ', 't', 'w',
                               'o', ' ', 'i', 'n', 'p', 'u', 't', ' ', 'a', 'r',
                               'g', 'u', 'm', 'e', 'n', 't', 's', '.' };
static mxArray * _mxarray1_;
static mxArray * _mxarray3_;
static mxArray * _mxarray4_;

static mxChar _array6_[47] = { 'R', 'e', 'q', 'u', 'i', 'r', 'e', 's', ' ', 'n',
                               'o', 'n', '-', 's', 'c', 'a', 'l', 'a', 'r', ' ',
                               'a', 'r', 'g', 'u', 'm', 'e', 'n', 't', 's', ' ',
                               't', 'o', ' ', 'm', 'a', 't', 'c', 'h', ' ', 'i',
                               'n', ' ', 's', 'i', 'z', 'e', '.' };
static mxArray * _mxarray5_;
static double _ieee_nan_;
static mxArray * _mxarray7_;
static mxArray * _mxarray8_;

void InitializeModule_gamcdf(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeString(38, _array2_);
    _mxarray3_ = mclInitializeDouble(3.0);
    _mxarray4_ = mclInitializeDouble(0.0);
    _mxarray5_ = mclInitializeString(47, _array6_);
    _ieee_nan_ = mclGetNaN();
    _mxarray7_ = mclInitializeDouble(_ieee_nan_);
    _mxarray8_ = mclInitializeDouble(1.7976931348623157e+308);
}

void TerminateModule_gamcdf(void) {
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mgamcdf(int nargout_, mxArray * x, mxArray * a, mxArray * b);

_mexLocalFunctionTable _local_function_table_gamcdf
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfGamcdf" contains the normal interface for the "gamcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gamcdf.m"
 * (lines 1-54). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfGamcdf(mxArray * x, mxArray * a, mxArray * b) {
    int nargout = 1;
    mxArray * p = NULL;
    mlfEnterNewContext(0, 3, x, a, b);
    p = Mgamcdf(nargout, x, a, b);
    mlfRestorePreviousContext(0, 3, x, a, b);
    return mlfReturnValue(p);
}

/*
 * The function "mlxGamcdf" contains the feval interface for the "gamcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gamcdf.m"
 * (lines 1-54). The feval function calls the implementation version of gamcdf
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxGamcdf(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: gamcdf Line: 1 Column: "
            "1 The function \"gamcdf\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: gamcdf Line: 1 Column: "
            "1 The function \"gamcdf\" was called with mor"
            "e than the declared number of inputs (3)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 3 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 3; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    mplhs[0] = Mgamcdf(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mgamcdf" is the implementation version of the "gamcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gamcdf.m"
 * (lines 1-54). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function p = gamcdf(x,a,b)
 */
static mxArray * Mgamcdf(int nargout_, mxArray * x, mxArray * a, mxArray * b) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_gamcdf);
    int nargin_ = mclNargin(3, x, a, b, NULL);
    mxArray * p = NULL;
    mxArray * k = NULL;
    mxArray * errorcode = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&x);
    mclCopyArray(&a);
    mclCopyArray(&b);
    /*
     * %GAMCDF Gamma cumulative distribution function.
     * %   P = GAMCDF(X,A,B) returns the gamma cumulative distribution
     * %   function with parameters A and B at the values in X.
     * %
     * %   The size of P is the common size of the input arguments. A scalar input  
     * %   functions as a constant matrix of the same size as the other inputs.    
     * %
     * %   Some references refer to the gamma distribution with a single
     * %   parameter. This corresponds to the default of B = 1. 
     * %
     * %   GAMMAINC does computational work.
     * 
     * %   References:
     * %      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
     * %      Springer-Verlag, 1986. p. 401.
     * %      [2]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
     * %      Functions", Government Printing Office, 1964, 26.1.32.
     * 
     * %   Copyright 1993-2002 The MathWorks, Inc. 
     * %   $Revision: 2.12 $  $Date: 2002/01/17 21:30:39 $
     * 
     * if nargin < 3, 
     */
    if (nargin_ < 3) {
        /*
         * b = 1; 
         */
        mlfAssign(&b, _mxarray0_);
    /*
     * end
     */
    }
    /*
     * 
     * if nargin < 2, 
     */
    if (nargin_ < 2) {
        /*
         * error('Requires at least two input arguments.'); 
         */
        mlfError(_mxarray1_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * [errorcode x a b] = distchck(3,x,a,b);
     */
    mlfAssign(
      &errorcode,
      mlfDistchck(
        &x,
        &a,
        &b,
        NULL,
        _mxarray3_,
        mclVa(x, "x"),
        mclVa(a, "a"),
        mclVa(b, "b"),
        NULL));
    /*
     * 
     * if errorcode > 0
     */
    if (mclGtBool(mclVv(errorcode, "errorcode"), _mxarray4_)) {
        /*
         * error('Requires non-scalar arguments to match in size.');
         */
        mlfError(_mxarray5_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * % Initialize P to zero.
     * p = zeros(size(x));
     */
    mlfAssign(
      &p, mlfZeros(mlfSize(mclValueVarargout(), mclVa(x, "x"), NULL), NULL));
    /*
     * 
     * %   Return NaN if the arguments are outside their respective limits.
     * p(a <= 0 | b <= 0) = NaN;
     */
    mclArrayAssign1(
      &p,
      _mxarray7_,
      mclOr(
        mclLe(mclVa(a, "a"), _mxarray4_), mclLe(mclVa(b, "b"), _mxarray4_)));
    /*
     * 
     * k = find(x > 0 & ~(a <= 0 | b <= 0));
     */
    mlfAssign(
      &k,
      mlfFind(
        NULL,
        NULL,
        mclAnd(
          mclGt(mclVa(x, "x"), _mxarray4_),
          mclNot(
            mclOr(
              mclLe(mclVa(a, "a"), _mxarray4_),
              mclLe(mclVa(b, "b"), _mxarray4_))))));
    /*
     * if any(k), 
     */
    if (mlfTobool(mlfAny(mclVv(k, "k"), NULL))) {
        /*
         * p(k) = gammainc(x(k) ./ b(k),a(k));
         */
        mclArrayAssign1(
          &p,
          mlfGammainc(
            mclRdivide(
              mclArrayRef1(mclVa(x, "x"), mclVv(k, "k")),
              mclArrayRef1(mclVa(b, "b"), mclVv(k, "k"))),
            mclArrayRef1(mclVa(a, "a"), mclVv(k, "k"))),
          mclVv(k, "k"));
    /*
     * end
     */
    }
    /*
     * 
     * % Make sure that round-off errors never make P greater than 1.
     * p(p > 1) = 1;
     */
    mclArrayAssign1(&p, _mxarray0_, mclGt(mclVv(p, "p"), _mxarray0_));
    /*
     * 
     * % If we have NaN or Inf, fix if possible
     * k = ~isfinite(p);
     */
    mlfAssign(&k, mclNot(mlfIsfinite(mclVv(p, "p"))));
    /*
     * if (any(k)), p(x>=sqrt(realmax)) = 1; end
     */
    if (mlfTobool(mlfAny(mclVv(k, "k"), NULL))) {
        mclArrayAssign1(
          &p, _mxarray0_, mclGe(mclVa(x, "x"), mlfSqrt(_mxarray8_)));
    }
    mclValidateOutput(p, 1, nargout_, "p", "gamcdf");
    mxDestroyArray(ans);
    mxDestroyArray(errorcode);
    mxDestroyArray(k);
    mxDestroyArray(b);
    mxDestroyArray(a);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return p;
}
