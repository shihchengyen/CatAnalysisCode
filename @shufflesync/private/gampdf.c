/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "gampdf.h"
#include "distchck.h"
#include "libmatlbm.h"
#include "libmmfile.h"
static mxArray * _mxarray0_;

static mxChar _array2_[37] = { 'R', 'e', 'q', 'u', 'i', 'r', 'e', 's', ' ', 'a',
                               't', ' ', 'l', 'e', 'a', 's', 't', ' ', 't', 'w',
                               'o', ' ', 'i', 'n', 'p', 'u', 't', ' ', 'a', 'r',
                               'g', 'u', 'm', 'e', 'n', 't', 's' };
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
static double _ieee_plusinf_;
static mxArray * _mxarray8_;

void InitializeModule_gampdf(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeString(37, _array2_);
    _mxarray3_ = mclInitializeDouble(3.0);
    _mxarray4_ = mclInitializeDouble(0.0);
    _mxarray5_ = mclInitializeString(47, _array6_);
    _ieee_nan_ = mclGetNaN();
    _mxarray7_ = mclInitializeDouble(_ieee_nan_);
    _ieee_plusinf_ = mclGetInf();
    _mxarray8_ = mclInitializeDouble(_ieee_plusinf_);
}

void TerminateModule_gampdf(void) {
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mgampdf(int nargout_, mxArray * x, mxArray * a, mxArray * b);

_mexLocalFunctionTable _local_function_table_gampdf
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfGampdf" contains the normal interface for the "gampdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gampdf.m"
 * (lines 1-49). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfGampdf(mxArray * x, mxArray * a, mxArray * b) {
    int nargout = 1;
    mxArray * y = NULL;
    mlfEnterNewContext(0, 3, x, a, b);
    y = Mgampdf(nargout, x, a, b);
    mlfRestorePreviousContext(0, 3, x, a, b);
    return mlfReturnValue(y);
}

/*
 * The function "mlxGampdf" contains the feval interface for the "gampdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gampdf.m"
 * (lines 1-49). The feval function calls the implementation version of gampdf
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxGampdf(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: gampdf Line: 1 Column: "
            "1 The function \"gampdf\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: gampdf Line: 1 Column: "
            "1 The function \"gampdf\" was called with mor"
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
    mplhs[0] = Mgampdf(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mgampdf" is the implementation version of the "gampdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gampdf.m"
 * (lines 1-49). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function y = gampdf(x,a,b)
 */
static mxArray * Mgampdf(int nargout_, mxArray * x, mxArray * a, mxArray * b) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_gampdf);
    int nargin_ = mclNargin(3, x, a, b, NULL);
    mxArray * y = NULL;
    mxArray * k2 = NULL;
    mxArray * k = NULL;
    mxArray * errorcode = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&x);
    mclCopyArray(&a);
    mclCopyArray(&b);
    /*
     * %GAMPDF Gamma probability density function.
     * %   Y = GAMPDF(X,A,B) returns the gamma probability density function 
     * %   with parameters A and B, at the values in X.
     * %
     * %   The size of Y is the common size of the input arguments. A scalar input  
     * %   functions as a constant matrix of the same size as the other inputs.    
     * %
     * %   Some references refer to the gamma distribution with a single
     * %   parameter. This corresponds to the default of B = 1.
     * 
     * %   References:
     * %      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
     * %      Springer-Verlag, 1986, pages 401-402.
     * 
     * %   Copyright 1993-2002 The MathWorks, Inc. 
     * %   $Revision: 2.10 $  $Date: 2002/01/17 21:30:41 $
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
         * error('Requires at least two input arguments'); 
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
     * % Initialize Y to zero.
     * y = zeros(size(x));
     */
    mlfAssign(
      &y, mlfZeros(mlfSize(mclValueVarargout(), mclVa(x, "x"), NULL), NULL));
    /*
     * 
     * %   Return NaN if the arguments are outside their respective limits.
     * y(a <= 0 | b <= 0) = NaN;     
     */
    mclArrayAssign1(
      &y,
      _mxarray7_,
      mclOr(
        mclLe(mclVa(a, "a"), _mxarray4_), mclLe(mclVa(b, "b"), _mxarray4_)));
    /*
     * 
     * k=find(x > 0 & ~(a <= 0 | b <= 0));
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
     * if any(k)
     */
    if (mlfTobool(mlfAny(mclVv(k, "k"), NULL))) {
        /*
         * y(k) = (a(k) - 1) .* log(x(k)) - (x(k) ./ b(k)) - gammaln(a(k)) - a(k) .* log(b(k));
         */
        mclArrayAssign1(
          &y,
          mclMinus(
            mclMinus(
              mclMinus(
                mclTimes(
                  mclMinus(
                    mclArrayRef1(mclVa(a, "a"), mclVv(k, "k")), _mxarray0_),
                  mlfLog(mclArrayRef1(mclVa(x, "x"), mclVv(k, "k")))),
                mclRdivide(
                  mclArrayRef1(mclVa(x, "x"), mclVv(k, "k")),
                  mclArrayRef1(mclVa(b, "b"), mclVv(k, "k")))),
              mlfNGammaln(
                0,
                mclValueVarargout(),
                mclArrayRef1(mclVa(a, "a"), mclVv(k, "k")),
                NULL)),
            mclTimes(
              mclArrayRef1(mclVa(a, "a"), mclVv(k, "k")),
              mlfLog(mclArrayRef1(mclVa(b, "b"), mclVv(k, "k"))))),
          mclVv(k, "k"));
        /*
         * y(k) = exp(y(k));
         */
        mclArrayAssign1(
          &y,
          mlfExp(mclArrayRef1(mclVv(y, "y"), mclVv(k, "k"))),
          mclVv(k, "k"));
    /*
     * end
     */
    }
    /*
     * y(x == 0 & a < 1) = Inf;
     */
    mclArrayAssign1(
      &y,
      _mxarray8_,
      mclAnd(
        mclEq(mclVa(x, "x"), _mxarray4_), mclLt(mclVa(a, "a"), _mxarray0_)));
    /*
     * k2 = find(x == 0 & a == 1);
     */
    mlfAssign(
      &k2,
      mlfFind(
        NULL,
        NULL,
        mclAnd(
          mclEq(mclVa(x, "x"), _mxarray4_), mclEq(mclVa(a, "a"), _mxarray0_))));
    /*
     * if any(k2)
     */
    if (mlfTobool(mlfAny(mclVv(k2, "k2"), NULL))) {
        /*
         * y(k2) = (1./b(k2));
         */
        mclArrayAssign1(
          &y,
          mclRdivide(_mxarray0_, mclArrayRef1(mclVa(b, "b"), mclVv(k2, "k2"))),
          mclVv(k2, "k2"));
    /*
     * end
     */
    }
    mclValidateOutput(y, 1, nargout_, "y", "gampdf");
    mxDestroyArray(ans);
    mxDestroyArray(errorcode);
    mxDestroyArray(k);
    mxDestroyArray(k2);
    mxDestroyArray(b);
    mxDestroyArray(a);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return y;
}
