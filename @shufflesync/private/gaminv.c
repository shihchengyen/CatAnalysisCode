/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "gaminv.h"
#include "distchck.h"
#include "gamcdf.h"
#include "gampdf.h"
#include "libmatlbm.h"
#include "norminv.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;
static mxArray * _mxarray2_;

static mxChar _array4_[47] = { 'R', 'e', 'q', 'u', 'i', 'r', 'e', 's', ' ', 'n',
                               'o', 'n', '-', 's', 'c', 'a', 'l', 'a', 'r', ' ',
                               'a', 'r', 'g', 'u', 'm', 'e', 'n', 't', 's', ' ',
                               't', 'o', ' ', 'm', 'a', 't', 'c', 'h', ' ', 'i',
                               'n', ' ', 's', 'i', 'z', 'e', '.' };
static mxArray * _mxarray3_;
static double _ieee_nan_;
static mxArray * _mxarray5_;
static double _ieee_plusinf_;
static mxArray * _mxarray6_;
static mxArray * _mxarray7_;
static mxArray * _mxarray8_;
static mxArray * _mxarray9_;
static mxArray * _mxarray10_;
static mxArray * _mxarray11_;
static mxArray * _mxarray12_;

static mxChar _array14_[37] = { 0x005c, 'n', 'W', 'a', 'r', 'n', 'i',
                                'n', 'g', ':', ' ', 'G', 'A', 'M', 'I',
                                'N', 'V', ' ', 'd', 'i', 'd', ' ', 'n',
                                'o', 't', ' ', 'c', 'o', 'n', 'v', 'e',
                                'r', 'g', 'e', '.', 0x005c, 'n' };
static mxArray * _mxarray13_;

static mxChar _array16_[20] = { 'T', 'h', 'e', ' ', 'l', 'a', 's',
                                't', ' ', 's', 't', 'e', 'p', ' ',
                                'w', 'a', 's', ':', ' ', ' ' };
static mxArray * _mxarray15_;

static mxChar _array18_[6] = { '%', '1', '3', '.', '8', 'f' };
static mxArray * _mxarray17_;

void InitializeModule_gaminv(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeDouble(3.0);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeString(47, _array4_);
    _ieee_nan_ = mclGetNaN();
    _mxarray5_ = mclInitializeDouble(_ieee_nan_);
    _ieee_plusinf_ = mclGetInf();
    _mxarray6_ = mclInitializeDouble(_ieee_plusinf_);
    _mxarray7_ = mclInitializeDouble(100.0);
    _mxarray8_ = mclInitializeDouble(2.0);
    _mxarray9_ = mclInitializeDouble(.5);
    _mxarray10_ = mclInitializeDouble(-2.0);
    _mxarray11_ = mclInitializeDouble(2.220446049250313e-16);
    _mxarray12_ = mclInitializeDouble(10.0);
    _mxarray13_ = mclInitializeString(37, _array14_);
    _mxarray15_ = mclInitializeString(20, _array16_);
    _mxarray17_ = mclInitializeString(6, _array18_);
}

void TerminateModule_gaminv(void) {
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mgaminv(int nargout_, mxArray * p, mxArray * a, mxArray * b);

_mexLocalFunctionTable _local_function_table_gaminv
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfGaminv" contains the normal interface for the "gaminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gaminv.m"
 * (lines 1-102). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
mxArray * mlfGaminv(mxArray * p, mxArray * a, mxArray * b) {
    int nargout = 1;
    mxArray * x = NULL;
    mlfEnterNewContext(0, 3, p, a, b);
    x = Mgaminv(nargout, p, a, b);
    mlfRestorePreviousContext(0, 3, p, a, b);
    return mlfReturnValue(x);
}

/*
 * The function "mlxGaminv" contains the feval interface for the "gaminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gaminv.m"
 * (lines 1-102). The feval function calls the implementation version of gaminv
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxGaminv(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: gaminv Line: 1 Column: "
            "1 The function \"gaminv\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: gaminv Line: 1 Column: "
            "1 The function \"gaminv\" was called with mor"
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
    mplhs[0] = Mgaminv(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mgaminv" is the implementation version of the "gaminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/gaminv.m"
 * (lines 1-102). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function x = gaminv(p,a,b);
 */
static mxArray * Mgaminv(int nargout_, mxArray * p, mxArray * a, mxArray * b) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_gaminv);
    int nargin_ = mclNargin(3, p, a, b, NULL);
    mxArray * x = NULL;
    mxArray * outstr = NULL;
    mxArray * str = NULL;
    mxArray * ksmall = NULL;
    mxArray * xnew = NULL;
    mxArray * h = NULL;
    mxArray * xk = NULL;
    mxArray * sigma = NULL;
    mxArray * mu = NULL;
    mxArray * temp = NULL;
    mxArray * v = NULL;
    mxArray * mn = NULL;
    mxArray * pk = NULL;
    mxArray * count = NULL;
    mxArray * count_limit = NULL;
    mxArray * k1 = NULL;
    mxArray * k0 = NULL;
    mxArray * tmp = NULL;
    mxArray * k = NULL;
    mxArray * ans = NULL;
    mxArray * errorcode = NULL;
    mclCopyArray(&p);
    mclCopyArray(&a);
    mclCopyArray(&b);
    /*
     * %GAMINV Inverse of the gamma cumulative distribution function (cdf).
     * %   X = GAMINV(P,A,B)  returns the inverse of the gamma cdf with  
     * %   parameters A and B, at the probabilities in P.
     * %
     * %   The size of X is the common size of the input arguments. A scalar input  
     * %   functions as a constant matrix of the same size as the other inputs.    
     * %
     * %   GAMINV uses Newton's method to converge to the solution.
     * 
     * %   References:
     * %      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
     * %      Functions", Government Printing Office, 1964, 6.5.
     * 
     * %   B.A. Jones 1-12-93
     * %   Copyright 1993-2002 The MathWorks, Inc. 
     * %   $Revision: 2.10 $  $Date: 2002/01/17 21:30:40 $
     * 
     * if nargin<3, 
     */
    if (nargin_ < 3) {
        /*
         * b=1;
         */
        mlfAssign(&b, _mxarray0_);
    /*
     * end
     */
    }
    /*
     * 
     * [errorcode p a b] = distchck(3,p,a,b);
     */
    mlfAssign(
      &errorcode,
      mlfDistchck(
        &p,
        &a,
        &b,
        NULL,
        _mxarray1_,
        mclVa(p, "p"),
        mclVa(a, "a"),
        mclVa(b, "b"),
        NULL));
    /*
     * 
     * if errorcode > 0
     */
    if (mclGtBool(mclVv(errorcode, "errorcode"), _mxarray2_)) {
        /*
         * error('Requires non-scalar arguments to match in size.');
         */
        mlfError(_mxarray3_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * %   Initialize X to zero.
     * x = zeros(size(p));
     */
    mlfAssign(
      &x, mlfZeros(mlfSize(mclValueVarargout(), mclVa(p, "p"), NULL), NULL));
    /*
     * 
     * k = find(p<0 | p>1 | a <= 0 | b <= 0);
     */
    mlfAssign(
      &k,
      mlfFind(
        NULL,
        NULL,
        mclOr(
          mclOr(
            mclOr(
              mclLt(mclVa(p, "p"), _mxarray2_),
              mclGt(mclVa(p, "p"), _mxarray0_)),
            mclLe(mclVa(a, "a"), _mxarray2_)),
          mclLe(mclVa(b, "b"), _mxarray2_))));
    /*
     * if any(k),
     */
    if (mlfTobool(mlfAny(mclVv(k, "k"), NULL))) {
        /*
         * tmp  = NaN;
         */
        mlfAssign(&tmp, _mxarray5_);
        /*
         * x(k) = tmp(ones(size(k)));
         */
        mclArrayAssign1(
          &x,
          mclArrayRef1(
            mclVv(tmp, "tmp"),
            mlfOnes(mlfSize(mclValueVarargout(), mclVv(k, "k"), NULL), NULL)),
          mclVv(k, "k"));
    /*
     * end
     */
    }
    /*
     * 
     * % The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.  
     * k0 = find(p == 0 & a > 0 & b > 0);
     */
    mlfAssign(
      &k0,
      mlfFind(
        NULL,
        NULL,
        mclAnd(
          mclAnd(
            mclEq(mclVa(p, "p"), _mxarray2_), mclGt(mclVa(a, "a"), _mxarray2_)),
          mclGt(mclVa(b, "b"), _mxarray2_))));
    /*
     * if any(k0),
     */
    if (mlfTobool(mlfAny(mclVv(k0, "k0"), NULL))) {
        /*
         * x(k0) = zeros(size(k0)); 
         */
        mclArrayAssign1(
          &x,
          mlfZeros(mlfSize(mclValueVarargout(), mclVv(k0, "k0"), NULL), NULL),
          mclVv(k0, "k0"));
    /*
     * end
     */
    }
    /*
     * 
     * k1 = find(p == 1 & a > 0 & b > 0);
     */
    mlfAssign(
      &k1,
      mlfFind(
        NULL,
        NULL,
        mclAnd(
          mclAnd(
            mclEq(mclVa(p, "p"), _mxarray0_), mclGt(mclVa(a, "a"), _mxarray2_)),
          mclGt(mclVa(b, "b"), _mxarray2_))));
    /*
     * if any(k1), 
     */
    if (mlfTobool(mlfAny(mclVv(k1, "k1"), NULL))) {
        /*
         * tmp = Inf;
         */
        mlfAssign(&tmp, _mxarray6_);
        /*
         * x(k1) = tmp(ones(size(k1))); 
         */
        mclArrayAssign1(
          &x,
          mclArrayRef1(
            mclVv(tmp, "tmp"),
            mlfOnes(mlfSize(mclValueVarargout(), mclVv(k1, "k1"), NULL), NULL)),
          mclVv(k1, "k1"));
    /*
     * end
     */
    }
    /*
     * 
     * % Newton's Method
     * % Permit no more than count_limit interations.
     * count_limit = 100;
     */
    mlfAssign(&count_limit, _mxarray7_);
    /*
     * count = 0;
     */
    mlfAssign(&count, _mxarray2_);
    /*
     * 
     * k = find(p > 0  &  p < 1 & a > 0 & b > 0);
     */
    mlfAssign(
      &k,
      mlfFind(
        NULL,
        NULL,
        mclAnd(
          mclAnd(
            mclAnd(
              mclGt(mclVa(p, "p"), _mxarray2_),
              mclLt(mclVa(p, "p"), _mxarray0_)),
            mclGt(mclVa(a, "a"), _mxarray2_)),
          mclGt(mclVa(b, "b"), _mxarray2_))));
    /*
     * if (~any(k(:))), return; end
     */
    if (mclNotBool(
          mlfAny(mclArrayRef1(mclVv(k, "k"), mlfCreateColonIndex()), NULL))) {
        goto return_;
    }
    /*
     * pk = p(k);
     */
    mlfAssign(&pk, mclArrayRef1(mclVa(p, "p"), mclVv(k, "k")));
    /*
     * 
     * % Supply a starting guess for the iteration.
     * %   Use a method of moments fit to the lognormal distribution. 
     * mn = a(k) .* b(k);
     */
    mlfAssign(
      &mn,
      mclTimes(
        mclArrayRef1(mclVa(a, "a"), mclVv(k, "k")),
        mclArrayRef1(mclVa(b, "b"), mclVv(k, "k"))));
    /*
     * v = mn .* b(k);
     */
    mlfAssign(
      &v,
      mclTimes(mclVv(mn, "mn"), mclArrayRef1(mclVa(b, "b"), mclVv(k, "k"))));
    /*
     * temp = log(v + mn .^ 2); 
     */
    mlfAssign(
      &temp,
      mlfLog(mclPlus(mclVv(v, "v"), mlfPower(mclVv(mn, "mn"), _mxarray8_))));
    /*
     * mu = 2 * log(mn) - 0.5 * temp;
     */
    mlfAssign(
      &mu,
      mclMinus(
        mclMtimes(_mxarray8_, mlfLog(mclVv(mn, "mn"))),
        mclMtimes(_mxarray9_, mclVv(temp, "temp"))));
    /*
     * sigma = -2 * log(mn) + temp;
     */
    mlfAssign(
      &sigma,
      mclPlus(
        mclMtimes(_mxarray10_, mlfLog(mclVv(mn, "mn"))), mclVv(temp, "temp")));
    /*
     * xk = exp(norminv(pk,mu,sigma));
     */
    mlfAssign(
      &xk,
      mlfExp(
        mlfNNorminv(
          1,
          NULL,
          NULL,
          mclVv(pk, "pk"),
          mclVv(mu, "mu"),
          mclVv(sigma, "sigma"),
          NULL,
          NULL)));
    /*
     * 
     * h = ones(size(pk)); 
     */
    mlfAssign(
      &h, mlfOnes(mlfSize(mclValueVarargout(), mclVv(pk, "pk"), NULL), NULL));
    /*
     * 
     * % Break out of the iteration loop for three reasons:
     * %  1) the last update is very small (compared to x)
     * %  2) the last update is very small (compared to sqrt(eps))
     * %  3) There are more than 100 iterations. This should NEVER happen. 
     * 
     * while(any(abs(h) > sqrt(eps)*abs(xk))  &  max(abs(h)) > sqrt(eps)    ...
     */
    for (;;) {
        mxArray * a_
          = mclInitialize(
              mlfAny(
                mclGt(
                  mlfAbs(mclVv(h, "h")),
                  mclMtimes(mlfSqrt(_mxarray11_), mlfAbs(mclVv(xk, "xk")))),
                NULL));
        if (mlfTobool(a_)) {
            mlfAssign(
              &a_,
              mclAnd(
                a_,
                mclGt(
                  mlfMax(NULL, mlfAbs(mclVv(h, "h")), NULL, NULL),
                  mlfSqrt(_mxarray11_))));
        } else {
            mlfAssign(&a_, mlfScalar(0));
        }
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(
                   a_,
                   mclLt(
                     mclVv(count, "count"),
                     mclVv(count_limit, "count_limit"))))) {
            mxDestroyArray(a_);
        } else {
            mxDestroyArray(a_);
            break;
        }
        /*
         * & count < count_limit), 
         * 
         * count = count + 1;
         */
        mlfAssign(&count, mclPlus(mclVv(count, "count"), _mxarray0_));
        /*
         * h = (gamcdf(xk,a(k),b(k)) - pk) ./ gampdf(xk,a(k),b(k));
         */
        mlfAssign(
          &h,
          mclRdivide(
            mclMinus(
              mlfGamcdf(
                mclVv(xk, "xk"),
                mclArrayRef1(mclVa(a, "a"), mclVv(k, "k")),
                mclArrayRef1(mclVa(b, "b"), mclVv(k, "k"))),
              mclVv(pk, "pk")),
            mlfGampdf(
              mclVv(xk, "xk"),
              mclArrayRef1(mclVa(a, "a"), mclVv(k, "k")),
              mclArrayRef1(mclVa(b, "b"), mclVv(k, "k")))));
        /*
         * xnew = xk - h;
         */
        mlfAssign(&xnew, mclMinus(mclVv(xk, "xk"), mclVv(h, "h")));
        /*
         * % Make sure that the current guess stays greater than zero.
         * % When Newton's Method suggests steps that lead to negative guesses
         * % take a step 9/10ths of the way to zero:
         * ksmall = find(xnew < 0);
         */
        mlfAssign(
          &ksmall, mlfFind(NULL, NULL, mclLt(mclVv(xnew, "xnew"), _mxarray2_)));
        /*
         * if any(ksmall),
         */
        if (mlfTobool(mlfAny(mclVv(ksmall, "ksmall"), NULL))) {
            /*
             * xnew(ksmall) = xk(ksmall) / 10;
             */
            mclArrayAssign1(
              &xnew,
              mclMrdivide(
                mclArrayRef1(mclVv(xk, "xk"), mclVv(ksmall, "ksmall")),
                _mxarray12_),
              mclVv(ksmall, "ksmall"));
            /*
             * h = xk-xnew;
             */
            mlfAssign(&h, mclMinus(mclVv(xk, "xk"), mclVv(xnew, "xnew")));
        /*
         * end
         */
        }
        /*
         * xk = xnew;
         */
        mlfAssign(&xk, mclVv(xnew, "xnew"));
    /*
     * end
     */
    }
    /*
     * 
     * 
     * % Store the converged value in the correct place
     * x(k) = xk;
     */
    mclArrayAssign1(&x, mclVv(xk, "xk"), mclVv(k, "k"));
    /*
     * 
     * if count == count_limit, 
     */
    if (mclEqBool(mclVv(count, "count"), mclVv(count_limit, "count_limit"))) {
        /*
         * fprintf('\nWarning: GAMINV did not converge.\n');
         */
        mclAssignAns(&ans, mlfNFprintf(0, _mxarray13_, NULL));
        /*
         * str = 'The last step was:  ';
         */
        mlfAssign(&str, _mxarray15_);
        /*
         * outstr = sprintf([str,'%13.8f'],h);
         */
        mlfAssign(
          &outstr,
          mlfSprintf(
            NULL,
            mlfHorzcat(mclVv(str, "str"), _mxarray17_, NULL),
            mclVv(h, "h"),
            NULL));
        /*
         * fprintf(outstr);
         */
        mclAssignAns(&ans, mlfNFprintf(0, mclVv(outstr, "outstr"), NULL));
    /*
     * end
     */
    }
    return_:
    mclValidateOutput(x, 1, nargout_, "x", "gaminv");
    mxDestroyArray(errorcode);
    mxDestroyArray(ans);
    mxDestroyArray(k);
    mxDestroyArray(tmp);
    mxDestroyArray(k0);
    mxDestroyArray(k1);
    mxDestroyArray(count_limit);
    mxDestroyArray(count);
    mxDestroyArray(pk);
    mxDestroyArray(mn);
    mxDestroyArray(v);
    mxDestroyArray(temp);
    mxDestroyArray(mu);
    mxDestroyArray(sigma);
    mxDestroyArray(xk);
    mxDestroyArray(h);
    mxDestroyArray(xnew);
    mxDestroyArray(ksmall);
    mxDestroyArray(str);
    mxDestroyArray(outstr);
    mxDestroyArray(b);
    mxDestroyArray(a);
    mxDestroyArray(p);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return x;
}
