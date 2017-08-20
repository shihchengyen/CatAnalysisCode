/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "norminv.h"
#include "libmatlbm.h"
#include "libmmfile.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;

static mxChar _array3_[60] = { 'M', 'u', 's', 't', ' ', 'p', 'r', 'o', 'v',
                               'i', 'd', 'e', ' ', 'c', 'o', 'v', 'a', 'r',
                               'i', 'a', 'n', 'c', 'e', ' ', 'm', 'a', 't',
                               'r', 'i', 'x', ' ', 't', 'o', ' ', 'c', 'o',
                               'm', 'p', 'u', 't', 'e', ' ', 'c', 'o', 'n',
                               'f', 'i', 'd', 'e', 'n', 'c', 'e', ' ', 'b',
                               'o', 'u', 'n', 'd', 's', '.' };
static mxArray * _mxarray2_;

static double _array5_[2] = { 2.0, 2.0 };
static mxArray * _mxarray4_;

static mxChar _array7_[47] = { 'C', 'o', 'v', 'a', 'r', 'i', 'a', 'n', 'c', 'e',
                               ' ', 'm', 'a', 't', 'r', 'i', 'x', ' ', 'm', 'u',
                               's', 't', ' ', 'h', 'a', 'v', 'e', ' ', '2', ' ',
                               'r', 'o', 'w', 's', ' ', 'a', 'n', 'd', ' ', 'c',
                               'o', 'l', 'u', 'm', 'n', 's', '.' };
static mxArray * _mxarray6_;
static mxArray * _mxarray8_;

static mxChar _array10_[39] = { 'A', 'L', 'P', 'H', 'A', ' ', 'm', 'u',
                                's', 't', ' ', 'b', 'e', ' ', 'a', ' ',
                                's', 'c', 'a', 'l', 'a', 'r', ' ', 'b',
                                'e', 't', 'w', 'e', 'e', 'n', ' ', '0',
                                ' ', 'a', 'n', 'd', ' ', '1', '.' };
static mxArray * _mxarray9_;
static double _ieee_nan_;
static mxArray * _mxarray11_;
static mxArray * _mxarray12_;

static mxChar _array14_[40] = { 'N', 'o', 'n', '-', 's', 'c', 'a', 'l',
                                'a', 'r', ' ', 'a', 'r', 'g', 'u', 'm',
                                'e', 'n', 't', 's', ' ', 'm', 'u', 's',
                                't', ' ', 'm', 'a', 't', 'c', 'h', ' ',
                                'i', 'n', ' ', 's', 'i', 'z', 'e', '.' };
static mxArray * _mxarray13_;

static mxChar _array16_[45] = { 'P', 'C', 'O', 'V', ' ', 'm', 'u', 's', 't',
                                ' ', 'b', 'e', ' ', 'a', ' ', 'p', 'o', 's',
                                'i', 't', 'i', 'v', 'e', ' ', 's', 'e', 'm',
                                'i', '-', 'd', 'e', 'f', 'i', 'n', 'i', 't',
                                'e', ' ', 'm', 'a', 't', 'r', 'i', 'x', '.' };
static mxArray * _mxarray15_;

void InitializeModule_norminv(void) {
    _mxarray0_ = mclInitializeDouble(0.0);
    _mxarray1_ = mclInitializeDouble(1.0);
    _mxarray2_ = mclInitializeString(60, _array3_);
    _mxarray4_ = mclInitializeDoubleVector(1, 2, _array5_);
    _mxarray6_ = mclInitializeString(47, _array7_);
    _mxarray8_ = mclInitializeDouble(.05);
    _mxarray9_ = mclInitializeString(39, _array10_);
    _ieee_nan_ = mclGetNaN();
    _mxarray11_ = mclInitializeDouble(_ieee_nan_);
    _mxarray12_ = mclInitializeDouble(2.0);
    _mxarray13_ = mclInitializeString(40, _array14_);
    _mxarray15_ = mclInitializeString(45, _array16_);
}

void TerminateModule_norminv(void) {
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mnorminv(mxArray * * xlo,
                          mxArray * * xup,
                          int nargout_,
                          mxArray * p,
                          mxArray * mu,
                          mxArray * sigma,
                          mxArray * pcov,
                          mxArray * alpha);

_mexLocalFunctionTable _local_function_table_norminv
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfNNorminv" contains the nargout interface for the "norminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/norminv.m"
 * (lines 1-77). This interface is only produced if the M-function uses the
 * special variable "nargout". The nargout interface allows the number of
 * requested outputs to be specified via the nargout argument, as opposed to
 * the normal interface which dynamically calculates the number of outputs
 * based on the number of non-NULL inputs it receives. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
mxArray * mlfNNorminv(int nargout,
                      mxArray * * xlo,
                      mxArray * * xup,
                      mxArray * p,
                      mxArray * mu,
                      mxArray * sigma,
                      mxArray * pcov,
                      mxArray * alpha) {
    mxArray * x = NULL;
    mxArray * xlo__ = NULL;
    mxArray * xup__ = NULL;
    mlfEnterNewContext(2, 5, xlo, xup, p, mu, sigma, pcov, alpha);
    x = Mnorminv(&xlo__, &xup__, nargout, p, mu, sigma, pcov, alpha);
    mlfRestorePreviousContext(2, 5, xlo, xup, p, mu, sigma, pcov, alpha);
    if (xlo != NULL) {
        mclCopyOutputArg(xlo, xlo__);
    } else {
        mxDestroyArray(xlo__);
    }
    if (xup != NULL) {
        mclCopyOutputArg(xup, xup__);
    } else {
        mxDestroyArray(xup__);
    }
    return mlfReturnValue(x);
}

/*
 * The function "mlfNorminv" contains the normal interface for the "norminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/norminv.m"
 * (lines 1-77). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfNorminv(mxArray * * xlo,
                     mxArray * * xup,
                     mxArray * p,
                     mxArray * mu,
                     mxArray * sigma,
                     mxArray * pcov,
                     mxArray * alpha) {
    int nargout = 1;
    mxArray * x = NULL;
    mxArray * xlo__ = NULL;
    mxArray * xup__ = NULL;
    mlfEnterNewContext(2, 5, xlo, xup, p, mu, sigma, pcov, alpha);
    if (xlo != NULL) {
        ++nargout;
    }
    if (xup != NULL) {
        ++nargout;
    }
    x = Mnorminv(&xlo__, &xup__, nargout, p, mu, sigma, pcov, alpha);
    mlfRestorePreviousContext(2, 5, xlo, xup, p, mu, sigma, pcov, alpha);
    if (xlo != NULL) {
        mclCopyOutputArg(xlo, xlo__);
    } else {
        mxDestroyArray(xlo__);
    }
    if (xup != NULL) {
        mclCopyOutputArg(xup, xup__);
    } else {
        mxDestroyArray(xup__);
    }
    return mlfReturnValue(x);
}

/*
 * The function "mlfVNorminv" contains the void interface for the "norminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/norminv.m"
 * (lines 1-77). The void interface is only produced if the M-function uses the
 * special variable "nargout", and has at least one output. The void interface
 * function specifies zero output arguments to the implementation version of
 * the function, and in the event that the implementation version still returns
 * an output (which, in MATLAB, would be assigned to the "ans" variable), it
 * deallocates the output. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlfVNorminv(mxArray * p,
                 mxArray * mu,
                 mxArray * sigma,
                 mxArray * pcov,
                 mxArray * alpha) {
    mxArray * x = NULL;
    mxArray * xlo = NULL;
    mxArray * xup = NULL;
    mlfEnterNewContext(0, 5, p, mu, sigma, pcov, alpha);
    x = Mnorminv(&xlo, &xup, 0, p, mu, sigma, pcov, alpha);
    mlfRestorePreviousContext(0, 5, p, mu, sigma, pcov, alpha);
    mxDestroyArray(x);
    mxDestroyArray(xlo);
    mxDestroyArray(xup);
}

/*
 * The function "mlxNorminv" contains the feval interface for the "norminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/norminv.m"
 * (lines 1-77). The feval function calls the implementation version of norminv
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxNorminv(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[5];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: norminv Line: 1 Column: "
            "1 The function \"norminv\" was called with mor"
            "e than the declared number of outputs (3)."),
          NULL);
    }
    if (nrhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: norminv Line: 1 Column:"
            " 1 The function \"norminv\" was called with m"
            "ore than the declared number of inputs (5)."),
          NULL);
    }
    for (i = 0; i < 3; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 5 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 5; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    mplhs[0]
      = Mnorminv(
          &mplhs[1],
          &mplhs[2],
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2],
          mprhs[3],
          mprhs[4]);
    mlfRestorePreviousContext(
      0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 3 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 3; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mnorminv" is the implementation version of the "norminv"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/norminv.m"
 * (lines 1-77). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [x,xlo,xup] = norminv(p,mu,sigma,pcov,alpha)
 */
static mxArray * Mnorminv(mxArray * * xlo,
                          mxArray * * xup,
                          int nargout_,
                          mxArray * p,
                          mxArray * mu,
                          mxArray * sigma,
                          mxArray * pcov,
                          mxArray * alpha) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_norminv);
    int nargin_ = mclNargin(5, p, mu, sigma, pcov, alpha, NULL);
    mxArray * x = NULL;
    mxArray * halfwidth = NULL;
    mxArray * normz = NULL;
    mxArray * xvar = NULL;
    mxArray * x0 = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&p);
    mclCopyArray(&mu);
    mclCopyArray(&sigma);
    mclCopyArray(&pcov);
    mclCopyArray(&alpha);
    /*
     * %NORMINV Inverse of the normal cumulative distribution function (cdf).
     * %   X = NORMINV(P,MU,SIGMA) returns the inverse cdf for the normal
     * %   distribution with mean MU and standard deviation SIGMA, evaluated at
     * %   the values in P.  The size of X is the common size of the input
     * %   arguments.  A scalar input functions as a constant matrix of the same
     * %   size as the other inputs.
     * %
     * %   Default values for MU and SIGMA are 0 and 1, respectively.
     * %
     * %   [X,XLO,XUP] = NORMINV(P,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
     * %   for X when the input parameters MU and SIGMA are estimates.  PCOV is a
     * %   2-by-2 matrix containing the covariance matrix of the estimated parameters.
     * %   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
     * %   bounds.  XLO and XUP are arrays of the same size as X containing the lower
     * %   and upper confidence bounds.
     * %
     * %   See also ERFINV, ERFCINV, NORMCDF, NORMFIT, NORMLIKE, NORMPDF,
     * %            NORMRND, NORMSTAT.
     * 
     * 
     * %   References:
     * %      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
     * %          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
     * %      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
     * %          Distributions, 2nd ed., Wiley, 170pp.
     * 
     * 
     * %   Copyright 1993-2003 The MathWorks, Inc. 
     * %   $Revision: 2.16 $  $Date: 2003/01/09 21:47:32 $
     * 
     * 
     * if nargin < 2
     */
    if (nargin_ < 2) {
        /*
         * mu = 0;
         */
        mlfAssign(&mu, _mxarray0_);
    /*
     * end
     */
    }
    /*
     * if nargin < 3
     */
    if (nargin_ < 3) {
        /*
         * sigma = 1;
         */
        mlfAssign(&sigma, _mxarray1_);
    /*
     * end
     */
    }
    /*
     * 
     * % More checking if we need to compute confidence bounds.
     * if nargout>2
     */
    if (nargout_ > 2) {
        /*
         * if nargin<4
         */
        if (nargin_ < 4) {
            /*
             * error('Must provide covariance matrix to compute confidence bounds.');
             */
            mlfError(_mxarray2_, NULL);
        /*
         * end
         */
        }
        /*
         * if ~isequal(size(pcov),[2 2])
         */
        if (mclNotBool(
              mlfIsequal(
                mlfSize(mclValueVarargout(), mclVa(pcov, "pcov"), NULL),
                _mxarray4_, NULL))) {
            /*
             * error('Covariance matrix must have 2 rows and columns.');
             */
            mlfError(_mxarray6_, NULL);
        /*
         * end
         */
        }
        /*
         * if nargin<5
         */
        if (nargin_ < 5) {
            /*
             * alpha = 0.05;
             */
            mlfAssign(&alpha, _mxarray8_);
        /*
         * elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
         */
        } else if (! mclScalarToBool(mlfIsnumeric(mclVa(alpha, "alpha")))
                   || mclScalarToBool(
                        mclNe(mlfNumel(mclVa(alpha, "alpha")), _mxarray1_))
                   || mclScalarToBool(mclLe(mclVa(alpha, "alpha"), _mxarray0_))
                   || mclScalarToBool(
                        mclGe(mclVa(alpha, "alpha"), _mxarray1_))) {
            /*
             * error('ALPHA must be a scalar between 0 and 1.');
             */
            mlfError(_mxarray9_, NULL);
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    /*
     * 
     * % Return NaN for out of range parameters or probabilities.
     * sigma(sigma <= 0) = NaN;
     */
    mclArrayAssign1(
      &sigma, _mxarray11_, mclLe(mclVa(sigma, "sigma"), _mxarray0_));
    /*
     * p(p < 0 | 1 < p) = NaN;
     */
    mclArrayAssign1(
      &p,
      _mxarray11_,
      mclOr(
        mclLt(mclVa(p, "p"), _mxarray0_), mclLt(_mxarray1_, mclVa(p, "p"))));
    /*
     * 
     * 
     * x0 = -sqrt(2).*erfcinv(2*p);
     */
    mlfAssign(
      &x0,
      mclTimes(
        mclUminus(mlfSqrt(_mxarray12_)),
        mlfErfcinv(mclMtimes(_mxarray12_, mclVa(p, "p")))));
    /*
     * try
     */
    mlfTry {
        /*
         * x = sigma.*x0 + mu;
         */
        mlfAssign(
          &x,
          mclPlus(
            mclTimes(mclVa(sigma, "sigma"), mclVv(x0, "x0")), mclVa(mu, "mu")));
    /*
     * catch
     */
    } mlfCatch {
        /*
         * error('Non-scalar arguments must match in size.');
         */
        mlfError(_mxarray13_, NULL);
    /*
     * end
     */
    } mlfEndCatch
    /*
     * 
     * % Compute confidence bounds if requested.
     * if nargout>=2
     */
    if (nargout_ >= 2) {
        /*
         * xvar = pcov(1,1) + 2*pcov(1,2)*x0 + pcov(2,2)*x0.^2;
         */
        mlfAssign(
          &xvar,
          mclPlus(
            mclPlus(
              mclIntArrayRef2(mclVa(pcov, "pcov"), 1, 1),
              mclMtimes(
                mclMtimes(
                  _mxarray12_, mclIntArrayRef2(mclVa(pcov, "pcov"), 1, 2)),
                mclVv(x0, "x0"))),
            mclMtimes(
              mclIntArrayRef2(mclVa(pcov, "pcov"), 2, 2),
              mlfPower(mclVv(x0, "x0"), _mxarray12_))));
        /*
         * if any(xvar<0)
         */
        if (mlfTobool(mlfAny(mclLt(mclVv(xvar, "xvar"), _mxarray0_), NULL))) {
            /*
             * error('PCOV must be a positive semi-definite matrix.');
             */
            mlfError(_mxarray15_, NULL);
        /*
         * end
         */
        }
        /*
         * normz = -norminv(alpha/2);
         */
        mlfAssign(
          &normz,
          mclUminus(
            mlfNNorminv(
              1,
              NULL,
              NULL,
              mclMrdivide(mclVa(alpha, "alpha"), _mxarray12_),
              NULL,
              NULL,
              NULL,
              NULL)));
        /*
         * halfwidth = normz * sqrt(xvar);
         */
        mlfAssign(
          &halfwidth,
          mclMtimes(mclVv(normz, "normz"), mlfSqrt(mclVv(xvar, "xvar"))));
        /*
         * xlo = x - halfwidth;
         */
        mlfAssign(xlo, mclMinus(mclVv(x, "x"), mclVv(halfwidth, "halfwidth")));
        /*
         * xup = x + halfwidth;
         */
        mlfAssign(xup, mclPlus(mclVv(x, "x"), mclVv(halfwidth, "halfwidth")));
    /*
     * end
     */
    }
    mclValidateOutput(x, 1, nargout_, "x", "norminv");
    mclValidateOutput(*xlo, 2, nargout_, "xlo", "norminv");
    mclValidateOutput(*xup, 3, nargout_, "xup", "norminv");
    mxDestroyArray(ans);
    mxDestroyArray(x0);
    mxDestroyArray(xvar);
    mxDestroyArray(normz);
    mxDestroyArray(halfwidth);
    mxDestroyArray(alpha);
    mxDestroyArray(pcov);
    mxDestroyArray(sigma);
    mxDestroyArray(mu);
    mxDestroyArray(p);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return x;
}
