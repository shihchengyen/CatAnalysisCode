/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "normcdf.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "norminv.h"
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

static mxChar _array13_[40] = { 'N', 'o', 'n', '-', 's', 'c', 'a', 'l',
                                'a', 'r', ' ', 'a', 'r', 'g', 'u', 'm',
                                'e', 'n', 't', 's', ' ', 'm', 'u', 's',
                                't', ' ', 'm', 'a', 't', 'c', 'h', ' ',
                                'i', 'n', ' ', 's', 'i', 'z', 'e', '.' };
static mxArray * _mxarray12_;
static mxArray * _mxarray14_;
static mxArray * _mxarray15_;

static mxChar _array17_[45] = { 'P', 'C', 'O', 'V', ' ', 'm', 'u', 's', 't',
                                ' ', 'b', 'e', ' ', 'a', ' ', 'p', 'o', 's',
                                'i', 't', 'i', 'v', 'e', ' ', 's', 'e', 'm',
                                'i', '-', 'd', 'e', 'f', 'i', 'n', 'i', 't',
                                'e', ' ', 'm', 'a', 't', 'r', 'i', 'x', '.' };
static mxArray * _mxarray16_;

void InitializeModule_normcdf(void) {
    _mxarray0_ = mclInitializeDouble(0.0);
    _mxarray1_ = mclInitializeDouble(1.0);
    _mxarray2_ = mclInitializeString(60, _array3_);
    _mxarray4_ = mclInitializeDoubleVector(1, 2, _array5_);
    _mxarray6_ = mclInitializeString(47, _array7_);
    _mxarray8_ = mclInitializeDouble(.05);
    _mxarray9_ = mclInitializeString(39, _array10_);
    _ieee_nan_ = mclGetNaN();
    _mxarray11_ = mclInitializeDouble(_ieee_nan_);
    _mxarray12_ = mclInitializeString(40, _array13_);
    _mxarray14_ = mclInitializeDouble(.5);
    _mxarray15_ = mclInitializeDouble(2.0);
    _mxarray16_ = mclInitializeString(45, _array17_);
}

void TerminateModule_normcdf(void) {
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray14_);
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

static mxArray * Mnormcdf(mxArray * * plo,
                          mxArray * * pup,
                          int nargout_,
                          mxArray * x,
                          mxArray * mu,
                          mxArray * sigma,
                          mxArray * pcov,
                          mxArray * alpha);

_mexLocalFunctionTable _local_function_table_normcdf
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfNNormcdf" contains the nargout interface for the "normcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/normcdf.m"
 * (lines 1-82). This interface is only produced if the M-function uses the
 * special variable "nargout". The nargout interface allows the number of
 * requested outputs to be specified via the nargout argument, as opposed to
 * the normal interface which dynamically calculates the number of outputs
 * based on the number of non-NULL inputs it receives. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
mxArray * mlfNNormcdf(int nargout,
                      mxArray * * plo,
                      mxArray * * pup,
                      mxArray * x,
                      mxArray * mu,
                      mxArray * sigma,
                      mxArray * pcov,
                      mxArray * alpha) {
    mxArray * p = NULL;
    mxArray * plo__ = NULL;
    mxArray * pup__ = NULL;
    mlfEnterNewContext(2, 5, plo, pup, x, mu, sigma, pcov, alpha);
    p = Mnormcdf(&plo__, &pup__, nargout, x, mu, sigma, pcov, alpha);
    mlfRestorePreviousContext(2, 5, plo, pup, x, mu, sigma, pcov, alpha);
    if (plo != NULL) {
        mclCopyOutputArg(plo, plo__);
    } else {
        mxDestroyArray(plo__);
    }
    if (pup != NULL) {
        mclCopyOutputArg(pup, pup__);
    } else {
        mxDestroyArray(pup__);
    }
    return mlfReturnValue(p);
}

/*
 * The function "mlfNormcdf" contains the normal interface for the "normcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/normcdf.m"
 * (lines 1-82). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfNormcdf(mxArray * * plo,
                     mxArray * * pup,
                     mxArray * x,
                     mxArray * mu,
                     mxArray * sigma,
                     mxArray * pcov,
                     mxArray * alpha) {
    int nargout = 1;
    mxArray * p = NULL;
    mxArray * plo__ = NULL;
    mxArray * pup__ = NULL;
    mlfEnterNewContext(2, 5, plo, pup, x, mu, sigma, pcov, alpha);
    if (plo != NULL) {
        ++nargout;
    }
    if (pup != NULL) {
        ++nargout;
    }
    p = Mnormcdf(&plo__, &pup__, nargout, x, mu, sigma, pcov, alpha);
    mlfRestorePreviousContext(2, 5, plo, pup, x, mu, sigma, pcov, alpha);
    if (plo != NULL) {
        mclCopyOutputArg(plo, plo__);
    } else {
        mxDestroyArray(plo__);
    }
    if (pup != NULL) {
        mclCopyOutputArg(pup, pup__);
    } else {
        mxDestroyArray(pup__);
    }
    return mlfReturnValue(p);
}

/*
 * The function "mlfVNormcdf" contains the void interface for the "normcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/normcdf.m"
 * (lines 1-82). The void interface is only produced if the M-function uses the
 * special variable "nargout", and has at least one output. The void interface
 * function specifies zero output arguments to the implementation version of
 * the function, and in the event that the implementation version still returns
 * an output (which, in MATLAB, would be assigned to the "ans" variable), it
 * deallocates the output. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlfVNormcdf(mxArray * x,
                 mxArray * mu,
                 mxArray * sigma,
                 mxArray * pcov,
                 mxArray * alpha) {
    mxArray * p = NULL;
    mxArray * plo = NULL;
    mxArray * pup = NULL;
    mlfEnterNewContext(0, 5, x, mu, sigma, pcov, alpha);
    p = Mnormcdf(&plo, &pup, 0, x, mu, sigma, pcov, alpha);
    mlfRestorePreviousContext(0, 5, x, mu, sigma, pcov, alpha);
    mxDestroyArray(p);
    mxDestroyArray(plo);
    mxDestroyArray(pup);
}

/*
 * The function "mlxNormcdf" contains the feval interface for the "normcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/normcdf.m"
 * (lines 1-82). The feval function calls the implementation version of normcdf
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxNormcdf(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[5];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: normcdf Line: 1 Column: "
            "1 The function \"normcdf\" was called with mor"
            "e than the declared number of outputs (3)."),
          NULL);
    }
    if (nrhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: normcdf Line: 1 Column:"
            " 1 The function \"normcdf\" was called with m"
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
      = Mnormcdf(
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
 * The function "Mnormcdf" is the implementation version of the "normcdf"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/normcdf.m"
 * (lines 1-82). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [p,plo,pup] = normcdf(x,mu,sigma,pcov,alpha)
 */
static mxArray * Mnormcdf(mxArray * * plo,
                          mxArray * * pup,
                          int nargout_,
                          mxArray * x,
                          mxArray * mu,
                          mxArray * sigma,
                          mxArray * pcov,
                          mxArray * alpha) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_normcdf);
    int nargin_ = mclNargin(5, x, mu, sigma, pcov, alpha, NULL);
    mxArray * p = NULL;
    mxArray * zup = NULL;
    mxArray * zlo = NULL;
    mxArray * halfwidth = NULL;
    mxArray * normz = NULL;
    mxArray * zvar = NULL;
    mxArray * z = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&x);
    mclCopyArray(&mu);
    mclCopyArray(&sigma);
    mclCopyArray(&pcov);
    mclCopyArray(&alpha);
    /*
     * %NORMCDF Normal cumulative distribution function (cdf).
     * %   P = NORMCDF(X,MU,SIGMA) returns the cdf of the normal distribution with
     * %   mean MU and standard deviation SIGMA, evaluated at the values in X.
     * %   The size of P is the common size of X, MU and SIGMA.  A scalar input
     * %   functions as a constant matrix of the same size as the other inputs.
     * %
     * %   Default values for MU and SIGMA are 0 and 1, respectively.
     * %
     * %   [P,PLO,PUP] = NORMCDF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
     * %   for P when the input parameters MU and SIGMA are estimates.  PCOV is a
     * %   2-by-2 matrix containing the covariance matrix of the estimated parameters.
     * %   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
     * %   bounds.  PLO and PUP are arrays of the same size as P containing the lower
     * %   and upper confidence bounds.
     * %
     * %   See also ERF, ERFC, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMRND, NORMSTAT.
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
     * %   $Revision: 2.16 $  $Date: 2003/01/09 21:47:31 $
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
     * if nargout>1
     */
    if (nargout_ > 1) {
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
     * % Return NaN for out of range parameters.
     * sigma(sigma <= 0) = NaN;
     */
    mclArrayAssign1(
      &sigma, _mxarray11_, mclLe(mclVa(sigma, "sigma"), _mxarray0_));
    /*
     * 
     * 
     * try
     */
    mlfTry {
        /*
         * z = (x-mu) ./ sigma;
         */
        mlfAssign(
          &z,
          mclRdivide(
            mclMinus(mclVa(x, "x"), mclVa(mu, "mu")), mclVa(sigma, "sigma")));
    /*
     * catch
     */
    } mlfCatch {
        /*
         * error('Non-scalar arguments must match in size.');
         */
        mlfError(_mxarray12_, NULL);
    /*
     * end
     */
    } mlfEndCatch
    /*
     * 
     * 
     * % Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
     * % to produce accurate near-zero results for large negative x.
     * p = 0.5 * erfc(-z ./ sqrt(2));
     */
    mlfAssign(
      &p,
      mclMtimes(
        _mxarray14_,
        mlfErfc(mclRdivide(mclUminus(mclVv(z, "z")), mlfSqrt(_mxarray15_)))));
    /*
     * 
     * % Compute confidence bounds if requested.
     * if nargout>=2
     */
    if (nargout_ >= 2) {
        /*
         * zvar = (pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2) ./ (sigma.^2);;
         */
        mlfAssign(
          &zvar,
          mclRdivide(
            mclPlus(
              mclPlus(
                mclIntArrayRef2(mclVa(pcov, "pcov"), 1, 1),
                mclMtimes(
                  mclMtimes(
                    _mxarray15_, mclIntArrayRef2(mclVa(pcov, "pcov"), 1, 2)),
                  mclVv(z, "z"))),
              mclMtimes(
                mclIntArrayRef2(mclVa(pcov, "pcov"), 2, 2),
                mlfPower(mclVv(z, "z"), _mxarray15_))),
            mlfPower(mclVa(sigma, "sigma"), _mxarray15_)));
        /*
         * if any(zvar<0)
         */
        if (mlfTobool(mlfAny(mclLt(mclVv(zvar, "zvar"), _mxarray0_), NULL))) {
            /*
             * error('PCOV must be a positive semi-definite matrix.');
             */
            mlfError(_mxarray16_, NULL);
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
              mclMrdivide(mclVa(alpha, "alpha"), _mxarray15_),
              NULL,
              NULL,
              NULL,
              NULL)));
        /*
         * halfwidth = normz * sqrt(zvar);
         */
        mlfAssign(
          &halfwidth,
          mclMtimes(mclVv(normz, "normz"), mlfSqrt(mclVv(zvar, "zvar"))));
        /*
         * zlo = z - halfwidth;
         */
        mlfAssign(&zlo, mclMinus(mclVv(z, "z"), mclVv(halfwidth, "halfwidth")));
        /*
         * zup = z + halfwidth;
         */
        mlfAssign(&zup, mclPlus(mclVv(z, "z"), mclVv(halfwidth, "halfwidth")));
        /*
         * 
         * plo = 0.5 * erfc(-zlo./sqrt(2));
         */
        mlfAssign(
          plo,
          mclMtimes(
            _mxarray14_,
            mlfErfc(
              mclRdivide(mclUminus(mclVv(zlo, "zlo")), mlfSqrt(_mxarray15_)))));
        /*
         * pup = 0.5 * erfc(-zup./sqrt(2));
         */
        mlfAssign(
          pup,
          mclMtimes(
            _mxarray14_,
            mlfErfc(
              mclRdivide(mclUminus(mclVv(zup, "zup")), mlfSqrt(_mxarray15_)))));
    /*
     * end
     */
    }
    mclValidateOutput(p, 1, nargout_, "p", "normcdf");
    mclValidateOutput(*plo, 2, nargout_, "plo", "normcdf");
    mclValidateOutput(*pup, 3, nargout_, "pup", "normcdf");
    mxDestroyArray(ans);
    mxDestroyArray(z);
    mxDestroyArray(zvar);
    mxDestroyArray(normz);
    mxDestroyArray(halfwidth);
    mxDestroyArray(zlo);
    mxDestroyArray(zup);
    mxDestroyArray(alpha);
    mxDestroyArray(pcov);
    mxDestroyArray(sigma);
    mxDestroyArray(mu);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return p;
}
