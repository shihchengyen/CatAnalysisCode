/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "jbtest.h"
#include "chi2cdf.h"
#include "chi2inv.h"
#include "libmatlbm.h"
#include "libmmfile.h"
static mxArray * _mxarray0_;

static mxChar _array2_[35] = { ' ', 'I', 'n', 'p', 'u', 't', ' ', 's', 'a',
                               'm', 'p', 'l', 'e', ' ', 0x0027, 'X', 0x0027,
                               ' ', 'm', 'u', 's', 't', ' ', 'b', 'e', ' ',
                               'a', ' ', 'v', 'e', 'c', 't', 'o', 'r', '.' };
static mxArray * _mxarray1_;

static mxChar _array4_[48] = { ' ', 'I', 'n', 'p', 'u', 't', ' ', 's', 'a',
                               'm', 'p', 'l', 'e', ' ', 0x0027, 'X',
                               0x0027, ' ', 'h', 'a', 's', ' ', 'n', 'o',
                               ' ', 'v', 'a', 'l', 'i', 'd', ' ', 'd', 'a',
                               't', 'a', ' ', '(', 'a', 'l', 'l', ' ', 'N',
                               'a', 'N', 0x0027, 's', ')', '.' };
static mxArray * _mxarray3_;

static mxChar _array6_[45] = { ' ', 'S', 'i', 'g', 'n', 'i', 'f', 'i', 'c', 'a',
                               'n', 'c', 'e', ' ', 'l', 'e', 'v', 'e', 'l', ' ',
                               0x0027, 'A', 'l', 'p', 'h', 'a', 0x0027, ' ',
                               'm', 'u', 's', 't', ' ', 'b', 'e', ' ', 'a', ' ',
                               's', 'c', 'a', 'l', 'a', 'r', '.' };
static mxArray * _mxarray5_;

static mxChar _array8_[52] = { ' ', 'S', 'i', 'g', 'n', 'i', 'f', 'i', 'c',
                               'a', 'n', 'c', 'e', ' ', 'l', 'e', 'v', 'e',
                               'l', ' ', 0x0027, 'A', 'l', 'p', 'h', 'a',
                               0x0027, ' ', 'm', 'u', 's', 't', ' ', 'b',
                               'e', ' ', 'b', 'e', 't', 'w', 'e', 'e', 'n',
                               ' ', '0', ' ', 'a', 'n', 'd', ' ', '1', '.' };
static mxArray * _mxarray7_;
static mxArray * _mxarray9_;
static mxArray * _mxarray10_;
static mxArray * _mxarray11_;
static mxArray * _mxarray12_;
static mxArray * _mxarray13_;
static mxArray * _mxarray14_;
static mxArray * _mxarray15_;

void InitializeModule_jbtest(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeString(35, _array2_);
    _mxarray3_ = mclInitializeString(48, _array4_);
    _mxarray5_ = mclInitializeString(45, _array6_);
    _mxarray7_ = mclInitializeString(52, _array8_);
    _mxarray9_ = mclInitializeDouble(0.0);
    _mxarray10_ = mclInitializeDouble(.05);
    _mxarray11_ = mclInitializeDouble(3.0);
    _mxarray12_ = mclInitializeDouble(4.0);
    _mxarray13_ = mclInitializeDouble(2.0);
    _mxarray14_ = mclInitializeDouble(6.0);
    _mxarray15_ = mclInitializeDouble(24.0);
}

void TerminateModule_jbtest(void) {
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mjbtest(mxArray * * pValue,
                         mxArray * * JBstatistic,
                         mxArray * * criticalValue,
                         int nargout_,
                         mxArray * x,
                         mxArray * alpha);

_mexLocalFunctionTable _local_function_table_jbtest
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfNJbtest" contains the nargout interface for the "jbtest"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/jbtest.m"
 * (lines 1-128). This interface is only produced if the M-function uses the
 * special variable "nargout". The nargout interface allows the number of
 * requested outputs to be specified via the nargout argument, as opposed to
 * the normal interface which dynamically calculates the number of outputs
 * based on the number of non-NULL inputs it receives. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
mxArray * mlfNJbtest(int nargout,
                     mxArray * * pValue,
                     mxArray * * JBstatistic,
                     mxArray * * criticalValue,
                     mxArray * x,
                     mxArray * alpha) {
    mxArray * H = NULL;
    mxArray * pValue__ = NULL;
    mxArray * JBstatistic__ = NULL;
    mxArray * criticalValue__ = NULL;
    mlfEnterNewContext(3, 2, pValue, JBstatistic, criticalValue, x, alpha);
    H = Mjbtest(&pValue__, &JBstatistic__, &criticalValue__, nargout, x, alpha);
    mlfRestorePreviousContext(
      3, 2, pValue, JBstatistic, criticalValue, x, alpha);
    if (pValue != NULL) {
        mclCopyOutputArg(pValue, pValue__);
    } else {
        mxDestroyArray(pValue__);
    }
    if (JBstatistic != NULL) {
        mclCopyOutputArg(JBstatistic, JBstatistic__);
    } else {
        mxDestroyArray(JBstatistic__);
    }
    if (criticalValue != NULL) {
        mclCopyOutputArg(criticalValue, criticalValue__);
    } else {
        mxDestroyArray(criticalValue__);
    }
    return mlfReturnValue(H);
}

/*
 * The function "mlfJbtest" contains the normal interface for the "jbtest"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/jbtest.m"
 * (lines 1-128). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
mxArray * mlfJbtest(mxArray * * pValue,
                    mxArray * * JBstatistic,
                    mxArray * * criticalValue,
                    mxArray * x,
                    mxArray * alpha) {
    int nargout = 1;
    mxArray * H = NULL;
    mxArray * pValue__ = NULL;
    mxArray * JBstatistic__ = NULL;
    mxArray * criticalValue__ = NULL;
    mlfEnterNewContext(3, 2, pValue, JBstatistic, criticalValue, x, alpha);
    if (pValue != NULL) {
        ++nargout;
    }
    if (JBstatistic != NULL) {
        ++nargout;
    }
    if (criticalValue != NULL) {
        ++nargout;
    }
    H = Mjbtest(&pValue__, &JBstatistic__, &criticalValue__, nargout, x, alpha);
    mlfRestorePreviousContext(
      3, 2, pValue, JBstatistic, criticalValue, x, alpha);
    if (pValue != NULL) {
        mclCopyOutputArg(pValue, pValue__);
    } else {
        mxDestroyArray(pValue__);
    }
    if (JBstatistic != NULL) {
        mclCopyOutputArg(JBstatistic, JBstatistic__);
    } else {
        mxDestroyArray(JBstatistic__);
    }
    if (criticalValue != NULL) {
        mclCopyOutputArg(criticalValue, criticalValue__);
    } else {
        mxDestroyArray(criticalValue__);
    }
    return mlfReturnValue(H);
}

/*
 * The function "mlfVJbtest" contains the void interface for the "jbtest"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/jbtest.m"
 * (lines 1-128). The void interface is only produced if the M-function uses
 * the special variable "nargout", and has at least one output. The void
 * interface function specifies zero output arguments to the implementation
 * version of the function, and in the event that the implementation version
 * still returns an output (which, in MATLAB, would be assigned to the "ans"
 * variable), it deallocates the output. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
void mlfVJbtest(mxArray * x, mxArray * alpha) {
    mxArray * H = NULL;
    mxArray * pValue = NULL;
    mxArray * JBstatistic = NULL;
    mxArray * criticalValue = NULL;
    mlfEnterNewContext(0, 2, x, alpha);
    H = Mjbtest(&pValue, &JBstatistic, &criticalValue, 0, x, alpha);
    mlfRestorePreviousContext(0, 2, x, alpha);
    mxDestroyArray(H);
    mxDestroyArray(pValue);
    mxDestroyArray(JBstatistic);
    mxDestroyArray(criticalValue);
}

/*
 * The function "mlxJbtest" contains the feval interface for the "jbtest"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/jbtest.m"
 * (lines 1-128). The feval function calls the implementation version of jbtest
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxJbtest(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[4];
    int i;
    if (nlhs > 4) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: jbtest Line: 1 Column: "
            "1 The function \"jbtest\" was called with mor"
            "e than the declared number of outputs (4)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: jbtest Line: 1 Column: "
            "1 The function \"jbtest\" was called with mor"
            "e than the declared number of inputs (2)."),
          NULL);
    }
    for (i = 0; i < 4; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0]
      = Mjbtest(&mplhs[1], &mplhs[2], &mplhs[3], nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 4 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 4; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mjbtest" is the implementation version of the "jbtest"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/jbtest.m"
 * (lines 1-128). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [H, pValue, JBstatistic, criticalValue] = jbtest(x , alpha)
 */
static mxArray * Mjbtest(mxArray * * pValue,
                         mxArray * * JBstatistic,
                         mxArray * * criticalValue,
                         int nargout_,
                         mxArray * x,
                         mxArray * alpha) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_jbtest);
    int nargin_ = mclNargin(2, x, alpha, NULL);
    mxArray * H = NULL;
    mxArray * B = NULL;
    mxArray * J = NULL;
    mxArray * n = NULL;
    mxArray * ans = NULL;
    mxArray * columns = NULL;
    mxArray * rows = NULL;
    mclCopyArray(&x);
    mclCopyArray(&alpha);
    /*
     * %JBTEST Bera-Jarque parametric hypothesis test of composite normality.
     * %   H = JBTEST(X,ALPHA) performs the Bera-Jarque test to determine if the null
     * %   hypothesis of composite normality is a reasonable assumption regarding the
     * %   population distribution of a random sample X. The desired significance 
     * %   level, ALPHA, is an optional scalar input (default = 0.05).
     * %
     * %   H indicates the Boolean result of the hypothesis test:
     * %     H = 0 => Do not reject the null hypothesis at significance level ALPHA.
     * %     H = 1 => Reject the null hypothesis at significance level ALPHA.
     * % 
     * %   The Bera-Jarque hypotheses are: 
     * %   Null Hypothesis:        X is normal with unspecified mean and variance.
     * %   Alternative Hypothesis: X is not normally distributed.
     * %
     * %   The Bera-Jarque test is a 2-sided test of composite normality with sample
     * %   mean and sample variance used as estimates of the population mean and 
     * %   variance, respectively. The test statistic is based on estimates of the 
     * %   sample skewness and kurtosis of the normalized data (the standardized 
     * %   Z-scores computed from X by subtracting the sample mean and normalizing by
     * %   the sample standard deviation). Under the null hypothesis, the standardized
     * %   3rd and 4th moments are asymptotically normal and independent, and the test
     * %   statistic has a Chi-square distribution with two degrees of freedom. Note 
     * %   that the Bera-Jarque test is an asymptotic test, and care should be taken
     * %   with small sample sizes.
     * %
     * %   X is a vector representing a random sample from some underlying 
     * %   distribution. Missing observations in X, indicated by NaN's (Not-a-Number),
     * %   are ignored.
     * %
     * %   [H,P] = JBTEST(...) also returns the P-value (significance level) at 
     * %   which the null hypothesis is rejected.
     * %
     * %   [H,P,JBSTAT] = JBTEST(...) also returns the test statistic JBSTAT.
     * %
     * %   [H,P,JBSTAT,CV] = JBTEST(...) also returns the critical value of the 
     * %   test CV.
     * %
     * %   See also LILLIETEST, KSTEST, KSTEST2.
     * %
     * 
     * % Author(s): R.A. Baker, 08/27/98
     * % Copyright 1993-2002 The MathWorks, Inc. 
     * % $Revision: 1.4 $   $ Date: 1998/01/30 13:45:34 $
     * 
     * %
     * % References:
     * %   Judge, G.G., Hill, R.C., et al, "Introduction to Theory & Practice of
     * %     Econometrics", 2nd edition, John Wiley & Sons, Inc., 1988, pages 890-892.
     * %
     * %   Spanos, A., "Statistical Foundations of Econometric Modelling",
     * %     Cambridge University Press, 1993, pages 453-455.
     * %
     * 
     * %
     * % Ensure the sample data is a VECTOR.
     * %
     * 
     * [rows , columns]  =  size(x);
     */
    mlfSize(mlfVarargout(&rows, &columns, NULL), mclVa(x, "x"), NULL);
    /*
     * 
     * if (rows ~= 1) & (columns ~= 1) 
     */
    {
        mxArray * a_ = mclInitialize(mclNe(mclVv(rows, "rows"), _mxarray0_));
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(a_, mclNe(mclVv(columns, "columns"), _mxarray0_)))) {
            mxDestroyArray(a_);
            /*
             * error(' Input sample ''X'' must be a vector.');
             */
            mlfError(_mxarray1_, NULL);
        } else {
            mxDestroyArray(a_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Remove missing observations indicated by NaN's.
     * %
     * 
     * x  =  x(~isnan(x));
     */
    mlfAssign(&x, mclArrayRef1(mclVa(x, "x"), mclNot(mlfIsnan(mclVa(x, "x")))));
    /*
     * 
     * if length(x) == 0
     */
    if (mclLengthInt(mclVa(x, "x")) == 0) {
        /*
         * error(' Input sample ''X'' has no valid data (all NaN''s).');
         */
        mlfError(_mxarray3_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * x  =  x(:);               % Ensure a column vector.
     */
    mlfAssign(&x, mclArrayRef1(mclVa(x, "x"), mlfCreateColonIndex()));
    /*
     * 
     * %
     * % Ensure the significance level, ALPHA, is a 
     * % scalar, and set default if necessary.
     * %
     * 
     * if (nargin >= 2) & ~isempty(alpha)
     */
    {
        mxArray * a_ = mclInitialize(mclBoolToArray(nargin_ >= 2));
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(a_, mclNot(mlfIsempty(mclVa(alpha, "alpha")))))) {
            mxDestroyArray(a_);
            /*
             * if prod(size(alpha)) > 1
             */
            if (mclGtBool(
                  mlfProd(
                    mlfSize(mclValueVarargout(), mclVa(alpha, "alpha"), NULL),
                    NULL),
                  _mxarray0_)) {
                /*
                 * error(' Significance level ''Alpha'' must be a scalar.');
                 */
                mlfError(_mxarray5_, NULL);
            /*
             * end
             */
            }
            /*
             * if (alpha <= 0 | alpha >= 1)
             */
            {
                mxArray * a_0
                  = mclInitialize(mclLe(mclVa(alpha, "alpha"), _mxarray9_));
                if (mlfTobool(a_0)
                    || mlfTobool(
                         mclOr(
                           a_0, mclGe(mclVa(alpha, "alpha"), _mxarray0_)))) {
                    mxDestroyArray(a_0);
                    /*
                     * error(' Significance level ''Alpha'' must be between 0 and 1.'); 
                     */
                    mlfError(_mxarray7_, NULL);
                } else {
                    mxDestroyArray(a_0);
                }
            /*
             * end
             */
            }
        /*
         * else
         */
        } else {
            mxDestroyArray(a_);
            /*
             * alpha  =  0.05;
             */
            mlfAssign(&alpha, _mxarray10_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Compute moments and test statistic.
     * %
     * 
     * n  =  length(x);              % Sample size.
     */
    mlfAssign(&n, mlfScalar(mclLengthInt(mclVa(x, "x"))));
    /*
     * x  =  (x - mean(x)) / std(x); % Standardized data (i.e., Z-scores).
     */
    mlfAssign(
      &x,
      mclMrdivide(
        mclMinus(mclVa(x, "x"), mlfMean(mclVa(x, "x"), NULL)),
        mlfStd(mclVa(x, "x"), NULL, NULL)));
    /*
     * J  =  sum(x.^3) / n;          % Moment 3.
     */
    mlfAssign(
      &J,
      mclMrdivide(
        mlfSum(mlfPower(mclVa(x, "x"), _mxarray11_), NULL), mclVv(n, "n")));
    /*
     * B  =  (sum(x.^4) / n) - 3;    % Moment 4 (excess kurtosis).
     */
    mlfAssign(
      &B,
      mclMinus(
        mclMrdivide(
          mlfSum(mlfPower(mclVa(x, "x"), _mxarray12_), NULL), mclVv(n, "n")),
        _mxarray11_));
    /*
     * 
     * JBstatistic =  n * ( (J^2)/6  +  (B^2)/24 );
     */
    mlfAssign(
      JBstatistic,
      mclMtimes(
        mclVv(n, "n"),
        mclPlus(
          mclMrdivide(mclMpower(mclVv(J, "J"), _mxarray13_), _mxarray14_),
          mclMrdivide(mclMpower(mclVv(B, "B"), _mxarray13_), _mxarray15_))));
    /*
     * 
     * %
     * % Compute the P-values (i.e., significance levels). Since the 
     * % CHI2INV function is a slow, iterative procedure, compute the 
     * % critical values ONLY if requested. Under the null hypothesis 
     * % of composite normality, the test statistic of the standardized 
     * % data is asymptotically Chi-Square distributed with 2 degrees 
     * % of freedom.
     * %
     * 
     * if nargout >= 4
     */
    if (nargout_ >= 4) {
        /*
         * criticalValue  =  chi2inv(1 - alpha , 2);
         */
        mlfAssign(
          criticalValue,
          mlfChi2inv(mclMinus(_mxarray0_, mclVa(alpha, "alpha")), _mxarray13_));
    /*
     * end
     */
    }
    /*
     * 
     * pValue =  1 - chi2cdf(JBstatistic , 2);
     */
    mlfAssign(
      pValue,
      mclMinus(
        _mxarray0_,
        mlfChi2cdf(mclVv(*JBstatistic, "JBstatistic"), _mxarray13_)));
    /*
     * 
     * %
     * % To maintain consistency with existing Statistics Toolbox hypothesis
     * % tests, returning 'H = 0' implies that we 'Do not reject the null 
     * % hypothesis at the significance level of alpha' and 'H = 1' implies 
     * % that we 'Reject the null hypothesis at significance level of alpha.'
     * %
     * 
     * H  = (alpha >= pValue);
     */
    mlfAssign(&H, mclGe(mclVa(alpha, "alpha"), mclVv(*pValue, "pValue")));
    mclValidateOutput(H, 1, nargout_, "H", "jbtest");
    mclValidateOutput(*pValue, 2, nargout_, "pValue", "jbtest");
    mclValidateOutput(*JBstatistic, 3, nargout_, "JBstatistic", "jbtest");
    mclValidateOutput(*criticalValue, 4, nargout_, "criticalValue", "jbtest");
    mxDestroyArray(rows);
    mxDestroyArray(columns);
    mxDestroyArray(ans);
    mxDestroyArray(n);
    mxDestroyArray(J);
    mxDestroyArray(B);
    mxDestroyArray(alpha);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return H;
    /*
     * 
     */
}
