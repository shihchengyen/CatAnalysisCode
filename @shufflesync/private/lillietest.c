/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "lillietest.h"
#include "cdfcalc.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "normcdf.h"
static mxArray * _mxarray0_;

static mxChar _array2_[35] = { ' ', 'I', 'n', 'p', 'u', 't', ' ', 's', 'a',
                               'm', 'p', 'l', 'e', ' ', 0x0027, 'X', 0x0027,
                               ' ', 'm', 'u', 's', 't', ' ', 'b', 'e', ' ',
                               'a', ' ', 'v', 'e', 'c', 't', 'o', 'r', '.' };
static mxArray * _mxarray1_;

static mxChar _array4_[59] = { ' ', 'S', 'a', 'm', 'p', 'l', 'e', ' ', 'v',
                               'e', 'c', 't', 'o', 'r', ' ', 0x0027, 'X',
                               0x0027, ' ', 'm', 'u', 's', 't', ' ', 'h',
                               'a', 'v', 'e', ' ', 'a', 't', ' ', 'l', 'e',
                               'a', 's', 't', ' ', '4', ' ', 'v', 'a', 'l',
                               'i', 'd', ' ', 'o', 'b', 's', 'e', 'r', 'v',
                               'a', 't', 'i', 'o', 'n', 's', '.' };
static mxArray * _mxarray3_;

static mxChar _array6_[45] = { ' ', 'S', 'i', 'g', 'n', 'i', 'f', 'i', 'c', 'a',
                               'n', 'c', 'e', ' ', 'l', 'e', 'v', 'e', 'l', ' ',
                               0x0027, 'A', 'l', 'p', 'h', 'a', 0x0027, ' ',
                               'm', 'u', 's', 't', ' ', 'b', 'e', ' ', 'a', ' ',
                               's', 'c', 'a', 'l', 'a', 'r', '.' };
static mxArray * _mxarray5_;
static mxArray * _mxarray7_;

static mxChar _array9_[47] = { ' ', 'S', 'i', 'g', 'n', 'i', 'f', 'i', 'c', 'a',
                               'n', 'c', 'e', ' ', 'l', 'e', 'v', 'e', 'l', ' ',
                               0x0027, 'A', 'l', 'p', 'h', 'a', 0x0027, ' ',
                               'm', 'u', 's', 't', ' ', 'i', 'n', ' ', '[', '0',
                               '.', '0', '1', ',', '0', '.', '2', ']', '.' };
static mxArray * _mxarray8_;
static mxArray * _mxarray10_;
static mxArray * _mxarray11_;

static double _array13_[5] = { .01, .05, .1, .15, .2 };
static mxArray * _mxarray12_;
static mxArray * _mxarray14_;

static double _array16_[19] = { 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                                11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
                                17.0, 18.0, 19.0, 20.0, 25.0, 30.0 };
static mxArray * _mxarray15_;

static double _array18_[95] = { .417, .405, .364, .348, .331, .311, .294, .284,
                                .275, .268, .261, .257, .25, .245, .239, .235,
                                .231, .2, .187, .381, .337, .319, .3, .285,
                                .271, .258, .249, .242, .234, .227, .22, .213,
                                .206, .2, .195, .19, .173, .161, .352, .315,
                                .294, .276, .261, .249, .239, .23, .223, .214,
                                .207, .201, .195, .189, .184, .179, .174, .158,
                                .144, .319, .299, .277, .258, .244, .233, .224,
                                .217, .212, .202, .194, .187, .182, .177, .173,
                                .169, .166, .147, .136, .3, .285, .265, .247,
                                .233, .223, .215, .206, .199, .19, .183, .177,
                                .173, .169, .166, .163, .16, .142, .131 };
static mxArray * _mxarray17_;

static mxChar _array20_[6] = { 'l', 'i', 'n', 'e', 'a', 'r' };
static mxArray * _mxarray19_;

static double _array22_[5] = { 1.1031, .886, .805, .768, .736 };
static mxArray * _mxarray21_;
static mxArray * _mxarray23_;
static mxArray * _mxarray24_;

void InitializeModule_lillietest(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeString(35, _array2_);
    _mxarray3_ = mclInitializeString(59, _array4_);
    _mxarray5_ = mclInitializeString(45, _array6_);
    _mxarray7_ = mclInitializeDouble(.2);
    _mxarray8_ = mclInitializeString(47, _array9_);
    _mxarray10_ = mclInitializeDouble(.01);
    _mxarray11_ = mclInitializeDouble(.05);
    _mxarray12_ = mclInitializeDoubleVector(5, 1, _array13_);
    _mxarray14_ = mclInitializeDouble(30.0);
    _mxarray15_ = mclInitializeDoubleVector(19, 1, _array16_);
    _mxarray17_ = mclInitializeDoubleVector(19, 5, _array18_);
    _mxarray19_ = mclInitializeString(6, _array20_);
    _mxarray21_ = mclInitializeDoubleVector(1, 5, _array22_);
    _mxarray23_ = mclInitializeDouble(0.0);
    _mxarray24_ = mclInitializeDouble(2.0);
}

void TerminateModule_lillietest(void) {
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray23_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mlillietest(mxArray * * pValue,
                             mxArray * * KSstatistic,
                             mxArray * * criticalValue,
                             int nargout_,
                             mxArray * x,
                             mxArray * alpha);

_mexLocalFunctionTable _local_function_table_lillietest
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfLillietest" contains the normal interface for the
 * "lillietest" M-function from file
 * "/Applications/MATLAB6p5p1/toolbox/stats/lillietest.m" (lines 1-192). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
mxArray * mlfLillietest(mxArray * * pValue,
                        mxArray * * KSstatistic,
                        mxArray * * criticalValue,
                        mxArray * x,
                        mxArray * alpha) {
    int nargout = 1;
    mxArray * H = NULL;
    mxArray * pValue__ = NULL;
    mxArray * KSstatistic__ = NULL;
    mxArray * criticalValue__ = NULL;
    mlfEnterNewContext(3, 2, pValue, KSstatistic, criticalValue, x, alpha);
    if (pValue != NULL) {
        ++nargout;
    }
    if (KSstatistic != NULL) {
        ++nargout;
    }
    if (criticalValue != NULL) {
        ++nargout;
    }
    H
      = Mlillietest(
          &pValue__, &KSstatistic__, &criticalValue__, nargout, x, alpha);
    mlfRestorePreviousContext(
      3, 2, pValue, KSstatistic, criticalValue, x, alpha);
    if (pValue != NULL) {
        mclCopyOutputArg(pValue, pValue__);
    } else {
        mxDestroyArray(pValue__);
    }
    if (KSstatistic != NULL) {
        mclCopyOutputArg(KSstatistic, KSstatistic__);
    } else {
        mxDestroyArray(KSstatistic__);
    }
    if (criticalValue != NULL) {
        mclCopyOutputArg(criticalValue, criticalValue__);
    } else {
        mxDestroyArray(criticalValue__);
    }
    return mlfReturnValue(H);
}

/*
 * The function "mlxLillietest" contains the feval interface for the
 * "lillietest" M-function from file
 * "/Applications/MATLAB6p5p1/toolbox/stats/lillietest.m" (lines 1-192). The
 * feval function calls the implementation version of lillietest through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxLillietest(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[4];
    int i;
    if (nlhs > 4) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: lillietest Line: 1 Column:"
            " 1 The function \"lillietest\" was called with m"
            "ore than the declared number of outputs (4)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: lillietest Line: 1 Column"
            ": 1 The function \"lillietest\" was called with"
            " more than the declared number of inputs (2)."),
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
      = Mlillietest(&mplhs[1], &mplhs[2], &mplhs[3], nlhs, mprhs[0], mprhs[1]);
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
 * The function "Mlillietest" is the implementation version of the "lillietest"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/lillietest.m"
 * (lines 1-192). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [H,pValue,KSstatistic,criticalValue] = lillietest(x , alpha)
 */
static mxArray * Mlillietest(mxArray * * pValue,
                             mxArray * * KSstatistic,
                             mxArray * * criticalValue,
                             int nargout_,
                             mxArray * x,
                             mxArray * alpha) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_lillietest);
    int nargin_ = mclNargin(2, x, alpha, NULL);
    mxArray * H = NULL;
    mxArray * deltaCDF = NULL;
    mxArray * delta2 = NULL;
    mxArray * delta1 = NULL;
    mxArray * nullCDF = NULL;
    mxArray * zScores = NULL;
    mxArray * Q = NULL;
    mxArray * quantiles = NULL;
    mxArray * sampleSize = NULL;
    mxArray * a = NULL;
    mxArray * emsg = NULL;
    mxArray * n = NULL;
    mxArray * xCDF = NULL;
    mxArray * sampleCDF = NULL;
    mxArray * ans = NULL;
    mxArray * columns = NULL;
    mxArray * rows = NULL;
    mclCopyArray(&x);
    mclCopyArray(&alpha);
    /*
     * %LILLIETEST Single sample Lilliefors hypothesis test of composite normality.
     * %   H = LILLIETEST(X,ALPHA) performs the Lilliefors modification of the 
     * %   Kolmogorov-Smirnov test to determine if the null hypothesis of composite
     * %   normality is a reasonable assumption regarding the population distribution
     * %   of a random sample X. The desired significance level, ALPHA, is an optional
     * %   scalar input (default = 0.05). The Lilliefors test is based on simulation, 
     * %   so the significance level is restricted to 0.01 <= ALPHA <= 0.20 (the 
     * %   region tabularized by Lilliefors).
     * %
     * %   H indicates the result of the hypothesis test:
     * %     H = 0 => Do not reject the null hypothesis at significance level ALPHA.
     * %     H = 1 => Reject the null hypothesis at significance level ALPHA.
     * % 
     * %   Let S(x) be the empirical c.d.f. estimated from the sample vector X, F(x) 
     * %   be the corresponding true (but unknown) population c.d.f., and CDF be a
     * %   normal c.d.f. with sample mean and standard deviation taken from X. The 
     * %   Lilliefors hypotheses and test statistic are:
     * % 
     * %   Null Hypothesis:        F(x) is normal with unspecified mean and variance.
     * %   Alternative Hypothesis: F(x) is not normally distributed.
     * %
     * %   Test Statistic:         T = max|S(x) - CDF|.
     * %
     * %   The decision to reject the null hypothesis occurs when the test statistic
     * %   exceeds the critical value. 
     * %
     * %   The Lilliefors test is a 2-sided test of composite normality with sample 
     * %   mean and sample variance used as estimates of the population mean and 
     * %   variance, respectively. The test statistic is based on the 'normalized'
     * %   samples (i.e., the Z-scores computed by subtracting the sample mean and
     * %   normalizing by the sample standard deviation).
     * %
     * %   X may be a row or column vector representing a random sample from some
     * %   underlying distribution. Missing observations in X, indicated by NaN's 
     * %   (Not-a-Number), are ignored.
     * %
     * %   [H,P] = LILLIETEST(...) also returns the approximate P-value P, computed
     * %   via interpolation into the Lilliefors simulation table. A NaN is 
     * %   returned when P is not found within the interval 0.01 <= P <= 0.20.
     * %
     * %   [H,P,LSTAT] = LILLIETEST(...) also returns the test statistic LSTAT, the 
     * %   largest absolute vertical deviation of S(x) from CDF.
     * %
     * %   [H,P,LSTAT,CV] = LILLIETEST(...) also returns the critical value of the 
     * %   test CV.
     * %
     * %   See also KSTEST, KSTEST2, CDFPLOT.
     * %
     * 
     * % Author(s): R.A. Baker, 08/27/98
     * % Copyright 1993-2002 The MathWorks, Inc. 
     * % $Revision: 1.4 $   $ Date: 1998/01/30 13:45:34 $
     * 
     * 
     * % References:
     * %   Conover, W.J., "Practical Nonparametric Statistics", 
     * %     John Wiley & Sons, Inc., 1980.
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
     * % Remove missing observations indicated by NaN's, and ensure that
     * % at least 4 valid observations remain.
     * %
     * 
     * x  =  x(~isnan(x));       % Remove missing observations indicated by NaN's.
     */
    mlfAssign(&x, mclArrayRef1(mclVa(x, "x"), mclNot(mlfIsnan(mclVa(x, "x")))));
    /*
     * x  =  x(:);               % Ensure a column vector.
     */
    mlfAssign(&x, mclArrayRef1(mclVa(x, "x"), mlfCreateColonIndex()));
    /*
     * 
     * if length(x) < 4
     */
    if (mclLengthInt(mclVa(x, "x")) < 4) {
        /*
         * error(' Sample vector ''X'' must have at least 4 valid observations.');
         */
        mlfError(_mxarray3_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Ensure the significance level, ALPHA, is a scalar between the significance 
     * % levels of Lilliefors' simulation, and set default if necessary.
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
             * if (alpha < 0.01) | (alpha > 0.20)
             */
            {
                mxArray * a_0
                  = mclInitialize(mclLt(mclVa(alpha, "alpha"), _mxarray10_));
                if (mlfTobool(a_0)
                    || mlfTobool(
                         mclOr(
                           a_0, mclGt(mclVa(alpha, "alpha"), _mxarray7_)))) {
                    mxDestroyArray(a_0);
                    /*
                     * error(' Significance level ''Alpha'' must in [0.01,0.2].');
                     */
                    mlfError(_mxarray8_, NULL);
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
            mlfAssign(&alpha, _mxarray11_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Calculate S(x), the sample CDF.
     * %
     * 
     * [sampleCDF,xCDF,n,emsg] = cdfcalc(x);
     */
    mlfAssign(&sampleCDF, mlfCdfcalc(&xCDF, &n, &emsg, mclVa(x, "x"), NULL));
    /*
     * error(emsg);
     */
    mlfError(mclVv(emsg, "emsg"), NULL);
    /*
     * 
     * %
     * % Compute the critical value upon which to 
     * % base the acceptance or rejection of the test.
     * %
     * 
     * a  =  [0.010  0.050  0.100  0.150  0.200]';  % Lilliefors' Alphas.
     */
    mlfAssign(&a, _mxarray12_);
    /*
     * 
     * if n <= 30
     */
    if (mclLeBool(mclVv(n, "n"), _mxarray14_)) {
        /*
         * 
         * sampleSize =  [4:20 25 30]'; % Sample sizes for each row of 'quantiles'.
         */
        mlfAssign(&sampleSize, _mxarray15_);
        /*
         * 
         * quantiles  =  [0.417  0.381  0.352  0.319  0.300
         */
        mlfAssign(&quantiles, _mxarray17_);
        /*
         * 0.405  0.337  0.315  0.299  0.285
         * 0.364  0.319  0.294  0.277  0.265
         * 0.348  0.300  0.276  0.258  0.247
         * 0.331  0.285  0.261  0.244  0.233
         * 0.311  0.271  0.249  0.233  0.223
         * 0.294  0.258  0.239  0.224  0.215
         * 0.284  0.249  0.230  0.217  0.206
         * 0.275  0.242  0.223  0.212  0.199
         * 0.268  0.234  0.214  0.202  0.190
         * 0.261  0.227  0.207  0.194  0.183
         * 0.257  0.220  0.201  0.187  0.177
         * 0.250  0.213  0.195  0.182  0.173
         * 0.245  0.206  0.189  0.177  0.169
         * 0.239  0.200  0.184  0.173  0.166
         * 0.235  0.195  0.179  0.169  0.163
         * 0.231  0.190  0.174  0.166  0.160
         * 0.200  0.173  0.158  0.147  0.142
         * 0.187  0.161  0.144  0.136  0.131];
         * %
         * %  Determine the interpolated 'row' of quantiles.
         * %
         * Q  =  interp2(a , sampleSize, quantiles , a , n , 'linear');
         */
        mlfAssign(
          &Q,
          mlfInterp2(
            mclVv(a, "a"),
            mclVv(sampleSize, "sampleSize"),
            mclVv(quantiles, "quantiles"),
            mclVv(a, "a"),
            mclVv(n, "n"),
            _mxarray19_,
            NULL));
        /*
         * 
         * %
         * %  Now compute the critical value.
         * %
         * criticalValue  =  interp1(a , Q , alpha , 'linear');
         */
        mlfAssign(
          criticalValue,
          mlfInterp1(
            mclVv(a, "a"),
            mclVv(Q, "Q"),
            mclVa(alpha, "alpha"),
            _mxarray19_,
            NULL));
    /*
     * 
     * else
     */
    } else {
        /*
         * 
         * %
         * %  Determine the asymptotic 'row' of quantiles.
         * %
         * Q  =  [1.1031 0.886 0.805 0.768 0.736] / sqrt(n);
         */
        mlfAssign(&Q, mclMrdivide(_mxarray21_, mlfSqrt(mclVv(n, "n"))));
        /*
         * 
         * %
         * %  1-D interpolation into Lilliefors' asymptotic expression.
         * %
         * criticalValue  =  interp1(a , Q , alpha , 'linear');
         */
        mlfAssign(
          criticalValue,
          mlfInterp1(
            mclVv(a, "a"),
            mclVv(Q, "Q"),
            mclVa(alpha, "alpha"),
            _mxarray19_,
            NULL));
    /*
     * 
     * end
     */
    }
    /*
     * 
     * %
     * % The theoretical CDF specified under the null hypothesis 
     * % is assumed to be normal (i.e., Gaussian).
     * %
     * 
     * zScores  =  (xCDF - mean(x))./std(x);
     */
    mlfAssign(
      &zScores,
      mclRdivide(
        mclMinus(mclVv(xCDF, "xCDF"), mlfMean(mclVa(x, "x"), NULL)),
        mlfStd(mclVa(x, "x"), NULL, NULL)));
    /*
     * nullCDF  =  normcdf(zScores , 0 , 1);
     */
    mlfAssign(
      &nullCDF,
      mlfNNormcdf(
        1,
        NULL,
        NULL,
        mclVv(zScores, "zScores"),
        _mxarray23_,
        _mxarray0_,
        NULL,
        NULL));
    /*
     * 
     * %
     * % Compute the test statistic of interest: T = max|S(x) - nullCDF(x)|.
     * %
     * 
     * delta1    =  sampleCDF(1:end-1) - nullCDF;   % Vertical difference at jumps approaching from the LEFT.
     */
    mlfAssign(
      &delta1,
      mclMinus(
        mclArrayRef1(
          mclVv(sampleCDF, "sampleCDF"),
          mlfColon(
            _mxarray0_,
            mclMinus(
              mlfEnd(mclVv(sampleCDF, "sampleCDF"), _mxarray0_, _mxarray0_),
              _mxarray0_),
            NULL)),
        mclVv(nullCDF, "nullCDF")));
    /*
     * delta2    =  sampleCDF(2:end)   - nullCDF;   % Vertical difference at jumps approaching from the RIGHT.
     */
    mlfAssign(
      &delta2,
      mclMinus(
        mclArrayRef1(
          mclVv(sampleCDF, "sampleCDF"),
          mlfColon(
            _mxarray24_,
            mlfEnd(mclVv(sampleCDF, "sampleCDF"), _mxarray0_, _mxarray0_),
            NULL)),
        mclVv(nullCDF, "nullCDF")));
    /*
     * deltaCDF  =  abs([delta1 ; delta2]);
     */
    mlfAssign(
      &deltaCDF,
      mlfAbs(
        mlfVertcat(mclVv(delta1, "delta1"), mclVv(delta2, "delta2"), NULL)));
    /*
     * 
     * KSstatistic =  max(deltaCDF);
     */
    mlfAssign(
      KSstatistic, mlfMax(NULL, mclVv(deltaCDF, "deltaCDF"), NULL, NULL));
    /*
     * 
     * %
     * % Compute the approximate P-value. A NaN is returned if the P-value 
     * % is not found within the available 'Alphas' of Lilliefors' table.
     * %
     * 
     * pValue  =  interp1(Q , a , KSstatistic , 'linear');
     */
    mlfAssign(
      pValue,
      mlfInterp1(
        mclVv(Q, "Q"),
        mclVv(a, "a"),
        mclVv(*KSstatistic, "KSstatistic"),
        _mxarray19_,
        NULL));
    /*
     * 
     * %
     * % To maintain consistency with existing Statistics Toolbox hypothesis
     * % tests, returning "H = 0" implies that we "Do not reject the null 
     * % hypothesis at the significance level of alpha" and "H = 1" implies 
     * % that we "Reject the null hypothesis at significance level of alpha."
     * %
     * 
     * H  =  (KSstatistic > criticalValue);
     */
    mlfAssign(
      &H,
      mclGt(
        mclVv(*KSstatistic, "KSstatistic"),
        mclVv(*criticalValue, "criticalValue")));
    mclValidateOutput(H, 1, nargout_, "H", "lillietest");
    mclValidateOutput(*pValue, 2, nargout_, "pValue", "lillietest");
    mclValidateOutput(*KSstatistic, 3, nargout_, "KSstatistic", "lillietest");
    mclValidateOutput(
      *criticalValue, 4, nargout_, "criticalValue", "lillietest");
    mxDestroyArray(rows);
    mxDestroyArray(columns);
    mxDestroyArray(ans);
    mxDestroyArray(sampleCDF);
    mxDestroyArray(xCDF);
    mxDestroyArray(n);
    mxDestroyArray(emsg);
    mxDestroyArray(a);
    mxDestroyArray(sampleSize);
    mxDestroyArray(quantiles);
    mxDestroyArray(Q);
    mxDestroyArray(zScores);
    mxDestroyArray(nullCDF);
    mxDestroyArray(delta1);
    mxDestroyArray(delta2);
    mxDestroyArray(deltaCDF);
    mxDestroyArray(alpha);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return H;
}
