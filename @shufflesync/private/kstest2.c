/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "kstest2.h"
#include "libmatlbm.h"
#include "libmmfile.h"

static mxChar _array1_[32] = { ' ', 'A', 't', ' ', 'l', 'e', 'a', 's',
                               't', ' ', '2', ' ', 'i', 'n', 'p', 'u',
                               't', 's', ' ', 'a', 'r', 'e', ' ', 'r',
                               'e', 'q', 'u', 'i', 'r', 'e', 'd', '.' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;

static mxChar _array4_[30] = { ' ', 'S', 'a', 'm', 'p', 'l', 'e', ' ',
                               0x0027, 'X', '1', 0x0027, ' ', 'm', 'u',
                               's', 't', ' ', 'b', 'e', ' ', 'a', ' ',
                               'v', 'e', 'c', 't', 'o', 'r', '.' };
static mxArray * _mxarray3_;

static mxChar _array6_[30] = { ' ', 'S', 'a', 'm', 'p', 'l', 'e', ' ',
                               0x0027, 'X', '2', 0x0027, ' ', 'm', 'u',
                               's', 't', ' ', 'b', 'e', ' ', 'a', ' ',
                               'v', 'e', 'c', 't', 'o', 'r', '.' };
static mxArray * _mxarray5_;

static mxChar _array8_[45] = { ' ', 'S', 'a', 'm', 'p', 'l', 'e', ' ', 'v', 'e',
                               'c', 't', 'o', 'r', ' ', 0x0027, 'X', '1',
                               0x0027, ' ', 'i', 's', ' ', 'c', 'o', 'm', 'p',
                               'o', 's', 'e', 'd', ' ', 'o', 'f', ' ', 'a', 'l',
                               'l', ' ', 'N', 'a', 'N', 0x0027, 's', '.' };
static mxArray * _mxarray7_;

static mxChar _array10_[45] = { ' ', 'S', 'a', 'm', 'p', 'l', 'e', ' ',
                                'v', 'e', 'c', 't', 'o', 'r', ' ',
                                0x0027, 'X', '2', 0x0027, ' ', 'i', 's',
                                ' ', 'c', 'o', 'm', 'p', 'o', 's', 'e',
                                'd', ' ', 'o', 'f', ' ', 'a', 'l', 'l',
                                ' ', 'N', 'a', 'N', 0x0027, 's', '.' };
static mxArray * _mxarray9_;

static mxChar _array12_[45] = { ' ', 'S', 'i', 'g', 'n', 'i', 'f', 'i',
                                'c', 'a', 'n', 'c', 'e', ' ', 'l', 'e',
                                'v', 'e', 'l', ' ', 0x0027, 'A', 'l',
                                'p', 'h', 'a', 0x0027, ' ', 'm', 'u',
                                's', 't', ' ', 'b', 'e', ' ', 'a', ' ',
                                's', 'c', 'a', 'l', 'a', 'r', '.' };
static mxArray * _mxarray11_;

static mxChar _array14_[52] = { ' ', 'S', 'i', 'g', 'n', 'i', 'f', 'i', 'c',
                                'a', 'n', 'c', 'e', ' ', 'l', 'e', 'v', 'e',
                                'l', ' ', 0x0027, 'A', 'l', 'p', 'h', 'a',
                                0x0027, ' ', 'm', 'u', 's', 't', ' ', 'b',
                                'e', ' ', 'b', 'e', 't', 'w', 'e', 'e', 'n',
                                ' ', '0', ' ', 'a', 'n', 'd', ' ', '1', '.' };
static mxArray * _mxarray13_;
static mxArray * _mxarray15_;
static mxArray * _mxarray16_;

static mxChar _array18_[48] = { ' ', 'T', 'y', 'p', 'e', '-', 'o', 'f', '-',
                                't', 'e', 's', 't', ' ', 'i', 'n', 'd', 'i',
                                'c', 'a', 't', 'o', 'r', ' ', 0x0027, 'T',
                                'a', 'i', 'l', 0x0027, ' ', 'm', 'u', 's',
                                't', ' ', 'b', 'e', ' ', 'a', ' ', 's', 'c',
                                'a', 'l', 'a', 'r', '.' };
static mxArray * _mxarray17_;

static mxChar _array20_[51] = { ' ', 'T', 'y', 'p', 'e', '-', 'o', 'f', '-',
                                't', 'e', 's', 't', ' ', 'i', 'n', 'd', 'i',
                                'c', 'a', 't', 'o', 'r', ' ', 0x0027, 'T',
                                'a', 'i', 'l', 0x0027, ' ', 'm', 'u', 's',
                                't', ' ', 'b', 'e', ' ', '-', '1', ',', ' ',
                                '0', ',', ' ', 'o', 'r', ' ', '1', '.' };
static mxArray * _mxarray19_;
static mxArray * _mxarray21_;
static double _ieee_minusinf_;
static mxArray * _mxarray22_;
static double _ieee_plusinf_;
static mxArray * _mxarray23_;
static mxArray * _mxarray24_;
static mxArray * _mxarray25_;
static mxArray * _mxarray26_;

static double _array28_[101] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                                 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
                                 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0,
                                 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0,
                                 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0,
                                 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
                                 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0,
                                 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0,
                                 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0,
                                 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0,
                                 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0,
                                 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0,
                                 98.0, 99.0, 100.0, 101.0 };
static mxArray * _mxarray27_;
static mxArray * _mxarray29_;

void InitializeModule_kstest2(void) {
    _mxarray0_ = mclInitializeString(32, _array1_);
    _mxarray2_ = mclInitializeDouble(1.0);
    _mxarray3_ = mclInitializeString(30, _array4_);
    _mxarray5_ = mclInitializeString(30, _array6_);
    _mxarray7_ = mclInitializeString(45, _array8_);
    _mxarray9_ = mclInitializeString(45, _array10_);
    _mxarray11_ = mclInitializeString(45, _array12_);
    _mxarray13_ = mclInitializeString(52, _array14_);
    _mxarray15_ = mclInitializeDouble(0.0);
    _mxarray16_ = mclInitializeDouble(.05);
    _mxarray17_ = mclInitializeString(48, _array18_);
    _mxarray19_ = mclInitializeString(51, _array20_);
    _mxarray21_ = mclInitializeDouble(-1.0);
    _ieee_minusinf_ = mclGetMinusInf();
    _mxarray22_ = mclInitializeDouble(_ieee_minusinf_);
    _ieee_plusinf_ = mclGetInf();
    _mxarray23_ = mclInitializeDouble(_ieee_plusinf_);
    _mxarray24_ = mclInitializeDouble(.12);
    _mxarray25_ = mclInitializeDouble(.11);
    _mxarray26_ = mclInitializeDouble(-2.0);
    _mxarray27_ = mclInitializeDoubleVector(101, 1, _array28_);
    _mxarray29_ = mclInitializeDouble(2.0);
}

void TerminateModule_kstest2(void) {
    mxDestroyArray(_mxarray29_);
    mxDestroyArray(_mxarray27_);
    mxDestroyArray(_mxarray26_);
    mxDestroyArray(_mxarray25_);
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray23_);
    mxDestroyArray(_mxarray22_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mkstest2(mxArray * * pValue,
                          mxArray * * KSstatistic,
                          int nargout_,
                          mxArray * x1,
                          mxArray * x2,
                          mxArray * alpha,
                          mxArray * tail);

_mexLocalFunctionTable _local_function_table_kstest2
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfKstest2" contains the normal interface for the "kstest2"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/kstest2.m"
 * (lines 1-186). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
mxArray * mlfKstest2(mxArray * * pValue,
                     mxArray * * KSstatistic,
                     mxArray * x1,
                     mxArray * x2,
                     mxArray * alpha,
                     mxArray * tail) {
    int nargout = 1;
    mxArray * H = NULL;
    mxArray * pValue__ = NULL;
    mxArray * KSstatistic__ = NULL;
    mlfEnterNewContext(2, 4, pValue, KSstatistic, x1, x2, alpha, tail);
    if (pValue != NULL) {
        ++nargout;
    }
    if (KSstatistic != NULL) {
        ++nargout;
    }
    H = Mkstest2(&pValue__, &KSstatistic__, nargout, x1, x2, alpha, tail);
    mlfRestorePreviousContext(2, 4, pValue, KSstatistic, x1, x2, alpha, tail);
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
    return mlfReturnValue(H);
}

/*
 * The function "mlxKstest2" contains the feval interface for the "kstest2"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/kstest2.m"
 * (lines 1-186). The feval function calls the implementation version of
 * kstest2 through this function. This function processes any input arguments
 * and passes them to the implementation version of the function, appearing
 * above.
 */
void mlxKstest2(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[4];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: kstest2 Line: 1 Column: "
            "1 The function \"kstest2\" was called with mor"
            "e than the declared number of outputs (3)."),
          NULL);
    }
    if (nrhs > 4) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: kstest2 Line: 1 Column:"
            " 1 The function \"kstest2\" was called with m"
            "ore than the declared number of inputs (4)."),
          NULL);
    }
    for (i = 0; i < 3; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 4 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 4; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 4, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);
    mplhs[0]
      = Mkstest2(
          &mplhs[1], &mplhs[2], nlhs, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);
    mlfRestorePreviousContext(0, 4, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 3 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 3; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mkstest2" is the implementation version of the "kstest2"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/kstest2.m"
 * (lines 1-186). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [H,pValue,KSstatistic] = kstest2(x1 , x2 , alpha , tail)
 */
static mxArray * Mkstest2(mxArray * * pValue,
                          mxArray * * KSstatistic,
                          int nargout_,
                          mxArray * x1,
                          mxArray * x2,
                          mxArray * alpha,
                          mxArray * tail) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_kstest2);
    int nargin_ = mclNargin(4, x1, x2, alpha, tail, NULL);
    mxArray * H = NULL;
    mxArray * j = NULL;
    mxArray * lambda = NULL;
    mxArray * n = NULL;
    mxArray * n2 = NULL;
    mxArray * n1 = NULL;
    mxArray * deltaCDF = NULL;
    mxArray * sampleCDF2 = NULL;
    mxArray * sampleCDF1 = NULL;
    mxArray * sumCounts2 = NULL;
    mxArray * sumCounts1 = NULL;
    mxArray * binCounts2 = NULL;
    mxArray * binCounts1 = NULL;
    mxArray * binEdges = NULL;
    mxArray * columns2 = NULL;
    mxArray * rows2 = NULL;
    mxArray * columns1 = NULL;
    mxArray * rows1 = NULL;
    mxArray * ans = NULL;
    mclCopyArray(&x1);
    mclCopyArray(&x2);
    mclCopyArray(&alpha);
    mclCopyArray(&tail);
    /*
     * %KSTEST2 Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
     * %   H = KSTEST2(X1,X2,ALPHA,TAIL) performs a Kolmogorov-Smirnov (K-S) test 
     * %   to determine if independent random samples, X1 and X2, are drawn from 
     * %   the same underlying continuous population. ALPHA and TAIL are optional
     * %   scalar inputs: ALPHA is the desired significance level (default = 0.05); 
     * %   TAIL indicates the type of test (default = 0). H indicates the result of
     * %   the hypothesis test:
     * %      H = 0 => Do not reject the null hypothesis at significance level ALPHA.
     * %      H = 1 => Reject the null hypothesis at significance level ALPHA.
     * % 
     * %   Let F1(x) and F2(x) be the empirical distribution functions from sample 
     * %   vectors X1 and X2, respectively. The 2-sample K-S test hypotheses and 
     * %   test statistic are:
     * %
     * %   Null Hypothesis: F1(x) = F2(x) for all x
     * %      For TAIL =  0 (2-sided test), alternative: F1(x) not equal to F2(x).
     * %      For TAIL =  1 (1-sided test), alternative: F1(x) > F2(x).
     * %      For TAIL = -1 (1-sided test), alternative: F1(x) < F2(x).
     * %
     * %   For TAIL = 0, 1, and -1, the test statistics are T = max|F1(x) - F2(x)|,
     * %   T = max[F1(x) - F2(x)], and T = max[F2(x) - F1(x)], respectively.
     * %
     * %   The decision to reject the null hypothesis occurs when the significance 
     * %   level, ALPHA, equals or exceeds the P-value.
     * %
     * %   X1 and X2 are row or column vectors of lengths N1 and N2, respectively, 
     * %   and represent random samples from some underlying distribution(s). 
     * %   Missing observations, indicated by NaN's (Not-a-Number), are ignored.
     * %
     * %   [H,P] = KSTEST2(...) also returns the asymptotic P-value P.
     * %
     * %   [H,P,KSSTAT] = KSTEST2(...) also returns the K-S test statistic KSSTAT
     * %   defined above for the test type indicated by TAIL.
     * %
     * %   The asymptotic P-value becomes very accurate for large sample sizes, and
     * %   is believed to be reasonably accurate for sample sizes N1 and N2 such 
     * %   that (N1*N2)/(N1 + N2) >= 4.
     * %
     * %   See also KSTEST, LILLIETEST, CDFPLOT.
     * %
     * 
     * % Author(s): R.A. Baker, 08/14/98
     * % Copyright 1993-2002 The MathWorks, Inc. 
     * % $Revision: 1.5 $   $ Date: 1998/01/30 13:45:34 $
     * 
     * %
     * % References:
     * %   (1) Massey, F.J., "The Kolmogorov-Smirnov Test for Goodness of Fit",
     * %         Journal of the American Statistical Association, 46 (March 1956), 68-77.
     * %   (2) Miller, L.H., "Table of Percentage Points of Kolmogorov Statistics",
     * %         Journal of the American Statistical Association, (March 1951), 111-121.
     * %   (3) Conover, W.J., "Practical Nonparametric Statistics", 
     * %         John Wiley & Sons, Inc., 1980.
     * %   (4) Press, W.H., et. al., "Numerical Recipes in C", 
     * %         Cambridge University Press, 1992.
     * 
     * if nargin < 2
     */
    if (nargin_ < 2) {
        /*
         * error(' At least 2 inputs are required.');
         */
        mlfError(_mxarray0_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Ensure each sample is a VECTOR.
     * %
     * 
     * [rows1 , columns1]  =  size(x1);
     */
    mlfSize(mlfVarargout(&rows1, &columns1, NULL), mclVa(x1, "x1"), NULL);
    /*
     * [rows2 , columns2]  =  size(x2);
     */
    mlfSize(mlfVarargout(&rows2, &columns2, NULL), mclVa(x2, "x2"), NULL);
    /*
     * 
     * if (rows1 ~= 1) & (columns1 ~= 1) 
     */
    {
        mxArray * a_ = mclInitialize(mclNe(mclVv(rows1, "rows1"), _mxarray2_));
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(a_, mclNe(mclVv(columns1, "columns1"), _mxarray2_)))) {
            mxDestroyArray(a_);
            /*
             * error(' Sample ''X1'' must be a vector.');
             */
            mlfError(_mxarray3_, NULL);
        } else {
            mxDestroyArray(a_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * if (rows2 ~= 1) & (columns2 ~= 1) 
     */
    {
        mxArray * a_ = mclInitialize(mclNe(mclVv(rows2, "rows2"), _mxarray2_));
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(a_, mclNe(mclVv(columns2, "columns2"), _mxarray2_)))) {
            mxDestroyArray(a_);
            /*
             * error(' Sample ''X2'' must be a vector.');
             */
            mlfError(_mxarray5_, NULL);
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
     * % Remove missing observations indicated by NaN's, and 
     * % ensure that valid observations remain.
     * %
     * 
     * x1  =  x1(~isnan(x1));
     */
    mlfAssign(
      &x1, mclArrayRef1(mclVa(x1, "x1"), mclNot(mlfIsnan(mclVa(x1, "x1")))));
    /*
     * x2  =  x2(~isnan(x2));
     */
    mlfAssign(
      &x2, mclArrayRef1(mclVa(x2, "x2"), mclNot(mlfIsnan(mclVa(x2, "x2")))));
    /*
     * x1  =  x1(:);
     */
    mlfAssign(&x1, mclArrayRef1(mclVa(x1, "x1"), mlfCreateColonIndex()));
    /*
     * x2  =  x2(:);
     */
    mlfAssign(&x2, mclArrayRef1(mclVa(x2, "x2"), mlfCreateColonIndex()));
    /*
     * 
     * if isempty(x1)
     */
    if (mlfTobool(mlfIsempty(mclVa(x1, "x1")))) {
        /*
         * error(' Sample vector ''X1'' is composed of all NaN''s.');
         */
        mlfError(_mxarray7_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * if isempty(x2)
     */
    if (mlfTobool(mlfIsempty(mclVa(x2, "x2")))) {
        /*
         * error(' Sample vector ''X2'' is composed of all NaN''s.');
         */
        mlfError(_mxarray9_, NULL);
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Ensure the significance level, ALPHA, is a scalar 
     * % between 0 and 1 and set default if necessary.
     * %
     * 
     * if (nargin >= 3) & ~isempty(alpha)
     */
    {
        mxArray * a_ = mclInitialize(mclBoolToArray(nargin_ >= 3));
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
                  _mxarray2_)) {
                /*
                 * error(' Significance level ''Alpha'' must be a scalar.');
                 */
                mlfError(_mxarray11_, NULL);
            /*
             * end
             */
            }
            /*
             * if (alpha <= 0 | alpha >= 1)
             */
            {
                mxArray * a_0
                  = mclInitialize(mclLe(mclVa(alpha, "alpha"), _mxarray15_));
                if (mlfTobool(a_0)
                    || mlfTobool(
                         mclOr(
                           a_0, mclGe(mclVa(alpha, "alpha"), _mxarray2_)))) {
                    mxDestroyArray(a_0);
                    /*
                     * error(' Significance level ''Alpha'' must be between 0 and 1.'); 
                     */
                    mlfError(_mxarray13_, NULL);
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
            mlfAssign(&alpha, _mxarray16_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Ensure the type-of-test indicator, TAIL, is a scalar integer from 
     * % the allowable set {-1 , 0 , 1}, and set default if necessary.
     * %
     * 
     * if (nargin >= 4) & ~isempty(tail)
     */
    {
        mxArray * a_ = mclInitialize(mclBoolToArray(nargin_ >= 4));
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(a_, mclNot(mlfIsempty(mclVa(tail, "tail")))))) {
            mxDestroyArray(a_);
            /*
             * if prod(size(tail)) > 1
             */
            if (mclGtBool(
                  mlfProd(
                    mlfSize(mclValueVarargout(), mclVa(tail, "tail"), NULL),
                    NULL),
                  _mxarray2_)) {
                /*
                 * error(' Type-of-test indicator ''Tail'' must be a scalar.');
                 */
                mlfError(_mxarray17_, NULL);
            /*
             * end
             */
            }
            /*
             * if (tail ~= -1) & (tail ~= 0) & (tail ~= 1)
             */
            {
                mxArray * a_1
                  = mclInitialize(mclNe(mclVa(tail, "tail"), _mxarray21_));
                if (mlfTobool(a_1)) {
                    mlfAssign(
                      &a_1,
                      mclAnd(a_1, mclNe(mclVa(tail, "tail"), _mxarray15_)));
                } else {
                    mlfAssign(&a_1, mlfScalar(0));
                }
                if (mlfTobool(a_1)
                    && mlfTobool(
                         mclAnd(
                           a_1, mclNe(mclVa(tail, "tail"), _mxarray2_)))) {
                    mxDestroyArray(a_1);
                    /*
                     * error(' Type-of-test indicator ''Tail'' must be -1, 0, or 1.');
                     */
                    mlfError(_mxarray19_, NULL);
                } else {
                    mxDestroyArray(a_1);
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
             * tail  =  0;
             */
            mlfAssign(&tail, _mxarray15_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
     * %
     * 
     * binEdges    =  [-inf ; sort([x1;x2]) ; inf];
     */
    mlfAssign(
      &binEdges,
      mlfVertcat(
        _mxarray22_,
        mlfSort(NULL, mlfVertcat(mclVa(x1, "x1"), mclVa(x2, "x2"), NULL), NULL),
        _mxarray23_,
        NULL));
    /*
     * 
     * binCounts1  =  histc (x1 , binEdges);
     */
    mlfAssign(
      &binCounts1,
      mlfNHistc(
        0,
        mclValueVarargout(),
        mclVa(x1, "x1"),
        mclVv(binEdges, "binEdges"),
        NULL));
    /*
     * binCounts2  =  histc (x2 , binEdges);
     */
    mlfAssign(
      &binCounts2,
      mlfNHistc(
        0,
        mclValueVarargout(),
        mclVa(x2, "x2"),
        mclVv(binEdges, "binEdges"),
        NULL));
    /*
     * 
     * sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
     */
    mlfAssign(
      &sumCounts1,
      mclRdivide(
        mlfCumsum(mclVv(binCounts1, "binCounts1"), NULL),
        mlfSum(mclVv(binCounts1, "binCounts1"), NULL)));
    /*
     * sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);
     */
    mlfAssign(
      &sumCounts2,
      mclRdivide(
        mlfCumsum(mclVv(binCounts2, "binCounts2"), NULL),
        mlfSum(mclVv(binCounts2, "binCounts2"), NULL)));
    /*
     * 
     * sampleCDF1  =  sumCounts1(1:end-1);
     */
    mlfAssign(
      &sampleCDF1,
      mclArrayRef1(
        mclVv(sumCounts1, "sumCounts1"),
        mlfColon(
          _mxarray2_,
          mclMinus(
            mlfEnd(mclVv(sumCounts1, "sumCounts1"), _mxarray2_, _mxarray2_),
            _mxarray2_),
          NULL)));
    /*
     * sampleCDF2  =  sumCounts2(1:end-1);
     */
    mlfAssign(
      &sampleCDF2,
      mclArrayRef1(
        mclVv(sumCounts2, "sumCounts2"),
        mlfColon(
          _mxarray2_,
          mclMinus(
            mlfEnd(mclVv(sumCounts2, "sumCounts2"), _mxarray2_, _mxarray2_),
            _mxarray2_),
          NULL)));
    /*
     * 
     * %
     * % Compute the test statistic of interest.
     * %
     * 
     * switch tail
     */
    {
        mxArray * v_ = mclInitialize(mclVa(tail, "tail"));
        if (mclSwitchCompare(v_, _mxarray15_)) {
            /*
             * case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
             * deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
             */
            mlfAssign(
              &deltaCDF,
              mlfAbs(
                mclMinus(
                  mclVv(sampleCDF1, "sampleCDF1"),
                  mclVv(sampleCDF2, "sampleCDF2"))));
        /*
         * 
         * case -1      %  1-sided test: T = max[F2(x) - F1(x)].
         */
        } else if (mclSwitchCompare(v_, _mxarray21_)) {
            /*
             * deltaCDF  =  sampleCDF2 - sampleCDF1;
             */
            mlfAssign(
              &deltaCDF,
              mclMinus(
                mclVv(sampleCDF2, "sampleCDF2"),
                mclVv(sampleCDF1, "sampleCDF1")));
        /*
         * 
         * case  1      %  1-sided test: T = max[F1(x) - F2(x)].
         */
        } else if (mclSwitchCompare(v_, _mxarray2_)) {
            /*
             * deltaCDF  =  sampleCDF1 - sampleCDF2;
             */
            mlfAssign(
              &deltaCDF,
              mclMinus(
                mclVv(sampleCDF1, "sampleCDF1"),
                mclVv(sampleCDF2, "sampleCDF2")));
        /*
         * end
         */
        }
        mxDestroyArray(v_);
    }
    /*
     * 
     * KSstatistic   =  max(deltaCDF);
     */
    mlfAssign(
      KSstatistic, mlfMax(NULL, mclVv(deltaCDF, "deltaCDF"), NULL, NULL));
    /*
     * 
     * %
     * % Compute the asymptotic P-value approximation and accept or
     * % reject the null hypothesis on the basis of the P-value.
     * %
     * 
     * n1     =  length(x1);
     */
    mlfAssign(&n1, mlfScalar(mclLengthInt(mclVa(x1, "x1"))));
    /*
     * n2     =  length(x2);
     */
    mlfAssign(&n2, mlfScalar(mclLengthInt(mclVa(x2, "x2"))));
    /*
     * n      =  n1 * n2 /(n1 + n2);
     */
    mlfAssign(
      &n,
      mclMrdivide(
        mclMtimes(mclVv(n1, "n1"), mclVv(n2, "n2")),
        mclPlus(mclVv(n1, "n1"), mclVv(n2, "n2"))));
    /*
     * lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);
     */
    mlfAssign(
      &lambda,
      mlfMax(
        NULL,
        mclMtimes(
          mclPlus(
            mclPlus(mlfSqrt(mclVv(n, "n")), _mxarray24_),
            mclMrdivide(_mxarray25_, mlfSqrt(mclVv(n, "n")))),
          mclVv(*KSstatistic, "KSstatistic")),
        _mxarray15_,
        NULL));
    /*
     * 
     * if tail ~= 0        % 1-sided test.
     */
    if (mclNeBool(mclVa(tail, "tail"), _mxarray15_)) {
        /*
         * 
         * pValue  =  exp(-2 * lambda * lambda);
         */
        mlfAssign(
          pValue,
          mlfExp(
            mclMtimes(
              mclMtimes(_mxarray26_, mclVv(lambda, "lambda")),
              mclVv(lambda, "lambda"))));
    /*
     * 
     * else                % 2-sided test (default).
     */
    } else {
        /*
         * %
         * %  Use the asymptotic Q-function to approximate the 2-sided P-value.
         * %
         * j       =  [1:101]';
         */
        mlfAssign(&j, _mxarray27_);
        /*
         * pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
         */
        mlfAssign(
          pValue,
          mclMtimes(
            _mxarray29_,
            mlfSum(
              mclTimes(
                mlfPower(_mxarray21_, mclMinus(mclVv(j, "j"), _mxarray2_)),
                mlfExp(
                  mclMtimes(
                    mclMtimes(
                      mclMtimes(_mxarray26_, mclVv(lambda, "lambda")),
                      mclVv(lambda, "lambda")),
                    mlfPower(mclVv(j, "j"), _mxarray29_)))),
              NULL)));
        /*
         * 
         * if pValue < 0 , pValue = 0; end
         */
        if (mclLtBool(mclVv(*pValue, "pValue"), _mxarray15_)) {
            mlfAssign(pValue, _mxarray15_);
        }
        /*
         * if pValue > 1 , pValue = 1; end
         */
        if (mclGtBool(mclVv(*pValue, "pValue"), _mxarray2_)) {
            mlfAssign(pValue, _mxarray2_);
        }
    /*
     * 
     * end
     */
    }
    /*
     * 
     * H  =  (alpha >= pValue);
     */
    mlfAssign(&H, mclGe(mclVa(alpha, "alpha"), mclVv(*pValue, "pValue")));
    mclValidateOutput(H, 1, nargout_, "H", "kstest2");
    mclValidateOutput(*pValue, 2, nargout_, "pValue", "kstest2");
    mclValidateOutput(*KSstatistic, 3, nargout_, "KSstatistic", "kstest2");
    mxDestroyArray(ans);
    mxDestroyArray(rows1);
    mxDestroyArray(columns1);
    mxDestroyArray(rows2);
    mxDestroyArray(columns2);
    mxDestroyArray(binEdges);
    mxDestroyArray(binCounts1);
    mxDestroyArray(binCounts2);
    mxDestroyArray(sumCounts1);
    mxDestroyArray(sumCounts2);
    mxDestroyArray(sampleCDF1);
    mxDestroyArray(sampleCDF2);
    mxDestroyArray(deltaCDF);
    mxDestroyArray(n1);
    mxDestroyArray(n2);
    mxDestroyArray(n);
    mxDestroyArray(lambda);
    mxDestroyArray(j);
    mxDestroyArray(tail);
    mxDestroyArray(alpha);
    mxDestroyArray(x2);
    mxDestroyArray(x1);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return H;
}
