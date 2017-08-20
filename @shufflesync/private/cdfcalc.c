/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "cdfcalc.h"
#include "libmatlbm.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;

static mxChar _array3_[1] = { 'X' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;

static mxChar _array6_[33] = { 'I', 'n', 'p', 'u', 't', ' ', 's', 'a', 'm',
                               'p', 'l', 'e', ' ', '%', 's', ' ', 'm', 'u',
                               's', 't', ' ', 'b', 'e', ' ', 'a', ' ', 'v',
                               'e', 'c', 't', 'o', 'r', '.' };
static mxArray * _mxarray5_;

static mxChar _array8_[46] = { 'I', 'n', 'p', 'u', 't', ' ', 's', 'a', 'm', 'p',
                               'l', 'e', ' ', '%', 's', ' ', 'h', 'a', 's', ' ',
                               'n', 'o', ' ', 'v', 'a', 'l', 'i', 'd', ' ', 'd',
                               'a', 't', 'a', ' ', '(', 'a', 'l', 'l', ' ', 'N',
                               'a', 'N', 0x0027, 's', ')', '.' };
static mxArray * _mxarray7_;
static mxArray * _mxarray9_;

void InitializeModule_cdfcalc(void) {
    _mxarray0_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray1_ = mclInitializeDouble(1.0);
    _mxarray2_ = mclInitializeString(1, _array3_);
    _mxarray4_ = mclInitializeDouble(0.0);
    _mxarray5_ = mclInitializeString(33, _array6_);
    _mxarray7_ = mclInitializeString(46, _array8_);
    _mxarray9_ = mclInitializeCharVector(0, 0, (mxChar *)NULL);
}

void TerminateModule_cdfcalc(void) {
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mcdfcalc(mxArray * * xCDF,
                          mxArray * * n,
                          mxArray * * emsg,
                          int nargout_,
                          mxArray * x,
                          mxArray * xname);

_mexLocalFunctionTable _local_function_table_cdfcalc
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfCdfcalc" contains the normal interface for the "cdfcalc"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/cdfcalc.m"
 * (lines 1-64). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfCdfcalc(mxArray * * xCDF,
                     mxArray * * n,
                     mxArray * * emsg,
                     mxArray * x,
                     mxArray * xname) {
    int nargout = 1;
    mxArray * yCDF = NULL;
    mxArray * xCDF__ = NULL;
    mxArray * n__ = NULL;
    mxArray * emsg__ = NULL;
    mlfEnterNewContext(3, 2, xCDF, n, emsg, x, xname);
    if (xCDF != NULL) {
        ++nargout;
    }
    if (n != NULL) {
        ++nargout;
    }
    if (emsg != NULL) {
        ++nargout;
    }
    yCDF = Mcdfcalc(&xCDF__, &n__, &emsg__, nargout, x, xname);
    mlfRestorePreviousContext(3, 2, xCDF, n, emsg, x, xname);
    if (xCDF != NULL) {
        mclCopyOutputArg(xCDF, xCDF__);
    } else {
        mxDestroyArray(xCDF__);
    }
    if (n != NULL) {
        mclCopyOutputArg(n, n__);
    } else {
        mxDestroyArray(n__);
    }
    if (emsg != NULL) {
        mclCopyOutputArg(emsg, emsg__);
    } else {
        mxDestroyArray(emsg__);
    }
    return mlfReturnValue(yCDF);
}

/*
 * The function "mlxCdfcalc" contains the feval interface for the "cdfcalc"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/cdfcalc.m"
 * (lines 1-64). The feval function calls the implementation version of cdfcalc
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxCdfcalc(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[4];
    int i;
    if (nlhs > 4) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: cdfcalc Line: 1 Column: "
            "1 The function \"cdfcalc\" was called with mor"
            "e than the declared number of outputs (4)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: cdfcalc Line: 1 Column:"
            " 1 The function \"cdfcalc\" was called with m"
            "ore than the declared number of inputs (2)."),
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
      = Mcdfcalc(&mplhs[1], &mplhs[2], &mplhs[3], nlhs, mprhs[0], mprhs[1]);
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
 * The function "Mcdfcalc" is the implementation version of the "cdfcalc"
 * M-function from file "/Applications/MATLAB6p5p1/toolbox/stats/cdfcalc.m"
 * (lines 1-64). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [yCDF,xCDF,n,emsg] = cdfcalc(x,xname)
 */
static mxArray * Mcdfcalc(mxArray * * xCDF,
                          mxArray * * n,
                          mxArray * * emsg,
                          int nargout_,
                          mxArray * x,
                          mxArray * xname) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_cdfcalc);
    int nargin_ = mclNargin(2, x, xname, NULL);
    mxArray * yCDF = NULL;
    mxArray * notdup = NULL;
    mclCopyArray(&x);
    mclCopyArray(&xname);
    /*
     * %CDFCALC Calculate an empirical cumulative distribution function.
     * %   [YCDF,XCDF] = CDFCALC(X) calculates an empirical cumulative
     * %   distribution function (CDF) of the observations in the data sample
     * %   vector X. X may be a row or column vector, and represents a random
     * %   sample of observations from some underlying distribution.  On
     * %   return XCDF is the set of X values at which the CDF increases.
     * %   At XCDF(i), the function increases from YCDF(i) to YCDF(i+1).
     * %
     * %   [YCDF,XCDF,N] = CDFCALC(X) also returns N, the sample size.
     * %
     * %   [YCDF,XCDF,N,EMSG] = CDFCALC(X) also returns an error message if X
     * %   is not a vector or if it contains no values other than NaN.
     * %
     * %   See also CDFPLOT.
     * 
     * %   Copyright 1993-2002 The MathWorks, Inc. 
     * %   $Revision: 1.5 $  $Date: 2002/01/17 21:30:11 $
     * 
     * % Ensure the data is a VECTOR.
     * yCDF = [];
     */
    mlfAssign(&yCDF, _mxarray0_);
    /*
     * xCDF = [];
     */
    mlfAssign(xCDF, _mxarray0_);
    /*
     * if (nargin < 2)
     */
    if (nargin_ < 2) {
        /*
         * if isempty(inputname(1))
         */
        if (mlfTobool(mlfIsempty(mlfInputname(_mxarray1_)))) {
            /*
             * xname = 'X';
             */
            mlfAssign(&xname, _mxarray2_);
        /*
         * else
         */
        } else {
            /*
             * xname = inputname(1);
             */
            mlfAssign(&xname, mlfInputname(_mxarray1_));
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    /*
     * n = 0;
     */
    mlfAssign(n, _mxarray4_);
    /*
     * if (min(size(x)) ~= 1)
     */
    if (mclNeBool(
          mlfMin(
            NULL,
            mlfSize(mclValueVarargout(), mclVa(x, "x"), NULL),
            NULL,
            NULL),
          _mxarray1_)) {
        /*
         * emsg = sprintf('Input sample %s must be a vector.', xname);
         */
        mlfAssign(
          emsg, mlfSprintf(NULL, _mxarray5_, mclVa(xname, "xname"), NULL));
        /*
         * return;
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * 
     * % Remove missing observations indicated by NaN's.
     * x = x(~isnan(x));
     */
    mlfAssign(&x, mclArrayRef1(mclVa(x, "x"), mclNot(mlfIsnan(mclVa(x, "x")))));
    /*
     * n = length(x);
     */
    mlfAssign(n, mlfScalar(mclLengthInt(mclVa(x, "x"))));
    /*
     * if n == 0
     */
    if (mclEqBool(mclVv(*n, "n"), _mxarray4_)) {
        /*
         * emsg = sprintf('Input sample %s has no valid data (all NaN''s).', ...
         */
        mlfAssign(
          emsg, mlfSprintf(NULL, _mxarray7_, mclVa(xname, "xname"), NULL));
        /*
         * xname);
         * return;
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * 
     * % Sort observation data in ascending order.
     * x = sort(x(:));
     */
    mlfAssign(
      &x,
      mlfSort(NULL, mclArrayRef1(mclVa(x, "x"), mlfCreateColonIndex()), NULL));
    /*
     * 
     * %
     * % Compute cumulative sum such that the sample CDF is
     * % F(x) = (number of observations <= x) / (total number of observations).
     * % Note that the bin edges are padded with +/- infinity for auto-scaling of
     * % the x-axis.
     * %
     * 
     * % Get cumulative sums
     * yCDF = (1:n)' / n;
     */
    mlfAssign(
      &yCDF,
      mclMrdivide(
        mlfCtranspose(mlfColon(_mxarray1_, mclVv(*n, "n"), NULL)),
        mclVv(*n, "n")));
    /*
     * 
     * % Remove duplicates; only need final one with total count
     * notdup = ([diff(x(:)); 1] > 0);
     */
    mlfAssign(
      &notdup,
      mclGt(
        mlfVertcat(
          mlfDiff(
            mclArrayRef1(mclVa(x, "x"), mlfCreateColonIndex()), NULL, NULL),
          _mxarray1_,
          NULL),
        _mxarray4_));
    /*
     * xCDF = x(notdup);
     */
    mlfAssign(xCDF, mclArrayRef1(mclVa(x, "x"), mclVv(notdup, "notdup")));
    /*
     * yCDF = [0; yCDF(notdup)];
     */
    mlfAssign(
      &yCDF,
      mlfVertcat(
        _mxarray4_,
        mclArrayRef1(mclVv(yCDF, "yCDF"), mclVv(notdup, "notdup")),
        NULL));
    /*
     * emsg = '';
     */
    mlfAssign(emsg, _mxarray9_);
    /*
     * 
     */
    return_:
    mclValidateOutput(yCDF, 1, nargout_, "yCDF", "cdfcalc");
    mclValidateOutput(*xCDF, 2, nargout_, "xCDF", "cdfcalc");
    mclValidateOutput(*n, 3, nargout_, "n", "cdfcalc");
    mclValidateOutput(*emsg, 4, nargout_, "emsg", "cdfcalc");
    mxDestroyArray(notdup);
    mxDestroyArray(xname);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return yCDF;
}
