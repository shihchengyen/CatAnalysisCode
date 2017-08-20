/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "histcie.h"
#include "concatenate.h"
#include "getOptArgs.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "vecc.h"

static mxChar _array1_[8] = { 'D', 'r', 'o', 'p', 'L', 'a', 's', 't' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;

static mxChar _array4_[8] = { 'D', 'a', 't', 'a', 'C', 'o', 'l', 's' };
static mxArray * _mxarray3_;

static mxChar _array6_[5] = { 'f', 'l', 'a', 'g', 's' };
static mxArray * _mxarray5_;

static mxArray * _array8_[2] = { NULL /*_mxarray0_*/, NULL /*_mxarray3_*/ };
static mxArray * _mxarray7_;
static mxArray * _mxarray9_;
static mxArray * _mxarray10_;
static double _ieee_nan_;
static mxArray * _mxarray11_;

static mxChar _array13_[5] = { 'h', 'i', 's', 't', 'c' };
static mxArray * _mxarray12_;
static mxArray * _mxarray14_;

void InitializeModule_histcie(void) {
    _mxarray0_ = mclInitializeString(8, _array1_);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeString(8, _array4_);
    _mxarray5_ = mclInitializeString(5, _array6_);
    _array8_[0] = _mxarray0_;
    _array8_[1] = _mxarray3_;
    _mxarray7_ = mclInitializeCellVector(1, 2, _array8_);
    _mxarray9_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray10_ = mclInitializeDouble(1.0);
    _ieee_nan_ = mclGetNaN();
    _mxarray11_ = mclInitializeDouble(_ieee_nan_);
    _mxarray12_ = mclInitializeString(5, _array13_);
    _mxarray14_ = mclInitializeDouble(2.0);
}

void TerminateModule_histcie(void) {
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mhistcie(mxArray * * bin,
                          int nargout_,
                          mxArray * x,
                          mxArray * edges,
                          mxArray * varargin);

_mexLocalFunctionTable _local_function_table_histcie
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfNHistcie" contains the nargout interface for the "histcie"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/histcie.m"
 * (lines 1-74). This interface is only produced if the M-function uses the
 * special variable "nargout". The nargout interface allows the number of
 * requested outputs to be specified via the nargout argument, as opposed to
 * the normal interface which dynamically calculates the number of outputs
 * based on the number of non-NULL inputs it receives. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
mxArray * mlfNHistcie(int nargout,
                      mxArray * * bin,
                      mxArray * x,
                      mxArray * edges,
                      ...) {
    mxArray * varargin = NULL;
    mxArray * n = NULL;
    mxArray * bin__ = NULL;
    mlfVarargin(&varargin, edges, 0);
    mlfEnterNewContext(1, -3, bin, x, edges, varargin);
    n = Mhistcie(&bin__, nargout, x, edges, varargin);
    mlfRestorePreviousContext(1, 2, bin, x, edges);
    mxDestroyArray(varargin);
    if (bin != NULL) {
        mclCopyOutputArg(bin, bin__);
    } else {
        mxDestroyArray(bin__);
    }
    return mlfReturnValue(n);
}

/*
 * The function "mlfHistcie" contains the normal interface for the "histcie"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/histcie.m"
 * (lines 1-74). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfHistcie(mxArray * * bin, mxArray * x, mxArray * edges, ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * n = NULL;
    mxArray * bin__ = NULL;
    mlfVarargin(&varargin, edges, 0);
    mlfEnterNewContext(1, -3, bin, x, edges, varargin);
    if (bin != NULL) {
        ++nargout;
    }
    n = Mhistcie(&bin__, nargout, x, edges, varargin);
    mlfRestorePreviousContext(1, 2, bin, x, edges);
    mxDestroyArray(varargin);
    if (bin != NULL) {
        mclCopyOutputArg(bin, bin__);
    } else {
        mxDestroyArray(bin__);
    }
    return mlfReturnValue(n);
}

/*
 * The function "mlfVHistcie" contains the void interface for the "histcie"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/histcie.m"
 * (lines 1-74). The void interface is only produced if the M-function uses the
 * special variable "nargout", and has at least one output. The void interface
 * function specifies zero output arguments to the implementation version of
 * the function, and in the event that the implementation version still returns
 * an output (which, in MATLAB, would be assigned to the "ans" variable), it
 * deallocates the output. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlfVHistcie(mxArray * x, mxArray * edges, ...) {
    mxArray * varargin = NULL;
    mxArray * n = NULL;
    mxArray * bin = NULL;
    mlfVarargin(&varargin, edges, 0);
    mlfEnterNewContext(0, -3, x, edges, varargin);
    n = Mhistcie(&bin, 0, x, edges, varargin);
    mlfRestorePreviousContext(0, 2, x, edges);
    mxDestroyArray(varargin);
    mxDestroyArray(n);
}

/*
 * The function "mlxHistcie" contains the feval interface for the "histcie"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/histcie.m"
 * (lines 1-74). The feval function calls the implementation version of histcie
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxHistcie(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: histcie Line: 1 Column: "
            "1 The function \"histcie\" was called with mor"
            "e than the declared number of outputs (2)."),
          NULL);
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mprhs[2] = NULL;
    mlfAssign(&mprhs[2], mclCreateVararginCell(nrhs - 2, prhs + 2));
    mplhs[0] = Mhistcie(&mplhs[1], nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
    mxDestroyArray(mprhs[2]);
}

/*
 * The function "Mhistcie" is the implementation version of the "histcie"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/histcie.m"
 * (lines 1-74). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [n,bin] = histcie(x,edges,varargin)
 */
static mxArray * Mhistcie(mxArray * * bin,
                          int nargout_,
                          mxArray * x,
                          mxArray * edges,
                          mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_histcie);
    mxArray * n = NULL;
    mxArray * ans = NULL;
    mxArray * bfi = NULL;
    mxArray * nh1 = NULL;
    mxArray * nhl = NULL;
    mxArray * xcols = NULL;
    mxArray * xrows = NULL;
    mxArray * binh = NULL;
    mxArray * nh = NULL;
    mxArray * Args = NULL;
    mclCopyArray(&x);
    mclCopyArray(&edges);
    mclCopyArray(&varargin);
    /*
     * % HISTCIE Histogram count including end point
     * %   [N,BIN] = HISTCIE(X,EDGES) is the same as the MATLAB function
     * %   HISTC except it adds the last value returned from HISTC (i.e.
     * %   number of values of X that match EDGES(end)) to the second 
     * %   last bin and sets the last bin to zero. The last value is
     * %   returned in N so that the function BAR(EDGES,N,'histc') can
     * %   still be used to create plots easily. BIN is also modified so
     * %   that any value in X matching EDGES(end) will be set to bin
     * %   length(EDGES)-1. HISTCIE will also return a zero vector in N 
     * %   as well as an empty matrix in BIN if X is empty.
     * %
     * %   [N,BIN] = HISTCIE(X,EDGES,'DropLast') drops the last value in
     * %   N, which is always 0.
     * %
     * %   [N,BIN] = HISTCIE(X,EDGES,'DataCols') forces the analysis to
     * %   use the row vector in X as separate data series with one data
     * %   point each.
     * 
     * Args = struct('DropLast',0,'DataCols',0);
     */
    mlfAssign(
      &Args, mlfStruct(_mxarray0_, _mxarray2_, _mxarray3_, _mxarray2_, NULL));
    /*
     * Args = getOptArgs(varargin,Args,'flags',{'DropLast','DataCols'});
     */
    mlfAssign(
      &Args,
      mlfGetOptArgs(
        NULL,
        mclVa(varargin, "varargin"),
        mclVv(Args, "Args"),
        _mxarray5_,
        _mxarray7_,
        NULL));
    /*
     * 
     * % check if x is empty
     * if(isempty(x))
     */
    if (mlfTobool(mlfIsempty(mclVa(x, "x")))) {
        /*
         * % return all zeros in n and emtpy matrix in bin
         * nh = zeros(size(edges));
         */
        mlfAssign(
          &nh,
          mlfZeros(
            mlfSize(mclValueVarargout(), mclVa(edges, "edges"), NULL), NULL));
        /*
         * % make sure nh is a column vector
         * nh = vecc(nh);
         */
        mlfAssign(&nh, mlfVecc(mclVv(nh, "nh")));
        /*
         * binh = [];
         */
        mlfAssign(&binh, _mxarray9_);
    /*
     * else
     */
    } else {
        /*
         * % get size of x
         * [xrows,xcols] = size(x);    
         */
        mlfSize(mlfVarargout(&xrows, &xcols, NULL), mclVa(x, "x"), NULL);
        /*
         * % check if it is a row vector
         * if( (xrows==1) && (xcols>1) )
         */
        if (mclScalarToBool(mclEq(mclVv(xrows, "xrows"), _mxarray10_))
            && mclScalarToBool(mclGt(mclVv(xcols, "xcols"), _mxarray10_))) {
            /*
             * if(Args.DataCols) 
             */
            if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".DataCols"))) {
                /*
                 * % add a row of NaN's to make sure the data is treated as
                 * % columns
                 * x = concatenate(x,NaN);
                 */
                mlfAssign(&x, mlfConcatenate(mclVa(x, "x"), _mxarray11_, NULL));
            /*
             * else
             */
            } else {
                /*
                 * % switch x to column vector
                 * x = vecc(x);
                 */
                mlfAssign(&x, mlfVecc(mclVa(x, "x")));
            /*
             * end
             */
            }
        /*
         * end
         */
        }
        /*
         * [nh,binh] = histc(x,edges);
         */
        mlfNHistc(
          0,
          mlfVarargout(&nh, &binh, NULL),
          mclVa(x, "x"),
          mclVa(edges, "edges"),
          NULL);
        /*
         * % make sure nh is a column vector, especially when x is 1x1
         * nh = vecc(nh);
         */
        mlfAssign(&nh, mlfVecc(mclVv(nh, "nh")));
        /*
         * % find length of n
         * nhl = size(nh,1);
         */
        mlfAssign(
          &nhl, mlfSize(mclValueVarargout(), mclVv(nh, "nh"), _mxarray10_));
        /*
         * % get second last index
         * nh1 = nhl - 1;
         */
        mlfAssign(&nh1, mclMinus(mclVv(nhl, "nhl"), _mxarray10_));
        /*
         * % add last value to second last value
         * nh(nh1,:) = nh(nh1,:) + nh(nhl,:);
         */
        mclArrayAssign2(
          &nh,
          mclPlus(
            mclArrayRef2(
              mclVv(nh, "nh"), mclVv(nh1, "nh1"), mlfCreateColonIndex()),
            mclArrayRef2(
              mclVv(nh, "nh"), mclVv(nhl, "nhl"), mlfCreateColonIndex())),
          mclVv(nh1, "nh1"),
          mlfCreateColonIndex());
        /*
         * % set last value to zero
         * nh(nhl,:) = 0;
         */
        mclArrayAssign2(
          &nh, _mxarray2_, mclVv(nhl, "nhl"), mlfCreateColonIndex());
        /*
         * % find binh == edges(end) and set them to nhl - 1 (i.e. nh1)
         * bfi = find(binh==nhl);
         */
        mlfAssign(
          &bfi,
          mlfFind(NULL, NULL, mclEq(mclVv(binh, "binh"), mclVv(nhl, "nhl"))));
        /*
         * if(~isempty(bfi))
         */
        if (mclNotBool(mlfIsempty(mclVv(bfi, "bfi")))) {
            /*
             * binh(bfi) = nh1;
             */
            mclArrayAssign1(&binh, mclVv(nh1, "nh1"), mclVv(bfi, "bfi"));
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
     * % if there were no output arguments, behave like hist and plot the histogram
     * if (nargout == 0)
     */
    if (nargout_ == 0) {
        /*
         * bar(edges,nh,'histc')
         */
        mclPrintAns(
          &ans,
          mlfNBar(
            0,
            NULL,
            mclVa(edges, "edges"),
            mclVv(nh, "nh"),
            _mxarray12_,
            NULL));
    /*
     * else
     */
    } else {
        /*
         * if(Args.DropLast)
         */
        if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".DropLast"))) {
            /*
             * % return output arguments with last value dropped
             * n = nh(1:(end-1),:);
             */
            mlfAssign(
              &n,
              mclArrayRef2(
                mclVv(nh, "nh"),
                mlfColon(
                  _mxarray10_,
                  mclMinus(
                    mlfEnd(mclVv(nh, "nh"), _mxarray10_, _mxarray14_),
                    _mxarray10_),
                  NULL),
                mlfCreateColonIndex()));
        /*
         * else
         */
        } else {
            /*
             * n = nh;
             */
            mlfAssign(&n, mclVv(nh, "nh"));
        /*
         * end
         */
        }
        /*
         * bin = binh;	
         */
        mlfAssign(bin, mclVv(binh, "binh"));
    /*
     * end
     */
    }
    mclValidateOutput(n, 1, nargout_, "n", "histcie");
    mclValidateOutput(*bin, 2, nargout_, "bin", "histcie");
    mxDestroyArray(Args);
    mxDestroyArray(nh);
    mxDestroyArray(binh);
    mxDestroyArray(xrows);
    mxDestroyArray(xcols);
    mxDestroyArray(nhl);
    mxDestroyArray(nh1);
    mxDestroyArray(bfi);
    mxDestroyArray(ans);
    mxDestroyArray(varargin);
    mxDestroyArray(edges);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return n;
}
