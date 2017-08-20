/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "convmtx.h"
#include "libmatlbm.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;

void InitializeModule_convmtx(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeDouble(-1.0);
    _mxarray2_ = mclInitializeDouble(2.0);
    _mxarray3_ = mclInitializeDouble(0.0);
}

void TerminateModule_convmtx(void) {
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mconvmtx(int nargout_, mxArray * v, mxArray * n);

_mexLocalFunctionTable _local_function_table_convmtx
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfConvmtx" contains the normal interface for the "convmtx"
 * M-function from file
 * "/Applications/MATLAB6p5p1/toolbox/signal/signal/convmtx.m" (lines 1-36).
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfConvmtx(mxArray * v, mxArray * n) {
    int nargout = 1;
    mxArray * t = NULL;
    mlfEnterNewContext(0, 2, v, n);
    t = Mconvmtx(nargout, v, n);
    mlfRestorePreviousContext(0, 2, v, n);
    return mlfReturnValue(t);
}

/*
 * The function "mlxConvmtx" contains the feval interface for the "convmtx"
 * M-function from file
 * "/Applications/MATLAB6p5p1/toolbox/signal/signal/convmtx.m" (lines 1-36).
 * The feval function calls the implementation version of convmtx through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxConvmtx(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: convmtx Line: 1 Column: "
            "1 The function \"convmtx\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: convmtx Line: 1 Column:"
            " 1 The function \"convmtx\" was called with m"
            "ore than the declared number of inputs (2)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = Mconvmtx(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mconvmtx" is the implementation version of the "convmtx"
 * M-function from file
 * "/Applications/MATLAB6p5p1/toolbox/signal/signal/convmtx.m" (lines 1-36). It
 * contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function t = convmtx(v,n)
 */
static mxArray * Mconvmtx(int nargout_, mxArray * v, mxArray * n) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_convmtx);
    mxArray * t = NULL;
    mxArray * ridx = NULL;
    mxArray * cidx = NULL;
    mxArray * x = NULL;
    mxArray * m = NULL;
    mxArray * r = NULL;
    mxArray * c = NULL;
    mxArray * mn = NULL;
    mxArray * lv = NULL;
    mxArray * nv = NULL;
    mxArray * mv = NULL;
    mclCopyArray(&v);
    mclCopyArray(&n);
    /*
     * %CONVMTX Convolution matrix.
     * %   CONVMTX(C,N) returns the convolution matrix for vector C.
     * %   If C is a column vector and X is a column vector of length N,
     * %   then CONVMTX(C,N)*X is the same as CONV(C,X).
     * %   If R is a row vector and X is a row vector of length N,
     * %   then X*CONVMTX(R,N) is the same as CONV(R,X).
     * %   See also CONV.
     * 
     * %   Author(s): L. Shure, 5-17-88
     * %   	   T. Krauss, 3-30-93, removed dependence on toeplitz
     * %   Copyright 1988-2002 The MathWorks, Inc.
     * %   $Revision: 1.6 $  $Date: 2002/03/28 17:27:27 $
     * 
     * [mv,nv] = size(v);
     */
    mlfSize(mlfVarargout(&mv, &nv, NULL), mclVa(v, "v"), NULL);
    /*
     * lv = max(mv,nv);
     */
    mlfAssign(&lv, mlfMax(NULL, mclVv(mv, "mv"), mclVv(nv, "nv"), NULL));
    /*
     * v = v(:);		% make v a column vector
     */
    mlfAssign(&v, mclArrayRef1(mclVa(v, "v"), mlfCreateColonIndex()));
    /*
     * mn = lv + n - 1;	% mn == number of rows of M; n == number of columns
     */
    mlfAssign(
      &mn, mclMinus(mclPlus(mclVv(lv, "lv"), mclVa(n, "n")), _mxarray0_));
    /*
     * 
     * %  t = toeplitz([v; zeros(n-1,1)],zeros(n,1));  put Toeplitz code inline
     * c = [v; zeros(n-1,1)];
     */
    mlfAssign(
      &c,
      mlfVertcat(
        mclVa(v, "v"),
        mlfZeros(mclMinus(mclVa(n, "n"), _mxarray0_), _mxarray0_, NULL),
        NULL));
    /*
     * r = zeros(n,1);
     */
    mlfAssign(&r, mlfZeros(mclVa(n, "n"), _mxarray0_, NULL));
    /*
     * m = length(c);
     */
    mlfAssign(&m, mlfScalar(mclLengthInt(mclVv(c, "c"))));
    /*
     * x = [r(n:-1:2) ; c(:)];                 % build vector of user data
     */
    mlfAssign(
      &x,
      mlfVertcat(
        mclArrayRef1(
          mclVv(r, "r"), mlfColon(mclVa(n, "n"), _mxarray1_, _mxarray2_)),
        mclArrayRef1(mclVv(c, "c"), mlfCreateColonIndex()),
        NULL));
    /*
     * %
     * cidx = (0:m-1)';
     */
    mlfAssign(
      &cidx,
      mlfCtranspose(
        mlfColon(_mxarray3_, mclMinus(mclVv(m, "m"), _mxarray0_), NULL)));
    /*
     * ridx = n:-1:1;
     */
    mlfAssign(&ridx, mlfColon(mclVa(n, "n"), _mxarray1_, _mxarray0_));
    /*
     * t = cidx(:,ones(n,1)) + ridx(ones(m,1),:);    % Toeplitz subscripts
     */
    mlfAssign(
      &t,
      mclPlus(
        mclArrayRef2(
          mclVv(cidx, "cidx"),
          mlfCreateColonIndex(),
          mlfOnes(mclVa(n, "n"), _mxarray0_, NULL)),
        mclArrayRef2(
          mclVv(ridx, "ridx"),
          mlfOnes(mclVv(m, "m"), _mxarray0_, NULL),
          mlfCreateColonIndex())));
    /*
     * t(:) = x(t);                            % actual data
     */
    mclArrayAssign1(
      &t, mclArrayRef1(mclVv(x, "x"), mclVv(t, "t")), mlfCreateColonIndex());
    /*
     * % end of toeplitz code
     * 
     * if mv < nv
     */
    if (mclLtBool(mclVv(mv, "mv"), mclVv(nv, "nv"))) {
        /*
         * t = t.';
         */
        mlfAssign(&t, mlfTranspose(mclVv(t, "t")));
    /*
     * end
     */
    }
    mclValidateOutput(t, 1, nargout_, "t", "convmtx");
    mxDestroyArray(mv);
    mxDestroyArray(nv);
    mxDestroyArray(lv);
    mxDestroyArray(mn);
    mxDestroyArray(c);
    mxDestroyArray(r);
    mxDestroyArray(m);
    mxDestroyArray(x);
    mxDestroyArray(cidx);
    mxDestroyArray(ridx);
    mxDestroyArray(n);
    mxDestroyArray(v);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return t;
    /*
     * 
     */
}
