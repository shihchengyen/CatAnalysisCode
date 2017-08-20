/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "concat.h"
#include "getOptArgs.h"
#include "libmatlbm.h"
#include "libmmfile.h"

static mxChar _array1_[10] = { 'C', 'o', 'l', 'u', 'm',
                               'n', 'w', 'i', 's', 'e' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;

static mxChar _array4_[13] = { 'D', 'i', 's', 'c', 'a', 'r', 'd',
                               'E', 'm', 'p', 't', 'y', 'A' };
static mxArray * _mxarray3_;

static mxChar _array6_[3] = { 'P', 'a', 'd' };
static mxArray * _mxarray5_;
static double _ieee_nan_;
static mxArray * _mxarray7_;

static mxArray * _array9_[2] = { NULL /*_mxarray0_*/, NULL /*_mxarray3_*/ };
static mxArray * _mxarray8_;
static mxArray * _mxarray10_;
static mxArray * _mxarray11_;
static mxArray * _mxarray12_;

void InitializeModule_concat(void) {
    _mxarray0_ = mclInitializeString(10, _array1_);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeString(13, _array4_);
    _mxarray5_ = mclInitializeString(3, _array6_);
    _ieee_nan_ = mclGetNaN();
    _mxarray7_ = mclInitializeDouble(_ieee_nan_);
    _array9_[0] = _mxarray0_;
    _array9_[1] = _mxarray3_;
    _mxarray8_ = mclInitializeCellVector(1, 2, _array9_);
    _mxarray10_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray11_ = mclInitializeDouble(1.0);
    _mxarray12_ = mclInitializeDouble(2.0);
}

void TerminateModule_concat(void) {
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mconcat(int nargout_, mxArray * a, mxArray * varargin);

_mexLocalFunctionTable _local_function_table_concat
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfConcat" contains the normal interface for the "concat"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/concat.m"
 * (lines 1-127). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
mxArray * mlfConcat(mxArray * a, ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * c = NULL;
    mlfVarargin(&varargin, a, 0);
    mlfEnterNewContext(0, -2, a, varargin);
    c = Mconcat(nargout, a, varargin);
    mlfRestorePreviousContext(0, 1, a);
    mxDestroyArray(varargin);
    return mlfReturnValue(c);
}

/*
 * The function "mlxConcat" contains the feval interface for the "concat"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/concat.m"
 * (lines 1-127). The feval function calls the implementation version of concat
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxConcat(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: concat Line: 1 Column: "
            "1 The function \"concat\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mprhs[1] = NULL;
    mlfAssign(&mprhs[1], mclCreateVararginCell(nrhs - 1, prhs + 1));
    mplhs[0] = Mconcat(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
    mxDestroyArray(mprhs[1]);
}

/*
 * The function "Mconcat" is the implementation version of the "concat"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/concat.m"
 * (lines 1-127). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function c = concat(a,varargin)
 */
static mxArray * Mconcat(int nargout_, mxArray * a, mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_concat);
    mxArray * c = NULL;
    mxArray * b0 = NULL;
    mxArray * a0 = NULL;
    mxArray * bcols = NULL;
    mxArray * brows = NULL;
    mxArray * acols = NULL;
    mxArray * arows = NULL;
    mxArray * b = NULL;
    mxArray * idx = NULL;
    mxArray * nNumArgs = NULL;
    mxArray * pad = NULL;
    mxArray * Args = NULL;
    mclCopyArray(&a);
    mclCopyArray(&varargin);
    /*
     * %CONCAT Concatenate two row vectors
     * %   C = CONCAT(A,B) concatenates two row vectors 
     * %   (i.e. C = [A; B]). The shorter vector is padded with NaNs.
     * %   Empty vectors are concatenated as well, i.e. if A and B 
     * %   are both of size [0 0], C will be of size [2 0].
     * %
     * %   C = CONCAT(A1,A2,A3,...) concatenates the first N numeric variables
     * %   together, padding each row to make sure they have the same number of
     * %   columns.
     * %
     * %   C = CONCAT(...,'Columnwise') concatenates column vectors.
     * %
     * %   C = CONCAT(...,'Pad',VALUE) pads the shorter vector with
     * %   VALUE instead. 
     * %   e.g. C = CONCAT(A,B,'Pad',0) pads the shorter vector with 0's.
     * %
     * %   C = CONCAT(...,'DiscardEmptyA') returns B if A is 
     * %   of size [0 0]. If B is also of size [0 0], an empty matrix
     * %   of size [1 0] will be returned.  
     * %
     * %   Dependencies: None.
     * 
     * Args = struct('Columnwise',0,'DiscardEmptyA',0,'Pad',nan);
     */
    mlfAssign(
      &Args,
      mlfStruct(
        _mxarray0_,
        _mxarray2_,
        _mxarray3_,
        _mxarray2_,
        _mxarray5_,
        _mxarray7_,
        NULL));
    /*
     * Args.flags = {'Columnwise','DiscardEmptyA'};
     */
    mlfIndexAssign(&Args, ".flags", _mxarray8_);
    /*
     * Args = getOptArgs(varargin,Args);
     */
    mlfAssign(
      &Args,
      mlfGetOptArgs(
        NULL, mclVa(varargin, "varargin"), mclVv(Args, "Args"), NULL));
    /*
     * 
     * pad = Args.Pad;
     */
    mlfAssign(&pad, mlfIndexRef(mclVv(Args, "Args"), ".Pad"));
    /*
     * 
     * % get number of numeric arguments
     * nNumArgs = length(Args.NumericArguments);
     */
    mlfAssign(
      &nNumArgs,
      mclFeval(
        mclValueVarargout(),
        mlxLength,
        mlfIndexRef(mclVv(Args, "Args"), ".NumericArguments"),
        NULL));
    /*
     * if(nNumArgs==0)
     */
    if (mclEqBool(mclVv(nNumArgs, "nNumArgs"), _mxarray2_)) {
        /*
         * c = a;
         */
        mlfAssign(&c, mclVa(a, "a"));
    /*
     * end
     */
    }
    /*
     * 
     * for idx = 1:nNumArgs
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclVv(nNumArgs, "nNumArgs"));
        if (v_ > e_) {
            mlfAssign(&idx, _mxarray10_);
        } else {
            /*
             * b = Args.NumericArguments{idx};
             * 
             * % get size of a and b
             * [arows acols] = size(a);
             * [brows bcols] = size(b);
             * % check if a or b are 0 by 0
             * if(sum([arows acols])==0)
             * a0 = 1;
             * else
             * a0 = 0;
             * end
             * if(sum([brows bcols])==0)
             * b0 = 1;
             * else
             * b0 = 0;
             * end
             * 
             * if(a0)
             * % a is 0 by 0
             * if(Args.DiscardEmptyA)
             * if(b0)
             * % both a and b are 0 by 0 but we are going to disregard a
             * % so return matrix that is 1 by 0 or 0 by 1
             * if(Args.Columnwise)
             * c = repmat(pad,0,1);
             * else
             * c = repmat(pad,1,0);
             * end
             * else % if(b0)
             * % a is 0 by 0 but b is not so just return b
             * c = b;
             * end % if(b0)
             * else % if(Args.DiscardEmptyA)
             * if(b0)
             * % both a and b are 0 by 0 so return matrix that is 2 by 0 or 0
             * % by 2
             * if(Args.Columnwise)
             * c = repmat(pad,0,2);
             * else
             * c = repmat(pad,2,0);
             * end
             * else % if(b0)
             * % a is 0 by 0 but b is not, so resize a to b and concatenate
             * if(Args.Columnwise)
             * c = [repmat(pad,brows,1) b];
             * else
             * c = [repmat(pad,1,bcols); b];
             * end
             * end % if(b0)
             * end % if(Args.DiscardEmptyA)
             * else % if(a0)
             * % a is not 0 by 0
             * if(b0)
             * % a is not 0 by 0 but b is, so resize b to a and concatenate
             * if(Args.Columnwise)
             * c = [a repmat(pad,arows,1)];
             * else
             * c = [a; repmat(pad,1,acols)];
             * end
             * else % if(b0)
             * % both a and b are not 0 by 0 so do normal concatenation which may
             * % include empty matrices that have one non-zero dimension
             * if(Args.Columnwise)
             * if(arows>brows)
             * % pad b before concatenating
             * c = [a [b; repmat(pad,arows-brows,bcols)]];
             * elseif(brows>arows)
             * % pad a before concatenating
             * c = [[a; repmat(pad,brows-arows,acols)] b];
             * else
             * % a and b are the same size so just concatenate
             * c = [a b];
             * end
             * else % if(Args.Columnwise)
             * if(acols>bcols)
             * % pad b before concatenating
             * c = [a; [b repmat(pad,brows,acols-bcols)]];
             * elseif(bcols>acols)
             * % pad a before concatenating
             * c = [[a repmat(pad,arows,bcols-acols)]; b];
             * else
             * % a and b are the same size so just concatenate
             * c = [a; b];
             * end
             * end % if(Args.Columnwise)
             * end % if(b0)
             * end % if(a0)
             * 
             * a = c;
             * end
             */
            for (; ; ) {
                mlfAssign(
                  &b,
                  mlfIndexRef(
                    mclVv(Args, "Args"),
                    ".NumericArguments{?}",
                    mlfScalar(v_)));
                mlfSize(
                  mlfVarargout(&arows, &acols, NULL), mclVa(a, "a"), NULL);
                mlfSize(
                  mlfVarargout(&brows, &bcols, NULL), mclVv(b, "b"), NULL);
                if (mclEqBool(
                      mlfSum(
                        mlfHorzcat(
                          mclVv(arows, "arows"), mclVv(acols, "acols"), NULL),
                        NULL),
                      _mxarray2_)) {
                    mlfAssign(&a0, _mxarray11_);
                } else {
                    mlfAssign(&a0, _mxarray2_);
                }
                if (mclEqBool(
                      mlfSum(
                        mlfHorzcat(
                          mclVv(brows, "brows"), mclVv(bcols, "bcols"), NULL),
                        NULL),
                      _mxarray2_)) {
                    mlfAssign(&b0, _mxarray11_);
                } else {
                    mlfAssign(&b0, _mxarray2_);
                }
                if (mlfTobool(mclVv(a0, "a0"))) {
                    if (mlfTobool(
                          mlfIndexRef(mclVv(Args, "Args"), ".DiscardEmptyA"))) {
                        if (mlfTobool(mclVv(b0, "b0"))) {
                            if (mlfTobool(
                                  mlfIndexRef(
                                    mclVv(Args, "Args"), ".Columnwise"))) {
                                mlfAssign(
                                  &c,
                                  mlfRepmat(
                                    mclVv(pad, "pad"),
                                    _mxarray2_,
                                    _mxarray11_));
                            } else {
                                mlfAssign(
                                  &c,
                                  mlfRepmat(
                                    mclVv(pad, "pad"),
                                    _mxarray11_,
                                    _mxarray2_));
                            }
                        } else {
                            mlfAssign(&c, mclVv(b, "b"));
                        }
                    } else {
                        if (mlfTobool(mclVv(b0, "b0"))) {
                            if (mlfTobool(
                                  mlfIndexRef(
                                    mclVv(Args, "Args"), ".Columnwise"))) {
                                mlfAssign(
                                  &c,
                                  mlfRepmat(
                                    mclVv(pad, "pad"),
                                    _mxarray2_,
                                    _mxarray12_));
                            } else {
                                mlfAssign(
                                  &c,
                                  mlfRepmat(
                                    mclVv(pad, "pad"),
                                    _mxarray12_,
                                    _mxarray2_));
                            }
                        } else {
                            if (mlfTobool(
                                  mlfIndexRef(
                                    mclVv(Args, "Args"), ".Columnwise"))) {
                                mlfAssign(
                                  &c,
                                  mlfHorzcat(
                                    mlfRepmat(
                                      mclVv(pad, "pad"),
                                      mclVv(brows, "brows"),
                                      _mxarray11_),
                                    mclVv(b, "b"),
                                    NULL));
                            } else {
                                mlfAssign(
                                  &c,
                                  mlfVertcat(
                                    mlfRepmat(
                                      mclVv(pad, "pad"),
                                      _mxarray11_,
                                      mclVv(bcols, "bcols")),
                                    mclVv(b, "b"),
                                    NULL));
                            }
                        }
                    }
                } else {
                    if (mlfTobool(mclVv(b0, "b0"))) {
                        if (mlfTobool(
                              mlfIndexRef(
                                mclVv(Args, "Args"), ".Columnwise"))) {
                            mlfAssign(
                              &c,
                              mlfHorzcat(
                                mclVa(a, "a"),
                                mlfRepmat(
                                  mclVv(pad, "pad"),
                                  mclVv(arows, "arows"),
                                  _mxarray11_),
                                NULL));
                        } else {
                            mlfAssign(
                              &c,
                              mlfVertcat(
                                mclVa(a, "a"),
                                mlfRepmat(
                                  mclVv(pad, "pad"),
                                  _mxarray11_,
                                  mclVv(acols, "acols")),
                                NULL));
                        }
                    } else {
                        if (mlfTobool(
                              mlfIndexRef(
                                mclVv(Args, "Args"), ".Columnwise"))) {
                            if (mclGtBool(
                                  mclVv(arows, "arows"),
                                  mclVv(brows, "brows"))) {
                                mlfAssign(
                                  &c,
                                  mlfHorzcat(
                                    mclVa(a, "a"),
                                    mlfVertcat(
                                      mclVv(b, "b"),
                                      mlfRepmat(
                                        mclVv(pad, "pad"),
                                        mclMinus(
                                          mclVv(arows, "arows"),
                                          mclVv(brows, "brows")),
                                        mclVv(bcols, "bcols")),
                                      NULL),
                                    NULL));
                            } else if (mclGtBool(
                                         mclVv(brows, "brows"),
                                         mclVv(arows, "arows"))) {
                                mlfAssign(
                                  &c,
                                  mlfHorzcat(
                                    mlfVertcat(
                                      mclVa(a, "a"),
                                      mlfRepmat(
                                        mclVv(pad, "pad"),
                                        mclMinus(
                                          mclVv(brows, "brows"),
                                          mclVv(arows, "arows")),
                                        mclVv(acols, "acols")),
                                      NULL),
                                    mclVv(b, "b"),
                                    NULL));
                            } else {
                                mlfAssign(
                                  &c,
                                  mlfHorzcat(
                                    mclVa(a, "a"), mclVv(b, "b"), NULL));
                            }
                        } else {
                            if (mclGtBool(
                                  mclVv(acols, "acols"),
                                  mclVv(bcols, "bcols"))) {
                                mlfAssign(
                                  &c,
                                  mlfVertcat(
                                    mclVa(a, "a"),
                                    mlfHorzcat(
                                      mclVv(b, "b"),
                                      mlfRepmat(
                                        mclVv(pad, "pad"),
                                        mclVv(brows, "brows"),
                                        mclMinus(
                                          mclVv(acols, "acols"),
                                          mclVv(bcols, "bcols"))),
                                      NULL),
                                    NULL));
                            } else if (mclGtBool(
                                         mclVv(bcols, "bcols"),
                                         mclVv(acols, "acols"))) {
                                mlfAssign(
                                  &c,
                                  mlfVertcat(
                                    mlfHorzcat(
                                      mclVa(a, "a"),
                                      mlfRepmat(
                                        mclVv(pad, "pad"),
                                        mclVv(arows, "arows"),
                                        mclMinus(
                                          mclVv(bcols, "bcols"),
                                          mclVv(acols, "acols"))),
                                      NULL),
                                    mclVv(b, "b"),
                                    NULL));
                            } else {
                                mlfAssign(
                                  &c,
                                  mlfVertcat(
                                    mclVa(a, "a"), mclVv(b, "b"), NULL));
                            }
                        }
                    }
                }
                mlfAssign(&a, mclVv(c, "c"));
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&idx, mlfScalar(v_));
        }
    }
    mclValidateOutput(c, 1, nargout_, "c", "concat");
    mxDestroyArray(Args);
    mxDestroyArray(pad);
    mxDestroyArray(nNumArgs);
    mxDestroyArray(idx);
    mxDestroyArray(b);
    mxDestroyArray(arows);
    mxDestroyArray(acols);
    mxDestroyArray(brows);
    mxDestroyArray(bcols);
    mxDestroyArray(a0);
    mxDestroyArray(b0);
    mxDestroyArray(varargin);
    mxDestroyArray(a);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return c;
}
