/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "concatenate.h"
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

static mxArray * _array6_[2] = { NULL /*_mxarray0_*/, NULL /*_mxarray3_*/ };
static mxArray * _mxarray5_;
static mxArray * _mxarray7_;
static double _ieee_nan_;
static mxArray * _mxarray8_;
static mxArray * _mxarray9_;

void InitializeModule_concatenate(void) {
    _mxarray0_ = mclInitializeString(10, _array1_);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeString(13, _array4_);
    _array6_[0] = _mxarray0_;
    _array6_[1] = _mxarray3_;
    _mxarray5_ = mclInitializeCellVector(1, 2, _array6_);
    _mxarray7_ = mclInitializeDouble(1.0);
    _ieee_nan_ = mclGetNaN();
    _mxarray8_ = mclInitializeDouble(_ieee_nan_);
    _mxarray9_ = mclInitializeDouble(2.0);
}

void TerminateModule_concatenate(void) {
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mconcatenate(int nargout_,
                              mxArray * a,
                              mxArray * b,
                              mxArray * varargin);

_mexLocalFunctionTable _local_function_table_concatenate
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfConcatenate" contains the normal interface for the
 * "concatenate" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/concatenate
 * .m" (lines 1-114). This function processes any input arguments and passes
 * them to the implementation version of the function, appearing above.
 */
mxArray * mlfConcatenate(mxArray * a, mxArray * b, ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * c = NULL;
    mlfVarargin(&varargin, b, 0);
    mlfEnterNewContext(0, -3, a, b, varargin);
    c = Mconcatenate(nargout, a, b, varargin);
    mlfRestorePreviousContext(0, 2, a, b);
    mxDestroyArray(varargin);
    return mlfReturnValue(c);
}

/*
 * The function "mlxConcatenate" contains the feval interface for the
 * "concatenate" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/concatenate
 * .m" (lines 1-114). The feval function calls the implementation version of
 * concatenate through this function. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
void mlxConcatenate(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: concatenate Line: 1 Column"
            ": 1 The function \"concatenate\" was called with"
            " more than the declared number of outputs (1)."),
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
    mprhs[2] = NULL;
    mlfAssign(&mprhs[2], mclCreateVararginCell(nrhs - 2, prhs + 2));
    mplhs[0] = Mconcatenate(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    mxDestroyArray(mprhs[2]);
}

/*
 * The function "Mconcatenate" is the implementation version of the
 * "concatenate" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/concatenate
 * .m" (lines 1-114). It contains the actual compiled code for that M-function.
 * It is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function c = concatenate (a,b,varargin)
 */
static mxArray * Mconcatenate(int nargout_,
                              mxArray * a,
                              mxArray * b,
                              mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_concatenate);
    mxArray * c = NULL;
    mxArray * b0 = NULL;
    mxArray * a0 = NULL;
    mxArray * bcols = NULL;
    mxArray * brows = NULL;
    mxArray * acols = NULL;
    mxArray * arows = NULL;
    mxArray * pad = NULL;
    mxArray * Args = NULL;
    mclCopyArray(&a);
    mclCopyArray(&b);
    mclCopyArray(&varargin);
    /*
     * %CONCATENATE Concatenate two row vectors
     * %   C = CONCATENATE(A,B) concatenates two row vectors 
     * %   (i.e. C = [A; B]). The shorter vector is padded with NaNs.
     * %   Empty vectors are concatenated as well, i.e. if A and B 
     * %   are both of size [0 0], C will be of size [2 0].
     * %
     * %   C = CONCATENATE(A,B,VALUE) pads the shorter vector with
     * %   VALUE instead. 
     * %   e.g. C = CONCATENATE(A,B,0) pads the shorter vector with 0's.
     * %
     * %   C = CONCATENATE(A,B,'Columnwise') concatenates two column vectors.
     * %
     * %   C = CONCATENATE(...,'DiscardEmptyA') returns B if A is 
     * %   of size [0 0]. If B is also of size [0 0], an empty matrix
     * %   of size [1 0] will be returned.  
     * %
     * %   Dependencies: None.
     * 
     * Args = struct('Columnwise',0,'DiscardEmptyA',0);
     */
    mlfAssign(
      &Args, mlfStruct(_mxarray0_, _mxarray2_, _mxarray3_, _mxarray2_, NULL));
    /*
     * Args.flags = {'Columnwise','DiscardEmptyA'};
     */
    mlfIndexAssign(&Args, ".flags", _mxarray5_);
    /*
     * Args = getOptArgs(varargin,Args);
     */
    mlfAssign(
      &Args,
      mlfGetOptArgs(
        NULL, mclVa(varargin, "varargin"), mclVv(Args, "Args"), NULL));
    /*
     * if(~isempty(Args.NumericArguments))
     */
    if (mclNotBool(
          mclFeval(
            mclValueVarargout(),
            mlxIsempty,
            mlfIndexRef(mclVv(Args, "Args"), ".NumericArguments"),
            NULL))) {
        /*
         * pad = Args.NumericArguments{1};
         */
        mlfAssign(
          &pad,
          mlfIndexRef(mclVv(Args, "Args"), ".NumericArguments{?}", _mxarray7_));
    /*
     * else
     */
    } else {
        /*
         * pad = NaN;
         */
        mlfAssign(&pad, _mxarray8_);
    /*
     * end
     */
    }
    /*
     * 
     * % get size of a and b
     * [arows acols] = size(a);
     */
    mlfSize(mlfVarargout(&arows, &acols, NULL), mclVa(a, "a"), NULL);
    /*
     * [brows bcols] = size(b);
     */
    mlfSize(mlfVarargout(&brows, &bcols, NULL), mclVa(b, "b"), NULL);
    /*
     * % check if a or b are 0 by 0
     * if(sum([arows acols])==0)
     */
    if (mclEqBool(
          mlfSum(
            mlfHorzcat(mclVv(arows, "arows"), mclVv(acols, "acols"), NULL),
            NULL),
          _mxarray2_)) {
        /*
         * a0 = 1;
         */
        mlfAssign(&a0, _mxarray7_);
    /*
     * else
     */
    } else {
        /*
         * a0 = 0;
         */
        mlfAssign(&a0, _mxarray2_);
    /*
     * end
     */
    }
    /*
     * if(sum([brows bcols])==0)
     */
    if (mclEqBool(
          mlfSum(
            mlfHorzcat(mclVv(brows, "brows"), mclVv(bcols, "bcols"), NULL),
            NULL),
          _mxarray2_)) {
        /*
         * b0 = 1;
         */
        mlfAssign(&b0, _mxarray7_);
    /*
     * else
     */
    } else {
        /*
         * b0 = 0;
         */
        mlfAssign(&b0, _mxarray2_);
    /*
     * end
     */
    }
    /*
     * 
     * if(a0)
     */
    if (mlfTobool(mclVv(a0, "a0"))) {
        /*
         * % a is 0 by 0
         * if(Args.DiscardEmptyA)
         */
        if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".DiscardEmptyA"))) {
            /*
             * if(b0)
             */
            if (mlfTobool(mclVv(b0, "b0"))) {
                /*
                 * % both a and b are 0 by 0 but we are going to disregard a
                 * % so return matrix that is 1 by 0 or 0 by 1
                 * if(Args.Columnwise)
                 */
                if (mlfTobool(
                      mlfIndexRef(mclVv(Args, "Args"), ".Columnwise"))) {
                    /*
                     * c = repmat(pad,0,1);
                     */
                    mlfAssign(
                      &c, mlfRepmat(mclVv(pad, "pad"), _mxarray2_, _mxarray7_));
                /*
                 * else
                 */
                } else {
                    /*
                     * c = repmat(pad,1,0);
                     */
                    mlfAssign(
                      &c, mlfRepmat(mclVv(pad, "pad"), _mxarray7_, _mxarray2_));
                /*
                 * end
                 */
                }
            /*
             * else % if(b0)
             */
            } else {
                /*
                 * % a is 0 by 0 but b is not so just return b
                 * c = b;
                 */
                mlfAssign(&c, mclVa(b, "b"));
            /*
             * end % if(b0)
             */
            }
        /*
         * else % if(Args.DiscardEmptyA)
         */
        } else {
            /*
             * if(b0)
             */
            if (mlfTobool(mclVv(b0, "b0"))) {
                /*
                 * % both a and b are 0 by 0 so return matrix that is 2 by 0 or 0
                 * % by 2
                 * if(Args.Columnwise)
                 */
                if (mlfTobool(
                      mlfIndexRef(mclVv(Args, "Args"), ".Columnwise"))) {
                    /*
                     * c = repmat(pad,0,2);
                     */
                    mlfAssign(
                      &c, mlfRepmat(mclVv(pad, "pad"), _mxarray2_, _mxarray9_));
                /*
                 * else
                 */
                } else {
                    /*
                     * c = repmat(pad,2,0);
                     */
                    mlfAssign(
                      &c, mlfRepmat(mclVv(pad, "pad"), _mxarray9_, _mxarray2_));
                /*
                 * end
                 */
                }
            /*
             * else % if(b0)
             */
            } else {
                /*
                 * % a is 0 by 0 but b is not, so resize a to b and concatenate
                 * if(Args.Columnwise)
                 */
                if (mlfTobool(
                      mlfIndexRef(mclVv(Args, "Args"), ".Columnwise"))) {
                    /*
                     * c = [repmat(pad,brows,1) b];
                     */
                    mlfAssign(
                      &c,
                      mlfHorzcat(
                        mlfRepmat(
                          mclVv(pad, "pad"), mclVv(brows, "brows"), _mxarray7_),
                        mclVa(b, "b"),
                        NULL));
                /*
                 * else
                 */
                } else {
                    /*
                     * c = [repmat(pad,1,bcols); b];
                     */
                    mlfAssign(
                      &c,
                      mlfVertcat(
                        mlfRepmat(
                          mclVv(pad, "pad"), _mxarray7_, mclVv(bcols, "bcols")),
                        mclVa(b, "b"),
                        NULL));
                /*
                 * end
                 */
                }
            /*
             * end % if(b0)
             */
            }
        /*
         * end % if(Args.DiscardEmptyA)
         */
        }
    /*
     * else % if(a0)
     */
    } else {
        /*
         * % a is not 0 by 0
         * if(b0)
         */
        if (mlfTobool(mclVv(b0, "b0"))) {
            /*
             * % a is not 0 by 0 but b is, so resize b to a and concatenate
             * if(Args.Columnwise)
             */
            if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".Columnwise"))) {
                /*
                 * c = [a repmat(pad,arows,1)];
                 */
                mlfAssign(
                  &c,
                  mlfHorzcat(
                    mclVa(a, "a"),
                    mlfRepmat(
                      mclVv(pad, "pad"), mclVv(arows, "arows"), _mxarray7_),
                    NULL));
            /*
             * else
             */
            } else {
                /*
                 * c = [a; repmat(pad,1,acols)];
                 */
                mlfAssign(
                  &c,
                  mlfVertcat(
                    mclVa(a, "a"),
                    mlfRepmat(
                      mclVv(pad, "pad"), _mxarray7_, mclVv(acols, "acols")),
                    NULL));
            /*
             * end
             */
            }
        /*
         * else % if(b0)
         */
        } else {
            /*
             * % both a and b are not 0 by 0 so do normal concatenation which may
             * % include empty matrices that have one non-zero dimension
             * if(Args.Columnwise)
             */
            if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".Columnwise"))) {
                /*
                 * if(arows>brows)
                 */
                if (mclGtBool(mclVv(arows, "arows"), mclVv(brows, "brows"))) {
                    /*
                     * % pad b before concatenating
                     * c = [a [b; repmat(pad,arows-brows,bcols)]];
                     */
                    mlfAssign(
                      &c,
                      mlfHorzcat(
                        mclVa(a, "a"),
                        mlfVertcat(
                          mclVa(b, "b"),
                          mlfRepmat(
                            mclVv(pad, "pad"),
                            mclMinus(
                              mclVv(arows, "arows"), mclVv(brows, "brows")),
                            mclVv(bcols, "bcols")),
                          NULL),
                        NULL));
                /*
                 * elseif(brows>arows)
                 */
                } else if (mclGtBool(
                             mclVv(brows, "brows"), mclVv(arows, "arows"))) {
                    /*
                     * % pad a before concatenating
                     * c = [[a; repmat(pad,brows-arows,acols)] b];
                     */
                    mlfAssign(
                      &c,
                      mlfHorzcat(
                        mlfVertcat(
                          mclVa(a, "a"),
                          mlfRepmat(
                            mclVv(pad, "pad"),
                            mclMinus(
                              mclVv(brows, "brows"), mclVv(arows, "arows")),
                            mclVv(acols, "acols")),
                          NULL),
                        mclVa(b, "b"),
                        NULL));
                /*
                 * else
                 */
                } else {
                    /*
                     * % a and b are the same size so just concatenate
                     * c = [a b];
                     */
                    mlfAssign(
                      &c, mlfHorzcat(mclVa(a, "a"), mclVa(b, "b"), NULL));
                /*
                 * end
                 */
                }
            /*
             * else % if(Args.Columnwise)
             */
            } else {
                /*
                 * if(acols>bcols)
                 */
                if (mclGtBool(mclVv(acols, "acols"), mclVv(bcols, "bcols"))) {
                    /*
                     * % pad b before concatenating
                     * c = [a; [b repmat(pad,brows,acols-bcols)]];
                     */
                    mlfAssign(
                      &c,
                      mlfVertcat(
                        mclVa(a, "a"),
                        mlfHorzcat(
                          mclVa(b, "b"),
                          mlfRepmat(
                            mclVv(pad, "pad"),
                            mclVv(brows, "brows"),
                            mclMinus(
                              mclVv(acols, "acols"), mclVv(bcols, "bcols"))),
                          NULL),
                        NULL));
                /*
                 * elseif(bcols>acols)
                 */
                } else if (mclGtBool(
                             mclVv(bcols, "bcols"), mclVv(acols, "acols"))) {
                    /*
                     * % pad a before concatenating
                     * c = [[a repmat(pad,arows,bcols-acols)]; b];
                     */
                    mlfAssign(
                      &c,
                      mlfVertcat(
                        mlfHorzcat(
                          mclVa(a, "a"),
                          mlfRepmat(
                            mclVv(pad, "pad"),
                            mclVv(arows, "arows"),
                            mclMinus(
                              mclVv(bcols, "bcols"), mclVv(acols, "acols"))),
                          NULL),
                        mclVa(b, "b"),
                        NULL));
                /*
                 * else
                 */
                } else {
                    /*
                     * % a and b are the same size so just concatenate
                     * c = [a; b];
                     */
                    mlfAssign(
                      &c, mlfVertcat(mclVa(a, "a"), mclVa(b, "b"), NULL));
                /*
                 * end
                 */
                }
            /*
             * end % if(Args.Columnwise)
             */
            }
        /*
         * end % if(b0)
         */
        }
    /*
     * end % if(a0)
     */
    }
    mclValidateOutput(c, 1, nargout_, "c", "concatenate");
    mxDestroyArray(Args);
    mxDestroyArray(pad);
    mxDestroyArray(arows);
    mxDestroyArray(acols);
    mxDestroyArray(brows);
    mxDestroyArray(bcols);
    mxDestroyArray(a0);
    mxDestroyArray(b0);
    mxDestroyArray(varargin);
    mxDestroyArray(b);
    mxDestroyArray(a);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return c;
}
