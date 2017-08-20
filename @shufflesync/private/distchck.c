/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "distchck.h"
#include "libmatlbm.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;

void InitializeModule_distchck(void) {
    _mxarray0_ = mclInitializeDouble(0.0);
    _mxarray1_ = mclInitializeDouble(1.0);
    _mxarray2_ = mclInitializeDouble(2.0);
    _mxarray3_ = mclInitializeDouble(3.0);
    _mxarray4_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray5_ = mclInitializeDouble(4.0);
}

void TerminateModule_distchck(void) {
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mdistchck(mxArray * * out1,
                           mxArray * * out2,
                           mxArray * * out3,
                           mxArray * * out4,
                           int nargout_,
                           mxArray * nparms,
                           mxArray * arg1,
                           mxArray * arg2,
                           mxArray * arg3,
                           mxArray * arg4);

_mexLocalFunctionTable _local_function_table_distchck
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfDistchck" contains the normal interface for the "distchck"
 * M-function from file "/Users/syen/Documents/ShihCheng/Matlab/distchck.m"
 * (lines 1-174). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
mxArray * mlfDistchck(mxArray * * out1,
                      mxArray * * out2,
                      mxArray * * out3,
                      mxArray * * out4,
                      mxArray * nparms,
                      mxArray * arg1,
                      mxArray * arg2,
                      mxArray * arg3,
                      mxArray * arg4) {
    int nargout = 1;
    mxArray * errorcode = NULL;
    mxArray * out1__ = NULL;
    mxArray * out2__ = NULL;
    mxArray * out3__ = NULL;
    mxArray * out4__ = NULL;
    mlfEnterNewContext(
      4, 5, out1, out2, out3, out4, nparms, arg1, arg2, arg3, arg4);
    if (out1 != NULL) {
        ++nargout;
    }
    if (out2 != NULL) {
        ++nargout;
    }
    if (out3 != NULL) {
        ++nargout;
    }
    if (out4 != NULL) {
        ++nargout;
    }
    errorcode
      = Mdistchck(
          &out1__,
          &out2__,
          &out3__,
          &out4__,
          nargout,
          nparms,
          arg1,
          arg2,
          arg3,
          arg4);
    mlfRestorePreviousContext(
      4, 5, out1, out2, out3, out4, nparms, arg1, arg2, arg3, arg4);
    if (out1 != NULL) {
        mclCopyOutputArg(out1, out1__);
    } else {
        mxDestroyArray(out1__);
    }
    if (out2 != NULL) {
        mclCopyOutputArg(out2, out2__);
    } else {
        mxDestroyArray(out2__);
    }
    if (out3 != NULL) {
        mclCopyOutputArg(out3, out3__);
    } else {
        mxDestroyArray(out3__);
    }
    if (out4 != NULL) {
        mclCopyOutputArg(out4, out4__);
    } else {
        mxDestroyArray(out4__);
    }
    return mlfReturnValue(errorcode);
}

/*
 * The function "mlxDistchck" contains the feval interface for the "distchck"
 * M-function from file "/Users/syen/Documents/ShihCheng/Matlab/distchck.m"
 * (lines 1-174). The feval function calls the implementation version of
 * distchck through this function. This function processes any input arguments
 * and passes them to the implementation version of the function, appearing
 * above.
 */
void mlxDistchck(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[5];
    mxArray * mplhs[5];
    int i;
    if (nlhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: distchck Line: 1 Column:"
            " 1 The function \"distchck\" was called with m"
            "ore than the declared number of outputs (5)."),
          NULL);
    }
    if (nrhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: distchck Line: 1 Column:"
            " 1 The function \"distchck\" was called with m"
            "ore than the declared number of inputs (5)."),
          NULL);
    }
    for (i = 0; i < 5; ++i) {
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
      = Mdistchck(
          &mplhs[1],
          &mplhs[2],
          &mplhs[3],
          &mplhs[4],
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2],
          mprhs[3],
          mprhs[4]);
    mlfRestorePreviousContext(
      0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 5 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 5; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mdistchck" is the implementation version of the "distchck"
 * M-function from file "/Users/syen/Documents/ShihCheng/Matlab/distchck.m"
 * (lines 1-174). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [errorcode,out1,out2,out3,out4] = distchck(nparms,arg1,arg2,arg3,arg4)
 */
static mxArray * Mdistchck(mxArray * * out1,
                           mxArray * * out2,
                           mxArray * * out3,
                           mxArray * * out4,
                           int nargout_,
                           mxArray * nparms,
                           mxArray * arg1,
                           mxArray * arg2,
                           mxArray * arg3,
                           mxArray * arg4) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_distchck);
    mxArray * errorcode = NULL;
    mxArray * scalararg4 = NULL;
    mxArray * c4 = NULL;
    mxArray * r4 = NULL;
    mxArray * columns = NULL;
    mxArray * rows = NULL;
    mxArray * scalararg3 = NULL;
    mxArray * c3 = NULL;
    mxArray * r3 = NULL;
    mxArray * scalararg2 = NULL;
    mxArray * scalararg1 = NULL;
    mxArray * c2 = NULL;
    mxArray * r2 = NULL;
    mxArray * c1 = NULL;
    mxArray * r1 = NULL;
    mclCopyArray(&nparms);
    mclCopyArray(&arg1);
    mclCopyArray(&arg2);
    mclCopyArray(&arg3);
    mclCopyArray(&arg4);
    /*
     * %DISTCHCK Checks the argument list for the probability functions.
     * 
     * %   B.A. Jones  1-22-93
     * %   Copyright (c) 1993-97 by The MathWorks, Inc.
     * %   $Revision: 2.5 $  $Date: 1997/04/08 15:02:06 $
     * 
     * errorcode = 0;
     */
    mlfAssign(&errorcode, _mxarray0_);
    /*
     * 
     * if nparms == 1
     */
    if (mclEqBool(mclVa(nparms, "nparms"), _mxarray1_)) {
        /*
         * out1 = arg1;
         */
        mlfAssign(out1, mclVa(arg1, "arg1"));
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
     * if nparms == 2
     */
    if (mclEqBool(mclVa(nparms, "nparms"), _mxarray2_)) {
        /*
         * [r1 c1] = size(arg1);
         */
        mlfSize(mlfVarargout(&r1, &c1, NULL), mclVa(arg1, "arg1"), NULL);
        /*
         * [r2 c2] = size(arg2);
         */
        mlfSize(mlfVarargout(&r2, &c2, NULL), mclVa(arg2, "arg2"), NULL);
        /*
         * scalararg1 = (prod(size(arg1)) == 1);
         */
        mlfAssign(
          &scalararg1,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg1, "arg1"), NULL), NULL),
            _mxarray1_));
        /*
         * scalararg2 = (prod(size(arg2)) == 1);
         */
        mlfAssign(
          &scalararg2,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg2, "arg2"), NULL), NULL),
            _mxarray1_));
        /*
         * if ~scalararg1 & ~scalararg2
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg1, "scalararg1")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg2, "scalararg2"))))) {
                mxDestroyArray(a_);
                /*
                 * if r1 ~= r2 | c1 ~= c2
                 */
                {
                    mxArray * a_0
                      = mclInitialize(mclNe(mclVv(r1, "r1"), mclVv(r2, "r2")));
                    if (mlfTobool(a_0)
                        || mlfTobool(
                             mclOr(
                               a_0,
                               mclNe(mclVv(c1, "c1"), mclVv(c2, "c2"))))) {
                        mxDestroyArray(a_0);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;         
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_0);
                    }
                /*
                 * end     
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * if scalararg1
         */
        if (mlfTobool(mclVv(scalararg1, "scalararg1"))) {
            /*
             * out1 = arg1(ones(r2,1),ones(c2,1));
             */
            mlfAssign(
              out1,
              mclArrayRef2(
                mclVa(arg1, "arg1"),
                mlfOnes(mclVv(r2, "r2"), _mxarray1_, NULL),
                mlfOnes(mclVv(c2, "c2"), _mxarray1_, NULL)));
        /*
         * else
         */
        } else {
            /*
             * out1 = arg1;
             */
            mlfAssign(out1, mclVa(arg1, "arg1"));
        /*
         * end
         */
        }
        /*
         * if scalararg2
         */
        if (mlfTobool(mclVv(scalararg2, "scalararg2"))) {
            /*
             * out2 = arg2(ones(r1,1),ones(c1,1));
             */
            mlfAssign(
              out2,
              mclArrayRef2(
                mclVa(arg2, "arg2"),
                mlfOnes(mclVv(r1, "r1"), _mxarray1_, NULL),
                mlfOnes(mclVv(c1, "c1"), _mxarray1_, NULL)));
        /*
         * else
         */
        } else {
            /*
             * out2 = arg2;
             */
            mlfAssign(out2, mclVa(arg2, "arg2"));
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
     * if nparms == 3
     */
    if (mclEqBool(mclVa(nparms, "nparms"), _mxarray3_)) {
        /*
         * [r1 c1] = size(arg1);
         */
        mlfSize(mlfVarargout(&r1, &c1, NULL), mclVa(arg1, "arg1"), NULL);
        /*
         * [r2 c2] = size(arg2);
         */
        mlfSize(mlfVarargout(&r2, &c2, NULL), mclVa(arg2, "arg2"), NULL);
        /*
         * [r3 c3] = size(arg3);
         */
        mlfSize(mlfVarargout(&r3, &c3, NULL), mclVa(arg3, "arg3"), NULL);
        /*
         * scalararg1 = (prod(size(arg1)) == 1);
         */
        mlfAssign(
          &scalararg1,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg1, "arg1"), NULL), NULL),
            _mxarray1_));
        /*
         * scalararg2 = (prod(size(arg2)) == 1);
         */
        mlfAssign(
          &scalararg2,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg2, "arg2"), NULL), NULL),
            _mxarray1_));
        /*
         * scalararg3 = (prod(size(arg3)) == 1);
         */
        mlfAssign(
          &scalararg3,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg3, "arg3"), NULL), NULL),
            _mxarray1_));
        /*
         * 
         * if ~scalararg1 & ~scalararg2
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg1, "scalararg1")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg2, "scalararg2"))))) {
                mxDestroyArray(a_);
                /*
                 * if r1 ~= r2 | c1 ~= c2
                 */
                {
                    mxArray * a_1
                      = mclInitialize(mclNe(mclVv(r1, "r1"), mclVv(r2, "r2")));
                    if (mlfTobool(a_1)
                        || mlfTobool(
                             mclOr(
                               a_1,
                               mclNe(mclVv(c1, "c1"), mclVv(c2, "c2"))))) {
                        mxDestroyArray(a_1);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;         
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_1);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg1 & ~scalararg3
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg1, "scalararg1")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg3, "scalararg3"))))) {
                mxDestroyArray(a_);
                /*
                 * if r1 ~= r3 | c1 ~= c3
                 */
                {
                    mxArray * a_2
                      = mclInitialize(mclNe(mclVv(r1, "r1"), mclVv(r3, "r3")));
                    if (mlfTobool(a_2)
                        || mlfTobool(
                             mclOr(
                               a_2,
                               mclNe(mclVv(c1, "c1"), mclVv(c3, "c3"))))) {
                        mxDestroyArray(a_2);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;                 
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_2);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg3 & ~scalararg2
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg3, "scalararg3")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg2, "scalararg2"))))) {
                mxDestroyArray(a_);
                /*
                 * if r3 ~= r2 | c3 ~= c2
                 */
                {
                    mxArray * a_3
                      = mclInitialize(mclNe(mclVv(r3, "r3"), mclVv(r2, "r2")));
                    if (mlfTobool(a_3)
                        || mlfTobool(
                             mclOr(
                               a_3,
                               mclNe(mclVv(c3, "c3"), mclVv(c2, "c2"))))) {
                        mxDestroyArray(a_3);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;         
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_3);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg1
         */
        if (mclNotBool(mclVv(scalararg1, "scalararg1"))) {
            /*
             * out1 = arg1;
             */
            mlfAssign(out1, mclVa(arg1, "arg1"));
        /*
         * end
         */
        }
        /*
         * if ~scalararg2
         */
        if (mclNotBool(mclVv(scalararg2, "scalararg2"))) {
            /*
             * out2 = arg2;
             */
            mlfAssign(out2, mclVa(arg2, "arg2"));
        /*
         * end
         */
        }
        /*
         * if ~scalararg3
         */
        if (mclNotBool(mclVv(scalararg3, "scalararg3"))) {
            /*
             * out3 = arg3;
             */
            mlfAssign(out3, mclVa(arg3, "arg3"));
        /*
         * end
         */
        }
        /*
         * rows = max([r1 r2 r3]);
         */
        mlfAssign(
          &rows,
          mlfMax(
            NULL,
            mlfHorzcat(mclVv(r1, "r1"), mclVv(r2, "r2"), mclVv(r3, "r3"), NULL),
            NULL,
            NULL));
        /*
         * columns = max([c1 c2 c3]);  
         */
        mlfAssign(
          &columns,
          mlfMax(
            NULL,
            mlfHorzcat(mclVv(c1, "c1"), mclVv(c2, "c2"), mclVv(c3, "c3"), NULL),
            NULL,
            NULL));
        /*
         * 
         * if scalararg1
         */
        if (mlfTobool(mclVv(scalararg1, "scalararg1"))) {
            /*
             * out1 = arg1(ones(rows,1),ones(columns,1));
             */
            mlfAssign(
              out1,
              mclArrayRef2(
                mclVa(arg1, "arg1"),
                mlfOnes(mclVv(rows, "rows"), _mxarray1_, NULL),
                mlfOnes(mclVv(columns, "columns"), _mxarray1_, NULL)));
        /*
         * end
         */
        }
        /*
         * if scalararg2
         */
        if (mlfTobool(mclVv(scalararg2, "scalararg2"))) {
            /*
             * out2 = arg2(ones(rows,1),ones(columns,1));
             */
            mlfAssign(
              out2,
              mclArrayRef2(
                mclVa(arg2, "arg2"),
                mlfOnes(mclVv(rows, "rows"), _mxarray1_, NULL),
                mlfOnes(mclVv(columns, "columns"), _mxarray1_, NULL)));
        /*
         * end
         */
        }
        /*
         * if scalararg3
         */
        if (mlfTobool(mclVv(scalararg3, "scalararg3"))) {
            /*
             * out3 = arg3(ones(rows,1),ones(columns,1));
             */
            mlfAssign(
              out3,
              mclArrayRef2(
                mclVa(arg3, "arg3"),
                mlfOnes(mclVv(rows, "rows"), _mxarray1_, NULL),
                mlfOnes(mclVv(columns, "columns"), _mxarray1_, NULL)));
        /*
         * end
         */
        }
        /*
         * out4 =[];
         */
        mlfAssign(out4, _mxarray4_);
    /*
     * 
     * end
     */
    }
    /*
     * 
     * if nparms == 4
     */
    if (mclEqBool(mclVa(nparms, "nparms"), _mxarray5_)) {
        /*
         * [r1 c1] = size(arg1);
         */
        mlfSize(mlfVarargout(&r1, &c1, NULL), mclVa(arg1, "arg1"), NULL);
        /*
         * [r2 c2] = size(arg2);
         */
        mlfSize(mlfVarargout(&r2, &c2, NULL), mclVa(arg2, "arg2"), NULL);
        /*
         * [r3 c3] = size(arg3);
         */
        mlfSize(mlfVarargout(&r3, &c3, NULL), mclVa(arg3, "arg3"), NULL);
        /*
         * [r4 c4] = size(arg4);
         */
        mlfSize(mlfVarargout(&r4, &c4, NULL), mclVa(arg4, "arg4"), NULL);
        /*
         * scalararg1 = (prod(size(arg1)) == 1);
         */
        mlfAssign(
          &scalararg1,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg1, "arg1"), NULL), NULL),
            _mxarray1_));
        /*
         * scalararg2 = (prod(size(arg2)) == 1);
         */
        mlfAssign(
          &scalararg2,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg2, "arg2"), NULL), NULL),
            _mxarray1_));
        /*
         * scalararg3 = (prod(size(arg3)) == 1);
         */
        mlfAssign(
          &scalararg3,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg3, "arg3"), NULL), NULL),
            _mxarray1_));
        /*
         * scalararg4 = (prod(size(arg4)) == 1);
         */
        mlfAssign(
          &scalararg4,
          mclEq(
            mlfProd(
              mlfSize(mclValueVarargout(), mclVa(arg4, "arg4"), NULL), NULL),
            _mxarray1_));
        /*
         * 
         * if ~scalararg1 & ~scalararg2
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg1, "scalararg1")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg2, "scalararg2"))))) {
                mxDestroyArray(a_);
                /*
                 * if r1 ~= r2 | c1 ~= c2
                 */
                {
                    mxArray * a_4
                      = mclInitialize(mclNe(mclVv(r1, "r1"), mclVv(r2, "r2")));
                    if (mlfTobool(a_4)
                        || mlfTobool(
                             mclOr(
                               a_4,
                               mclNe(mclVv(c1, "c1"), mclVv(c2, "c2"))))) {
                        mxDestroyArray(a_4);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;         
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_4);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg1 & ~scalararg3
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg1, "scalararg1")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg3, "scalararg3"))))) {
                mxDestroyArray(a_);
                /*
                 * if r1 ~= r3 | c1 ~= c3
                 */
                {
                    mxArray * a_5
                      = mclInitialize(mclNe(mclVv(r1, "r1"), mclVv(r3, "r3")));
                    if (mlfTobool(a_5)
                        || mlfTobool(
                             mclOr(
                               a_5,
                               mclNe(mclVv(c1, "c1"), mclVv(c3, "c3"))))) {
                        mxDestroyArray(a_5);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;                 
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_5);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg1 & ~scalararg4
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg1, "scalararg1")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg4, "scalararg4"))))) {
                mxDestroyArray(a_);
                /*
                 * if r1 ~= r4 | c1 ~= c4
                 */
                {
                    mxArray * a_6
                      = mclInitialize(mclNe(mclVv(r1, "r1"), mclVv(r4, "r4")));
                    if (mlfTobool(a_6)
                        || mlfTobool(
                             mclOr(
                               a_6,
                               mclNe(mclVv(c1, "c1"), mclVv(c4, "c4"))))) {
                        mxDestroyArray(a_6);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;                 
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_6);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg3 & ~scalararg2
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg3, "scalararg3")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg2, "scalararg2"))))) {
                mxDestroyArray(a_);
                /*
                 * if r3 ~= r2 | c3 ~= c2
                 */
                {
                    mxArray * a_7
                      = mclInitialize(mclNe(mclVv(r3, "r3"), mclVv(r2, "r2")));
                    if (mlfTobool(a_7)
                        || mlfTobool(
                             mclOr(
                               a_7,
                               mclNe(mclVv(c3, "c3"), mclVv(c2, "c2"))))) {
                        mxDestroyArray(a_7);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;         
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_7);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg4 & ~scalararg2
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg4, "scalararg4")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg2, "scalararg2"))))) {
                mxDestroyArray(a_);
                /*
                 * if r4 ~= r2 | c4 ~= c2
                 */
                {
                    mxArray * a_8
                      = mclInitialize(mclNe(mclVv(r4, "r4"), mclVv(r2, "r2")));
                    if (mlfTobool(a_8)
                        || mlfTobool(
                             mclOr(
                               a_8,
                               mclNe(mclVv(c4, "c4"), mclVv(c2, "c2"))))) {
                        mxDestroyArray(a_8);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;         
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_8);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * if ~scalararg3 & ~scalararg4
         */
        {
            mxArray * a_
              = mclInitialize(mclNot(mclVv(scalararg3, "scalararg3")));
            if (mlfTobool(a_)
                && mlfTobool(
                     mclAnd(a_, mclNot(mclVv(scalararg4, "scalararg4"))))) {
                mxDestroyArray(a_);
                /*
                 * if r3 ~= r4 | c3 ~= c4
                 */
                {
                    mxArray * a_9
                      = mclInitialize(mclNe(mclVv(r3, "r3"), mclVv(r4, "r4")));
                    if (mlfTobool(a_9)
                        || mlfTobool(
                             mclOr(
                               a_9,
                               mclNe(mclVv(c3, "c3"), mclVv(c4, "c4"))))) {
                        mxDestroyArray(a_9);
                        /*
                         * errorcode = 1;
                         */
                        mlfAssign(&errorcode, _mxarray1_);
                        /*
                         * return;         
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_9);
                    }
                /*
                 * end
                 */
                }
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * 
         * 
         * if ~scalararg1
         */
        if (mclNotBool(mclVv(scalararg1, "scalararg1"))) {
            /*
             * out1 = arg1;
             */
            mlfAssign(out1, mclVa(arg1, "arg1"));
        /*
         * end
         */
        }
        /*
         * if ~scalararg2
         */
        if (mclNotBool(mclVv(scalararg2, "scalararg2"))) {
            /*
             * out2 = arg2;
             */
            mlfAssign(out2, mclVa(arg2, "arg2"));
        /*
         * end
         */
        }
        /*
         * if ~scalararg3
         */
        if (mclNotBool(mclVv(scalararg3, "scalararg3"))) {
            /*
             * out3 = arg3;
             */
            mlfAssign(out3, mclVa(arg3, "arg3"));
        /*
         * end
         */
        }
        /*
         * if ~scalararg4
         */
        if (mclNotBool(mclVv(scalararg4, "scalararg4"))) {
            /*
             * out4 = arg4;
             */
            mlfAssign(out4, mclVa(arg4, "arg4"));
        /*
         * end
         */
        }
        /*
         * 
         * rows = max([r1 r2 r3 r4]);
         */
        mlfAssign(
          &rows,
          mlfMax(
            NULL,
            mlfHorzcat(
              mclVv(r1, "r1"),
              mclVv(r2, "r2"),
              mclVv(r3, "r3"),
              mclVv(r4, "r4"),
              NULL),
            NULL,
            NULL));
        /*
         * columns = max([c1 c2 c3 c4]);     
         */
        mlfAssign(
          &columns,
          mlfMax(
            NULL,
            mlfHorzcat(
              mclVv(c1, "c1"),
              mclVv(c2, "c2"),
              mclVv(c3, "c3"),
              mclVv(c4, "c4"),
              NULL),
            NULL,
            NULL));
        /*
         * if scalararg1
         */
        if (mlfTobool(mclVv(scalararg1, "scalararg1"))) {
            /*
             * out1 = arg1(ones(rows,1),ones(columns,1));
             */
            mlfAssign(
              out1,
              mclArrayRef2(
                mclVa(arg1, "arg1"),
                mlfOnes(mclVv(rows, "rows"), _mxarray1_, NULL),
                mlfOnes(mclVv(columns, "columns"), _mxarray1_, NULL)));
        /*
         * end
         */
        }
        /*
         * if scalararg2
         */
        if (mlfTobool(mclVv(scalararg2, "scalararg2"))) {
            /*
             * out2 = arg2(ones(rows,1),ones(columns,1));
             */
            mlfAssign(
              out2,
              mclArrayRef2(
                mclVa(arg2, "arg2"),
                mlfOnes(mclVv(rows, "rows"), _mxarray1_, NULL),
                mlfOnes(mclVv(columns, "columns"), _mxarray1_, NULL)));
        /*
         * end
         */
        }
        /*
         * if scalararg3
         */
        if (mlfTobool(mclVv(scalararg3, "scalararg3"))) {
            /*
             * out3 = arg3(ones(rows,1),ones(columns,1));
             */
            mlfAssign(
              out3,
              mclArrayRef2(
                mclVa(arg3, "arg3"),
                mlfOnes(mclVv(rows, "rows"), _mxarray1_, NULL),
                mlfOnes(mclVv(columns, "columns"), _mxarray1_, NULL)));
        /*
         * end
         */
        }
        /*
         * if scalararg4
         */
        if (mlfTobool(mclVv(scalararg4, "scalararg4"))) {
            /*
             * out4 = arg4(ones(rows,1),ones(columns,1));
             */
            mlfAssign(
              out4,
              mclArrayRef2(
                mclVa(arg4, "arg4"),
                mlfOnes(mclVv(rows, "rows"), _mxarray1_, NULL),
                mlfOnes(mclVv(columns, "columns"), _mxarray1_, NULL)));
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
     */
    return_:
    mclValidateOutput(errorcode, 1, nargout_, "errorcode", "distchck");
    mclValidateOutput(*out1, 2, nargout_, "out1", "distchck");
    mclValidateOutput(*out2, 3, nargout_, "out2", "distchck");
    mclValidateOutput(*out3, 4, nargout_, "out3", "distchck");
    mclValidateOutput(*out4, 5, nargout_, "out4", "distchck");
    mxDestroyArray(r1);
    mxDestroyArray(c1);
    mxDestroyArray(r2);
    mxDestroyArray(c2);
    mxDestroyArray(scalararg1);
    mxDestroyArray(scalararg2);
    mxDestroyArray(r3);
    mxDestroyArray(c3);
    mxDestroyArray(scalararg3);
    mxDestroyArray(rows);
    mxDestroyArray(columns);
    mxDestroyArray(r4);
    mxDestroyArray(c4);
    mxDestroyArray(scalararg4);
    mxDestroyArray(arg4);
    mxDestroyArray(arg3);
    mxDestroyArray(arg2);
    mxDestroyArray(arg1);
    mxDestroyArray(nparms);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return errorcode;
}
