/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "removeargs.h"
#include "libmatlbm.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;

void InitializeModule_removeargs(void) {
    _mxarray0_ = mclInitializeDouble(2.0);
    _mxarray1_ = mclInitializeDouble(1.0);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeCellVector(0, 0, (mxArray * *)NULL);
}

void TerminateModule_removeargs(void) {
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mremoveargs(mxArray * * num_args,
                             int nargout_,
                             mxArray * args,
                             mxArray * i,
                             mxArray * number);

_mexLocalFunctionTable _local_function_table_removeargs
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfRemoveargs" contains the normal interface for the
 * "removeargs" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/removeargs.
 * m" (lines 1-21). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
mxArray * mlfRemoveargs(mxArray * * num_args,
                        mxArray * args,
                        mxArray * i,
                        mxArray * number) {
    int nargout = 1;
    mxArray * rargs = NULL;
    mxArray * num_args__ = NULL;
    mlfEnterNewContext(1, 3, num_args, args, i, number);
    if (num_args != NULL) {
        ++nargout;
    }
    rargs = Mremoveargs(&num_args__, nargout, args, i, number);
    mlfRestorePreviousContext(1, 3, num_args, args, i, number);
    if (num_args != NULL) {
        mclCopyOutputArg(num_args, num_args__);
    } else {
        mxDestroyArray(num_args__);
    }
    return mlfReturnValue(rargs);
}

/*
 * The function "mlxRemoveargs" contains the feval interface for the
 * "removeargs" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/removeargs.
 * m" (lines 1-21). The feval function calls the implementation version of
 * removeargs through this function. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
void mlxRemoveargs(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: removeargs Line: 1 Column:"
            " 1 The function \"removeargs\" was called with m"
            "ore than the declared number of outputs (2)."),
          NULL);
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: removeargs Line: 1 Column"
            ": 1 The function \"removeargs\" was called with"
            " more than the declared number of inputs (3)."),
          NULL);
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 3 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 3; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    mplhs[0] = Mremoveargs(&mplhs[1], nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mremoveargs" is the implementation version of the "removeargs"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/removeargs.
 * m" (lines 1-21). It contains the actual compiled code for that M-function.
 * It is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [rargs,num_args] = removeargs(args,i,number)
 */
static mxArray * Mremoveargs(mxArray * * num_args,
                             int nargout_,
                             mxArray * args,
                             mxArray * i,
                             mxArray * number) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_removeargs);
    mxArray * rargs = NULL;
    mxArray * v2 = NULL;
    mxArray * v1 = NULL;
    mxArray * vmax = NULL;
    mxArray * vmin = NULL;
    mxArray * nargs = NULL;
    mclCopyArray(&args);
    mclCopyArray(&i);
    mclCopyArray(&number);
    /*
     * %REMOVEARGS Remove arguments from list
     * 
     * nargs = size(args,2);
     */
    mlfAssign(
      &nargs, mlfSize(mclValueVarargout(), mclVa(args, "args"), _mxarray0_));
    /*
     * 
     * vmin = i - 1;
     */
    mlfAssign(&vmin, mclMinus(mclVa(i, "i"), _mxarray1_));
    /*
     * vmax = i + number;
     */
    mlfAssign(&vmax, mclPlus(mclVa(i, "i"), mclVa(number, "number")));
    /*
     * if vmin > 0
     */
    if (mclGtBool(mclVv(vmin, "vmin"), _mxarray2_)) {
        /*
         * v1 = {args{1:vmin}};
         */
        mlfAssign(
          &v1,
          mlfCellhcat(
            mlfIndexRef(
              mclVa(args, "args"),
              "{?}",
              mlfColon(_mxarray1_, mclVv(vmin, "vmin"), NULL)),
            NULL));
    /*
     * else
     */
    } else {
        /*
         * v1 = {};
         */
        mlfAssign(&v1, _mxarray3_);
    /*
     * end
     */
    }
    /*
     * if vmax <= nargs
     */
    if (mclLeBool(mclVv(vmax, "vmax"), mclVv(nargs, "nargs"))) {
        /*
         * v2 = {args{vmax:nargs}};
         */
        mlfAssign(
          &v2,
          mlfCellhcat(
            mlfIndexRef(
              mclVa(args, "args"),
              "{?}",
              mlfColon(mclVv(vmax, "vmax"), mclVv(nargs, "nargs"), NULL)),
            NULL));
    /*
     * else
     */
    } else {
        /*
         * v2 = {};
         */
        mlfAssign(&v2, _mxarray3_);
    /*
     * end
     */
    }
    /*
     * 
     * rargs = {v1{:},v2{:}};
     */
    mlfAssign(
      &rargs,
      mlfCellhcat(
        mlfIndexRef(mclVv(v1, "v1"), "{?}", mlfCreateColonIndex()),
        mlfIndexRef(mclVv(v2, "v2"), "{?}", mlfCreateColonIndex()),
        NULL));
    /*
     * num_args = nargs - number;
     */
    mlfAssign(
      num_args, mclMinus(mclVv(nargs, "nargs"), mclVa(number, "number")));
    mclValidateOutput(rargs, 1, nargout_, "rargs", "removeargs");
    mclValidateOutput(*num_args, 2, nargout_, "num_args", "removeargs");
    mxDestroyArray(nargs);
    mxDestroyArray(vmin);
    mxDestroyArray(vmax);
    mxDestroyArray(v1);
    mxDestroyArray(v2);
    mxDestroyArray(number);
    mxDestroyArray(i);
    mxDestroyArray(args);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return rargs;
}
