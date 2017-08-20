/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:39 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "nptFileParts.h"
#include "libmatlbm.h"
#include "libmmfile.h"

static mxChar _array1_[4] = { 'M', 'A', 'C', '2' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;

void InitializeModule_nptFileParts(void) {
    _mxarray0_ = mclInitializeString(4, _array1_);
    _mxarray2_ = mclInitializeDouble(2.0);
    _mxarray3_ = mclInitializeDouble(1.0);
}

void TerminateModule_nptFileParts(void) {
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * MnptFileParts(mxArray * * name,
                               mxArray * * ext,
                               mxArray * * ver,
                               int nargout_,
                               mxArray * file);

_mexLocalFunctionTable _local_function_table_nptFileParts
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfNptFileParts" contains the normal interface for the
 * "nptFileParts" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptFilePart
 * s.m" (lines 1-19). This function processes any input arguments and passes
 * them to the implementation version of the function, appearing above.
 */
mxArray * mlfNptFileParts(mxArray * * name,
                          mxArray * * ext,
                          mxArray * * ver,
                          mxArray * file) {
    int nargout = 1;
    mxArray * path = NULL;
    mxArray * name__ = NULL;
    mxArray * ext__ = NULL;
    mxArray * ver__ = NULL;
    mlfEnterNewContext(3, 1, name, ext, ver, file);
    if (name != NULL) {
        ++nargout;
    }
    if (ext != NULL) {
        ++nargout;
    }
    if (ver != NULL) {
        ++nargout;
    }
    path = MnptFileParts(&name__, &ext__, &ver__, nargout, file);
    mlfRestorePreviousContext(3, 1, name, ext, ver, file);
    if (name != NULL) {
        mclCopyOutputArg(name, name__);
    } else {
        mxDestroyArray(name__);
    }
    if (ext != NULL) {
        mclCopyOutputArg(ext, ext__);
    } else {
        mxDestroyArray(ext__);
    }
    if (ver != NULL) {
        mclCopyOutputArg(ver, ver__);
    } else {
        mxDestroyArray(ver__);
    }
    return mlfReturnValue(path);
}

/*
 * The function "mlxNptFileParts" contains the feval interface for the
 * "nptFileParts" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptFilePart
 * s.m" (lines 1-19). The feval function calls the implementation version of
 * nptFileParts through this function. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
void mlxNptFileParts(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[4];
    int i;
    if (nlhs > 4) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: nptFileParts Line: 1 Column"
            ": 1 The function \"nptFileParts\" was called with"
            " more than the declared number of outputs (4)."),
          NULL);
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: nptFileParts Line: 1 Column"
            ": 1 The function \"nptFileParts\" was called with"
            " more than the declared number of inputs (1)."),
          NULL);
    }
    for (i = 0; i < 4; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = MnptFileParts(&mplhs[1], &mplhs[2], &mplhs[3], nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 4 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 4; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "MnptFileParts" is the implementation version of the
 * "nptFileParts" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptFilePart
 * s.m" (lines 1-19). It contains the actual compiled code for that M-function.
 * It is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [path,name,ext,ver] = nptFileParts(file)
 */
static mxArray * MnptFileParts(mxArray * * name,
                               mxArray * * ext,
                               mxArray * * ver,
                               int nargout_,
                               mxArray * file) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_nptFileParts);
    mxArray * path = NULL;
    mxArray * plength = NULL;
    mclCopyArray(&file);
    /*
     * %nptFileParts Platform independent version of FILEPARTS.
     * %	[PATH,NAME,EXT,VER] = nptFileParts(FILE) returns the path, filename, 
     * %	extension and version for the specified file.  VER will be non-empty
     * %	only on VMS.  The only difference from the FILEPARTS function is that
     * %	it removes the path separator from the end of the path argument when
     * %	run on Classic Mac OS.
     * %
     * %	You can reconstruct the file from the parts using
     * %	fullfile(path,[name ext ver])
     * 
     * [path,name,ext,ver] = fileparts(file);
     */
    mlfAssign(&path, mlfFileparts(name, ext, ver, mclVa(file, "file")));
    /*
     * if strcmp(computer,'MAC2')
     */
    if (mlfTobool(mlfStrcmp(mlfComputer(NULL, NULL), _mxarray0_))) {
        /*
         * if ~isempty(path)
         */
        if (mclNotBool(mlfIsempty(mclVv(path, "path")))) {
            /*
             * plength = size(path,2);
             */
            mlfAssign(
              &plength,
              mlfSize(mclValueVarargout(), mclVv(path, "path"), _mxarray2_));
            /*
             * path = path(1:plength-1);
             */
            mlfAssign(
              &path,
              mclArrayRef1(
                mclVv(path, "path"),
                mlfColon(
                  _mxarray3_,
                  mclMinus(mclVv(plength, "plength"), _mxarray3_),
                  NULL)));
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    mclValidateOutput(path, 1, nargout_, "path", "nptFileParts");
    mclValidateOutput(*name, 2, nargout_, "name", "nptFileParts");
    mclValidateOutput(*ext, 3, nargout_, "ext", "nptFileParts");
    mclValidateOutput(*ver, 4, nargout_, "ver", "nptFileParts");
    mxDestroyArray(plength);
    mxDestroyArray(file);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return path;
}
