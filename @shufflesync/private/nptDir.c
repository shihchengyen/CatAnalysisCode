/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "nptDir.h"
#include "getOptArgs.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "nptFileParts.h"

static mxChar _array1_[15] = { 'C', 'a', 's', 'e', 'I', 'n', 's', 'e',
                               'n', 's', 'i', 't', 'i', 'v', 'e' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;

static mxChar _array5_[5] = { 'P', 'C', 'W', 'I', 'N' };
static mxArray * _mxarray4_;

static mxChar _array7_[3] = { 'M', 'A', 'C' };
static mxArray * _mxarray6_;
static mxArray * _mxarray8_;
static mxArray * _mxarray9_;

static mxArray * _array11_[2] = { NULL /*_mxarray8_*/, NULL /*_mxarray9_*/ };
static mxArray * _mxarray10_;

static mxChar _array13_[24] = { 'T', 'o', 'o', ' ', 'm', 'a', 'n', 'y',
                                ' ', 'i', 'n', 'p', 'u', 't', ' ', 'a',
                                'r', 'g', 'u', 'm', 'e', 'n', 't', 's' };
static mxArray * _mxarray12_;
static mxArray * _mxarray14_;

static mxChar _array16_[1] = { '.' };
static mxArray * _mxarray15_;

static mxChar _array18_[2] = { '.', '.' };
static mxArray * _mxarray17_;
static mxArray * _mxarray19_;

static mxChar _array21_[1] = { '*' };
static mxArray * _mxarray20_;

static mxChar _array23_[1] = { '^' };
static mxArray * _mxarray22_;

static mxChar _array25_[1] = { '$' };
static mxArray * _mxarray24_;

static mxChar _array27_[2] = { 0x005c, '.' };
static mxArray * _mxarray26_;

static mxChar _array29_[2] = { '.', '*' };
static mxArray * _mxarray28_;

static mxChar _array31_[7] = { 'i', 's', 'e', 'm', 'p', 't', 'y' };
static mxArray * _mxarray30_;

void InitializeModule_nptDir(void) {
    _mxarray0_ = mclInitializeString(15, _array1_);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeCell(_mxarray0_);
    _mxarray4_ = mclInitializeString(5, _array5_);
    _mxarray6_ = mclInitializeString(3, _array7_);
    _mxarray8_ = mclInitializeDouble(1.0);
    _mxarray9_ = mclInitializeDouble(2.0);
    _array11_[0] = _mxarray8_;
    _array11_[1] = _mxarray9_;
    _mxarray10_ = mclInitializeCellVector(1, 2, _array11_);
    _mxarray12_ = mclInitializeString(24, _array13_);
    _mxarray14_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray15_ = mclInitializeString(1, _array16_);
    _mxarray17_ = mclInitializeString(2, _array18_);
    _mxarray19_ = mclInitializeCharVector(0, 0, (mxChar *)NULL);
    _mxarray20_ = mclInitializeString(1, _array21_);
    _mxarray22_ = mclInitializeString(1, _array23_);
    _mxarray24_ = mclInitializeString(1, _array25_);
    _mxarray26_ = mclInitializeString(2, _array27_);
    _mxarray28_ = mclInitializeString(2, _array29_);
    _mxarray30_ = mclInitializeString(7, _array31_);
}

void TerminateModule_nptDir(void) {
    mxDestroyArray(_mxarray30_);
    mxDestroyArray(_mxarray28_);
    mxDestroyArray(_mxarray26_);
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray22_);
    mxDestroyArray(_mxarray20_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * mlfNptDir_doDir(mxArray * matchstr, mxArray * Args);
static void mlxNptDir_doDir(int nlhs,
                            mxArray * plhs[],
                            int nrhs,
                            mxArray * prhs[]);
static mxArray * MnptDir(int nargout_, mxArray * varargin);
static mxArray * MnptDir_doDir(int nargout_,
                               mxArray * matchstr,
                               mxArray * Args);

static mexFunctionTableEntry local_function_table_[1]
  = { { "doDir", mlxNptDir_doDir, 2, 1, NULL } };

_mexLocalFunctionTable _local_function_table_nptDir
  = { 1, local_function_table_ };

/*
 * The function "mlfNptDir" contains the normal interface for the "nptDir"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptDir.m"
 * (lines 1-74). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfNptDir(mxArray * synthetic_varargin_argument, ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * a = NULL;
    mlfVarargin(&varargin, synthetic_varargin_argument, 1);
    mlfEnterNewContext(0, -1, varargin);
    a = MnptDir(nargout, varargin);
    mlfRestorePreviousContext(0, 0);
    mxDestroyArray(varargin);
    return mlfReturnValue(a);
}

/*
 * The function "mlxNptDir" contains the feval interface for the "nptDir"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptDir.m"
 * (lines 1-74). The feval function calls the implementation version of nptDir
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxNptDir(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: nptDir Line: 1 Column: "
            "1 The function \"nptDir\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    mlfEnterNewContext(0, 0);
    mprhs[0] = NULL;
    mlfAssign(&mprhs[0], mclCreateVararginCell(nrhs, prhs));
    mplhs[0] = MnptDir(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 0);
    plhs[0] = mplhs[0];
    mxDestroyArray(mprhs[0]);
}

/*
 * The function "mlfNptDir_doDir" contains the normal interface for the
 * "nptDir/doDir" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptDir.m"
 * (lines 74-144). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfNptDir_doDir(mxArray * matchstr, mxArray * Args) {
    int nargout = 1;
    mxArray * dirlist = NULL;
    mlfEnterNewContext(0, 2, matchstr, Args);
    dirlist = MnptDir_doDir(nargout, matchstr, Args);
    mlfRestorePreviousContext(0, 2, matchstr, Args);
    return mlfReturnValue(dirlist);
}

/*
 * The function "mlxNptDir_doDir" contains the feval interface for the
 * "nptDir/doDir" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptDir.m"
 * (lines 74-144). The feval function calls the implementation version of
 * nptDir/doDir through this function. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
static void mlxNptDir_doDir(int nlhs,
                            mxArray * plhs[],
                            int nrhs,
                            mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: nptDir/doDir Line: 74 Colum"
            "n: 1 The function \"nptDir/doDir\" was called wit"
            "h more than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: nptDir/doDir Line: 74 Colum"
            "n: 1 The function \"nptDir/doDir\" was called wit"
            "h more than the declared number of inputs (2)."),
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
    mplhs[0] = MnptDir_doDir(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "MnptDir" is the implementation version of the "nptDir"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptDir.m"
 * (lines 1-74). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function a = nptDir(varargin)
 */
static mxArray * MnptDir(int nargout_, mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_nptDir);
    int nargin_ = mclNargin(-1, varargin, NULL);
    mxArray * a = NULL;
    mxArray * i = NULL;
    mxArray * dirsize = NULL;
    mxArray * ans = NULL;
    mxArray * dirlist = NULL;
    mxArray * osnum = NULL;
    mxArray * platform = NULL;
    mxArray * ignorecase = NULL;
    mxArray * Args = NULL;
    mclCopyArray(&varargin);
    /*
     * %nptDir Platform independent version of DIR.
     * %   nptDir directory_name lists the files in the current directory. 
     * %   When used with no arguments, nptDir acts like the DIR function in 
     * %   Matlab, except it will remove entries that begin with a period, 
     * %   e.g. '.', '..', and any hidden files created on Unix operating 
     * %   systems. 
     * %
     * %   nptDir('directory_name') lists the files in a directory. Pathnames 
     * %   and wildcards may be used. By default, this function is case
     * %   sensitive on non-Windows systems. 
     * %      e.g. nptDir('*.ini')
     * %   will return a list consisting of files different from 
     * %   nptDir('*.INI').
     * %
     * %   nptDir('directory_name','CaseInsensitive') will return a list of 
     * %   case insensitive files that match 'directory_name'. This will have
     * %   no effect on Windows systems. 
     * %
     * %   D = nptDIR('directory_name') returns the results in an M-by-1
     * %   structure with the fields: 
     * %      name  -- filename
     * %      date  -- modification date
     * %      bytes -- number of bytes allocated to the file
     * %      isdir -- 1 if name is a directory and 0 if not
     * 
     * Args = struct('CaseInsensitive',0);
     */
    mlfAssign(&Args, mlfStruct(_mxarray0_, _mxarray2_, NULL));
    /*
     * Args.flags = {'CaseInsensitive'};
     */
    mlfIndexAssign(&Args, ".flags", _mxarray3_);
    /*
     * Args = getOptArgs(varargin,Args);
     */
    mlfAssign(
      &Args,
      mlfGetOptArgs(
        NULL, mclVa(varargin, "varargin"), mclVv(Args, "Args"), NULL));
    /*
     * 
     * ignorecase = 0;
     */
    mlfAssign(&ignorecase, _mxarray2_);
    /*
     * % get platform
     * platform = computer;
     */
    mlfAssign(&platform, mlfComputer(NULL, NULL));
    /*
     * if(strcmp(platform,'PCWIN'))
     */
    if (mlfTobool(mlfStrcmp(mclVv(platform, "platform"), _mxarray4_))) {
        /*
         * osnum = 0;
         */
        mlfAssign(&osnum, _mxarray2_);
    /*
     * elseif(strcmp(platform,'MAC'))
     */
    } else if (mlfTobool(mlfStrcmp(mclVv(platform, "platform"), _mxarray6_))) {
        /*
         * osnum = 1;
         */
        mlfAssign(&osnum, _mxarray8_);
    /*
     * else
     */
    } else {
        /*
         * osnum = 2;
         */
        mlfAssign(&osnum, _mxarray9_);
    /*
     * end
     */
    }
    /*
     * 
     * switch nargin
     */
    {
        mxArray * v_ = mclInitialize(mlfScalar(nargin_));
        if (mclSwitchCompare(v_, _mxarray2_)) {
            /*
             * case 0
             * % no argument so we are going to list all files
             * % directory listings without arguments work fine on Windows share
             * % directories mounted on Mac OS X so we don't need to do anything
             * % special
             * dirlist = dir;
             */
            mlfAssign(&dirlist, mlfNDir(1, NULL));
        /*
         * case {1,2}
         */
        } else if (mclSwitchCompare(v_, _mxarray10_)) {
            /*
             * if(osnum == 0)
             */
            if (mclEqBool(mclVv(osnum, "osnum"), _mxarray2_)) {
                /*
                 * % Windows machine so ignore the second argument
                 * dirlist = dir(varargin{1});
                 */
                mlfAssign(
                  &dirlist,
                  mclFeval(
                    mclValueVarargout(),
                    mlxDir,
                    mlfIndexRef(mclVa(varargin, "varargin"), "{?}", _mxarray8_),
                    NULL));
            /*
             * else
             */
            } else {
                /*
                 * % non-Windows machine so use workaround to do listing
                 * % tecnically dir works on Windows shares mounted on Linux
                 * % but we will use the workaround since it simplifies dealing
                 * % with the case-sensitivity issue
                 * dirlist = doDir(varargin{1},Args);
                 */
                mlfAssign(
                  &dirlist,
                  mclFeval(
                    mclValueVarargout(),
                    mlxNptDir_doDir,
                    mlfIndexRef(mclVa(varargin, "varargin"), "{?}", _mxarray8_),
                    mclVv(Args, "Args"),
                    NULL));
            /*
             * end
             */
            }
        /*
         * otherwise
         */
        } else {
            /*
             * error('Too many input arguments');
             */
            mlfError(_mxarray12_, NULL);
        /*
         * end
         */
        }
        mxDestroyArray(v_);
    }
    /*
     * 
     * dirsize = size(dirlist,1);
     */
    mlfAssign(
      &dirsize,
      mlfSize(mclValueVarargout(), mclVv(dirlist, "dirlist"), _mxarray8_));
    /*
     * a = [];
     */
    mlfAssign(&a, _mxarray14_);
    /*
     * for i = 1:dirsize
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclVv(dirsize, "dirsize"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray14_);
        } else {
            /*
             * if ~strcmp(dirlist(i).name(1),'.')
             * % first entry is a '.', which means that we should remove
             * % it from our list
             * a = [a; dirlist(i)];
             * end
             * end
             */
            for (; ; ) {
                if (mclNotBool(
                      mclFeval(
                        mclValueVarargout(),
                        mlxStrcmp,
                        mlfIndexRef(
                          mclVv(dirlist, "dirlist"),
                          "(?).name(?)",
                          mlfScalar(v_),
                          _mxarray8_),
                        _mxarray15_,
                        NULL))) {
                    mlfAssign(
                      &a,
                      mlfVertcat(
                        mclVv(a, "a"),
                        mclIntArrayRef1(mclVv(dirlist, "dirlist"), v_),
                        NULL));
                }
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    mclValidateOutput(a, 1, nargout_, "a", "nptDir");
    mxDestroyArray(Args);
    mxDestroyArray(ignorecase);
    mxDestroyArray(platform);
    mxDestroyArray(osnum);
    mxDestroyArray(dirlist);
    mxDestroyArray(ans);
    mxDestroyArray(dirsize);
    mxDestroyArray(i);
    mxDestroyArray(varargin);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return a;
    /*
     * 
     */
}

/*
 * The function "MnptDir_doDir" is the implementation version of the
 * "nptDir/doDir" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/nptDir.m"
 * (lines 74-144). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function dirlist = doDir(matchstr,Args)
 */
static mxArray * MnptDir_doDir(int nargout_,
                               mxArray * matchstr,
                               mxArray * Args) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_nptDir);
    mxArray * dirlist = NULL;
    mxArray * goodcells = NULL;
    mxArray * names = NULL;
    mxArray * list = NULL;
    mxArray * ans = NULL;
    mxArray * cwd = NULL;
    mxArray * e = NULL;
    mxArray * n = NULL;
    mxArray * p = NULL;
    mclCopyArray(&matchstr);
    mclCopyArray(&Args);
    /*
     * 
     * % check if matchstr consists of a path
     * [p,n,e] = nptFileParts(matchstr);
     */
    mlfAssign(&p, mlfNptFileParts(&n, &e, NULL, mclVa(matchstr, "matchstr")));
    /*
     * % check if we were supposed to just list all the files in another
     * % directory. If matchstr was something like '..' or '../..' then the
     * % matlab function fileparts used in nptFileParts returns n as '.' and 
     * % e as '.'. So we need to fix that before we continue.
     * if(strcmp(n,'.') && strcmp(e,'.'))
     */
    if (mclScalarToBool(mlfStrcmp(mclVv(n, "n"), _mxarray15_))
        && mclScalarToBool(mlfStrcmp(mclVv(e, "e"), _mxarray15_))) {
        /*
         * if(isempty(p))
         */
        if (mlfTobool(mlfIsempty(mclVv(p, "p")))) {
            /*
             * p = '..';
             */
            mlfAssign(&p, _mxarray17_);
        /*
         * else
         */
        } else {
            /*
             * p = [p filesep '..'];
             */
            mlfAssign(
              &p, mlfHorzcat(mclVv(p, "p"), mlfFilesep(), _mxarray17_, NULL));
        /*
         * end
         */
        }
        /*
         * % set n and e to empty string so we know to just return the entire dir 
         * % listing
         * n = '';
         */
        mlfAssign(&n, _mxarray19_);
        /*
         * e = '';
         */
        mlfAssign(&e, _mxarray19_);
    /*
     * end
     */
    }
    /*
     * % temporary fix for bug in MATLAB when listing directories on a
     * % mounted Windows share on MAC OS X
     * cwd = '';
     */
    mlfAssign(&cwd, _mxarray19_);
    /*
     * if(~isempty(p))
     */
    if (mclNotBool(mlfIsempty(mclVv(p, "p")))) {
        /*
         * % save current directory
         * cwd = pwd;
         */
        mlfAssign(&cwd, mlfPwd());
        /*
         * % change to path
         * cd(p)
         */
        mclPrintAns(&ans, mlfNCd(0, mclVv(p, "p")));
        /*
         * % get directory listing
         * list = dir;
         */
        mlfAssign(&list, mlfNDir(1, NULL));
        /*
         * % change back to previous directory so we don't change
         * % directories on the user without them knowing
         * cd(cwd)
         */
        mclPrintAns(&ans, mlfNCd(0, mclVv(cwd, "cwd")));
        /*
         * % remove the path prefix from matchstr
         * matchstr = [n e];
         */
        mlfAssign(&matchstr, mlfHorzcat(mclVv(n, "n"), mclVv(e, "e"), NULL));
    /*
     * else
     */
    } else {
        /*
         * % get directory listing
         * list = dir;
         */
        mlfAssign(&list, mlfNDir(1, NULL));
    /*
     * end
     */
    }
    /*
     * % If matchstr was something like '../' then n would be empty so we will
     * % just return the dir listing
     * if(isempty(n))
     */
    if (mlfTobool(mlfIsempty(mclVv(n, "n")))) {
        /*
         * dirlist = list;
         */
        mlfAssign(&dirlist, mclVv(list, "list"));
        /*
         * return
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * % put the names in a cell array;
     * names = {list.name};
     */
    mlfAssign(
      &names, mlfCellhcat(mlfIndexRef(mclVv(list, "list"), ".name"), NULL));
    /*
     * % parse matchstr to convert it to pattern that regexp
     * % recognizes
     * if(matchstr(1)~='*')
     */
    if (mclNeBool(
          mclIntArrayRef1(mclVa(matchstr, "matchstr"), 1), _mxarray20_)) {
        /*
         * % add ^ to beginning of matchstr
         * matchstr = ['^' matchstr];
         */
        mlfAssign(
          &matchstr,
          mlfHorzcat(_mxarray22_, mclVa(matchstr, "matchstr"), NULL));
    /*
     * end
     */
    }
    /*
     * if(matchstr(end)~='*')
     */
    if (mclNeBool(
          mclArrayRef1(
            mclVa(matchstr, "matchstr"),
            mlfEnd(mclVa(matchstr, "matchstr"), _mxarray8_, _mxarray8_)),
          _mxarray20_)) {
        /*
         * % add $ to end of matchstr
         * matchstr = [matchstr '$'];
         */
        mlfAssign(
          &matchstr,
          mlfHorzcat(mclVa(matchstr, "matchstr"), _mxarray24_, NULL));
    /*
     * end
     */
    }
    /*
     * % replace . with \.
     * matchstr = strrep(matchstr,'.','\.');
     */
    mlfAssign(
      &matchstr,
      mlfStrrep(mclVa(matchstr, "matchstr"), _mxarray15_, _mxarray26_));
    /*
     * % replace * with .*
     * matchstr = strrep(matchstr,'*','.*');
     */
    mlfAssign(
      &matchstr,
      mlfStrrep(mclVa(matchstr, "matchstr"), _mxarray20_, _mxarray28_));
    /*
     * if(Args.CaseInsensitive)
     */
    if (mlfTobool(mlfIndexRef(mclVa(Args, "Args"), ".CaseInsensitive"))) {
        /*
         * % find which cells have contents that match matchstr regardless
         * % of case
         * goodcells = regexpi(names,matchstr);
         */
        mlfAssign(
          &goodcells,
          mlfRegexpi(
            NULL,
            NULL,
            mclVv(names, "names"),
            mclVa(matchstr, "matchstr"),
            NULL,
            NULL,
            NULL,
            NULL));
    /*
     * else
     */
    } else {
        /*
         * % find which cells have contents that match matchstr
         * goodcells = regexp(names,matchstr);
         */
        mlfAssign(
          &goodcells,
          mlfRegexp(
            NULL,
            NULL,
            mclVv(names, "names"),
            mclVa(matchstr, "matchstr"),
            NULL,
            NULL,
            NULL,
            NULL));
    /*
     * end
     */
    }
    /*
     * % extract those cell indices from the original list
     * dirlist = list(~cellfun('isempty',goodcells));
     */
    mlfAssign(
      &dirlist,
      mclArrayRef1(
        mclVv(list, "list"),
        mclNot(
          mlfNCellfun(
            0,
            mclValueVarargout(),
            _mxarray30_,
            mclVv(goodcells, "goodcells"),
            NULL))));
    return_:
    mclValidateOutput(dirlist, 1, nargout_, "dirlist", "nptDir/doDir");
    mxDestroyArray(p);
    mxDestroyArray(n);
    mxDestroyArray(e);
    mxDestroyArray(cwd);
    mxDestroyArray(ans);
    mxDestroyArray(list);
    mxDestroyArray(names);
    mxDestroyArray(goodcells);
    mxDestroyArray(Args);
    mxDestroyArray(matchstr);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return dirlist;
}
