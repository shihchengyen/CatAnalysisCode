/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "getOptArgs.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "removeargs.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;

static mxChar _array5_[5] = { 'f', 'l', 'a', 'g', 's' };
static mxArray * _mxarray4_;
static mxArray * _mxarray6_;

static mxChar _array8_[7] = { 'a', 'l', 'i', 'a', 's', 'e', 's' };
static mxArray * _mxarray7_;

static mxChar _array10_[6] = { 'r', 'e', 'm', 'o', 'v', 'e' };
static mxArray * _mxarray9_;

static mxChar _array12_[8] = { 's', 'u', 'b', 't', 'r', 'a', 'c', 't' };
static mxArray * _mxarray11_;

static mxChar _array14_[9] = { 's', 'h', 'o', 'r', 't', 'c', 'u', 't', 's' };
static mxArray * _mxarray13_;

static mxChar _array16_[11] = { 's', 't', 'o', 'p', 'o', 'n',
                                'e', 'r', 'r', 'o', 'r' };
static mxArray * _mxarray15_;
static mxArray * _mxarray17_;
static mxArray * _mxarray18_;

static mxChar _array20_[1] = { ' ' };
static mxArray * _mxarray19_;
static mxArray * _mxarray21_;

static mxChar _array23_[5] = { 'e', 'x', 'a', 'c', 't' };
static mxArray * _mxarray22_;

static mxChar _array25_[25] = { 'U', 'n', 'k', 'n', 'o', 'w', 'n', ' ', 'n',
                                'a', 'm', 'e', 'd', ' ', 'p', 'a', 'r', 'a',
                                'm', 'e', 't', 'e', 'r', ':', ' ' };
static mxArray * _mxarray24_;

static mxChar _array27_[28] = { 'E', 'x', 'p', 'e', 'c', 't', 'e',
                                'd', ' ', 'a', ' ', 'n', 'a', 'm',
                                'e', 'd', ' ', 'p', 'a', 'r', 'a',
                                'm', 'e', 't', 'e', 'r', ':', ' ' };
static mxArray * _mxarray26_;

void InitializeModule_getOptArgs(void) {
    _mxarray0_ = mclInitializeCellVector(0, 0, (mxArray * *)NULL);
    _mxarray1_ = mclInitializeCharVector(0, 0, (mxChar *)NULL);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeDouble(1.0);
    _mxarray4_ = mclInitializeString(5, _array5_);
    _mxarray6_ = mclInitializeDouble(2.0);
    _mxarray7_ = mclInitializeString(7, _array8_);
    _mxarray9_ = mclInitializeString(6, _array10_);
    _mxarray11_ = mclInitializeString(8, _array12_);
    _mxarray13_ = mclInitializeString(9, _array14_);
    _mxarray15_ = mclInitializeString(11, _array16_);
    _mxarray17_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray18_ = mclInitializeDouble(3.0);
    _mxarray19_ = mclInitializeString(1, _array20_);
    _mxarray21_ = mclInitializeDouble(4.0);
    _mxarray22_ = mclInitializeString(5, _array23_);
    _mxarray24_ = mclInitializeString(25, _array25_);
    _mxarray26_ = mclInitializeString(28, _array27_);
}

void TerminateModule_getOptArgs(void) {
    mxDestroyArray(_mxarray26_);
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray22_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * MgetOptArgs(mxArray * * args,
                             int nargout_,
                             mxArray * args_in,
                             mxArray * ArgStruct_in,
                             mxArray * varargin);

_mexLocalFunctionTable _local_function_table_getOptArgs
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfGetOptArgs" contains the normal interface for the
 * "getOptArgs" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/getOptArgs.
 * m" (lines 1-261). This function processes any input arguments and passes
 * them to the implementation version of the function, appearing above.
 */
mxArray * mlfGetOptArgs(mxArray * * args,
                        mxArray * args_in,
                        mxArray * ArgStruct_in,
                        ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * ArgStruct = NULL;
    mxArray * args__ = NULL;
    mlfVarargin(&varargin, ArgStruct_in, 0);
    mlfEnterNewContext(1, -3, args, args_in, ArgStruct_in, varargin);
    if (args != NULL) {
        ++nargout;
    }
    ArgStruct = MgetOptArgs(&args__, nargout, args_in, ArgStruct_in, varargin);
    mlfRestorePreviousContext(1, 2, args, args_in, ArgStruct_in);
    mxDestroyArray(varargin);
    if (args != NULL) {
        mclCopyOutputArg(args, args__);
    } else {
        mxDestroyArray(args__);
    }
    return mlfReturnValue(ArgStruct);
}

/*
 * The function "mlxGetOptArgs" contains the feval interface for the
 * "getOptArgs" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/getOptArgs.
 * m" (lines 1-261). The feval function calls the implementation version of
 * getOptArgs through this function. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
void mlxGetOptArgs(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: getOptArgs Line: 1 Column:"
            " 1 The function \"getOptArgs\" was called with m"
            "ore than the declared number of outputs (2)."),
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
    mplhs[0] = MgetOptArgs(&mplhs[1], nlhs, mprhs[0], mprhs[1], mprhs[2]);
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
 * The function "MgetOptArgs" is the implementation version of the "getOptArgs"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/getOptArgs.
 * m" (lines 1-261). It contains the actual compiled code for that M-function.
 * It is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function [ArgStruct,args]=getOptArgs(args,ArgStruct,varargin)
 */
static mxArray * MgetOptArgs(mxArray * * args,
                             int nargout_,
                             mxArray * args_in,
                             mxArray * ArgStruct_in,
                             mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_getOptArgs);
    int nargin_ = mclNargin(-3, args_in, ArgStruct_in, varargin, NULL);
    mxArray * ArgStruct = NULL;
    mxArray * val = NULL;
    mxArray * curField = NULL;
    mxArray * ans = NULL;
    mxArray * ShortCutIdx = NULL;
    mxArray * a = NULL;
    mxArray * al = NULL;
    mxArray * l = NULL;
    mxArray * idx = NULL;
    mxArray * j = NULL;
    mxArray * FieldIdx = NULL;
    mxArray * naliases = NULL;
    mxArray * nAliases = NULL;
    mxArray * FnamesAbbr = NULL;
    mxArray * FnamesFull = NULL;
    mxArray * AbbrevIdx = NULL;
    mxArray * name = NULL;
    mxArray * Fnames = NULL;
    mxArray * NumArgCount = NULL;
    mxArray * numargs = NULL;
    mxArray * ShortCutParams = NULL;
    mxArray * arg = NULL;
    mxArray * i = NULL;
    mxArray * nArgs = NULL;
    mxArray * bError = NULL;
    mxArray * ShortCuts = NULL;
    mxArray * ShortcutParams = NULL;
    mxArray * SubtractParams = NULL;
    mxArray * RemoveParams = NULL;
    mxArray * FlagTypeParams = NULL;
    mxArray * Aliases = NULL;
    mclCopyInputArg(args, args_in);
    mclCopyInputArg(&ArgStruct, ArgStruct_in);
    mclCopyArray(&varargin);
    /*
     * % getOptArgs Helper function for parsing varargin. 
     * %    [ArgStruct,ARGS] = getOptArgs(ARGS,ArgStruct,VARARGIN) parses 
     * %    ARGS for parameter-value pairs and enters them into ArgStruct,
     * %    which is a structure containing named arguments with default
     * %    values. Any numeric arguments at the beginning of ARGS that are
     * %    not matched to a parameter string are stored in the field 
     * %    NumericArguments of ArgStruct. The ARGS cell array, which may 
     * %    be modified via optional input arguments, is also returned.
     * %
     * %    The optional input arguments are:
     * %    'flags' - followed by a cell array indicating which arguments are 
     * %              flagtype arguments, i.e. arguments that don't require a 
     * %              value (the value will be set to 1 if it is present). If
     * %              this argument is not specified, but there is a field in
     * %              ArgStruct named 'flags', that will be used instead. If
     * %              this argument is present, it will be added to ArgStruct
     * %              if it is not already present, and will replace the flags
     * %              field if present.
     * %    'aliases' - followed by a cell array indicating aliases, which
     * %                can be used to map one argument-name to several 
     * %                argstruct fields.
     * %    'remove' - followed by a cell array indicating which arguments 
     * %               in ARGS are to be removed if present.
     * %    'subtract' - following by a cell array indicating which arguments
     * %                 in ARGS will have their values subtracted by 1.
     * %    'shortcuts' - followed by a cell array with the name of the 
     * %                  shortcut in the first column and a cell array in the
     * %                  second column with the full arguments.
     * %    'stopOnError' - if present, indicates that the function should
     * %                    stop execution when it finds unknown arguments.
     * %    example usage: 
     * %    --------------
     * %    function parseargtest(varargin)
     * %
     * %    %define the acceptable named arguments and assign default values
     * %    Args=struct('Holdaxis',0, ...
     * %           'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
     * %           'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0, ...
     * %           'PaddingBottom',0,'MarginLeft',.1,'MarginRight',.1, ...
     * %           'MarginTop',.1,'MarginBottom',.1,'rows',[],'cols',[], ...
     * %           'RedoLevels',0,'SaveLevels',0,'Dir','','flags',{'Holdaxis'}); 
     * %
     * %    % The capital letters define abrreviations.  
     * %    % Eg. parseargtest('spacingvertical',0) is equivalent to  
     * %    % parseargtest('sv',0) 
     * %    % fill the arg-struct with values entered by the user
     * %    [Args,varargin]=getOptArgs(varargin,Args,'flags',{'Holdaxis'}, ... 
     * %        'aliases',{'Spacing' {'sh','sv'}; ...
     * %                   'Padding' {'pl','pr','pt','pb'}; ...
     * %                   'Margin' {'ml','mr','mt','mb'}}, ...
     * %        'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
     * %        'subtract',{'RedoLevels','SaveLevels'},...
     * %        'remove',{'Dir'}, ...
     * %        'stopOnError');
     * %
     * %    disp(Args)
     * 
     * % Based on parseArgs.m from Matlab Central by Aslak Grinsted 2003
     * 
     * Aliases={};
     */
    mlfAssign(&Aliases, _mxarray0_);
    /*
     * FlagTypeParams='';
     */
    mlfAssign(&FlagTypeParams, _mxarray1_);
    /*
     * RemoveParams='';
     */
    mlfAssign(&RemoveParams, _mxarray1_);
    /*
     * SubtractParams='';
     */
    mlfAssign(&SubtractParams, _mxarray1_);
    /*
     * ShortcutParams='';
     */
    mlfAssign(&ShortcutParams, _mxarray1_);
    /*
     * ShortCuts='';
     */
    mlfAssign(&ShortCuts, _mxarray1_);
    /*
     * % flag to indicate if we want to stop execution with an error or 
     * % continue
     * bError = 0;
     */
    mlfAssign(&bError, _mxarray2_);
    /*
     * 
     * % get length of varargin
     * nArgs = nargin - 2;
     */
    mlfAssign(&nArgs, mlfScalar(nargin_ - 2));
    /*
     * i = 1;
     */
    mlfAssign(&i, _mxarray3_);
    /*
     * 
     * % check if ArgStruct has a field called flags
     * if(isfield(ArgStruct,'flags'))
     */
    if (mlfTobool(mlfIsfield(mclVa(ArgStruct, "ArgStruct"), _mxarray4_))) {
        /*
         * FlagTypeParams = lower(ArgStruct.flags);
         */
        mlfAssign(
          &FlagTypeParams,
          mclFeval(
            mclValueVarargout(),
            mlxLower,
            mlfIndexRef(mclVa(ArgStruct, "ArgStruct"), ".flags"),
            NULL));
    /*
     * end
     */
    }
    /*
     * 
     * % look for optional arguments
     * while(i <= nArgs)
     */
    while (mclLeBool(mclVv(i, "i"), mclVv(nArgs, "nArgs"))) {
        /*
         * arg = varargin{i};	
         */
        mlfAssign(
          &arg, mlfIndexRef(mclVa(varargin, "varargin"), "{?}", mclVv(i, "i")));
        /*
         * if ischar(arg)
         */
        if (mlfTobool(mlfIschar(mclVv(arg, "arg")))) {
            /*
             * arg = lower(arg);
             */
            mlfAssign(&arg, mlfLower(mclVv(arg, "arg")));
            /*
             * switch arg
             */
            {
                mxArray * v_ = mclInitialize(mclVv(arg, "arg"));
                if (mclSwitchCompare(v_, _mxarray4_)) {
                    /*
                     * case('flags')
                     * % make argument case insensitive
                     * FlagTypeParams = lower(varargin{i+1});
                     */
                    mlfAssign(
                      &FlagTypeParams,
                      mclFeval(
                        mclValueVarargout(),
                        mlxLower,
                        mlfIndexRef(
                          mclVa(varargin, "varargin"),
                          "{?}",
                          mclPlus(mclVv(i, "i"), _mxarray3_)),
                        NULL));
                    /*
                     * i = i + 2;
                     */
                    mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray6_));
                    /*
                     * % set flags field in ArgStruct
                     * ArgStruct.flags = FlagTypeParams;
                     */
                    mlfIndexAssign(
                      &ArgStruct,
                      ".flags",
                      mclVv(FlagTypeParams, "FlagTypeParams"));
                /*
                 * case('aliases')
                 */
                } else if (mclSwitchCompare(v_, _mxarray7_)) {
                    /*
                     * Aliases = varargin{i+1};
                     */
                    mlfAssign(
                      &Aliases,
                      mlfIndexRef(
                        mclVa(varargin, "varargin"),
                        "{?}",
                        mclPlus(mclVv(i, "i"), _mxarray3_)));
                    /*
                     * i = i + 2;
                     */
                    mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray6_));
                /*
                 * case('remove')
                 */
                } else if (mclSwitchCompare(v_, _mxarray9_)) {
                    /*
                     * RemoveParams = lower(varargin{i+1});
                     */
                    mlfAssign(
                      &RemoveParams,
                      mclFeval(
                        mclValueVarargout(),
                        mlxLower,
                        mlfIndexRef(
                          mclVa(varargin, "varargin"),
                          "{?}",
                          mclPlus(mclVv(i, "i"), _mxarray3_)),
                        NULL));
                    /*
                     * i = i + 2;
                     */
                    mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray6_));
                /*
                 * case('subtract')
                 */
                } else if (mclSwitchCompare(v_, _mxarray11_)) {
                    /*
                     * SubtractParams = lower(varargin{i+1});
                     */
                    mlfAssign(
                      &SubtractParams,
                      mclFeval(
                        mclValueVarargout(),
                        mlxLower,
                        mlfIndexRef(
                          mclVa(varargin, "varargin"),
                          "{?}",
                          mclPlus(mclVv(i, "i"), _mxarray3_)),
                        NULL));
                    /*
                     * i = i + 2;
                     */
                    mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray6_));
                /*
                 * case('shortcuts')
                 */
                } else if (mclSwitchCompare(v_, _mxarray13_)) {
                    /*
                     * ShortCutParams = varargin{i+1};
                     */
                    mlfAssign(
                      &ShortCutParams,
                      mlfIndexRef(
                        mclVa(varargin, "varargin"),
                        "{?}",
                        mclPlus(mclVv(i, "i"), _mxarray3_)));
                    /*
                     * ShortCuts = lower(strvcat(ShortCutParams{:,1}));
                     */
                    mlfAssign(
                      &ShortCuts,
                      mlfLower(
                        mlfStrvcat(
                          mlfIndexRef(
                            mclVv(ShortCutParams, "ShortCutParams"),
                            "{?,?}",
                            mlfCreateColonIndex(),
                            _mxarray3_),
                          NULL)));
                    /*
                     * i = i + 2;
                     */
                    mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray6_));
                /*
                 * case('stoponerror')
                 */
                } else if (mclSwitchCompare(v_, _mxarray15_)) {
                    /*
                     * bError = 1;
                     */
                    mlfAssign(&bError, _mxarray3_);
                    /*
                     * i = i + 1;
                     */
                    mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray3_));
                /*
                 * otherwise
                 */
                } else {
                    /*
                     * i = i + 1;
                     */
                    mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray3_));
                /*
                 * end
                 */
                }
                mxDestroyArray(v_);
            }
        /*
         * else
         */
        } else {
            /*
             * i = i + 1;
             */
            mlfAssign(&i, mclPlus(mclVv(i, "i"), _mxarray3_));
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
     * %
     * % Skip "numeric" arguments preceeding first param,value pair
     * %
     * % get number of arguments in args
     * numargs = size(args,2);
     */
    mlfAssign(
      &numargs, mlfSize(mclValueVarargout(), mclVa(*args, "args"), _mxarray6_));
    /*
     * % start at 1 instead of 0 and subtract 1 later so we can actually 
     * % index into args
     * NumArgCount=1;
     */
    mlfAssign(&NumArgCount, _mxarray3_);
    /*
     * while (NumArgCount<=numargs)&(~ischar(args{NumArgCount}))
     */
    for (;;) {
        mxArray * a_
          = mclInitialize(
              mclLe(
                mclVv(NumArgCount, "NumArgCount"), mclVv(numargs, "numargs")));
        if (mlfTobool(a_)
            && mlfTobool(
                 mclAnd(
                   a_,
                   mclNot(
                     mclFeval(
                       mclValueVarargout(),
                       mlxIschar,
                       mlfIndexRef(
                         mclVa(*args, "args"),
                         "{?}",
                         mclVv(NumArgCount, "NumArgCount")),
                       NULL))))) {
            mxDestroyArray(a_);
        } else {
            mxDestroyArray(a_);
            break;
        }
        /*
         * NumArgCount=NumArgCount+1;
         */
        mlfAssign(
          &NumArgCount, mclPlus(mclVv(NumArgCount, "NumArgCount"), _mxarray3_));
    /*
     * end
     */
    }
    /*
     * if(NumArgCount>1)
     */
    if (mclGtBool(mclVv(NumArgCount, "NumArgCount"), _mxarray3_)) {
        /*
         * ArgStruct.NumericArguments = {args{1:(NumArgCount-1)}};
         */
        mlfIndexAssign(
          &ArgStruct,
          ".NumericArguments",
          mlfCellhcat(
            mlfIndexRef(
              mclVa(*args, "args"),
              "{?}",
              mlfColon(
                _mxarray3_,
                mclMinus(mclVv(NumArgCount, "NumArgCount"), _mxarray3_),
                NULL)),
            NULL));
    /*
     * else
     */
    } else {
        /*
         * ArgStruct.NumericArguments = [];
         */
        mlfIndexAssign(&ArgStruct, ".NumericArguments", _mxarray17_);
    /*
     * end
     */
    }
    /*
     * 
     * %
     * % Make an accepted fieldname matrix (case insensitive)
     * %
     * Fnames=fieldnames(ArgStruct);
     */
    mlfAssign(&Fnames, mlfFieldnames(mclVa(ArgStruct, "ArgStruct")));
    /*
     * for i=1:length(Fnames)
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclLengthInt(mclVv(Fnames, "Fnames"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray17_);
        } else {
            /*
             * name=lower(Fnames{i,1});
             * Fnames{i,2}=name; %col2=lower
             * % find characters that are uppercase
             * AbbrevIdx=find(Fnames{i,1}~=name);
             * % store uppercase letters as abbreviations
             * %col3=abreviation letters (those that are uppercase in the ArgStruct) 
             * % e.g. SpacingHoriz->sh
             * Fnames{i,3}=[name(AbbrevIdx) ' '];
             * %the space prevents strvcat from removing empty lines
             * %Does this parameter have a value? (e.g. not flagtype)
             * Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams,'exact'));
             * end
             */
            for (; ; ) {
                mlfAssign(
                  &name,
                  mclFeval(
                    mclValueVarargout(),
                    mlxLower,
                    mlfIndexRef(
                      mclVv(Fnames, "Fnames"),
                      "{?,?}",
                      mlfScalar(v_),
                      _mxarray3_),
                    NULL));
                mlfIndexAssign(
                  &Fnames,
                  "{?,?}",
                  mlfScalar(v_),
                  _mxarray6_,
                  mclVv(name, "name"));
                mlfAssign(
                  &AbbrevIdx,
                  mlfFind(
                    NULL,
                    NULL,
                    mclFeval(
                      mclValueVarargout(),
                      mlxNe,
                      mlfIndexRef(
                        mclVv(Fnames, "Fnames"),
                        "{?,?}",
                        mlfScalar(v_),
                        _mxarray3_),
                      mclVv(name, "name"),
                      NULL)));
                mlfIndexAssign(
                  &Fnames,
                  "{?,?}",
                  mlfScalar(v_),
                  _mxarray18_,
                  mlfHorzcat(
                    mclArrayRef1(
                      mclVv(name, "name"), mclVv(AbbrevIdx, "AbbrevIdx")),
                    _mxarray19_,
                    NULL));
                mlfIndexAssign(
                  &Fnames,
                  "{?,?}",
                  mlfScalar(v_),
                  _mxarray21_,
                  mlfIsempty(
                    mclFeval(
                      mclValueVarargout(),
                      mlxStrmatch,
                      mlfIndexRef(
                        mclVv(Fnames, "Fnames"),
                        "{?,?}",
                        mlfScalar(v_),
                        _mxarray6_),
                      mclVv(FlagTypeParams, "FlagTypeParams"),
                      _mxarray22_,
                      NULL)));
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    /*
     * % convert to character matrices with 1 string in each row
     * FnamesFull=strvcat(Fnames{:,2});
     */
    mlfAssign(
      &FnamesFull,
      mlfStrvcat(
        mlfIndexRef(
          mclVv(Fnames, "Fnames"), "{?,?}", mlfCreateColonIndex(), _mxarray6_),
        NULL));
    /*
     * FnamesAbbr=strvcat(Fnames{:,3});
     */
    mlfAssign(
      &FnamesAbbr,
      mlfStrvcat(
        mlfIndexRef(
          mclVv(Fnames, "Fnames"), "{?,?}", mlfCreateColonIndex(), _mxarray18_),
        NULL));
    /*
     * 
     * nAliases = size(Aliases,1);
     */
    mlfAssign(
      &nAliases,
      mlfSize(mclValueVarargout(), mclVv(Aliases, "Aliases"), _mxarray3_));
    /*
     * if nAliases>0  
     */
    if (mclGtBool(mclVv(nAliases, "nAliases"), _mxarray2_)) {
        /*
         * for i=1:nAliases
         */
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclVv(nAliases, "nAliases"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray17_);
        } else {
            /*
             * name=lower(Aliases{i,2});
             * % get length of aliases
             * naliases = size(Aliases{i,2},2);
             * % find indices corresponding to aliases
             * FieldIdx = [];
             * for j = 1:naliases
             * idx = strmatch(name{j},FnamesAbbr,'exact');
             * if isempty(idx)
             * idx = strmatch(name{j},FnamesFull,'exact');
             * end
             * if ~isempty(idx)
             * FieldIdx = [FieldIdx; idx];
             * end
             * end
             * Aliases{i,2}=FieldIdx;
             * name = lower(Aliases{i,1});
             * % find characters that are uppercase
             * AbbrevIdx=find(Aliases{i,1}~=name);
             * % the space prevents strvcat from removing empty lines
             * Aliases{i,3}=[name(AbbrevIdx) ' '];
             * %dont need the name in uppercase anymore for aliases
             * Aliases{i,1}=name;
             * end
             */
            for (; ; ) {
                mlfAssign(
                  &name,
                  mclFeval(
                    mclValueVarargout(),
                    mlxLower,
                    mlfIndexRef(
                      mclVv(Aliases, "Aliases"),
                      "{?,?}",
                      mlfScalar(v_),
                      _mxarray6_),
                    NULL));
                mlfAssign(
                  &naliases,
                  mclFeval(
                    mclValueVarargout(),
                    mlxSize,
                    mlfIndexRef(
                      mclVv(Aliases, "Aliases"),
                      "{?,?}",
                      mlfScalar(v_),
                      _mxarray6_),
                    _mxarray6_,
                    NULL));
                mlfAssign(&FieldIdx, _mxarray17_);
                {
                    int v_0 = mclForIntStart(1);
                    int e_0 = mclForIntEnd(mclVv(naliases, "naliases"));
                    if (v_0 > e_0) {
                        mlfAssign(&j, _mxarray17_);
                    } else {
                        for (; ; ) {
                            mlfAssign(
                              &idx,
                              mclFeval(
                                mclValueVarargout(),
                                mlxStrmatch,
                                mlfIndexRef(
                                  mclVv(name, "name"), "{?}", mlfScalar(v_0)),
                                mclVv(FnamesAbbr, "FnamesAbbr"),
                                _mxarray22_,
                                NULL));
                            if (mlfTobool(mlfIsempty(mclVv(idx, "idx")))) {
                                mlfAssign(
                                  &idx,
                                  mclFeval(
                                    mclValueVarargout(),
                                    mlxStrmatch,
                                    mlfIndexRef(
                                      mclVv(name, "name"),
                                      "{?}",
                                      mlfScalar(v_0)),
                                    mclVv(FnamesFull, "FnamesFull"),
                                    _mxarray22_,
                                    NULL));
                            }
                            if (mclNotBool(mlfIsempty(mclVv(idx, "idx")))) {
                                mlfAssign(
                                  &FieldIdx,
                                  mlfVertcat(
                                    mclVv(FieldIdx, "FieldIdx"),
                                    mclVv(idx, "idx"),
                                    NULL));
                            }
                            if (v_0 == e_0) {
                                break;
                            }
                            ++v_0;
                        }
                        mlfAssign(&j, mlfScalar(v_0));
                    }
                }
                mlfIndexAssign(
                  &Aliases,
                  "{?,?}",
                  mlfScalar(v_),
                  _mxarray6_,
                  mclVv(FieldIdx, "FieldIdx"));
                mlfAssign(
                  &name,
                  mclFeval(
                    mclValueVarargout(),
                    mlxLower,
                    mlfIndexRef(
                      mclVv(Aliases, "Aliases"),
                      "{?,?}",
                      mlfScalar(v_),
                      _mxarray3_),
                    NULL));
                mlfAssign(
                  &AbbrevIdx,
                  mlfFind(
                    NULL,
                    NULL,
                    mclFeval(
                      mclValueVarargout(),
                      mlxNe,
                      mlfIndexRef(
                        mclVv(Aliases, "Aliases"),
                        "{?,?}",
                        mlfScalar(v_),
                        _mxarray3_),
                      mclVv(name, "name"),
                      NULL)));
                mlfIndexAssign(
                  &Aliases,
                  "{?,?}",
                  mlfScalar(v_),
                  _mxarray18_,
                  mlfHorzcat(
                    mclArrayRef1(
                      mclVv(name, "name"), mclVv(AbbrevIdx, "AbbrevIdx")),
                    _mxarray19_,
                    NULL));
                mlfIndexAssign(
                  &Aliases,
                  "{?,?}",
                  mlfScalar(v_),
                  _mxarray3_,
                  mclVv(name, "name"));
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
        /*
         * %Append aliases to the end of FnamesFull and FnamesAbbr
         * FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); 
         */
        mlfAssign(
          &FnamesFull,
          mlfStrvcat(
            mclVv(FnamesFull, "FnamesFull"),
            mlfStrvcat(
              mlfIndexRef(
                mclVv(Aliases, "Aliases"),
                "{?,?}",
                mlfCreateColonIndex(),
                _mxarray3_),
              NULL),
            NULL));
        /*
         * FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3}));
         */
        mlfAssign(
          &FnamesAbbr,
          mlfStrvcat(
            mclVv(FnamesAbbr, "FnamesAbbr"),
            mlfStrvcat(
              mlfIndexRef(
                mclVv(Aliases, "Aliases"),
                "{?,?}",
                mlfCreateColonIndex(),
                _mxarray18_),
              NULL),
            NULL));
    /*
     * end
     */
    }
    /*
     * 
     * %--------------get parameters--------------------
     * l = NumArgCount;
     */
    mlfAssign(&l, mclVv(NumArgCount, "NumArgCount"));
    /*
     * al = length(args);
     */
    mlfAssign(&al, mlfScalar(mclLengthInt(mclVa(*args, "args"))));
    /*
     * while (l<=al)
     */
    while (mclLeBool(mclVv(l, "l"), mclVv(al, "al"))) {
        /*
         * a=args{l};
         */
        mlfAssign(&a, mlfIndexRef(mclVa(*args, "args"), "{?}", mclVv(l, "l")));
        /*
         * if ischar(a) 
         */
        if (mlfTobool(mlfIschar(mclVv(a, "a")))) {
            /*
             * if ~isempty(a)
             */
            if (mclNotBool(mlfIsempty(mclVv(a, "a")))) {
                /*
                 * a=lower(a);
                 */
                mlfAssign(&a, mlfLower(mclVv(a, "a")));
                /*
                 * %try abbreviations (must be exact)
                 * FieldIdx=strmatch(a,FnamesAbbr,'exact'); 
                 */
                mlfAssign(
                  &FieldIdx,
                  mlfStrmatch(
                    mclVv(a, "a"),
                    mclVv(FnamesAbbr, "FnamesAbbr"),
                    _mxarray22_));
                /*
                 * if isempty(FieldIdx) 
                 */
                if (mlfTobool(mlfIsempty(mclVv(FieldIdx, "FieldIdx")))) {
                    /*
                     * FieldIdx=strmatch(a,FnamesFull,'exact'); 
                     */
                    mlfAssign(
                      &FieldIdx,
                      mlfStrmatch(
                        mclVv(a, "a"),
                        mclVv(FnamesFull, "FnamesFull"),
                        _mxarray22_));
                /*
                 * end
                 */
                }
                /*
                 * if FieldIdx>length(Fnames) %then it's an alias type.
                 */
                if (mclGtBool(
                      mclVv(FieldIdx, "FieldIdx"),
                      mlfScalar(mclLengthInt(mclVv(Fnames, "Fnames"))))) {
                    /*
                     * FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
                     */
                    mlfAssign(
                      &FieldIdx,
                      mlfIndexRef(
                        mclVv(Aliases, "Aliases"),
                        "{?,?}",
                        mclMinus(
                          mclVv(FieldIdx, "FieldIdx"),
                          mlfScalar(mclLengthInt(mclVv(Fnames, "Fnames")))),
                        _mxarray6_));
                /*
                 * end
                 */
                }
                /*
                 * 
                 * if isempty(FieldIdx) 
                 */
                if (mlfTobool(mlfIsempty(mclVv(FieldIdx, "FieldIdx")))) {
                    /*
                     * % check if it is a shortcut
                     * ShortCutIdx = strmatch(a,ShortCuts,'exact');
                     */
                    mlfAssign(
                      &ShortCutIdx,
                      mlfStrmatch(
                        mclVv(a, "a"),
                        mclVv(ShortCuts, "ShortCuts"),
                        _mxarray22_));
                    /*
                     * if(~isempty(ShortCutIdx))
                     */
                    if (mclNotBool(
                          mlfIsempty(mclVv(ShortCutIdx, "ShortCutIdx")))) {
                        /*
                         * % replace shortcut with full string in args
                         * if(l==al)
                         */
                        if (mclEqBool(mclVv(l, "l"), mclVv(al, "al"))) {
                            /*
                             * % if we are already at the end of the cell array
                             * args = {args{1:(l-1)} ShortCutParams{ShortCutIdx,2}{:}};
                             */
                            mlfAssign(
                              args,
                              mlfCellhcat(
                                mlfIndexRef(
                                  mclVa(*args, "args"),
                                  "{?}",
                                  mlfColon(
                                    _mxarray3_,
                                    mclMinus(mclVv(l, "l"), _mxarray3_),
                                    NULL)),
                                mlfIndexRef(
                                  mclVv(ShortCutParams, "ShortCutParams"),
                                  "{?,?}{?}",
                                  mclVv(ShortCutIdx, "ShortCutIdx"),
                                  _mxarray6_,
                                  mlfCreateColonIndex()),
                                NULL));
                        /*
                         * elseif(l==1)
                         */
                        } else if (mclEqBool(mclVv(l, "l"), _mxarray3_)) {
                            /*
                             * % if we are at the beginning of the cell array
                             * args = {ShortCutParams{ShortCutIdx,2}{:} args{l+1:end}};
                             */
                            mlfAssign(
                              args,
                              mlfCellhcat(
                                mlfIndexRef(
                                  mclVv(ShortCutParams, "ShortCutParams"),
                                  "{?,?}{?}",
                                  mclVv(ShortCutIdx, "ShortCutIdx"),
                                  _mxarray6_,
                                  mlfCreateColonIndex()),
                                mlfIndexRef(
                                  mclVa(*args, "args"),
                                  "{?}",
                                  mlfColon(
                                    mclPlus(mclVv(l, "l"), _mxarray3_),
                                    mlfEnd(
                                      mclVa(*args, "args"),
                                      _mxarray3_,
                                      _mxarray3_),
                                    NULL)),
                                NULL));
                        /*
                         * else
                         */
                        } else {
                            /*
                             * args = {args{1:(l-1)} ShortCutParams{ShortCutIdx,2}{:} args{l+1:end}};
                             */
                            mlfAssign(
                              args,
                              mlfCellhcat(
                                mlfIndexRef(
                                  mclVa(*args, "args"),
                                  "{?}",
                                  mlfColon(
                                    _mxarray3_,
                                    mclMinus(mclVv(l, "l"), _mxarray3_),
                                    NULL)),
                                mlfIndexRef(
                                  mclVv(ShortCutParams, "ShortCutParams"),
                                  "{?,?}{?}",
                                  mclVv(ShortCutIdx, "ShortCutIdx"),
                                  _mxarray6_,
                                  mlfCreateColonIndex()),
                                mlfIndexRef(
                                  mclVa(*args, "args"),
                                  "{?}",
                                  mlfColon(
                                    mclPlus(mclVv(l, "l"), _mxarray3_),
                                    mlfEnd(
                                      mclVa(*args, "args"),
                                      _mxarray3_,
                                      _mxarray3_),
                                    NULL)),
                                NULL));
                        /*
                         * end
                         */
                        }
                        /*
                         * al = length(args);
                         */
                        mlfAssign(
                          &al, mlfScalar(mclLengthInt(mclVa(*args, "args"))));
                    /*
                     * elseif bError
                     */
                    } else if (mlfTobool(mclVv(bError, "bError"))) {
                        /*
                         * error(['Unknown named parameter: ' a])
                         */
                        mlfError(
                          mlfHorzcat(_mxarray24_, mclVv(a, "a"), NULL), NULL);
                    /*
                     * else
                     */
                    } else {
                        /*
                         * l = l + 1;
                         */
                        mlfAssign(&l, mclPlus(mclVv(l, "l"), _mxarray3_));
                    /*
                     * end
                     */
                    }
                /*
                 * else
                 */
                } else {
                    mclForLoopIterator viter__;
                    /*
                     * for curField=FieldIdx' %if it is an alias it could be more than one.
                     */
                    for (mclForStart(
                           &viter__,
                           mlfCtranspose(mclVv(FieldIdx, "FieldIdx")),
                           NULL,
                           NULL);
                         mclForNext(&viter__, &curField);
                         ) {
                        /*
                         * if (Fnames{curField,4})
                         */
                        if (mlfTobool(
                              mlfIndexRef(
                                mclVv(Fnames, "Fnames"),
                                "{?,?}",
                                mclVv(curField, "curField"),
                                _mxarray21_))) {
                            /*
                             * val=args{l+1};
                             */
                            mlfAssign(
                              &val,
                              mlfIndexRef(
                                mclVa(*args, "args"),
                                "{?}",
                                mclPlus(mclVv(l, "l"), _mxarray3_)));
                        /*
                         * else
                         */
                        } else {
                            /*
                             * %parameter is of flag type and is set (1=true)....
                             * val=1; 
                             */
                            mlfAssign(&val, _mxarray3_);
                        /*
                         * end
                         */
                        }
                        /*
                         * ArgStruct.(Fnames{curField,1})=val;
                         */
                        mlfIndexAssign(
                          &ArgStruct,
                          ".(?)",
                          mlfIndexRef(
                            mclVv(Fnames, "Fnames"),
                            "{?,?}",
                            mclVv(curField, "curField"),
                            _mxarray3_),
                          mclVv(val, "val"));
                        /*
                         * if(strmatch(a,SubtractParams,'exact'))
                         */
                        if (mlfTobool(
                              mlfStrmatch(
                                mclVv(a, "a"),
                                mclVv(SubtractParams, "SubtractParams"),
                                _mxarray22_))) {
                            /*
                             * args{l+1} = args{l+1} - 1;
                             */
                            mlfIndexAssign(
                              args,
                              "{?}",
                              mclPlus(mclVv(l, "l"), _mxarray3_),
                              mclFeval(
                                mclValueVarargout(),
                                mlxMinus,
                                mlfIndexRef(
                                  mclVa(*args, "args"),
                                  "{?}",
                                  mclPlus(mclVv(l, "l"), _mxarray3_)),
                                _mxarray3_,
                                NULL));
                        /*
                         * end
                         */
                        }
                    /*
                     * end
                     */
                    }
                    mclDestroyForLoopIterator(viter__);
                    /*
                     * if(strmatch(a,RemoveParams,'exact'))
                     */
                    if (mlfTobool(
                          mlfStrmatch(
                            mclVv(a, "a"),
                            mclVv(RemoveParams, "RemoveParams"),
                            _mxarray22_))) {
                        /*
                         * if(Fnames{FieldIdx(1),4})
                         */
                        if (mlfTobool(
                              mlfIndexRef(
                                mclVv(Fnames, "Fnames"),
                                "{?,?}",
                                mclIntArrayRef1(mclVv(FieldIdx, "FieldIdx"), 1),
                                _mxarray21_))) {
                            /*
                             * % parameter with value so remove both argument and 
                             * % value
                             * [args,al] = removeargs(args,l,2);
                             */
                            mlfAssign(
                              args,
                              mlfRemoveargs(
                                &al,
                                mclVa(*args, "args"),
                                mclVv(l, "l"),
                                _mxarray6_));
                        /*
                         * else
                         */
                        } else {
                            /*
                             * % flag type so just remove argument
                             * [args,al] = removeargs(args,l,1);
                             */
                            mlfAssign(
                              args,
                              mlfRemoveargs(
                                &al,
                                mclVa(*args, "args"),
                                mclVv(l, "l"),
                                _mxarray3_));
                        /*
                         * end
                         */
                        }
                    /*
                     * else
                     */
                    } else {
                        /*
                         * % if flag type then go to next argument but if 
                         * % param/value, then skip next argument
                         * l=l+1+Fnames{FieldIdx(1),4}; 
                         */
                        mlfAssign(
                          &l,
                          mclFeval(
                            mclValueVarargout(),
                            mlxPlus,
                            mclPlus(mclVv(l, "l"), _mxarray3_),
                            mlfIndexRef(
                              mclVv(Fnames, "Fnames"),
                              "{?,?}",
                              mclIntArrayRef1(mclVv(FieldIdx, "FieldIdx"), 1),
                              _mxarray21_),
                            NULL));
                    /*
                     * end
                     */
                    }
                /*
                 * end
                 */
                }
            /*
             * else 
             */
            } else {
                /*
                 * l=l+1;
                 */
                mlfAssign(&l, mclPlus(mclVv(l, "l"), _mxarray3_));
            /*
             * end
             */
            }
        /*
         * else
         */
        } else {
            /*
             * if bError
             */
            if (mlfTobool(mclVv(bError, "bError"))) {
                /*
                 * error(['Expected a named parameter: ' num2str(a)])
                 */
                mlfError(
                  mlfHorzcat(
                    _mxarray26_, mlfNum2str(mclVv(a, "a"), NULL), NULL),
                  NULL);
            /*
             * else
             */
            } else {
                /*
                 * l = l + 1;
                 */
                mlfAssign(&l, mclPlus(mclVv(l, "l"), _mxarray3_));
            /*
             * end
             */
            }
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    mclValidateOutput(ArgStruct, 1, nargout_, "ArgStruct", "getOptArgs");
    mclValidateOutput(*args, 2, nargout_, "args", "getOptArgs");
    mxDestroyArray(Aliases);
    mxDestroyArray(FlagTypeParams);
    mxDestroyArray(RemoveParams);
    mxDestroyArray(SubtractParams);
    mxDestroyArray(ShortcutParams);
    mxDestroyArray(ShortCuts);
    mxDestroyArray(bError);
    mxDestroyArray(nArgs);
    mxDestroyArray(i);
    mxDestroyArray(arg);
    mxDestroyArray(ShortCutParams);
    mxDestroyArray(numargs);
    mxDestroyArray(NumArgCount);
    mxDestroyArray(Fnames);
    mxDestroyArray(name);
    mxDestroyArray(AbbrevIdx);
    mxDestroyArray(FnamesFull);
    mxDestroyArray(FnamesAbbr);
    mxDestroyArray(nAliases);
    mxDestroyArray(naliases);
    mxDestroyArray(FieldIdx);
    mxDestroyArray(j);
    mxDestroyArray(idx);
    mxDestroyArray(l);
    mxDestroyArray(al);
    mxDestroyArray(a);
    mxDestroyArray(ShortCutIdx);
    mxDestroyArray(ans);
    mxDestroyArray(curField);
    mxDestroyArray(val);
    mxDestroyArray(varargin);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return ArgStruct;
}
