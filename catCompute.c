/*
 * MATLAB Compiler: 3.0
 * Date: Thu Apr 22 15:58:01 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "catCompute.m" 
 */
#include "catCompute.h"
#include "getSurrogateSpikes.h"
#include "libmatlbm.h"
#include "libmmfile.h"

static mxChar _array1_[14] = { 'r', 'e', 'f', 'r', 'a', 'c', 't',
                               'o', 'r', 'y', '.', 'm', 'a', 't' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;
static mxArray * _mxarray4_;

static mxChar _array6_[8] = { 'S', 'e', 't', ' ', '%', 'd', 0x005c, 'n' };
static mxArray * _mxarray5_;

static mxChar _array8_[8] = { 'D', 'r', 'o', 'p', 'L', 'a', 's', 't' };
static mxArray * _mxarray7_;

static mxChar _array10_[4] = { '.', 'b', 'i', 'n' };
static mxArray * _mxarray9_;

static mxChar _array12_[1] = { 'w' };
static mxArray * _mxarray11_;

static mxChar _array14_[7] = { 'i', 'e', 'e', 'e', '-', 'l', 'e' };
static mxArray * _mxarray13_;

static mxChar _array16_[5] = { 'i', 'n', 't', '3', '2' };
static mxArray * _mxarray15_;
static mxArray * _mxarray17_;

static mxChar _array19_[2] = { 'F', 'F' };
static mxArray * _mxarray18_;

void InitializeModule_catCompute(void) {
    _mxarray0_ = mclInitializeString(14, _array1_);
    _mxarray2_ = mclInitializeDouble(1.0);
    _mxarray3_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray4_ = mclInitializeDouble(2.0);
    _mxarray5_ = mclInitializeString(8, _array6_);
    _mxarray7_ = mclInitializeString(8, _array8_);
    _mxarray9_ = mclInitializeString(4, _array10_);
    _mxarray11_ = mclInitializeString(1, _array12_);
    _mxarray13_ = mclInitializeString(7, _array14_);
    _mxarray15_ = mclInitializeString(5, _array16_);
    _mxarray17_ = mclInitializeDouble(10.0);
    _mxarray18_ = mclInitializeString(2, _array19_);
}

void TerminateModule_catCompute(void) {
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static void McatCompute(mxArray * sets, mxArray * filename);

_mexLocalFunctionTable _local_function_table_catCompute
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfCatCompute" contains the normal interface for the
 * "catCompute" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/catCompute.m" (lines 1-60). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlfCatCompute(mxArray * sets, mxArray * filename) {
    mlfEnterNewContext(0, 2, sets, filename);
    McatCompute(sets, filename);
    mlfRestorePreviousContext(0, 2, sets, filename);
}

/*
 * The function "mlxCatCompute" contains the feval interface for the
 * "catCompute" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/catCompute.m" (lines 1-60). The
 * feval function calls the implementation version of catCompute through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxCatCompute(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    int i;
    if (nlhs > 0) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: catCompute Line: 1 Column:"
            " 1 The function \"catCompute\" was called with m"
            "ore than the declared number of outputs (0)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: catCompute Line: 1 Column"
            ": 1 The function \"catCompute\" was called with"
            " more than the declared number of inputs (2)."),
          NULL);
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    McatCompute(mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
}

/*
 * The function "McatCompute" is the implementation version of the "catCompute"
 * M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/catCompute.m" (lines 1-60). It
 * contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function catCompute(sets,filename)
 */
static void McatCompute(mxArray * sets, mxArray * filename) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_catCompute);
    mxArray * st = NULL;
    mxArray * j = NULL;
    mxArray * fid = NULL;
    mxArray * sc1 = NULL;
    mxArray * sc = NULL;
    mxArray * cell2array = NULL;
    mxArray * histcie = NULL;
    mxArray * ans = NULL;
    mxArray * i = NULL;
    mxArray * sg = NULL;
    mxArray * sptrain = NULL;
    mxArray * reps = NULL;
    mxArray * numframes = NULL;
    mxArray * fbins = NULL;
    mxArray * rt1 = NULL;
    mxArray * rtlength = NULL;
    mxArray * rf = NULL;
    mxArray * mat = NULL;
    mclCopyArray(&sets);
    mclCopyArray(&filename);
    /*
     * %catCompute Computes various statistics for cat data
     * %   catCompute(sets,filename)
     * 
     * if(ischar(sets))
     */
    if (mlfTobool(mlfIschar(mclVa(sets, "sets")))) {
        /*
         * sets = str2num(sets);
         */
        mlfAssign(&sets, mlfStr2num(NULL, mclVa(sets, "sets")));
    /*
     * end
     */
    }
    /*
     * 
     * % load refractory data
     * mat = load('refractory.mat');
     */
    mlfAssign(&mat, mlfLoadStruct(_mxarray0_, NULL));
    /*
     * rf = struct(mat.rf);
     */
    mlfAssign(&rf, mlfStruct(mlfIndexRef(mclVv(mat, "mat"), ".rf"), NULL));
    /*
     * 
     * % get frame boundaries
     * rtlength = length(rf.data.rtEdges{1});
     */
    mlfAssign(
      &rtlength,
      mclFeval(
        mclValueVarargout(),
        mlxLength,
        mlfIndexRef(mclVv(rf, "rf"), ".data.rtEdges{?}", _mxarray2_),
        NULL));
    /*
     * % take out last value in rtEdges so there will be an even number of points
     * rt1 = reshape(rf.data.rtEdges{1}(1:(rtlength-1)),rf.data.qtframebins,[]);
     */
    mlfAssign(
      &rt1,
      mclFeval(
        mclValueVarargout(),
        mlxReshape,
        mlfIndexRef(
          mclVv(rf, "rf"),
          ".data.rtEdges{?}(?)",
          _mxarray2_,
          mlfColon(
            _mxarray2_,
            mclMinus(mclVv(rtlength, "rtlength"), _mxarray2_),
            NULL)),
        mlfIndexRef(mclVv(rf, "rf"), ".data.qtframebins"),
        _mxarray3_,
        NULL));
    /*
     * % grab the first row which should be the frame limits and add the last
     * % point in rtEdges
     * fbins = [rt1(1,:)'; rf.data.rtEdges{1}(rtlength)];
     */
    mlfAssign(
      &fbins,
      mlfVertcat(
        mlfCtranspose(
          mclArrayRef2(mclVv(rt1, "rt1"), _mxarray2_, mlfCreateColonIndex())),
        mlfIndexRef(
          mclVv(rf, "rf"),
          ".data.rtEdges{?}(?)",
          _mxarray2_,
          mclVv(rtlength, "rtlength")),
        NULL));
    /*
     * % number of frames should be number of columns in rt1
     * numframes = size(rt1,2);
     */
    mlfAssign(
      &numframes, mlfSize(mclValueVarargout(), mclVv(rt1, "rt1"), _mxarray4_));
    /*
     * % get number of repetitions
     * reps = rf.data.repetitions;
     */
    mlfAssign(&reps, mlfIndexRef(mclVv(rf, "rf"), ".data.repetitions"));
    /*
     * 
     * % allocate memory
     * sptrain = cell(1,sets);
     */
    mlfAssign(&sptrain, mlfCell(_mxarray2_, mclVa(sets, "sets"), NULL));
    /*
     * sg.scmean = zeros(numframes,sets);
     */
    mlfIndexAssign(
      &sg,
      ".scmean",
      mlfZeros(mclVv(numframes, "numframes"), mclVa(sets, "sets"), NULL));
    /*
     * sg.scstd = zeros(numframes,sets);
     */
    mlfIndexAssign(
      &sg,
      ".scstd",
      mlfZeros(mclVv(numframes, "numframes"), mclVa(sets, "sets"), NULL));
    /*
     * 
     * % generate surrogates
     * for i = 1:sets
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclVa(sets, "sets"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray3_);
        } else {
            /*
             * fprintf('Set %d\n',i);
             * % generate surrogate spike trains, each repetition is in a column
             * sptrain{i} = getSurrogateSpikes(rf.data,1);
             * % get spike counts for each frame of sptrain, one column for each rep
             * sc = histcie(cell2array(sptrain{i}),fbins,'DropLast');
             * % reshape into repetitions and then transpose so we can take the mean
             * % and std easily
             * sc1 = reshape(sc,numframes,[])';
             * % store mean and std in columns
             * sg.scmean(:,i) = mean(sc1)';
             * sg.scstd(:,i) = std(sc1)';
             * end
             */
            for (; ; ) {
                mclAssignAns(
                  &ans, mlfNFprintf(0, _mxarray5_, mlfScalar(v_), NULL));
                mlfIndexAssign(
                  &sptrain,
                  "{?}",
                  mlfScalar(v_),
                  mclFeval(
                    mclValueVarargout(),
                    mlxGetSurrogateSpikes,
                    mlfIndexRef(mclVv(rf, "rf"), ".data"),
                    _mxarray2_,
                    NULL));
                mlfAssign(
                  &sc,
                  mlfIndexRef(
                    mclVv(histcie, "histcie"),
                    "(?,?,?)",
                    mclArrayRef1(
                      mclVv(cell2array, "cell2array"),
                      mlfIndexRef(
                        mclVv(sptrain, "sptrain"), "{?}", mlfScalar(v_))),
                    mclVv(fbins, "fbins"),
                    _mxarray7_));
                mlfAssign(
                  &sc1,
                  mlfCtranspose(
                    mlfReshape(
                      mclVv(sc, "sc"),
                      mclVv(numframes, "numframes"), _mxarray3_, NULL)));
                mlfIndexAssign(
                  &sg,
                  ".scmean(?,?)",
                  mlfCreateColonIndex(),
                  mlfScalar(v_),
                  mlfCtranspose(mlfMean(mclVv(sc1, "sc1"), NULL)));
                mlfIndexAssign(
                  &sg,
                  ".scstd(?,?)",
                  mlfCreateColonIndex(),
                  mlfScalar(v_),
                  mlfCtranspose(mlfStd(mclVv(sc1, "sc1"), NULL, NULL)));
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    /*
     * 
     * % save surrogates and statistics for surrogates
     * fid = fopen([filename '.bin'],'w','ieee-le');
     */
    mlfAssign(
      &fid,
      mlfFopen(
        NULL,
        NULL,
        mlfHorzcat(mclVa(filename, "filename"), _mxarray9_, NULL),
        _mxarray11_,
        _mxarray13_));
    /*
     * fwrite(fid,[sets reps],'int32');
     */
    mclAssignAns(
      &ans,
      mlfFwrite(
        mclVv(fid, "fid"),
        mlfHorzcat(mclVa(sets, "sets"), mclVv(reps, "reps"), NULL),
        _mxarray15_,
        NULL));
    /*
     * for i = 1:sets
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclVa(sets, "sets"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray3_);
        } else {
            /*
             * for j = 1:reps
             * st = sptrain{i}{j};
             * % write number of spikes for set i, repetition j
             * fwrite(fid,size(st,1),'int32');
             * % write spike times for set i, repetition j
             * % multiply by 10 so we can convert tenth of a ms to integer
             * fwrite(fid,st*10,'int32');
             * end
             * end
             */
            for (; ; ) {
                int v_0 = mclForIntStart(1);
                int e_0 = mclForIntEnd(mclVv(reps, "reps"));
                if (v_0 > e_0) {
                    mlfAssign(&j, _mxarray3_);
                } else {
                    for (; ; ) {
                        mlfAssign(
                          &st,
                          mlfIndexRef(
                            mclVv(sptrain, "sptrain"),
                            "{?}{?}",
                            mlfScalar(v_),
                            mlfScalar(v_0)));
                        mclAssignAns(
                          &ans,
                          mlfFwrite(
                            mclVv(fid, "fid"),
                            mlfSize(
                              mclValueVarargout(), mclVv(st, "st"), _mxarray2_),
                            _mxarray15_,
                            NULL));
                        mclAssignAns(
                          &ans,
                          mlfFwrite(
                            mclVv(fid, "fid"),
                            mclMtimes(mclVv(st, "st"), _mxarray17_),
                            _mxarray15_,
                            NULL));
                        if (v_0 == e_0) {
                            break;
                        }
                        ++v_0;
                    }
                    mlfAssign(&j, mlfScalar(v_0));
                }
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    /*
     * fclose(fid);
     */
    mclAssignAns(&ans, mlfFclose(mclVv(fid, "fid")));
    /*
     * save([filename 'FF'],'sg');
     */
    mlfSave(
      mlfHorzcat(mclVa(filename, "filename"), _mxarray18_, NULL),
      "w",
      "sg",
      sg,
      NULL);
    mxDestroyArray(mat);
    mxDestroyArray(rf);
    mxDestroyArray(rtlength);
    mxDestroyArray(rt1);
    mxDestroyArray(fbins);
    mxDestroyArray(numframes);
    mxDestroyArray(reps);
    mxDestroyArray(sptrain);
    mxDestroyArray(sg);
    mxDestroyArray(i);
    mxDestroyArray(ans);
    mxDestroyArray(histcie);
    mxDestroyArray(cell2array);
    mxDestroyArray(sc);
    mxDestroyArray(sc1);
    mxDestroyArray(fid);
    mxDestroyArray(j);
    mxDestroyArray(st);
    mxDestroyArray(filename);
    mxDestroyArray(sets);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
}
