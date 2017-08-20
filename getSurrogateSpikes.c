/*
 * MATLAB Compiler: 3.0
 * Date: Thu Apr 22 15:58:01 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "catCompute.m" 
 */
#include "getSurrogateSpikes.h"
#include "libmatlbm.h"
#include "libmmfile.h"

static mxChar _array1_[6] = { 'd', 'b', 'f', 'l', 'a', 'g' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;

static mxChar _array4_[5] = { 'c', 's', 't', 'e', 'p' };
static mxArray * _mxarray3_;

static mxChar _array6_[11] = { 'd', 'i', 's', 'p', 'l', 'a',
                               'y', 'f', 'l', 'a', 'g' };
static mxArray * _mxarray5_;

static mxChar _array8_[5] = { 'f', 'l', 'a', 'g', 's' };
static mxArray * _mxarray7_;

static mxArray * _array10_[2] = { NULL /*_mxarray0_*/, NULL /*_mxarray5_*/ };
static mxArray * _mxarray9_;

static mxChar _array12_[5] = { 's', 't', 'a', 't', 'e' };
static mxArray * _mxarray11_;
static mxArray * _mxarray13_;
static mxArray * _mxarray14_;
static mxArray * _mxarray15_;

static mxChar _array17_[30] = { 'R', 'e', 'p', 's', ':', ' ', '%', 'i',
                                ' ', 'n', 'b', 'i', 'n', 's', ':', ' ',
                                '%', 'i', ' ', 'c', 's', 't', 'e', 'p',
                                ':', ' ', '%', 'i', 0x005c, 'n' };
static mxArray * _mxarray16_;
static mxArray * _mxarray18_;

static mxChar _array20_[2] = { 'o', 'n' };
static mxArray * _mxarray19_;

static mxChar _array22_[17] = { 'F', 'i', 'n', 'i', 's', 'h', 'e', 'd', ' ',
                                'r', 'e', 'p', ' ', '%', 'i', 0x005c, 'n' };
static mxArray * _mxarray21_;

static mxChar _array24_[40] = { 'c', 's', 't', 'a', 'r', 't', ':', ' ', '%',
                                'i', ' ', 'c', 'e', 'n', 'd', ':', ' ', '%',
                                'i', ' ', 'e', 'V', 'a', 'l', 'u', 'e', ':',
                                ' ', '%', 'i', ' ', 'r', 'l', 'n', ':', ' ',
                                '%', 'f', 0x005c, 'n' };
static mxArray * _mxarray23_;

static mxChar _array26_[17] = { 'q', 't', ':', ' ', '%', 'f', ' ', 'c', 'S',
                                'u', 'm', ':', ' ', '%', 'f', 0x005c, 'n' };
static mxArray * _mxarray25_;

static mxChar _array28_[2] = { '.', '-' };
static mxArray * _mxarray27_;

static mxChar _array30_[3] = { 'r', '.', '-' };
static mxArray * _mxarray29_;

static mxChar _array32_[3] = { 'g', '.', '-' };
static mxArray * _mxarray31_;

static mxChar _array34_[11] = { 'x', 'r', 'l', 'n', 'i', ':',
                                ' ', '%', 'i', 0x005c, 'n' };
static mxArray * _mxarray33_;

static mxChar _array36_[28] = { 'x', 'r', 'l', 'n', 'i', ':', ' ', '%',
                                'i', ' ', 's', 'p', 't', 'i', ':', ' ',
                                '%', 'i', ' ', 's', 'p', 't', ':', ' ',
                                '%', 'f', 0x005c, 'n' };
static mxArray * _mxarray35_;

static mxChar _array38_[29] = { 'c', 's', 't', 'a', 'r', 't', ':', ' ',
                                '%', 'i', ' ', 'c', 'e', 'n', 'd', ':',
                                ' ', '%', 'i', ' ', 'r', 'l', 'n', ':',
                                ' ', '%', 'f', 0x005c, 'n' };
static mxArray * _mxarray37_;

static mxChar _array40_[24] = { 'q', 't', ':', ' ', '%', 'f', ' ', 'w', 't',
                                ':', ' ', '%', 'f', ' ', 'c', 's', 'u', 'm',
                                ':', ' ', '%', 'f', 0x005c, 'n' };
static mxArray * _mxarray39_;

void InitializeModule_getSurrogateSpikes(void) {
    _mxarray0_ = mclInitializeString(6, _array1_);
    _mxarray2_ = mclInitializeDouble(0.0);
    _mxarray3_ = mclInitializeString(5, _array4_);
    _mxarray5_ = mclInitializeString(11, _array6_);
    _mxarray7_ = mclInitializeString(5, _array8_);
    _array10_[0] = _mxarray0_;
    _array10_[1] = _mxarray5_;
    _mxarray9_ = mclInitializeCellVector(1, 2, _array10_);
    _mxarray11_ = mclInitializeString(5, _array12_);
    _mxarray13_ = mclInitializeDouble(100.0);
    _mxarray14_ = mclInitializeDouble(1000.0);
    _mxarray15_ = mclInitializeDouble(1.0);
    _mxarray16_ = mclInitializeString(30, _array17_);
    _mxarray18_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray19_ = mclInitializeString(2, _array20_);
    _mxarray21_ = mclInitializeString(17, _array22_);
    _mxarray23_ = mclInitializeString(40, _array24_);
    _mxarray25_ = mclInitializeString(17, _array26_);
    _mxarray27_ = mclInitializeString(2, _array28_);
    _mxarray29_ = mclInitializeString(3, _array30_);
    _mxarray31_ = mclInitializeString(3, _array32_);
    _mxarray33_ = mclInitializeString(11, _array34_);
    _mxarray35_ = mclInitializeString(28, _array36_);
    _mxarray37_ = mclInitializeString(29, _array38_);
    _mxarray39_ = mclInitializeString(24, _array40_);
}

void TerminateModule_getSurrogateSpikes(void) {
    mxDestroyArray(_mxarray39_);
    mxDestroyArray(_mxarray37_);
    mxDestroyArray(_mxarray35_);
    mxDestroyArray(_mxarray33_);
    mxDestroyArray(_mxarray31_);
    mxDestroyArray(_mxarray29_);
    mxDestroyArray(_mxarray27_);
    mxDestroyArray(_mxarray25_);
    mxDestroyArray(_mxarray23_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * MgetSurrogateSpikes(int nargout_,
                                     mxArray * rfd,
                                     mxArray * n,
                                     mxArray * varargin);

_mexLocalFunctionTable _local_function_table_getSurrogateSpikes
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfGetSurrogateSpikes" contains the normal interface for the
 * "getSurrogateSpikes" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/getSurrogateSpikes.m" (lines
 * 1-191). This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfGetSurrogateSpikes(mxArray * rfd, mxArray * n, ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * sptrain = NULL;
    mlfVarargin(&varargin, n, 0);
    mlfEnterNewContext(0, -3, rfd, n, varargin);
    sptrain = MgetSurrogateSpikes(nargout, rfd, n, varargin);
    mlfRestorePreviousContext(0, 2, rfd, n);
    mxDestroyArray(varargin);
    return mlfReturnValue(sptrain);
}

/*
 * The function "mlxGetSurrogateSpikes" contains the feval interface for the
 * "getSurrogateSpikes" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/getSurrogateSpikes.m" (lines
 * 1-191). The feval function calls the implementation version of
 * getSurrogateSpikes through this function. This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
void mlxGetSurrogateSpikes(int nlhs,
                           mxArray * plhs[],
                           int nrhs,
                           mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: getSurrogateSpikes Line: 1 Colu"
            "mn: 1 The function \"getSurrogateSpikes\" was called "
            "with more than the declared number of outputs (1)."),
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
    mplhs[0] = MgetSurrogateSpikes(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    mxDestroyArray(mprhs[2]);
}

/*
 * The function "MgetSurrogateSpikes" is the implementation version of the
 * "getSurrogateSpikes" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/getSurrogateSpikes.m" (lines
 * 1-191). It contains the actual compiled code for that M-function. It is a
 * static function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function sptrain = getSurrogateSpikes(rfd,n,varargin)
 */
static mxArray * MgetSurrogateSpikes(int nargout_,
                                     mxArray * rfd,
                                     mxArray * n,
                                     mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(
          &_local_function_table_getSurrogateSpikes);
    mxArray * sptrain = NULL;
    mxArray * spt = NULL;
    mxArray * xrlni = NULL;
    mxArray * pxvals = NULL;
    mxArray * cSum = NULL;
    mxArray * cstart = NULL;
    mxArray * rln = NULL;
    mxArray * r = NULL;
    mxArray * cla = NULL;
    mxArray * bstop = NULL;
    mxArray * cSi = NULL;
    mxArray * eValue = NULL;
    mxArray * cend = NULL;
    mxArray * spti = NULL;
    mxArray * rep = NULL;
    mxArray * matSum = NULL;
    mxArray * qt = NULL;
    mxArray * nbins = NULL;
    mxArray * reps = NULL;
    mxArray * rqt = NULL;
    mxArray * ans = NULL;
    mxArray * cstep = NULL;
    mxArray * getOptArgs = NULL;
    mxArray * Args = NULL;
    mxArray * lwt = NULL;
    mxArray * wt = NULL;
    mclCopyArray(&rfd);
    mclCopyArray(&n);
    mclCopyArray(&varargin);
    /*
     * %getSurrogateSpikes Generate surrogate spikes using q(t)
     * %   SPTRAIN = getSurrogateSpikes(RFD,CELLN,VARARGIN) generates
     * %   surrogate spikes for CELLN in the RFD structure. The following
     * %   fields in RFD are required:
     * %      rfd.qt(n) - the free firing rate, q(t) for cell n.
     * %      rfd.repetitions(n) - the number of repetitions for cell n.
     * %      rfd.wt{n} - the recovery function for cell n.
     * %      rfd.duration - the duration of q(t).
     * %      rfd.qtbinsize - the bin size used to compute q(t) in ms.
     * %      rfd.rtEdges - the time vector for the q(t) function.
     * %
     * %   These are the optional input arguments:
     * %      dbflag - prints out debug information while running.
     * %      cstep - number of data points to step through at a time
     * %              when creating spikes (default is the length of 
     * %              rfd.wt{n}).
     * %      displayflag - plots the intermediate calculations that were
     * %                    used to create surrogate spikes.
     * %
     * %   sptrain = getSurrogateSpikes(rfd,celln,'dbflag','cstep', ...
     * %                length(rfd.wt{n}),'displayflag')
     * 
     * wt = rfd.wt{n};
     */
    mlfAssign(&wt, mlfIndexRef(mclVa(rfd, "rfd"), ".wt{?}", mclVa(n, "n")));
    /*
     * % use length of wt to step through qt
     * lwt = length(wt);
     */
    mlfAssign(&lwt, mlfScalar(mclLengthInt(mclVv(wt, "wt"))));
    /*
     * Args = struct('dbflag',0,'cstep',lwt,'displayflag',0);
     */
    mlfAssign(
      &Args,
      mlfStruct(
        _mxarray0_,
        _mxarray2_,
        _mxarray3_,
        mclVv(lwt, "lwt"),
        _mxarray5_,
        _mxarray2_,
        NULL));
    /*
     * 
     * Args = getOptArgs(varargin,Args,'flags',{'dbflag','displayflag'});
     */
    mlfAssign(
      &Args,
      mlfIndexRef(
        mclVv(getOptArgs, "getOptArgs"),
        "(?,?,?,?)",
        mclVa(varargin, "varargin"),
        mclVv(Args, "Args"),
        _mxarray7_,
        _mxarray9_));
    /*
     * cstep = Args.cstep;
     */
    mlfAssign(&cstep, mlfIndexRef(mclVv(Args, "Args"), ".cstep"));
    /*
     * 
     * % initialize rand seed
     * rand('state',sum(100*clock)); 
     */
    mclAssignAns(
      &ans,
      mlfNRand(
        0,
        _mxarray11_,
        mlfSum(mclMtimes(_mxarray13_, mlfClock()), NULL),
        NULL));
    /*
     * 
     * % get free firing rate, q(t)
     * rqt = rfd.qt{n} * rfd.qtbinsize(n) / 1000;
     */
    mlfAssign(
      &rqt,
      mclMrdivide(
        mclFeval(
          mclValueVarargout(),
          mlxMtimes,
          mlfIndexRef(mclVa(rfd, "rfd"), ".qt{?}", mclVa(n, "n")),
          mlfIndexRef(mclVa(rfd, "rfd"), ".qtbinsize(?)", mclVa(n, "n")),
          NULL),
        _mxarray14_));
    /*
     * % get number of repetitions
     * reps = rfd.repetitions(n);
     */
    mlfAssign(
      &reps, mlfIndexRef(mclVa(rfd, "rfd"), ".repetitions(?)", mclVa(n, "n")));
    /*
     * % get recovery function, w(t) and pad with 1's
     * wt = [wt; ones(cstep-lwt,1)];
     */
    mlfAssign(
      &wt,
      mlfVertcat(
        mclVv(wt, "wt"),
        mlfOnes(
          mclMinus(mclVv(cstep, "cstep"), mclVv(lwt, "lwt")),
          _mxarray15_,
          NULL),
        NULL));
    /*
     * 
     * % get number of bins in qt
     * nbins = length(rqt);
     */
    mlfAssign(&nbins, mlfScalar(mclLengthInt(mclVv(rqt, "rqt"))));
    /*
     * % wrap q(t) function around but since we are only moving in windows of
     * % cstep size, we just need to pad with cstep values
     * qt = [rqt; rqt(1:cstep)];
     */
    mlfAssign(
      &qt,
      mlfVertcat(
        mclVv(rqt, "rqt"),
        mclArrayRef1(
          mclVv(rqt, "rqt"),
          mlfColon(_mxarray15_, mclVv(cstep, "cstep"), NULL)),
        NULL));
    /*
     * % calculate running sum for cstep points
     * matSum = tril(ones(cstep,cstep));
     */
    mlfAssign(
      &matSum,
      mlfTril(
        mlfOnes(mclVv(cstep, "cstep"), mclVv(cstep, "cstep"), NULL), NULL));
    /*
     * 
     * if(Args.dbflag)
     */
    if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".dbflag"))) {
        /*
         * fprintf('Reps: %i nbins: %i cstep: %i\n',reps,nbins,cstep);
         */
        mclAssignAns(
          &ans,
          mlfNFprintf(
            0,
            _mxarray16_,
            mclVv(reps, "reps"),
            mclVv(nbins, "nbins"),
            mclVv(cstep, "cstep"),
            NULL));
    /*
     * end
     */
    }
    /*
     * 
     * % get first spike by using w(t) = 1
     * % create sptrain so we don't have to keep changing memory size
     * % assume you can't have more than 1 spike in 1 ms so make sptrain
     * % equal to the duration of a repetition in ms
     * rep = 1;
     */
    mlfAssign(&rep, _mxarray15_);
    /*
     * sptrain{rep} = zeros(rfd.duration,1);
     */
    mlfIndexAssign(
      &sptrain,
      "{?}",
      mclVv(rep, "rep"),
      mlfZeros(mlfIndexRef(mclVa(rfd, "rfd"), ".duration"), _mxarray15_, NULL));
    /*
     * spti = 0;
     */
    mlfAssign(&spti, _mxarray2_);
    /*
     * % initialize cend for loop 
     * cend = 0;
     */
    mlfAssign(&cend, _mxarray2_);
    /*
     * eValue = 0;
     */
    mlfAssign(&eValue, _mxarray2_);
    /*
     * cSi = [];
     */
    mlfAssign(&cSi, _mxarray18_);
    /*
     * bstop = 0;
     */
    mlfAssign(&bstop, _mxarray2_);
    /*
     * 
     * if(Args.displayflag)
     */
    if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".displayflag"))) {
        /*
         * % clear the figure
         * cla
         */
        mclPrintAns(&ans, mclVv(cla, "cla"));
        /*
         * % set hold to on
         * hold on
         */
        mlfHold(_mxarray19_);
    /*
     * end
     */
    }
    /*
     * 
     * % get random number uniformly distributed in the range {0,1}
     * r = rand;
     */
    mlfAssign(&r, mlfNRand(1, NULL));
    /*
     * rln = -log(r);
     */
    mlfAssign(&rln, mclUminus(mlfLog(mclVv(r, "r"))));
    /*
     * while(1)
     */
    for (;;) {
        /*
         * while( isempty(cSi) )
         */
        while (mlfTobool(mlfIsempty(mclVv(cSi, "cSi")))) {
            /*
             * % restart the indicies if they exceed nbins
             * if(cend>nbins)
             */
            if (mclGtBool(mclVv(cend, "cend"), mclVv(nbins, "nbins"))) {
                /*
                 * % finalize old repetition by removing extraneous zeros
                 * sptrain{rep} = sptrain{rep}(1:spti);
                 */
                mlfIndexAssign(
                  &sptrain,
                  "{?}",
                  mclVv(rep, "rep"),
                  mlfIndexRef(
                    mclVv(sptrain, "sptrain"),
                    "{?}(?)",
                    mclVv(rep, "rep"),
                    mlfColon(_mxarray15_, mclVv(spti, "spti"), NULL)));
                /*
                 * fprintf('Finished rep %i\n',rep);
                 */
                mclAssignAns(
                  &ans, mlfNFprintf(0, _mxarray21_, mclVv(rep, "rep"), NULL));
                /*
                 * % start a new repetition
                 * rep = rep + 1;
                 */
                mlfAssign(&rep, mclPlus(mclVv(rep, "rep"), _mxarray15_));
                /*
                 * if(rep>reps)
                 */
                if (mclGtBool(mclVv(rep, "rep"), mclVv(reps, "reps"))) {
                    /*
                     * % set flag so outer while loop will break as well
                     * bstop = 1;
                     */
                    mlfAssign(&bstop, _mxarray15_);
                    /*
                     * break;
                     */
                    break;
                /*
                 * end
                 */
                }
                /*
                 * sptrain{rep} = zeros(rfd.duration,1);
                 */
                mlfIndexAssign(
                  &sptrain,
                  "{?}",
                  mclVv(rep, "rep"),
                  mlfZeros(
                    mlfIndexRef(mclVa(rfd, "rfd"), ".duration"),
                    _mxarray15_,
                    NULL));
                /*
                 * spti = 0;
                 */
                mlfAssign(&spti, _mxarray2_);
                /*
                 * cend = cend - nbins;
                 */
                mlfAssign(
                  &cend, mclMinus(mclVv(cend, "cend"), mclVv(nbins, "nbins")));
            /*
             * end
             */
            }
            /*
             * % get the next start and end values
             * cstart = cend + 1;
             */
            mlfAssign(&cstart, mclPlus(mclVv(cend, "cend"), _mxarray15_));
            /*
             * cend = cend + cstep;
             */
            mlfAssign(
              &cend, mclPlus(mclVv(cend, "cend"), mclVv(cstep, "cstep")));
            /*
             * % compute the next cstep cummulative sums
             * cSum = eValue + matSum * qt(cstart:cend);
             */
            mlfAssign(
              &cSum,
              mclPlus(
                mclVv(eValue, "eValue"),
                mclMtimes(
                  mclVv(matSum, "matSum"),
                  mclArrayRef1(
                    mclVv(qt, "qt"),
                    mlfColon(
                      mclVv(cstart, "cstart"), mclVv(cend, "cend"), NULL)))));
            /*
             * if(Args.dbflag)
             */
            if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".dbflag"))) {
                /*
                 * fprintf('cstart: %i cend: %i eValue: %i rln: %f\n',cstart,cend,eValue,rln);
                 */
                mclAssignAns(
                  &ans,
                  mlfNFprintf(
                    0,
                    _mxarray23_,
                    mclVv(cstart, "cstart"),
                    mclVv(cend, "cend"),
                    mclVv(eValue, "eValue"),
                    mclVv(rln, "rln"),
                    NULL));
                /*
                 * % take the transpose since fprintf grabs values by columns
                 * fprintf('qt: %f cSum: %f\n',[qt(cstart:cend)';cSum']);
                 */
                mclAssignAns(
                  &ans,
                  mlfNFprintf(
                    0,
                    _mxarray25_,
                    mlfVertcat(
                      mlfCtranspose(
                        mclArrayRef1(
                          mclVv(qt, "qt"),
                          mlfColon(
                            mclVv(cstart, "cstart"),
                            mclVv(cend, "cend"),
                            NULL))),
                      mlfCtranspose(mclVv(cSum, "cSum")),
                      NULL),
                    NULL));
            /*
             * end
             */
            }
            /*
             * % find if there is an index greater than rln
             * cSi = find(cSum>rln);
             */
            mlfAssign(
              &cSi,
              mlfFind(
                NULL, NULL, mclGt(mclVv(cSum, "cSum"), mclVv(rln, "rln"))));
            /*
             * if(Args.displayflag & isempty(cSi))
             */
            {
                mxArray * a_
                  = mclInitialize(
                      mlfIndexRef(mclVv(Args, "Args"), ".displayflag"));
                if (mlfTobool(a_)
                    && mlfTobool(mclAnd(a_, mlfIsempty(mclVv(cSi, "cSi"))))) {
                    mxDestroyArray(a_);
                    /*
                     * pxvals = rfd.rtEdges{n}(cstart:cend);
                     */
                    mlfAssign(
                      &pxvals,
                      mlfIndexRef(
                        mclVa(rfd, "rfd"),
                        ".rtEdges{?}(?)",
                        mclVa(n, "n"),
                        mlfColon(
                          mclVv(cstart, "cstart"), mclVv(cend, "cend"), NULL)));
                    /*
                     * % plot qt from cstart to cend
                     * plot(pxvals,qt(cstart:cend),'.-')
                     */
                    mclPrintAns(
                      &ans,
                      mlfNPlot(
                        0,
                        mclVv(pxvals, "pxvals"),
                        mclArrayRef1(
                          mclVv(qt, "qt"),
                          mlfColon(
                            mclVv(cstart, "cstart"),
                            mclVv(cend, "cend"),
                            NULL)),
                        _mxarray27_,
                        NULL));
                    /*
                     * % plot cSum from cstart to cend
                     * plot(pxvals,cSum,'r.-')
                     */
                    mclPrintAns(
                      &ans,
                      mlfNPlot(
                        0,
                        mclVv(pxvals, "pxvals"),
                        mclVv(cSum, "cSum"),
                        _mxarray29_,
                        NULL));
                    /*
                     * % plot rln from cstart to cend
                     * plot(pxvals,rln,'g.-')
                     */
                    mclPrintAns(
                      &ans,
                      mlfNPlot(
                        0,
                        mclVv(pxvals, "pxvals"),
                        mclVv(rln, "rln"),
                        _mxarray31_,
                        NULL));
                } else {
                    mxDestroyArray(a_);
                }
            /*
             * end
             */
            }
            /*
             * % get the last value from the cummulative sum for next
             * % calculation
             * eValue = cSum(cstep);
             */
            mlfAssign(
              &eValue,
              mclArrayRef1(mclVv(cSum, "cSum"), mclVv(cstep, "cstep")));
        /*
         * end
         */
        }
        /*
         * if(bstop)
         */
        if (mlfTobool(mclVv(bstop, "bstop"))) {
            /*
             * % inner loop exceeded reps so stop
             * break;
             */
            break;
        /*
         * end
         */
        }
        /*
         * % value found so figure out where to put the spike
         * % get index that first exceeds rln
         * % if it was index 22, cstart will be 21, cSi(1) will be 2 so in 
         * % order to get back 22, we subtract 1 from cstart + cSi(1)
         * xrlni = cstart+cSi(1)-1;
         */
        mlfAssign(
          &xrlni,
          mclMinus(
            mclPlus(
              mclVv(cstart, "cstart"), mclIntArrayRef1(mclVv(cSi, "cSi"), 1)),
            _mxarray15_));
        /*
         * if(Args.dbflag)
         */
        if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".dbflag"))) {
            /*
             * fprintf('xrlni: %i\n',xrlni);
             */
            mclAssignAns(
              &ans, mlfNFprintf(0, _mxarray33_, mclVv(xrlni, "xrlni"), NULL));
        /*
         * end
         */
        }
        /*
         * if(Args.displayflag)
         */
        if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".displayflag"))) {
            /*
             * pxvals = rfd.rtEdges{n}(cstart:xrlni);
             */
            mlfAssign(
              &pxvals,
              mlfIndexRef(
                mclVa(rfd, "rfd"),
                ".rtEdges{?}(?)",
                mclVa(n, "n"),
                mlfColon(
                  mclVv(cstart, "cstart"), mclVv(xrlni, "xrlni"), NULL)));
            /*
             * plot(pxvals,qt(cstart:xrlni),'.-')
             */
            mclPrintAns(
              &ans,
              mlfNPlot(
                0,
                mclVv(pxvals, "pxvals"),
                mclArrayRef1(
                  mclVv(qt, "qt"),
                  mlfColon(
                    mclVv(cstart, "cstart"), mclVv(xrlni, "xrlni"), NULL)),
                _mxarray27_,
                NULL));
            /*
             * plot(pxvals,cSum(1:cSi(1)),'r.-')
             */
            mclPrintAns(
              &ans,
              mlfNPlot(
                0,
                mclVv(pxvals, "pxvals"),
                mclArrayRef1(
                  mclVv(cSum, "cSum"),
                  mlfColon(
                    _mxarray15_, mclIntArrayRef1(mclVv(cSi, "cSi"), 1), NULL)),
                _mxarray29_,
                NULL));
            /*
             * plot(pxvals,rln,'g.-')
             */
            mclPrintAns(
              &ans,
              mlfNPlot(
                0,
                mclVv(pxvals, "pxvals"),
                mclVv(rln, "rln"),
                _mxarray31_,
                NULL));
        /*
         * end
         */
        }
        /*
         * % check to see if we need to go to the next repetition
         * if(xrlni>nbins)
         */
        if (mclGtBool(mclVv(xrlni, "xrlni"), mclVv(nbins, "nbins"))) {
            /*
             * % finalize old repetition by removing extraneous zeros
             * sptrain{rep} = sptrain{rep}(1:spti);
             */
            mlfIndexAssign(
              &sptrain,
              "{?}",
              mclVv(rep, "rep"),
              mlfIndexRef(
                mclVv(sptrain, "sptrain"),
                "{?}(?)",
                mclVv(rep, "rep"),
                mlfColon(_mxarray15_, mclVv(spti, "spti"), NULL)));
            /*
             * fprintf('Finished rep %i\n',rep);
             */
            mclAssignAns(
              &ans, mlfNFprintf(0, _mxarray21_, mclVv(rep, "rep"), NULL));
            /*
             * % create new repetition
             * rep = rep + 1;
             */
            mlfAssign(&rep, mclPlus(mclVv(rep, "rep"), _mxarray15_));
            /*
             * if(rep>reps)
             */
            if (mclGtBool(mclVv(rep, "rep"), mclVv(reps, "reps"))) {
                /*
                 * % break out of while loop
                 * break;
                 */
                break;
            /*
             * end
             */
            }
            /*
             * sptrain{rep} = zeros(rfd.duration,1);
             */
            mlfIndexAssign(
              &sptrain,
              "{?}",
              mclVv(rep, "rep"),
              mlfZeros(
                mlfIndexRef(mclVa(rfd, "rfd"), ".duration"),
                _mxarray15_,
                NULL));
            /*
             * spti = 0;
             */
            mlfAssign(&spti, _mxarray2_);
            /*
             * % adjust xrlni to the new repetition's time scale
             * xrlni = xrlni - nbins;
             */
            mlfAssign(
              &xrlni, mclMinus(mclVv(xrlni, "xrlni"), mclVv(nbins, "nbins")));
        /*
         * end
         */
        }
        /*
         * % convert index to spike time
         * spt = rfd.rtEdges{n}(xrlni);
         */
        mlfAssign(
          &spt,
          mlfIndexRef(
            mclVa(rfd, "rfd"),
            ".rtEdges{?}(?)",
            mclVa(n, "n"),
            mclVv(xrlni, "xrlni")));
        /*
         * % increment spike count
         * spti = spti + 1;
         */
        mlfAssign(&spti, mclPlus(mclVv(spti, "spti"), _mxarray15_));
        /*
         * % add value to spike train
         * sptrain{rep}(spti) = spt;
         */
        mlfIndexAssign(
          &sptrain,
          "{?}(?)",
          mclVv(rep, "rep"),
          mclVv(spti, "spti"),
          mclVv(spt, "spt"));
        /*
         * % get new random number
         * r = rand;
         */
        mlfAssign(&r, mlfNRand(1, NULL));
        /*
         * rln =  -log(r);	
         */
        mlfAssign(&rln, mclUminus(mlfLog(mclVv(r, "r"))));
        /*
         * % reset cstart, cend, eValue and cSi
         * cstart = xrlni;
         */
        mlfAssign(&cstart, mclVv(xrlni, "xrlni"));
        /*
         * % if cstep is 10 and cstart is 22 then cend = 22 + 10 - 1
         * % will be 10 values
         * cend = cstart + cstep - 1;
         */
        mlfAssign(
          &cend,
          mclMinus(
            mclPlus(mclVv(cstart, "cstart"), mclVv(cstep, "cstep")),
            _mxarray15_));
        /*
         * % find cummulative sum again taking into account the
         * % relative refractory period
         * cSum = matSum * (qt(cstart:cend) .* wt);
         */
        mlfAssign(
          &cSum,
          mclMtimes(
            mclVv(matSum, "matSum"),
            mclTimes(
              mclArrayRef1(
                mclVv(qt, "qt"),
                mlfColon(mclVv(cstart, "cstart"), mclVv(cend, "cend"), NULL)),
              mclVv(wt, "wt"))));
        /*
         * if(Args.dbflag)
         */
        if (mlfTobool(mlfIndexRef(mclVv(Args, "Args"), ".dbflag"))) {
            /*
             * fprintf('xrlni: %i spti: %i spt: %f\n',xrlni,spti,spt);
             */
            mclAssignAns(
              &ans,
              mlfNFprintf(
                0,
                _mxarray35_,
                mclVv(xrlni, "xrlni"),
                mclVv(spti, "spti"),
                mclVv(spt, "spt"),
                NULL));
            /*
             * fprintf('cstart: %i cend: %i rln: %f\n',cstart,cend,rln);
             */
            mclAssignAns(
              &ans,
              mlfNFprintf(
                0,
                _mxarray37_,
                mclVv(cstart, "cstart"),
                mclVv(cend, "cend"),
                mclVv(rln, "rln"),
                NULL));
            /*
             * fprintf('qt: %f wt: %f csum: %f\n',[qt(cstart:cend)'; wt'; cSum']);
             */
            mclAssignAns(
              &ans,
              mlfNFprintf(
                0,
                _mxarray39_,
                mlfVertcat(
                  mlfCtranspose(
                    mclArrayRef1(
                      mclVv(qt, "qt"),
                      mlfColon(
                        mclVv(cstart, "cstart"), mclVv(cend, "cend"), NULL))),
                  mlfCtranspose(mclVv(wt, "wt")),
                  mlfCtranspose(mclVv(cSum, "cSum")),
                  NULL),
                NULL));
        /*
         * end
         */
        }
        /*
         * % find if there is an index greater than rln
         * cSi = find(cSum>rln);
         */
        mlfAssign(
          &cSi,
          mlfFind(NULL, NULL, mclGt(mclVv(cSum, "cSum"), mclVv(rln, "rln"))));
        /*
         * if(Args.displayflag & isempty(cSi))
         */
        {
            mxArray * a_
              = mclInitialize(
                  mlfIndexRef(mclVv(Args, "Args"), ".displayflag"));
            if (mlfTobool(a_)
                && mlfTobool(mclAnd(a_, mlfIsempty(mclVv(cSi, "cSi"))))) {
                mxDestroyArray(a_);
                /*
                 * pxvals = rfd.rtEdges{n}(cstart:cend);
                 */
                mlfAssign(
                  &pxvals,
                  mlfIndexRef(
                    mclVa(rfd, "rfd"),
                    ".rtEdges{?}(?)",
                    mclVa(n, "n"),
                    mlfColon(
                      mclVv(cstart, "cstart"), mclVv(cend, "cend"), NULL)));
                /*
                 * % plot qt from cstart to cend
                 * plot(pxvals,qt(cstart:cend),'.-')
                 */
                mclPrintAns(
                  &ans,
                  mlfNPlot(
                    0,
                    mclVv(pxvals, "pxvals"),
                    mclArrayRef1(
                      mclVv(qt, "qt"),
                      mlfColon(
                        mclVv(cstart, "cstart"), mclVv(cend, "cend"), NULL)),
                    _mxarray27_,
                    NULL));
                /*
                 * % plot cSum from cstart to cend
                 * plot(pxvals,cSum,'r.-')
                 */
                mclPrintAns(
                  &ans,
                  mlfNPlot(
                    0,
                    mclVv(pxvals, "pxvals"),
                    mclVv(cSum, "cSum"),
                    _mxarray29_,
                    NULL));
                /*
                 * % plot rln from cstart to cend
                 * plot(pxvals,rln,'g.-')
                 */
                mclPrintAns(
                  &ans,
                  mlfNPlot(
                    0,
                    mclVv(pxvals, "pxvals"),
                    mclVv(rln, "rln"),
                    _mxarray31_,
                    NULL));
            } else {
                mxDestroyArray(a_);
            }
        /*
         * end
         */
        }
        /*
         * % get the last value from the cummulative sum for next
         * % calculation
         * eValue = cSum(cstep);
         */
        mlfAssign(
          &eValue, mclArrayRef1(mclVv(cSum, "cSum"), mclVv(cstep, "cstep")));
    /*
     * end	% end of while(1) loop
     */
    }
    mclValidateOutput(sptrain, 1, nargout_, "sptrain", "getSurrogateSpikes");
    mxDestroyArray(wt);
    mxDestroyArray(lwt);
    mxDestroyArray(Args);
    mxDestroyArray(getOptArgs);
    mxDestroyArray(cstep);
    mxDestroyArray(ans);
    mxDestroyArray(rqt);
    mxDestroyArray(reps);
    mxDestroyArray(nbins);
    mxDestroyArray(qt);
    mxDestroyArray(matSum);
    mxDestroyArray(rep);
    mxDestroyArray(spti);
    mxDestroyArray(cend);
    mxDestroyArray(eValue);
    mxDestroyArray(cSi);
    mxDestroyArray(bstop);
    mxDestroyArray(cla);
    mxDestroyArray(r);
    mxDestroyArray(rln);
    mxDestroyArray(cstart);
    mxDestroyArray(cSum);
    mxDestroyArray(pxvals);
    mxDestroyArray(xrlni);
    mxDestroyArray(spt);
    mxDestroyArray(varargin);
    mxDestroyArray(n);
    mxDestroyArray(rfd);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return sptrain;
}
