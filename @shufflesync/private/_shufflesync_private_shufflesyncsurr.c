/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jun 26 12:13:21 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "shufflesyncsurr" 
 */
#include "_shufflesync_private_shufflesyncsurr.h"
#include "concatenate.h"
#include "histcie.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "vecc.h"
#include "vecr.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;
static mxArray * _mxarray6_;

static mxChar _array8_[5] = { 's', 't', 'a', 't', 'e' };
static mxArray * _mxarray7_;
static mxArray * _mxarray9_;
static double _ieee_nan_;
static mxArray * _mxarray10_;

static mxChar _array12_[8] = { 'D', 'r', 'o', 'p', 'L', 'a', 's', 't' };
static mxArray * _mxarray11_;

static mxChar _array14_[4] = { '%', '0', '4', 'd' };
static mxArray * _mxarray13_;

void InitializeModule__shufflesync_private_shufflesyncsurr(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeDouble(2.0);
    _mxarray2_ = mclInitializeDouble(3.0);
    _mxarray3_ = mclInitializeDouble(4.0);
    _mxarray4_ = mclInitializeDouble(0.0);
    _mxarray5_ = mclInitializeCellVector(0, 0, (mxArray * *)NULL);
    _mxarray6_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray7_ = mclInitializeString(5, _array8_);
    _mxarray9_ = mclInitializeDouble(-1.0);
    _ieee_nan_ = mclGetNaN();
    _mxarray10_ = mclInitializeDouble(_ieee_nan_);
    _mxarray11_ = mclInitializeString(8, _array12_);
    _mxarray13_ = mclInitializeString(4, _array14_);
}

void TerminateModule__shufflesync_private_shufflesyncsurr(void) {
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static void M_shufflesync_private_shufflesyncsurr(mxArray * sdatafile,
                                                  mxArray * varargin);

_mexLocalFunctionTable _local_function_table__shufflesync_private_shufflesyncsurr
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlf_shufflesync_private_shufflesyncsurr" contains the normal
 * interface for the "@shufflesync/private/shufflesyncsurr" M-function from
 * file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/@shufflesync/private/shufflesyncs
 * urr.m" (lines 1-532). This function processes any input arguments and passes
 * them to the implementation version of the function, appearing above.
 */
void mlf_shufflesync_private_shufflesyncsurr(mxArray * sdatafile, ...) {
    mxArray * varargin = NULL;
    mlfVarargin(&varargin, sdatafile, 0);
    mlfEnterNewContext(0, -2, sdatafile, varargin);
    M_shufflesync_private_shufflesyncsurr(sdatafile, varargin);
    mlfRestorePreviousContext(0, 1, sdatafile);
    mxDestroyArray(varargin);
}

/*
 * The function "mlx_shufflesync_private_shufflesyncsurr" contains the feval
 * interface for the "@shufflesync/private/shufflesyncsurr" M-function from
 * file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/@shufflesync/private/shufflesyncs
 * urr.m" (lines 1-532). The feval function calls the implementation version of
 * @shufflesync/private/shufflesyncsurr through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
void mlx_shufflesync_private_shufflesyncsurr(int nlhs,
                                             mxArray * plhs[],
                                             int nrhs,
                                             mxArray * prhs[]) {
    mxArray * mprhs[2];
    int i;
    if (nlhs > 0) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: @shufflesync/private/shufflesyncsurr Line: "
            "1 Column: 1 The function \"@shufflesync/private/shufflesyncsurr\""
            " was called with more than the declared number of outputs (0)."),
          NULL);
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
    M_shufflesync_private_shufflesyncsurr(mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    mxDestroyArray(mprhs[1]);
}

/*
 * The function "M_shufflesync_private_shufflesyncsurr" is the implementation
 * version of the "@shufflesync/private/shufflesyncsurr" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/@shufflesync/private/shufflesyncs
 * urr.m" (lines 1-532). It contains the actual compiled code for that
 * M-function. It is a static function and must only be called from one of the
 * interface functions, appearing below.
 */
/*
 * function shufflesyncsurr(sdatafile,varargin)
 */
static void M_shufflesync_private_shufflesyncsurr(mxArray * sdatafile,
                                                  mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(
          &_local_function_table__shufflesync_private_shufflesyncsurr);
    int nargin_ = mclNargin(-2, sdatafile, varargin, NULL);
    mxArray * surrprefix = NULL;
    mxArray * excsptime = NULL;
    mxArray * enspikes = NULL;
    mxArray * esptime = NULL;
    mxArray * eclindices = NULL;
    mxArray * esrepi = NULL;
    mxArray * nspikes = NULL;
    mxArray * sptime = NULL;
    mxArray * clindices = NULL;
    mxArray * srepi = NULL;
    mxArray * repi2 = NULL;
    mxArray * a = NULL;
    mxArray * esrinopair = NULL;
    mxArray * discard = NULL;
    mxArray * srinopair = NULL;
    mxArray * _T0_ = NULL;
    mxArray * repsnopair = NULL;
    mxArray * erepnumscell = NULL;
    mxArray * repnumscell = NULL;
    mxArray * ecl1cell = NULL;
    mxArray * eclcell = NULL;
    mxArray * cl1cell = NULL;
    mxArray * clcell = NULL;
    mxArray * excsptimes2 = NULL;
    mxArray * excsptimes1 = NULL;
    mxArray * sptimes2 = NULL;
    mxArray * sptimes1 = NULL;
    mxArray * sprepnl = NULL;
    mxArray * sprepnvec = NULL;
    mxArray * sprepsboth = NULL;
    mxArray * spsria2 = NULL;
    mxArray * spsrib1 = NULL;
    mxArray * spintrepnums2 = NULL;
    mxArray * spsrib2 = NULL;
    mxArray * spsria1 = NULL;
    mxArray * spintrepnums1 = NULL;
    mxArray * sprepnums2 = NULL;
    mxArray * sperepnums2 = NULL;
    mxArray * lagbins2 = NULL;
    mxArray * xcmat = NULL;
    mxArray * spt1 = NULL;
    mxArray * randnums = NULL;
    mxArray * shuffle = NULL;
    mxArray * nisi0 = NULL;
    mxArray * surri2 = NULL;
    mxArray * surrvec2 = NULL;
    mxArray * surri1 = NULL;
    mxArray * surrvec1 = NULL;
    mxArray * surrvec = NULL;
    mxArray * surri = NULL;
    mxArray * nref = NULL;
    mxArray * bisi = NULL;
    mxArray * shufflevecl = NULL;
    mxArray * bspikes = NULL;
    mxArray * celli = NULL;
    mxArray * nfc = NULL;
    mxArray * needrefcheck = NULL;
    mxArray * noskip2 = NULL;
    mxArray * numSurrogates1 = NULL;
    mxArray * noskip1 = NULL;
    mxArray * cellvecr = NULL;
    mxArray * cellvec = NULL;
    mxArray * xc1 = NULL;
    mxArray * numSurrogates2 = NULL;
    mxArray * nexcvals2 = NULL;
    mxArray * excvals2 = NULL;
    mxArray * emat1 = NULL;
    mxArray * nexcvals1 = NULL;
    mxArray * excvals1 = NULL;
    mxArray * emat2 = NULL;
    mxArray * nxcvals = NULL;
    mxArray * xcvals = NULL;
    mxArray * mat2 = NULL;
    mxArray * mat1 = NULL;
    mxArray * excnspikescells = NULL;
    mxArray * excnsc2 = NULL;
    mxArray * excnsc1 = NULL;
    mxArray * excsptimes = NULL;
    mxArray * excspt2 = NULL;
    mxArray * excspt1 = NULL;
    mxArray * enspikescells = NULL;
    mxArray * ensc2 = NULL;
    mxArray * ensc1 = NULL;
    mxArray * esptimes = NULL;
    mxArray * espt2 = NULL;
    mxArray * espt1 = NULL;
    mxArray * ec2lindices = NULL;
    mxArray * ec1lindices = NULL;
    mxArray * nspikescells = NULL;
    mxArray * nsc2 = NULL;
    mxArray * nsc1 = NULL;
    mxArray * sptimes = NULL;
    mxArray * c2lindices = NULL;
    mxArray * c1lindices = NULL;
    mxArray * repn = NULL;
    mxArray * repi = NULL;
    mxArray * endlimits = NULL;
    mxArray * tmpendlimits = NULL;
    mxArray * startlimits = NULL;
    mxArray * tmpstartlimits = NULL;
    mxArray * c2l1 = NULL;
    mxArray * c1l1 = NULL;
    mxArray * blr1 = NULL;
    mxArray * cell2limits = NULL;
    mxArray * cell1limits = NULL;
    mxArray * ec2l1 = NULL;
    mxArray * ec1l1 = NULL;
    mxArray * ecell2limits = NULL;
    mxArray * ecell1limits = NULL;
    mxArray * esptrains = NULL;
    mxArray * elrepn2 = NULL;
    mxArray * elrepn1 = NULL;
    mxArray * spiketrains = NULL;
    mxArray * lrepn2 = NULL;
    mxArray * lrepn1 = NULL;
    mxArray * xcreps = NULL;
    mxArray * repswspike2 = NULL;
    mxArray * repswspike1 = NULL;
    mxArray * sptrains = NULL;
    mxArray * xc = NULL;
    mxArray * xctemp = NULL;
    mxArray * repnl = NULL;
    mxArray * repnvec = NULL;
    mxArray * repsboth = NULL;
    mxArray * sria2 = NULL;
    mxArray * srib1 = NULL;
    mxArray * intrepnums2 = NULL;
    mxArray * srib2 = NULL;
    mxArray * sria1 = NULL;
    mxArray * intrepnums1 = NULL;
    mxArray * repnums2 = NULL;
    mxArray * repnums1 = NULL;
    mxArray * cell1i = NULL;
    mxArray * cell2i = NULL;
    mxArray * urepspikeb = NULL;
    mxArray * urepspikea = NULL;
    mxArray * urepspike = NULL;
    mxArray * crepspike = NULL;
    mxArray * ctrialspikebin = NULL;
    mxArray * ctsbi = NULL;
    mxArray * fmax2 = NULL;
    mxArray * endbins = NULL;
    mxArray * fmin2 = NULL;
    mxArray * startbins = NULL;
    mxArray * erepnums2 = NULL;
    mxArray * erepnums1 = NULL;
    mxArray * ecell1i = NULL;
    mxArray * ecell2i = NULL;
    mxArray * reps = NULL;
    mxArray * eurepspikeb = NULL;
    mxArray * eurepspikea = NULL;
    mxArray * eurepspike = NULL;
    mxArray * trialspikebin = NULL;
    mxArray * binidxsize = NULL;
    mxArray * emptytsn = NULL;
    mxArray * repspike = NULL;
    mxArray * trialspiken = NULL;
    mxArray * binidx = NULL;
    mxArray * fmax = NULL;
    mxArray * endbins2 = NULL;
    mxArray * fmin = NULL;
    mxArray * startbins2 = NULL;
    mxArray * windown = NULL;
    mxArray * numwindows = NULL;
    mxArray * endwindow = NULL;
    mxArray * startwindow = NULL;
    mxArray * cellpairdatafile = NULL;
    mxArray * ans = NULL;
    mxArray * CEN2COL = NULL;
    mxArray * CEN1COL = NULL;
    mxArray * EXT2COL = NULL;
    mxArray * EXT1COL = NULL;
    mxArray * SHIFTCOL = NULL;
    mxArray * DATACOL = NULL;
    mclCopyArray(&sdatafile);
    mclCopyArray(&varargin);
    /*
     * 
     * % define some constants
     * % the values for the data are stored in the first column of the xcorr 
     * % matrix as well as the spike times matrix
     * DATACOL = 1;
     */
    mlfAssign(&DATACOL, _mxarray0_);
    /*
     * % data for the shift predictor is stored in the second column of the xcorr 
     * % matrix
     * SHIFTCOL = 2;
     */
    mlfAssign(&SHIFTCOL, _mxarray1_);
    /*
     * % column numbers for startlimits and endlimits
     * EXT1COL = 1;
     */
    mlfAssign(&EXT1COL, _mxarray0_);
    /*
     * EXT2COL = 2;
     */
    mlfAssign(&EXT2COL, _mxarray1_);
    /*
     * CEN1COL = 3;
     */
    mlfAssign(&CEN1COL, _mxarray2_);
    /*
     * CEN2COL = 4;
     */
    mlfAssign(&CEN2COL, _mxarray3_);
    /*
     * 
     * % load surrogate data file
     * load(sdatafile);
     */
    {
        mxArray * name_ = mclInitialize(mclVa(sdatafile, "sdatafile"));
        mclLoadConditional(
          name_,
          "CEN1COL",
          &CEN1COL,
          "CEN2COL",
          &CEN2COL,
          "DATACOL",
          &DATACOL,
          "EXT1COL",
          &EXT1COL,
          "EXT2COL",
          &EXT2COL,
          "SHIFTCOL",
          &SHIFTCOL,
          "a",
          &a,
          "ans",
          &ans,
          "binidx",
          &binidx,
          "binidxsize",
          &binidxsize,
          "bisi",
          &bisi,
          "blr1",
          &blr1,
          "bspikes",
          &bspikes,
          "c1l1",
          &c1l1,
          "c1lindices",
          &c1lindices,
          NULL);
        {
            mxArray * name_0 = mclInitialize(name_);
            mclLoadConditional(
              name_0,
              "c2l1",
              &c2l1,
              "c2lindices",
              &c2lindices,
              "cell1i",
              &cell1i,
              "cell1limits",
              &cell1limits,
              "cell2i",
              &cell2i,
              "cell2limits",
              &cell2limits,
              "celli",
              &celli,
              "cellpairdatafile",
              &cellpairdatafile,
              "cellvec",
              &cellvec,
              "cellvecr",
              &cellvecr,
              "cl1cell",
              &cl1cell,
              "clcell",
              &clcell,
              "clindices",
              &clindices,
              "crepspike",
              &crepspike,
              "ctrialspikebin",
              &ctrialspikebin,
              NULL);
            {
                mxArray * name_1 = mclInitialize(name_0);
                mclLoadConditional(
                  name_1,
                  "ctsbi",
                  &ctsbi,
                  "discard",
                  &discard,
                  "ec1l1",
                  &ec1l1,
                  "ec1lindices",
                  &ec1lindices,
                  "ec2l1",
                  &ec2l1,
                  "ec2lindices",
                  &ec2lindices,
                  "ecell1i",
                  &ecell1i,
                  "ecell1limits",
                  &ecell1limits,
                  "ecell2i",
                  &ecell2i,
                  "ecell2limits",
                  &ecell2limits,
                  "ecl1cell",
                  &ecl1cell,
                  "eclcell",
                  &eclcell,
                  "eclindices",
                  &eclindices,
                  "elrepn1",
                  &elrepn1,
                  "elrepn2",
                  &elrepn2,
                  NULL);
                {
                    mxArray * name_2 = mclInitialize(name_1);
                    mclLoadConditional(
                      name_2,
                      "emat1",
                      &emat1,
                      "emat2",
                      &emat2,
                      "emptytsn",
                      &emptytsn,
                      "endbins",
                      &endbins,
                      "endbins2",
                      &endbins2,
                      "endlimits",
                      &endlimits,
                      "endwindow",
                      &endwindow,
                      "ensc1",
                      &ensc1,
                      "ensc2",
                      &ensc2,
                      "enspikes",
                      &enspikes,
                      "enspikescells",
                      &enspikescells,
                      "erepnums1",
                      &erepnums1,
                      "erepnums2",
                      &erepnums2,
                      "erepnumscell",
                      &erepnumscell,
                      "espt1",
                      &espt1,
                      NULL);
                    {
                        mxArray * name_3 = mclInitialize(name_2);
                        mclLoadConditional(
                          name_3,
                          "espt2",
                          &espt2,
                          "esptime",
                          &esptime,
                          "esptimes",
                          &esptimes,
                          "esptrains",
                          &esptrains,
                          "esrepi",
                          &esrepi,
                          "esrinopair",
                          &esrinopair,
                          "eurepspike",
                          &eurepspike,
                          "eurepspikea",
                          &eurepspikea,
                          "eurepspikeb",
                          &eurepspikeb,
                          "excnsc1",
                          &excnsc1,
                          "excnsc2",
                          &excnsc2,
                          "excnspikescells",
                          &excnspikescells,
                          "excspt1",
                          &excspt1,
                          "excspt2",
                          &excspt2,
                          "excsptime",
                          &excsptime,
                          NULL);
                        {
                            mxArray * name_4 = mclInitialize(name_3);
                            mclLoadConditional(
                              name_4,
                              "excsptimes",
                              &excsptimes,
                              "excsptimes1",
                              &excsptimes1,
                              "excsptimes2",
                              &excsptimes2,
                              "excvals1",
                              &excvals1,
                              "excvals2",
                              &excvals2,
                              "fmax",
                              &fmax,
                              "fmax2",
                              &fmax2,
                              "fmin",
                              &fmin,
                              "fmin2",
                              &fmin2,
                              "intrepnums1",
                              &intrepnums1,
                              "intrepnums2",
                              &intrepnums2,
                              "lagbins2",
                              &lagbins2,
                              "lrepn1",
                              &lrepn1,
                              "lrepn2",
                              &lrepn2,
                              "mat1",
                              &mat1,
                              NULL);
                            {
                                mxArray * name_5 = mclInitialize(name_4);
                                mclLoadConditional(
                                  name_5,
                                  "mat2",
                                  &mat2,
                                  "needrefcheck",
                                  &needrefcheck,
                                  "nexcvals1",
                                  &nexcvals1,
                                  "nexcvals2",
                                  &nexcvals2,
                                  "nfc",
                                  &nfc,
                                  "nisi0",
                                  &nisi0,
                                  "noskip1",
                                  &noskip1,
                                  "noskip2",
                                  &noskip2,
                                  "nref",
                                  &nref,
                                  "nsc1",
                                  &nsc1,
                                  "nsc2",
                                  &nsc2,
                                  "nspikes",
                                  &nspikes,
                                  "nspikescells",
                                  &nspikescells,
                                  "numSurrogates1",
                                  &numSurrogates1,
                                  "numSurrogates2",
                                  &numSurrogates2,
                                  NULL);
                                {
                                    mxArray * name_6 = mclInitialize(name_5);
                                    mclLoadConditional(
                                      name_6,
                                      "numwindows",
                                      &numwindows,
                                      "nxcvals",
                                      &nxcvals,
                                      "randnums",
                                      &randnums,
                                      "repi",
                                      &repi,
                                      "repi2",
                                      &repi2,
                                      "repn",
                                      &repn,
                                      "repnl",
                                      &repnl,
                                      "repnums1",
                                      &repnums1,
                                      "repnums2",
                                      &repnums2,
                                      "repnumscell",
                                      &repnumscell,
                                      "repnvec",
                                      &repnvec,
                                      "reps",
                                      &reps,
                                      "repsboth",
                                      &repsboth,
                                      "repsnopair",
                                      &repsnopair,
                                      "repspike",
                                      &repspike,
                                      NULL);
                                    {
                                        mxArray * name_7
                                          = mclInitialize(name_6);
                                        mclLoadConditional(
                                          name_7,
                                          "repswspike1",
                                          &repswspike1,
                                          "repswspike2",
                                          &repswspike2,
                                          "sdatafile",
                                          &sdatafile,
                                          "shuffle",
                                          &shuffle,
                                          "shufflevecl",
                                          &shufflevecl,
                                          "sperepnums2",
                                          &sperepnums2,
                                          "spiketrains",
                                          &spiketrains,
                                          "spintrepnums1",
                                          &spintrepnums1,
                                          "spintrepnums2",
                                          &spintrepnums2,
                                          "sprepnl",
                                          &sprepnl,
                                          "sprepnums2",
                                          &sprepnums2,
                                          "sprepnvec",
                                          &sprepnvec,
                                          "sprepsboth",
                                          &sprepsboth,
                                          "spsria1",
                                          &spsria1,
                                          "spsria2",
                                          &spsria2,
                                          NULL);
                                        {
                                            mxArray * name_8
                                              = mclInitialize(name_7);
                                            mclLoadConditional(
                                              name_8,
                                              "spsrib1",
                                              &spsrib1,
                                              "spsrib2",
                                              &spsrib2,
                                              "spt1",
                                              &spt1,
                                              "sptime",
                                              &sptime,
                                              "sptimes",
                                              &sptimes,
                                              "sptimes1",
                                              &sptimes1,
                                              "sptimes2",
                                              &sptimes2,
                                              "sptrains",
                                              &sptrains,
                                              "srepi",
                                              &srepi,
                                              "sria1",
                                              &sria1,
                                              "sria2",
                                              &sria2,
                                              "srib1",
                                              &srib1,
                                              "srib2",
                                              &srib2,
                                              "srinopair",
                                              &srinopair,
                                              "startbins",
                                              &startbins,
                                              NULL);
                                            {
                                                mxArray * name_9
                                                  = mclInitialize(name_8);
                                                mclLoadConditional(
                                                  name_9,
                                                  "startbins2",
                                                  &startbins2,
                                                  "startlimits",
                                                  &startlimits,
                                                  "startwindow",
                                                  &startwindow,
                                                  "surri",
                                                  &surri,
                                                  "surri1",
                                                  &surri1,
                                                  "surri2",
                                                  &surri2,
                                                  "surrprefix",
                                                  &surrprefix,
                                                  "surrvec",
                                                  &surrvec,
                                                  "surrvec1",
                                                  &surrvec1,
                                                  "surrvec2",
                                                  &surrvec2,
                                                  "tmpendlimits",
                                                  &tmpendlimits,
                                                  "tmpstartlimits",
                                                  &tmpstartlimits,
                                                  "trialspikebin",
                                                  &trialspikebin,
                                                  "trialspiken",
                                                  &trialspiken,
                                                  "urepspike",
                                                  &urepspike,
                                                  NULL);
                                                mclLoadConditional(
                                                  name_9,
                                                  "urepspikea",
                                                  &urepspikea,
                                                  "urepspikeb",
                                                  &urepspikeb,
                                                  "varargin",
                                                  &varargin,
                                                  "windown",
                                                  &windown,
                                                  "xc",
                                                  &xc,
                                                  "xc1",
                                                  &xc1,
                                                  "xcmat",
                                                  &xcmat,
                                                  "xcreps",
                                                  &xcreps,
                                                  "xctemp",
                                                  &xctemp,
                                                  "xcvals",
                                                  &xcvals,
                                                  NULL);
                                                mxDestroyArray(name_9);
                                            }
                                            mxDestroyArray(name_8);
                                        }
                                        mxDestroyArray(name_7);
                                    }
                                    mxDestroyArray(name_6);
                                }
                                mxDestroyArray(name_5);
                            }
                            mxDestroyArray(name_4);
                        }
                        mxDestroyArray(name_3);
                    }
                    mxDestroyArray(name_2);
                }
                mxDestroyArray(name_1);
            }
            mxDestroyArray(name_0);
        }
        mxDestroyArray(name_);
    }
    /*
     * % load data file
     * load(cellpairdatafile);
     */
    {
        mxArray * name_
          = mclInitialize(mclVv(cellpairdatafile, "cellpairdatafile"));
        mclLoadConditional(
          name_,
          "CEN1COL",
          &CEN1COL,
          "CEN2COL",
          &CEN2COL,
          "DATACOL",
          &DATACOL,
          "EXT1COL",
          &EXT1COL,
          "EXT2COL",
          &EXT2COL,
          "SHIFTCOL",
          &SHIFTCOL,
          "a",
          &a,
          "ans",
          &ans,
          "binidx",
          &binidx,
          "binidxsize",
          &binidxsize,
          "bisi",
          &bisi,
          "blr1",
          &blr1,
          "bspikes",
          &bspikes,
          "c1l1",
          &c1l1,
          "c1lindices",
          &c1lindices,
          NULL);
        {
            mxArray * name_00 = mclInitialize(name_);
            mclLoadConditional(
              name_00,
              "c2l1",
              &c2l1,
              "c2lindices",
              &c2lindices,
              "cell1i",
              &cell1i,
              "cell1limits",
              &cell1limits,
              "cell2i",
              &cell2i,
              "cell2limits",
              &cell2limits,
              "celli",
              &celli,
              "cellpairdatafile",
              &cellpairdatafile,
              "cellvec",
              &cellvec,
              "cellvecr",
              &cellvecr,
              "cl1cell",
              &cl1cell,
              "clcell",
              &clcell,
              "clindices",
              &clindices,
              "crepspike",
              &crepspike,
              "ctrialspikebin",
              &ctrialspikebin,
              NULL);
            {
                mxArray * name_01 = mclInitialize(name_00);
                mclLoadConditional(
                  name_01,
                  "ctsbi",
                  &ctsbi,
                  "discard",
                  &discard,
                  "ec1l1",
                  &ec1l1,
                  "ec1lindices",
                  &ec1lindices,
                  "ec2l1",
                  &ec2l1,
                  "ec2lindices",
                  &ec2lindices,
                  "ecell1i",
                  &ecell1i,
                  "ecell1limits",
                  &ecell1limits,
                  "ecell2i",
                  &ecell2i,
                  "ecell2limits",
                  &ecell2limits,
                  "ecl1cell",
                  &ecl1cell,
                  "eclcell",
                  &eclcell,
                  "eclindices",
                  &eclindices,
                  "elrepn1",
                  &elrepn1,
                  "elrepn2",
                  &elrepn2,
                  NULL);
                {
                    mxArray * name_02 = mclInitialize(name_01);
                    mclLoadConditional(
                      name_02,
                      "emat1",
                      &emat1,
                      "emat2",
                      &emat2,
                      "emptytsn",
                      &emptytsn,
                      "endbins",
                      &endbins,
                      "endbins2",
                      &endbins2,
                      "endlimits",
                      &endlimits,
                      "endwindow",
                      &endwindow,
                      "ensc1",
                      &ensc1,
                      "ensc2",
                      &ensc2,
                      "enspikes",
                      &enspikes,
                      "enspikescells",
                      &enspikescells,
                      "erepnums1",
                      &erepnums1,
                      "erepnums2",
                      &erepnums2,
                      "erepnumscell",
                      &erepnumscell,
                      "espt1",
                      &espt1,
                      NULL);
                    {
                        mxArray * name_03 = mclInitialize(name_02);
                        mclLoadConditional(
                          name_03,
                          "espt2",
                          &espt2,
                          "esptime",
                          &esptime,
                          "esptimes",
                          &esptimes,
                          "esptrains",
                          &esptrains,
                          "esrepi",
                          &esrepi,
                          "esrinopair",
                          &esrinopair,
                          "eurepspike",
                          &eurepspike,
                          "eurepspikea",
                          &eurepspikea,
                          "eurepspikeb",
                          &eurepspikeb,
                          "excnsc1",
                          &excnsc1,
                          "excnsc2",
                          &excnsc2,
                          "excnspikescells",
                          &excnspikescells,
                          "excspt1",
                          &excspt1,
                          "excspt2",
                          &excspt2,
                          "excsptime",
                          &excsptime,
                          NULL);
                        {
                            mxArray * name_04 = mclInitialize(name_03);
                            mclLoadConditional(
                              name_04,
                              "excsptimes",
                              &excsptimes,
                              "excsptimes1",
                              &excsptimes1,
                              "excsptimes2",
                              &excsptimes2,
                              "excvals1",
                              &excvals1,
                              "excvals2",
                              &excvals2,
                              "fmax",
                              &fmax,
                              "fmax2",
                              &fmax2,
                              "fmin",
                              &fmin,
                              "fmin2",
                              &fmin2,
                              "intrepnums1",
                              &intrepnums1,
                              "intrepnums2",
                              &intrepnums2,
                              "lagbins2",
                              &lagbins2,
                              "lrepn1",
                              &lrepn1,
                              "lrepn2",
                              &lrepn2,
                              "mat1",
                              &mat1,
                              NULL);
                            {
                                mxArray * name_05 = mclInitialize(name_04);
                                mclLoadConditional(
                                  name_05,
                                  "mat2",
                                  &mat2,
                                  "needrefcheck",
                                  &needrefcheck,
                                  "nexcvals1",
                                  &nexcvals1,
                                  "nexcvals2",
                                  &nexcvals2,
                                  "nfc",
                                  &nfc,
                                  "nisi0",
                                  &nisi0,
                                  "noskip1",
                                  &noskip1,
                                  "noskip2",
                                  &noskip2,
                                  "nref",
                                  &nref,
                                  "nsc1",
                                  &nsc1,
                                  "nsc2",
                                  &nsc2,
                                  "nspikes",
                                  &nspikes,
                                  "nspikescells",
                                  &nspikescells,
                                  "numSurrogates1",
                                  &numSurrogates1,
                                  "numSurrogates2",
                                  &numSurrogates2,
                                  NULL);
                                {
                                    mxArray * name_06 = mclInitialize(name_05);
                                    mclLoadConditional(
                                      name_06,
                                      "numwindows",
                                      &numwindows,
                                      "nxcvals",
                                      &nxcvals,
                                      "randnums",
                                      &randnums,
                                      "repi",
                                      &repi,
                                      "repi2",
                                      &repi2,
                                      "repn",
                                      &repn,
                                      "repnl",
                                      &repnl,
                                      "repnums1",
                                      &repnums1,
                                      "repnums2",
                                      &repnums2,
                                      "repnumscell",
                                      &repnumscell,
                                      "repnvec",
                                      &repnvec,
                                      "reps",
                                      &reps,
                                      "repsboth",
                                      &repsboth,
                                      "repsnopair",
                                      &repsnopair,
                                      "repspike",
                                      &repspike,
                                      NULL);
                                    {
                                        mxArray * name_07
                                          = mclInitialize(name_06);
                                        mclLoadConditional(
                                          name_07,
                                          "repswspike1",
                                          &repswspike1,
                                          "repswspike2",
                                          &repswspike2,
                                          "sdatafile",
                                          &sdatafile,
                                          "shuffle",
                                          &shuffle,
                                          "shufflevecl",
                                          &shufflevecl,
                                          "sperepnums2",
                                          &sperepnums2,
                                          "spiketrains",
                                          &spiketrains,
                                          "spintrepnums1",
                                          &spintrepnums1,
                                          "spintrepnums2",
                                          &spintrepnums2,
                                          "sprepnl",
                                          &sprepnl,
                                          "sprepnums2",
                                          &sprepnums2,
                                          "sprepnvec",
                                          &sprepnvec,
                                          "sprepsboth",
                                          &sprepsboth,
                                          "spsria1",
                                          &spsria1,
                                          "spsria2",
                                          &spsria2,
                                          NULL);
                                        {
                                            mxArray * name_08
                                              = mclInitialize(name_07);
                                            mclLoadConditional(
                                              name_08,
                                              "spsrib1",
                                              &spsrib1,
                                              "spsrib2",
                                              &spsrib2,
                                              "spt1",
                                              &spt1,
                                              "sptime",
                                              &sptime,
                                              "sptimes",
                                              &sptimes,
                                              "sptimes1",
                                              &sptimes1,
                                              "sptimes2",
                                              &sptimes2,
                                              "sptrains",
                                              &sptrains,
                                              "srepi",
                                              &srepi,
                                              "sria1",
                                              &sria1,
                                              "sria2",
                                              &sria2,
                                              "srib1",
                                              &srib1,
                                              "srib2",
                                              &srib2,
                                              "srinopair",
                                              &srinopair,
                                              "startbins",
                                              &startbins,
                                              NULL);
                                            {
                                                mxArray * name_09
                                                  = mclInitialize(name_08);
                                                mclLoadConditional(
                                                  name_09,
                                                  "startbins2",
                                                  &startbins2,
                                                  "startlimits",
                                                  &startlimits,
                                                  "startwindow",
                                                  &startwindow,
                                                  "surri",
                                                  &surri,
                                                  "surri1",
                                                  &surri1,
                                                  "surri2",
                                                  &surri2,
                                                  "surrprefix",
                                                  &surrprefix,
                                                  "surrvec",
                                                  &surrvec,
                                                  "surrvec1",
                                                  &surrvec1,
                                                  "surrvec2",
                                                  &surrvec2,
                                                  "tmpendlimits",
                                                  &tmpendlimits,
                                                  "tmpstartlimits",
                                                  &tmpstartlimits,
                                                  "trialspikebin",
                                                  &trialspikebin,
                                                  "trialspiken",
                                                  &trialspiken,
                                                  "urepspike",
                                                  &urepspike,
                                                  NULL);
                                                mclLoadConditional(
                                                  name_09,
                                                  "urepspikea",
                                                  &urepspikea,
                                                  "urepspikeb",
                                                  &urepspikeb,
                                                  "varargin",
                                                  &varargin,
                                                  "windown",
                                                  &windown,
                                                  "xc",
                                                  &xc,
                                                  "xc1",
                                                  &xc1,
                                                  "xcmat",
                                                  &xcmat,
                                                  "xcreps",
                                                  &xcreps,
                                                  "xctemp",
                                                  &xctemp,
                                                  "xcvals",
                                                  &xcvals,
                                                  NULL);
                                                mxDestroyArray(name_09);
                                            }
                                            mxDestroyArray(name_08);
                                        }
                                        mxDestroyArray(name_07);
                                    }
                                    mxDestroyArray(name_06);
                                }
                                mxDestroyArray(name_05);
                            }
                            mxDestroyArray(name_04);
                        }
                        mxDestroyArray(name_03);
                    }
                    mxDestroyArray(name_02);
                }
                mxDestroyArray(name_01);
            }
            mxDestroyArray(name_00);
        }
        mxDestroyArray(name_);
    }
    /*
     * 
     * % needed to compile standalone executable
     * if(nargin>1)
     */
    if (nargin_ > 1) {
        /*
         * startwindow = varargin{1};
         */
        mlfAssign(
          &startwindow,
          mlfIndexRef(mclVa(varargin, "varargin"), "{?}", _mxarray0_));
        /*
         * if(ischar(startwindow))
         */
        if (mlfTobool(mlfIschar(mclVv(startwindow, "startwindow")))) {
            /*
             * startwindow = str2num(startwindow);
             */
            mlfAssign(
              &startwindow,
              mlfStr2num(NULL, mclVv(startwindow, "startwindow")));
        /*
         * end
         */
        }
        /*
         * endwindow = varargin{2};
         */
        mlfAssign(
          &endwindow,
          mlfIndexRef(mclVa(varargin, "varargin"), "{?}", _mxarray1_));
        /*
         * if(ischar(endwindow))
         */
        if (mlfTobool(mlfIschar(mclVv(endwindow, "endwindow")))) {
            /*
             * endwindow = str2num(endwindow);
             */
            mlfAssign(
              &endwindow, mlfStr2num(NULL, mclVv(endwindow, "endwindow")));
        /*
         * end
         */
        }
    /*
     * else
     */
    } else {
        /*
         * startwindow = 1;
         */
        mlfAssign(&startwindow, _mxarray0_);
        /*
         * endwindow = numwindows;
         */
        mlfAssign(&endwindow, mclVv(numwindows, "numwindows"));
    /*
     * end
     */
    }
    /*
     * 
     * for windown = startwindow:endwindow
     */
    {
        mclForLoopIterator viter__;
        for (mclForStart(
               &viter__,
               mclVv(startwindow, "startwindow"),
               mclVv(endwindow, "endwindow"),
               NULL);
             mclForNext(&viter__, &windown);
             ) {
            /*
             * % get bins corresponding to window n
             * % use startbins2 and endbins2 so we can eliminate boundary effects
             * fmin = startbins2(windown) - 1;
             */
            mlfAssign(
              &fmin,
              mclMinus(
                mclArrayRef1(
                  mclVv(startbins2, "startbins2"), mclVv(windown, "windown")),
                _mxarray0_));
            /*
             * fmax = endbins2(windown) + 1;
             */
            mlfAssign(
              &fmax,
              mclPlus(
                mclArrayRef1(
                  mclVv(endbins2, "endbins2"), mclVv(windown, "windown")),
                _mxarray0_));
            /*
             * % find indices corresponding to window n
             * % subtract and add 1 to fmin and fmax so we don't have to use <= and >=
             * [trialspiken,repspike] = find(binidx>fmin & binidx<fmax);
             */
            mlfAssign(
              &trialspiken,
              mlfFind(
                &repspike,
                NULL,
                mclAnd(
                  mclGt(mclVv(binidx, "binidx"), mclVv(fmin, "fmin")),
                  mclLt(mclVv(binidx, "binidx"), mclVv(fmax, "fmax")))));
            /*
             * % check for empty trialspiken so we don't try to use sub2ind which will
             * % fail with empty matrices
             * emptytsn = isempty(trialspiken);
             */
            mlfAssign(&emptytsn, mlfIsempty(mclVv(trialspiken, "trialspiken")));
            /*
             * if(~emptytsn)
             */
            if (mclNotBool(mclVv(emptytsn, "emptytsn"))) {
                /*
                 * % grab the actual binidx values. Can't use the values returned by find 
                 * % since it will be all 1's since the matrix passed to find is binary		
                 * trialspikebin = binidx(sub2ind(binidxsize,trialspiken,repspike));
                 */
                mlfAssign(
                  &trialspikebin,
                  mclArrayRef1(
                    mclVv(binidx, "binidx"),
                    mlfSub2ind(
                      mclVv(binidxsize, "binidxsize"),
                      mclVv(trialspiken, "trialspiken"),
                      mclVv(repspike, "repspike"),
                      NULL)));
                /*
                 * % check to see if no indices were found, which means that there were no
                 * % spikes in either cell in this frame, or if there were no repetitions where
                 * % there were spikes in both cells in this frame
                 * % we first do unique to identify repetitions with spikes
                 * % then we subtract 1 and take the mod using the number of repetitions, e.g. if
                 * % there were 100 reps, columns 1 to 100 become 0 to 99 and columns 101 to 200 
                 * % become 0 to 99
                 * % then we sort the output so that if there were at least one spike in both cells
                 * % we will have a repeated repetition number, which will show up as a 0 in the
                 * % diff
                 * % find the unique cols (i.e. repetition)
                 * [eurepspike,eurepspikea,eurepspikeb] = unique(repspike);    
                 */
                mlfAssign(
                  &eurepspike,
                  mlfNUnique(
                    3,
                    &eurepspikea,
                    &eurepspikeb,
                    mclVv(repspike, "repspike"),
                    NULL));
                /*
                 * % find the reps belonging to the second cell
                 * ecell2i = eurepspike>reps;
                 */
                mlfAssign(
                  &ecell2i,
                  mclGt(mclVv(eurepspike, "eurepspike"), mclVv(reps, "reps")));
                /*
                 * ecell1i = ~ecell2i;
                 */
                mlfAssign(&ecell1i, mclNot(mclVv(ecell2i, "ecell2i")));
                /*
                 * erepnums1 = eurepspike(ecell1i);
                 */
                mlfAssign(
                  &erepnums1,
                  mclArrayRef1(
                    mclVv(eurepspike, "eurepspike"),
                    mclVv(ecell1i, "ecell1i")));
                /*
                 * erepnums2 = eurepspike(ecell2i) - reps;
                 */
                mlfAssign(
                  &erepnums2,
                  mclMinus(
                    mclArrayRef1(
                      mclVv(eurepspike, "eurepspike"),
                      mclVv(ecell2i, "ecell2i")),
                    mclVv(reps, "reps")));
                /*
                 * 
                 * % find spikes in central window
                 * fmin2 = startbins(windown) - 1;
                 */
                mlfAssign(
                  &fmin2,
                  mclMinus(
                    mclArrayRef1(
                      mclVv(startbins, "startbins"), mclVv(windown, "windown")),
                    _mxarray0_));
                /*
                 * fmax2 = endbins(windown) + 1;
                 */
                mlfAssign(
                  &fmax2,
                  mclPlus(
                    mclArrayRef1(
                      mclVv(endbins, "endbins"), mclVv(windown, "windown")),
                    _mxarray0_));
                /*
                 * ctsbi = (trialspikebin>fmin2 & trialspikebin<fmax2);
                 */
                mlfAssign(
                  &ctsbi,
                  mclAnd(
                    mclGt(
                      mclVv(trialspikebin, "trialspikebin"),
                      mclVv(fmin2, "fmin2")),
                    mclLt(
                      mclVv(trialspikebin, "trialspikebin"),
                      mclVv(fmax2, "fmax2"))));
                /*
                 * % grab the actual binidx values. Can't use the values returned by find 
                 * % since it will be all 1's since the matrix passed to find is binary		
                 * ctrialspikebin = trialspikebin(ctsbi);
                 */
                mlfAssign(
                  &ctrialspikebin,
                  mclArrayRef1(
                    mclVv(trialspikebin, "trialspikebin"),
                    mclVv(ctsbi, "ctsbi")));
                /*
                 * crepspike = repspike(ctsbi);
                 */
                mlfAssign(
                  &crepspike,
                  mclArrayRef1(
                    mclVv(repspike, "repspike"), mclVv(ctsbi, "ctsbi")));
                /*
                 * % use the above indices to get bin numbers
                 * % can't use the values returned by find since the array passed to
                 * % find is a logical array
                 * [urepspike,urepspikea,urepspikeb] = unique(crepspike);
                 */
                mlfAssign(
                  &urepspike,
                  mlfNUnique(
                    3,
                    &urepspikea,
                    &urepspikeb,
                    mclVv(crepspike, "crepspike"),
                    NULL));
                /*
                 * cell2i = urepspike>reps;
                 */
                mlfAssign(
                  &cell2i,
                  mclGt(mclVv(urepspike, "urepspike"), mclVv(reps, "reps")));
                /*
                 * cell1i = ~cell2i;
                 */
                mlfAssign(&cell1i, mclNot(mclVv(cell2i, "cell2i")));
                /*
                 * repnums1 = urepspike(cell1i);
                 */
                mlfAssign(
                  &repnums1,
                  mclArrayRef1(
                    mclVv(urepspike, "urepspike"), mclVv(cell1i, "cell1i")));
                /*
                 * repnums2 = urepspike(cell2i) - reps;
                 */
                mlfAssign(
                  &repnums2,
                  mclMinus(
                    mclArrayRef1(
                      mclVv(urepspike, "urepspike"), mclVv(cell2i, "cell2i")),
                    mclVv(reps, "reps")));
                /*
                 * 
                 * % find the reps that have spikes in the central window of at least 1 cell
                 * % can't use union since it will include reps with spikes in the central
                 * % window for 1 cell and no spikes in the extended window for the other
                 * [intrepnums1,sria1,srib2] = intersect(erepnums1,repnums2);
                 */
                mlfAssign(
                  &intrepnums1,
                  mlfNIntersect(
                    3,
                    &sria1,
                    &srib2,
                    mclVv(erepnums1, "erepnums1"),
                    mclVv(repnums2, "repnums2"),
                    NULL));
                /*
                 * [intrepnums2,srib1,sria2] = intersect(erepnums2,repnums1);
                 */
                mlfAssign(
                  &intrepnums2,
                  mlfNIntersect(
                    3,
                    &srib1,
                    &sria2,
                    mclVv(erepnums2, "erepnums2"),
                    mclVv(repnums1, "repnums1"),
                    NULL));
                /*
                 * repsboth = unique([intrepnums1; intrepnums2]);
                 */
                mlfAssign(
                  &repsboth,
                  mlfNUnique(
                    1,
                    NULL,
                    NULL,
                    mlfVertcat(
                      mclVv(intrepnums1, "intrepnums1"),
                      mclVv(intrepnums2, "intrepnums2"),
                      NULL),
                    NULL));
                /*
                 * % save rep numbers so we can compute synchrony only on these reps
                 * repnvec = vecr(repsboth);
                 */
                mlfAssign(&repnvec, mlfVecr(mclVv(repsboth, "repsboth")));
                /*
                 * repnl = length(repnvec);
                 */
                mlfAssign(
                  &repnl, mlfScalar(mclLengthInt(mclVv(repnvec, "repnvec"))));
            /*
             * end
             */
            }
            /*
             * if(emptytsn || (repnl==0))
             */
            if (mclScalarToBool(mclVv(emptytsn, "emptytsn"))
                || mclScalarToBool(mclEq(mclVv(repnl, "repnl"), _mxarray4_))) {
                /*
                 * % no pairs of spikes found in any rep so just save matrix of all zeros
                 * % don't use matrix of nan's since a matrix of zeros is more accurate
                 * xc = xctemp;
                 */
                mlfAssign(&xc, mclVv(xctemp, "xctemp"));
                /*
                 * sptrains = {};
                 */
                mlfAssign(&sptrains, _mxarray5_);
                /*
                 * repswspike1 = [];
                 */
                mlfAssign(&repswspike1, _mxarray6_);
                /*
                 * repswspike2 = [];
                 */
                mlfAssign(&repswspike2, _mxarray6_);
            /*
             * else % if(~isempty(trialspiken))
             */
            } else {
                /*
                 * xcreps = cell(repnl,1);
                 */
                mlfAssign(
                  &xcreps, mlfCell(mclVv(repnl, "repnl"), _mxarray0_, NULL));
                /*
                 * % get number of repetitions with spikes in 1st cell
                 * lrepn1 = length(repnums1);
                 */
                mlfAssign(
                  &lrepn1,
                  mlfScalar(mclLengthInt(mclVv(repnums1, "repnums1"))));
                /*
                 * % get number of repetitions with spikes in 2nd cell
                 * lrepn2 = length(repnums2);
                 */
                mlfAssign(
                  &lrepn2,
                  mlfScalar(mclLengthInt(mclVv(repnums2, "repnums2"))));
                /*
                 * % allocate memory for the required number of repetitions
                 * sptrains = spiketrains;
                 */
                mlfAssign(&sptrains, mclVv(spiketrains, "spiketrains"));
                /*
                 * sptrains{1} = cell(lrepn1,1);
                 */
                mlfIndexAssign(
                  &sptrains,
                  "{?}",
                  _mxarray0_,
                  mlfCell(mclVv(lrepn1, "lrepn1"), _mxarray0_, NULL));
                /*
                 * sptrains{2} = cell(lrepn2,1);
                 */
                mlfIndexAssign(
                  &sptrains,
                  "{?}",
                  _mxarray1_,
                  mlfCell(mclVv(lrepn2, "lrepn2"), _mxarray0_, NULL));
                /*
                 * % get number of repetitions with spikes in the extended window of 
                 * % the 1st cell
                 * elrepn1 = length(erepnums1);
                 */
                mlfAssign(
                  &elrepn1,
                  mlfScalar(mclLengthInt(mclVv(erepnums1, "erepnums1"))));
                /*
                 * % get number of repetitions with spikes in 2nd cell
                 * elrepn2 = length(erepnums2);
                 */
                mlfAssign(
                  &elrepn2,
                  mlfScalar(mclLengthInt(mclVv(erepnums2, "erepnums2"))));
                /*
                 * % allocate memory for the required number of repetitions
                 * esptrains = spiketrains;
                 */
                mlfAssign(&esptrains, mclVv(spiketrains, "spiketrains"));
                /*
                 * esptrains{1} = cell(elrepn1,1);
                 */
                mlfIndexAssign(
                  &esptrains,
                  "{?}",
                  _mxarray0_,
                  mlfCell(mclVv(elrepn1, "elrepn1"), _mxarray0_, NULL));
                /*
                 * esptrains{2} = cell(elrepn2,1);
                 */
                mlfIndexAssign(
                  &esptrains,
                  "{?}",
                  _mxarray1_,
                  mlfCell(mclVv(elrepn2, "elrepn2"), _mxarray0_, NULL));
                /*
                 * 
                 * % grab the limits for each group from the output of unique
                 * ecell1limits = eurepspikea(ecell1i);
                 */
                mlfAssign(
                  &ecell1limits,
                  mclArrayRef1(
                    mclVv(eurepspikea, "eurepspikea"),
                    mclVv(ecell1i, "ecell1i")));
                /*
                 * ecell2limits = eurepspikea(ecell2i);
                 */
                mlfAssign(
                  &ecell2limits,
                  mclArrayRef1(
                    mclVv(eurepspikea, "eurepspikea"),
                    mclVv(ecell2i, "ecell2i")));
                /*
                 * % add 1 to cell1limits and cell2limits so we don't have to add 1
                 * % inside the for-loop
                 * ec1l1 = [0; ecell1limits(1:(elrepn1-1))] + 1;
                 */
                mlfAssign(
                  &ec1l1,
                  mclPlus(
                    mlfVertcat(
                      _mxarray4_,
                      mclArrayRef1(
                        mclVv(ecell1limits, "ecell1limits"),
                        mlfColon(
                          _mxarray0_,
                          mclMinus(mclVv(elrepn1, "elrepn1"), _mxarray0_),
                          NULL)),
                      NULL),
                    _mxarray0_));
                /*
                 * ec2l1 = [ecell1limits(elrepn1); ecell2limits(1:(elrepn2-1))] + 1;
                 */
                mlfAssign(
                  &ec2l1,
                  mclPlus(
                    mlfVertcat(
                      mclArrayRef1(
                        mclVv(ecell1limits, "ecell1limits"),
                        mclVv(elrepn1, "elrepn1")),
                      mclArrayRef1(
                        mclVv(ecell2limits, "ecell2limits"),
                        mlfColon(
                          _mxarray0_,
                          mclMinus(mclVv(elrepn2, "elrepn2"), _mxarray0_),
                          NULL)),
                      NULL),
                    _mxarray0_));
                /*
                 * % grab the limits for each group from the output of unique
                 * cell1limits = urepspikea(cell1i);
                 */
                mlfAssign(
                  &cell1limits,
                  mclArrayRef1(
                    mclVv(urepspikea, "urepspikea"), mclVv(cell1i, "cell1i")));
                /*
                 * cell2limits = urepspikea(cell2i);
                 */
                mlfAssign(
                  &cell2limits,
                  mclArrayRef1(
                    mclVv(urepspikea, "urepspikea"), mclVv(cell2i, "cell2i")));
                /*
                 * % add 1 to cell1limits and cell2limits so we don't have to add 1
                 * % inside the for-loop
                 * blr1 = lrepn1 > 0;
                 */
                mlfAssign(&blr1, mclGt(mclVv(lrepn1, "lrepn1"), _mxarray4_));
                /*
                 * if(blr1 && lrepn2>0)
                 */
                if (mclScalarToBool(mclVv(blr1, "blr1"))
                    && mclScalarToBool(
                         mclGt(mclVv(lrepn2, "lrepn2"), _mxarray4_))) {
                    /*
                     * c1l1 = [0; cell1limits(1:(lrepn1-1))] + 1;
                     */
                    mlfAssign(
                      &c1l1,
                      mclPlus(
                        mlfVertcat(
                          _mxarray4_,
                          mclArrayRef1(
                            mclVv(cell1limits, "cell1limits"),
                            mlfColon(
                              _mxarray0_,
                              mclMinus(mclVv(lrepn1, "lrepn1"), _mxarray0_),
                              NULL)),
                          NULL),
                        _mxarray0_));
                    /*
                     * c2l1 = [cell1limits(lrepn1); cell2limits(1:(lrepn2-1))] + 1;
                     */
                    mlfAssign(
                      &c2l1,
                      mclPlus(
                        mlfVertcat(
                          mclArrayRef1(
                            mclVv(cell1limits, "cell1limits"),
                            mclVv(lrepn1, "lrepn1")),
                          mclArrayRef1(
                            mclVv(cell2limits, "cell2limits"),
                            mlfColon(
                              _mxarray0_,
                              mclMinus(mclVv(lrepn2, "lrepn2"), _mxarray0_),
                              NULL)),
                          NULL),
                        _mxarray0_));
                /*
                 * elseif(blr1)
                 */
                } else if (mlfTobool(mclVv(blr1, "blr1"))) {
                    /*
                     * c1l1 = [0; cell1limits(1:(lrepn1-1))] + 1;
                     */
                    mlfAssign(
                      &c1l1,
                      mclPlus(
                        mlfVertcat(
                          _mxarray4_,
                          mclArrayRef1(
                            mclVv(cell1limits, "cell1limits"),
                            mlfColon(
                              _mxarray0_,
                              mclMinus(mclVv(lrepn1, "lrepn1"), _mxarray0_),
                              NULL)),
                          NULL),
                        _mxarray0_));
                    /*
                     * c2l1 = [];
                     */
                    mlfAssign(&c2l1, _mxarray6_);
                /*
                 * else
                 */
                } else {
                    /*
                     * c1l1 = [];
                     */
                    mlfAssign(&c1l1, _mxarray6_);
                    /*
                     * c2l1 = [0; cell2limits(1:(lrepn2-1))] + 1;
                     */
                    mlfAssign(
                      &c2l1,
                      mclPlus(
                        mlfVertcat(
                          _mxarray4_,
                          mclArrayRef1(
                            mclVv(cell2limits, "cell2limits"),
                            mlfColon(
                              _mxarray0_,
                              mclMinus(mclVv(lrepn2, "lrepn2"), _mxarray0_),
                              NULL)),
                          NULL),
                        _mxarray0_));
                /*
                 * end
                 */
                }
                /*
                 * 
                 * % set the limits for the indices in trialspikebin
                 * % this array will keep track of which windows are empty so we don't
                 * % have to add if statements inside the loop to make sure things work
                 * startlimits = tmpstartlimits;
                 */
                mlfAssign(
                  &startlimits, mclVv(tmpstartlimits, "tmpstartlimits"));
                /*
                 * endlimits = tmpendlimits;
                 */
                mlfAssign(&endlimits, mclVv(tmpendlimits, "tmpendlimits"));
                /*
                 * % set start and end limits for extended window in cell 1
                 * startlimits(erepnums1,EXT1COL) = ec1l1;
                 */
                mclArrayAssign2(
                  &startlimits,
                  mclVv(ec1l1, "ec1l1"),
                  mclVv(erepnums1, "erepnums1"),
                  mclVv(EXT1COL, "EXT1COL"));
                /*
                 * endlimits(erepnums1,EXT1COL) = ecell1limits;
                 */
                mclArrayAssign2(
                  &endlimits,
                  mclVv(ecell1limits, "ecell1limits"),
                  mclVv(erepnums1, "erepnums1"),
                  mclVv(EXT1COL, "EXT1COL"));
                /*
                 * % set start and end limits for extended window in cell 2
                 * startlimits(erepnums2,EXT2COL) = ec2l1;
                 */
                mclArrayAssign2(
                  &startlimits,
                  mclVv(ec2l1, "ec2l1"),
                  mclVv(erepnums2, "erepnums2"),
                  mclVv(EXT2COL, "EXT2COL"));
                /*
                 * endlimits(erepnums2,EXT2COL) = ecell2limits;
                 */
                mclArrayAssign2(
                  &endlimits,
                  mclVv(ecell2limits, "ecell2limits"),
                  mclVv(erepnums2, "erepnums2"),
                  mclVv(EXT2COL, "EXT2COL"));
                /*
                 * % set start and end limits for central window in cell 1
                 * startlimits(repnums1,CEN1COL) = c1l1;
                 */
                mclArrayAssign2(
                  &startlimits,
                  mclVv(c1l1, "c1l1"),
                  mclVv(repnums1, "repnums1"),
                  mclVv(CEN1COL, "CEN1COL"));
                /*
                 * endlimits(repnums1,CEN1COL) = cell1limits;
                 */
                mclArrayAssign2(
                  &endlimits,
                  mclVv(cell1limits, "cell1limits"),
                  mclVv(repnums1, "repnums1"),
                  mclVv(CEN1COL, "CEN1COL"));
                /*
                 * % set start and end limits for central window in cell 2
                 * startlimits(repnums2,CEN2COL) = c2l1;
                 */
                mclArrayAssign2(
                  &startlimits,
                  mclVv(c2l1, "c2l1"),
                  mclVv(repnums2, "repnums2"),
                  mclVv(CEN2COL, "CEN2COL"));
                /*
                 * endlimits(repnums2,CEN2COL) = cell2limits;
                 */
                mclArrayAssign2(
                  &endlimits,
                  mclVv(cell2limits, "cell2limits"),
                  mclVv(repnums2, "repnums2"),
                  mclVv(CEN2COL, "CEN2COL"));
                /*
                 * 
                 * % initialize random number generator which is capable of generating 
                 * % 2^1492 values so a typical run uses 12x200x3x1000=7200000 random
                 * % numbers so we won't have to reset the state between calls to randn
                 * rand('state',sum(windown*clock));
                 */
                mclAssignAns(
                  &ans,
                  mlfNRand(
                    0,
                    _mxarray7_,
                    mlfSum(
                      mclMtimes(mclVv(windown, "windown"), mlfClock()), NULL),
                    NULL));
                /*
                 * 
                 * for repi = 1:repnl
                 */
                {
                    int v_ = mclForIntStart(1);
                    int e_ = mclForIntEnd(mclVv(repnl, "repnl"));
                    if (v_ > e_) {
                        mlfAssign(&repi, _mxarray6_);
                    } else {
                        /*
                         * % get rep number
                         * repn = repsboth(repi);
                         * % get the spike times for the central window
                         * c1lindices = startlimits(repn,CEN1COL):endlimits(repn,CEN1COL);
                         * c2lindices = startlimits(repn,CEN2COL):endlimits(repn,CEN2COL);
                         * sptimes{2} = ctrialspikebin(c2lindices);
                         * sptimes{1} = ctrialspikebin(c1lindices);
                         * nsc1 = length(c1lindices);
                         * nsc2 = length(c2lindices);
                         * nspikescells = [nsc1 nsc2];
                         * 
                         * % get the spike times for the extended window
                         * ec1lindices = startlimits(repn,EXT1COL):endlimits(repn,EXT1COL);
                         * ec2lindices = startlimits(repn,EXT2COL):endlimits(repn,EXT2COL);
                         * espt1 = trialspikebin(ec1lindices);
                         * espt2 = trialspikebin(ec2lindices);
                         * esptimes{2} = espt2;
                         * esptimes{1} = espt1;
                         * ensc1 = length(espt1);
                         * ensc2 = length(espt2);
                         * enspikescells = [ensc1 ensc2];
                         * 
                         * % get the spike times for the extended window after removing the
                         * % spikes in the central window. This makes computing the xcorr 
                         * % easier since we don't have to subtract out the double-counting 
                         * % of the central spikes
                         * excspt1 = setdiff(esptimes{1},sptimes{1});
                         * excspt2 = setdiff(esptimes{2},sptimes{2});
                         * excsptimes{2} = excspt2;
                         * excsptimes{1} = excspt1;
                         * excnsc1 = length(excspt1);
                         * excnsc2 = length(excspt2);
                         * excnspikescells = [excnsc1 excnsc2];
                         * 
                         * % check spike times
                         * %             concat([full([trialspikebin(ec1lindices) trialspiken(ec1lindices) repspike(ec1lindices)])], ...
                         * %             		[full([trialspikebin(ec2lindices) trialspiken(ec2lindices) repspike(ec2lindices)])], ...
                         * %             		[full([ctrialspikebin(c1lindices) crepspike(c1lindices)])], ...
                         * %             		[full([ctrialspikebin(c2lindices) crepspike(c2lindices)])],'Columnwise')
                         * 
                         * % create matrix to compute cross-correlation between central 
                         * % windows of cells 1 and 2
                         * mat1 = [sptimes{1}(:) repmat(-1,nsc1,1)];
                         * mat2 = [ones(1,nsc2); sptimes{2}(:)'];
                         * % compute xcorr for the data
                         * xcvals = mat1 * mat2;
                         * nxcvals = nsc1 * nsc2;
                         * 
                         * % compute xcorr between central window of cell 1 and extended
                         * % window of cell 2
                         * emat2 = [ones(1,excnsc2); excsptimes{2}(:)'];
                         * excvals1 = mat1 * emat2;
                         * nexcvals1 = nsc1 * excnsc2;
                         * 
                         * % compute xcorr between central window of cell 2 and extended
                         * % window of cell 1
                         * emat1 = [excsptimes{1}(:) repmat(-1,excnsc1,1)];
                         * excvals2 = emat1 * mat2;
                         * nexcvals2 = excnsc1 * nsc2;
                         * 
                         * % initialize with nan since 0's mean 0 time difference or
                         * % synchronous spikes
                         * xc1 = repmat(nan,nxcvals+nexcvals1+nexcvals2,numSurrogates2);
                         * xc1(:,DATACOL) = [xcvals(:); excvals1(:); excvals2(:)];
                         * % initialize cellvecr to cellvec
                         * cellvecr = cellvec;
                         * % create array for storing spike times
                         * if(nsc1==0)
                         * % sptrains{1}{repi} = [];
                         * cellvecr = 2;
                         * spiketrains{1} = [];
                         * noskip1 = 0;
                         * else
                         * sptrains{1}{repi} = zeros(nsc1,numSurrogates1);
                         * % save the spike times of the data for the psth analysis
                         * sptrains{1}{repi}(:,DATACOL) = sptimes{1}(:);
                         * noskip1 = 1;
                         * end
                         * if(nsc2==0)
                         * % sptrains{2}{repi} = [];
                         * cellvecr = 1;
                         * spiketrains{2} = [];
                         * noskip2 = 0;
                         * else
                         * sptrains{2}{repi} = zeros(nsc2,numSurrogates1);
                         * % save the spike times of the data for the psth analysis
                         * sptrains{2}{repi}(:,DATACOL) = sptimes{2}(:);
                         * noskip2 = 1;
                         * end
                         * % set flags to see if we can skip the while loop
                         * % do this outside the for surri = surrvec loop so we don't do
                         * % it over and over again
                         * nfc = needrefcheck;
                         * for celli = cellvecr
                         * % only do while loop if there is more than 1 spike
                         * % and the intervals between the spikes is such that 
                         * % there will be overlaps
                         * bspikes = enspikescells(celli) > 1;
                         * bisi = ~isempty(find(diff(esptimes{celli})<shufflevecl));
                         * % check if there are isi's smaller than 1 in the extended
                         * % window
                         * nref = length(find(diff(excsptimes{celli})<1));
                         * if( bspikes && bisi)
                         * nfc(celli) = 1;
                         * end
                         * end
                         * 
                         * % loop over number of surrogates
                         * for surri = surrvec
                         * surri1 = surrvec1(surri);
                         * surri2 = surrvec2(surri);
                         * for celli = cellvecr
                         * % check if we need to check for refractory period 
                         * % violations
                         * if(nfc(celli))
                         * % make sure we go into the while loop
                         * nisi0 = 1;
                         * % check for refractory period violations
                         * while(nisi0)
                         * % get random numbers for cell i
                         * randnums = floor(rand(nspikescells(celli),1) ...
                         * * shufflevecl) - shuffle;
                         * % add random numbers to spike times for cell i and
                         * % sort the spike times
                         * spt1 = sort(sptimes{celli}+randnums);
                         * % take the diff and find 0's and then subtract
                         * % the pre-existing 0's to see if they were 
                         * % introduced by the surrogates
                         * % no need to worry about negative numbers since
                         * % the first number is always going to be larger
                         * % or equal to nref
                         * nisi0 = length(find(diff(spt1)==0)) - nref;
                         * end % while(nisi0)
                         * spiketrains{celli} = spt1;
                         * else % if(nfc(celli))
                         * % don't have to worry about overlap between the 
                         * % shuffled spikes so just add the shuffle to the
                         * % spike times
                         * spiketrains{celli} = sort(sptimes{celli} + ...
                         * floor(rand(nspikescells(celli),1) ...
                         * * shufflevecl) - shuffle);
                         * end % if(nfc(celli))
                         * end % for celli = cellvec
                         * % now set the first cell's spiketrains
                         * if(noskip1)
                         * mat1(:,1) = spiketrains{1}(:);
                         * end
                         * % now set the second cell's spiketrains
                         * if(noskip2)
                         * mat2(2,:) = spiketrains{2}(:)';
                         * end
                         * % now compute the cross-correlation between the first and
                         * % second cell
                         * xcvals = mat1 * mat2;                
                         * % compute xcorr between central window of cell 1 and extended
                         * % window of cell 2
                         * excvals1 = mat1 * emat2;
                         * % compute xcorr between central window of cell 2 and extended
                         * % window of cell 1
                         * excvals2 = emat1 * mat2;
                         * 
                         * % save the differences so we can do a histogram at the end
                         * % instead of after every surrogate
                         * xc1(:,surri2) = [xcvals(:); excvals1(:); excvals2(:)];
                         * % store the spike times for psth analysis
                         * if(noskip1)
                         * sptrains{1}{repi}(:,surri1) = spiketrains{1}(:);
                         * end
                         * if(noskip2)
                         * sptrains{2}{repi}(:,surri1) = spiketrains{2}(:);
                         * end
                         * end % loop over number of surrogates
                         * xcreps{repi} = xc1;
                         * end % for repi = 1:repnl
                         */
                        for (; ; ) {
                            mlfAssign(
                              &repn,
                              mclIntArrayRef1(mclVv(repsboth, "repsboth"), v_));
                            mlfAssign(
                              &c1lindices,
                              mlfColon(
                                mclArrayRef2(
                                  mclVv(startlimits, "startlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(CEN1COL, "CEN1COL")),
                                mclArrayRef2(
                                  mclVv(endlimits, "endlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(CEN1COL, "CEN1COL")),
                                NULL));
                            mlfAssign(
                              &c2lindices,
                              mlfColon(
                                mclArrayRef2(
                                  mclVv(startlimits, "startlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(CEN2COL, "CEN2COL")),
                                mclArrayRef2(
                                  mclVv(endlimits, "endlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(CEN2COL, "CEN2COL")),
                                NULL));
                            mlfIndexAssign(
                              &sptimes,
                              "{?}",
                              _mxarray1_,
                              mclArrayRef1(
                                mclVv(ctrialspikebin, "ctrialspikebin"),
                                mclVv(c2lindices, "c2lindices")));
                            mlfIndexAssign(
                              &sptimes,
                              "{?}",
                              _mxarray0_,
                              mclArrayRef1(
                                mclVv(ctrialspikebin, "ctrialspikebin"),
                                mclVv(c1lindices, "c1lindices")));
                            mlfAssign(
                              &nsc1,
                              mlfScalar(
                                mclLengthInt(mclVv(c1lindices, "c1lindices"))));
                            mlfAssign(
                              &nsc2,
                              mlfScalar(
                                mclLengthInt(mclVv(c2lindices, "c2lindices"))));
                            mlfAssign(
                              &nspikescells,
                              mlfHorzcat(
                                mclVv(nsc1, "nsc1"),
                                mclVv(nsc2, "nsc2"),
                                NULL));
                            mlfAssign(
                              &ec1lindices,
                              mlfColon(
                                mclArrayRef2(
                                  mclVv(startlimits, "startlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(EXT1COL, "EXT1COL")),
                                mclArrayRef2(
                                  mclVv(endlimits, "endlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(EXT1COL, "EXT1COL")),
                                NULL));
                            mlfAssign(
                              &ec2lindices,
                              mlfColon(
                                mclArrayRef2(
                                  mclVv(startlimits, "startlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(EXT2COL, "EXT2COL")),
                                mclArrayRef2(
                                  mclVv(endlimits, "endlimits"),
                                  mclVv(repn, "repn"),
                                  mclVv(EXT2COL, "EXT2COL")),
                                NULL));
                            mlfAssign(
                              &espt1,
                              mclArrayRef1(
                                mclVv(trialspikebin, "trialspikebin"),
                                mclVv(ec1lindices, "ec1lindices")));
                            mlfAssign(
                              &espt2,
                              mclArrayRef1(
                                mclVv(trialspikebin, "trialspikebin"),
                                mclVv(ec2lindices, "ec2lindices")));
                            mlfIndexAssign(
                              &esptimes,
                              "{?}",
                              _mxarray1_,
                              mclVv(espt2, "espt2"));
                            mlfIndexAssign(
                              &esptimes,
                              "{?}",
                              _mxarray0_,
                              mclVv(espt1, "espt1"));
                            mlfAssign(
                              &ensc1,
                              mlfScalar(mclLengthInt(mclVv(espt1, "espt1"))));
                            mlfAssign(
                              &ensc2,
                              mlfScalar(mclLengthInt(mclVv(espt2, "espt2"))));
                            mlfAssign(
                              &enspikescells,
                              mlfHorzcat(
                                mclVv(ensc1, "ensc1"),
                                mclVv(ensc2, "ensc2"),
                                NULL));
                            mlfAssign(
                              &excspt1,
                              mclFeval(
                                mclValueVarargout(),
                                mlxSetdiff,
                                mlfIndexRef(
                                  mclVv(esptimes, "esptimes"),
                                  "{?}",
                                  _mxarray0_),
                                mlfIndexRef(
                                  mclVv(sptimes, "sptimes"), "{?}", _mxarray0_),
                                NULL));
                            mlfAssign(
                              &excspt2,
                              mclFeval(
                                mclValueVarargout(),
                                mlxSetdiff,
                                mlfIndexRef(
                                  mclVv(esptimes, "esptimes"),
                                  "{?}",
                                  _mxarray1_),
                                mlfIndexRef(
                                  mclVv(sptimes, "sptimes"), "{?}", _mxarray1_),
                                NULL));
                            mlfIndexAssign(
                              &excsptimes,
                              "{?}",
                              _mxarray1_,
                              mclVv(excspt2, "excspt2"));
                            mlfIndexAssign(
                              &excsptimes,
                              "{?}",
                              _mxarray0_,
                              mclVv(excspt1, "excspt1"));
                            mlfAssign(
                              &excnsc1,
                              mlfScalar(
                                mclLengthInt(mclVv(excspt1, "excspt1"))));
                            mlfAssign(
                              &excnsc2,
                              mlfScalar(
                                mclLengthInt(mclVv(excspt2, "excspt2"))));
                            mlfAssign(
                              &excnspikescells,
                              mlfHorzcat(
                                mclVv(excnsc1, "excnsc1"),
                                mclVv(excnsc2, "excnsc2"),
                                NULL));
                            mlfAssign(
                              &mat1,
                              mlfHorzcat(
                                mlfIndexRef(
                                  mclVv(sptimes, "sptimes"),
                                  "{?}(?)",
                                  _mxarray0_,
                                  mlfCreateColonIndex()),
                                mlfRepmat(
                                  _mxarray9_, mclVv(nsc1, "nsc1"), _mxarray0_),
                                NULL));
                            mlfAssign(
                              &mat2,
                              mlfVertcat(
                                mlfOnes(_mxarray0_, mclVv(nsc2, "nsc2"), NULL),
                                mclFeval(
                                  mclValueVarargout(),
                                  mlxCtranspose,
                                  mlfIndexRef(
                                    mclVv(sptimes, "sptimes"),
                                    "{?}(?)",
                                    _mxarray1_,
                                    mlfCreateColonIndex()),
                                  NULL),
                                NULL));
                            mlfAssign(
                              &xcvals,
                              mclMtimes(
                                mclVv(mat1, "mat1"), mclVv(mat2, "mat2")));
                            mlfAssign(
                              &nxcvals,
                              mclMtimes(
                                mclVv(nsc1, "nsc1"), mclVv(nsc2, "nsc2")));
                            mlfAssign(
                              &emat2,
                              mlfVertcat(
                                mlfOnes(
                                  _mxarray0_, mclVv(excnsc2, "excnsc2"), NULL),
                                mclFeval(
                                  mclValueVarargout(),
                                  mlxCtranspose,
                                  mlfIndexRef(
                                    mclVv(excsptimes, "excsptimes"),
                                    "{?}(?)",
                                    _mxarray1_,
                                    mlfCreateColonIndex()),
                                  NULL),
                                NULL));
                            mlfAssign(
                              &excvals1,
                              mclMtimes(
                                mclVv(mat1, "mat1"), mclVv(emat2, "emat2")));
                            mlfAssign(
                              &nexcvals1,
                              mclMtimes(
                                mclVv(nsc1, "nsc1"),
                                mclVv(excnsc2, "excnsc2")));
                            mlfAssign(
                              &emat1,
                              mlfHorzcat(
                                mlfIndexRef(
                                  mclVv(excsptimes, "excsptimes"),
                                  "{?}(?)",
                                  _mxarray0_,
                                  mlfCreateColonIndex()),
                                mlfRepmat(
                                  _mxarray9_,
                                  mclVv(excnsc1, "excnsc1"),
                                  _mxarray0_),
                                NULL));
                            mlfAssign(
                              &excvals2,
                              mclMtimes(
                                mclVv(emat1, "emat1"), mclVv(mat2, "mat2")));
                            mlfAssign(
                              &nexcvals2,
                              mclMtimes(
                                mclVv(excnsc1, "excnsc1"),
                                mclVv(nsc2, "nsc2")));
                            mlfAssign(
                              &xc1,
                              mlfRepmat(
                                _mxarray10_,
                                mclPlus(
                                  mclPlus(
                                    mclVv(nxcvals, "nxcvals"),
                                    mclVv(nexcvals1, "nexcvals1")),
                                  mclVv(nexcvals2, "nexcvals2")),
                                mclVv(numSurrogates2, "numSurrogates2")));
                            mclArrayAssign2(
                              &xc1,
                              mlfVertcat(
                                mclArrayRef1(
                                  mclVv(xcvals, "xcvals"),
                                  mlfCreateColonIndex()),
                                mclArrayRef1(
                                  mclVv(excvals1, "excvals1"),
                                  mlfCreateColonIndex()),
                                mclArrayRef1(
                                  mclVv(excvals2, "excvals2"),
                                  mlfCreateColonIndex()),
                                NULL),
                              mlfCreateColonIndex(),
                              mclVv(DATACOL, "DATACOL"));
                            mlfAssign(&cellvecr, mclVv(cellvec, "cellvec"));
                            if (mclEqBool(mclVv(nsc1, "nsc1"), _mxarray4_)) {
                                mlfAssign(&cellvecr, _mxarray1_);
                                mlfIndexAssign(
                                  &spiketrains, "{?}", _mxarray0_, _mxarray6_);
                                mlfAssign(&noskip1, _mxarray4_);
                            } else {
                                mlfIndexAssign(
                                  &sptrains,
                                  "{?}{?}",
                                  _mxarray0_,
                                  mlfScalar(v_),
                                  mlfZeros(
                                    mclVv(nsc1, "nsc1"),
                                    mclVv(numSurrogates1, "numSurrogates1"),
                                    NULL));
                                mlfIndexAssign(
                                  &sptrains,
                                  "{?}{?}(?,?)",
                                  _mxarray0_,
                                  mlfScalar(v_),
                                  mlfCreateColonIndex(),
                                  mclVv(DATACOL, "DATACOL"),
                                  mlfIndexRef(
                                    mclVv(sptimes, "sptimes"),
                                    "{?}(?)",
                                    _mxarray0_,
                                    mlfCreateColonIndex()));
                                mlfAssign(&noskip1, _mxarray0_);
                            }
                            if (mclEqBool(mclVv(nsc2, "nsc2"), _mxarray4_)) {
                                mlfAssign(&cellvecr, _mxarray0_);
                                mlfIndexAssign(
                                  &spiketrains, "{?}", _mxarray1_, _mxarray6_);
                                mlfAssign(&noskip2, _mxarray4_);
                            } else {
                                mlfIndexAssign(
                                  &sptrains,
                                  "{?}{?}",
                                  _mxarray1_,
                                  mlfScalar(v_),
                                  mlfZeros(
                                    mclVv(nsc2, "nsc2"),
                                    mclVv(numSurrogates1, "numSurrogates1"),
                                    NULL));
                                mlfIndexAssign(
                                  &sptrains,
                                  "{?}{?}(?,?)",
                                  _mxarray1_,
                                  mlfScalar(v_),
                                  mlfCreateColonIndex(),
                                  mclVv(DATACOL, "DATACOL"),
                                  mlfIndexRef(
                                    mclVv(sptimes, "sptimes"),
                                    "{?}(?)",
                                    _mxarray1_,
                                    mlfCreateColonIndex()));
                                mlfAssign(&noskip2, _mxarray0_);
                            }
                            mlfAssign(
                              &nfc, mclVv(needrefcheck, "needrefcheck"));
                            {
                                mclForLoopIterator viter__0;
                                for (mclForStart(
                                       &viter__0,
                                       mclVv(cellvecr, "cellvecr"),
                                       NULL,
                                       NULL);
                                     mclForNext(&viter__0, &celli);
                                     ) {
                                    mlfAssign(
                                      &bspikes,
                                      mclGt(
                                        mclArrayRef1(
                                          mclVv(enspikescells, "enspikescells"),
                                          mclVv(celli, "celli")),
                                        _mxarray0_));
                                    mlfAssign(
                                      &bisi,
                                      mclNot(
                                        mlfIsempty(
                                          mlfFind(
                                            NULL,
                                            NULL,
                                            mclLt(
                                              mclFeval(
                                                mclValueVarargout(),
                                                mlxDiff,
                                                mlfIndexRef(
                                                  mclVv(esptimes, "esptimes"),
                                                  "{?}",
                                                  mclVv(celli, "celli")),
                                                NULL),
                                              mclVv(
                                                shufflevecl,
                                                "shufflevecl"))))));
                                    mlfAssign(
                                      &nref,
                                      mlfScalar(
                                        mclLengthInt(
                                          mlfFind(
                                            NULL,
                                            NULL,
                                            mclLt(
                                              mclFeval(
                                                mclValueVarargout(),
                                                mlxDiff,
                                                mlfIndexRef(
                                                  mclVv(
                                                    excsptimes, "excsptimes"),
                                                  "{?}",
                                                  mclVv(celli, "celli")),
                                                NULL),
                                              _mxarray0_)))));
                                    if (mclScalarToBool(
                                          mclVv(bspikes, "bspikes"))
                                        && mclScalarToBool(
                                             mclVv(bisi, "bisi"))) {
                                        mclArrayAssign1(
                                          &nfc,
                                          _mxarray0_,
                                          mclVv(celli, "celli"));
                                    }
                                }
                                mclDestroyForLoopIterator(viter__0);
                            }
                            {
                                mclForLoopIterator viter__1;
                                for (mclForStart(
                                       &viter__1,
                                       mclVv(surrvec, "surrvec"),
                                       NULL,
                                       NULL);
                                     mclForNext(&viter__1, &surri);
                                     ) {
                                    mlfAssign(
                                      &surri1,
                                      mclArrayRef1(
                                        mclVv(surrvec1, "surrvec1"),
                                        mclVv(surri, "surri")));
                                    mlfAssign(
                                      &surri2,
                                      mclArrayRef1(
                                        mclVv(surrvec2, "surrvec2"),
                                        mclVv(surri, "surri")));
                                    {
                                        mclForLoopIterator viter__2;
                                        for (mclForStart(
                                               &viter__2,
                                               mclVv(cellvecr, "cellvecr"),
                                               NULL,
                                               NULL);
                                             mclForNext(&viter__2, &celli);
                                             ) {
                                            if (mlfTobool(
                                                  mclArrayRef1(
                                                    mclVv(nfc, "nfc"),
                                                    mclVv(celli, "celli")))) {
                                                mlfAssign(&nisi0, _mxarray0_);
                                                while (mlfTobool(
                                                         mclVv(
                                                           nisi0, "nisi0"))) {
                                                    mlfAssign(
                                                      &randnums,
                                                      mclMinus(
                                                        mlfFloor(
                                                          mclMtimes(
                                                            mlfNRand(
                                                              1,
                                                              mclArrayRef1(
                                                                mclVv(
                                                                  nspikescells,
                                                                  "nspikes"
                                                                  "cells"),
                                                                mclVv(
                                                                  celli,
                                                                  "celli")),
                                                              _mxarray0_,
                                                              NULL),
                                                            mclVv(
                                                              shufflevecl,
                                                              "shufflevecl"))),
                                                        mclVv(
                                                          shuffle, "shuffle")));
                                                    mlfAssign(
                                                      &spt1,
                                                      mlfSort(
                                                        NULL,
                                                        mclFeval(
                                                          mclValueVarargout(),
                                                          mlxPlus,
                                                          mlfIndexRef(
                                                            mclVv(
                                                              sptimes,
                                                              "sptimes"),
                                                            "{?}",
                                                            mclVv(
                                                              celli, "celli")),
                                                          mclVv(
                                                            randnums,
                                                            "randnums"),
                                                          NULL),
                                                        NULL));
                                                    mlfAssign(
                                                      &nisi0,
                                                      mclMinus(
                                                        mlfScalar(
                                                          mclLengthInt(
                                                            mlfFind(
                                                              NULL,
                                                              NULL,
                                                              mclEq(
                                                                mlfDiff(
                                                                  mclVv(
                                                                    spt1,
                                                                    "spt1"),
                                                                  NULL,
                                                                  NULL),
                                                                _mxarray4_)))),
                                                        mclVv(nref, "nref")));
                                                }
                                                mlfIndexAssign(
                                                  &spiketrains,
                                                  "{?}",
                                                  mclVv(celli, "celli"),
                                                  mclVv(spt1, "spt1"));
                                            } else {
                                                mlfIndexAssign(
                                                  &spiketrains,
                                                  "{?}",
                                                  mclVv(celli, "celli"),
                                                  mlfSort(
                                                    NULL,
                                                    mclMinus(
                                                      mclFeval(
                                                        mclValueVarargout(),
                                                        mlxPlus,
                                                        mlfIndexRef(
                                                          mclVv(
                                                            sptimes, "sptimes"),
                                                          "{?}",
                                                          mclVv(
                                                            celli, "celli")),
                                                        mlfFloor(
                                                          mclMtimes(
                                                            mlfNRand(
                                                              1,
                                                              mclArrayRef1(
                                                                mclVv(
                                                                  nspikescells,
                                                                  "nspikes"
                                                                  "cells"),
                                                                mclVv(
                                                                  celli,
                                                                  "celli")),
                                                              _mxarray0_,
                                                              NULL),
                                                            mclVv(
                                                              shufflevecl,
                                                              "shufflevecl"))),
                                                        NULL),
                                                      mclVv(
                                                        shuffle, "shuffle")),
                                                    NULL));
                                            }
                                        }
                                        mclDestroyForLoopIterator(viter__2);
                                    }
                                    if (mlfTobool(mclVv(noskip1, "noskip1"))) {
                                        mclArrayAssign2(
                                          &mat1,
                                          mlfIndexRef(
                                            mclVv(spiketrains, "spiketrains"),
                                            "{?}(?)",
                                            _mxarray0_,
                                            mlfCreateColonIndex()),
                                          mlfCreateColonIndex(),
                                          _mxarray0_);
                                    }
                                    if (mlfTobool(mclVv(noskip2, "noskip2"))) {
                                        mclArrayAssign2(
                                          &mat2,
                                          mclFeval(
                                            mclValueVarargout(),
                                            mlxCtranspose,
                                            mlfIndexRef(
                                              mclVv(spiketrains, "spiketrains"),
                                              "{?}(?)",
                                              _mxarray1_,
                                              mlfCreateColonIndex()),
                                            NULL),
                                          _mxarray1_,
                                          mlfCreateColonIndex());
                                    }
                                    mlfAssign(
                                      &xcvals,
                                      mclMtimes(
                                        mclVv(mat1, "mat1"),
                                        mclVv(mat2, "mat2")));
                                    mlfAssign(
                                      &excvals1,
                                      mclMtimes(
                                        mclVv(mat1, "mat1"),
                                        mclVv(emat2, "emat2")));
                                    mlfAssign(
                                      &excvals2,
                                      mclMtimes(
                                        mclVv(emat1, "emat1"),
                                        mclVv(mat2, "mat2")));
                                    mclArrayAssign2(
                                      &xc1,
                                      mlfVertcat(
                                        mclArrayRef1(
                                          mclVv(xcvals, "xcvals"),
                                          mlfCreateColonIndex()),
                                        mclArrayRef1(
                                          mclVv(excvals1, "excvals1"),
                                          mlfCreateColonIndex()),
                                        mclArrayRef1(
                                          mclVv(excvals2, "excvals2"),
                                          mlfCreateColonIndex()),
                                        NULL),
                                      mlfCreateColonIndex(),
                                      mclVv(surri2, "surri2"));
                                    if (mlfTobool(mclVv(noskip1, "noskip1"))) {
                                        mlfIndexAssign(
                                          &sptrains,
                                          "{?}{?}(?,?)",
                                          _mxarray0_,
                                          mlfScalar(v_),
                                          mlfCreateColonIndex(),
                                          mclVv(surri1, "surri1"),
                                          mlfIndexRef(
                                            mclVv(spiketrains, "spiketrains"),
                                            "{?}(?)",
                                            _mxarray0_,
                                            mlfCreateColonIndex()));
                                    }
                                    if (mlfTobool(mclVv(noskip2, "noskip2"))) {
                                        mlfIndexAssign(
                                          &sptrains,
                                          "{?}{?}(?,?)",
                                          _mxarray1_,
                                          mlfScalar(v_),
                                          mlfCreateColonIndex(),
                                          mclVv(surri1, "surri1"),
                                          mlfIndexRef(
                                            mclVv(spiketrains, "spiketrains"),
                                            "{?}(?)",
                                            _mxarray1_,
                                            mlfCreateColonIndex()));
                                    }
                                }
                                mclDestroyForLoopIterator(viter__1);
                            }
                            mlfIndexAssign(
                              &xcreps, "{?}", mlfScalar(v_), mclVv(xc1, "xc1"));
                            if (v_ == e_) {
                                break;
                            }
                            ++v_;
                        }
                        mlfAssign(&repi, mlfScalar(v_));
                    }
                }
                /*
                 * % convert cell array xcreps to matrix
                 * xcmat = cell2mat(xcreps);
                 */
                mlfAssign(&xcmat, mlfCell2mat(mclVv(xcreps, "xcreps")));
                /*
                 * % check if xcmat is a single row (i.e. only 1 rep had exactly 1
                 * % spike in each cell)
                 * if(size(xcmat,1)==1)
                 */
                if (mclEqBool(
                      mlfSize(
                        mclValueVarargout(), mclVv(xcmat, "xcmat"), _mxarray0_),
                      _mxarray0_)) {
                    /*
                     * % add a second row of nan so the histogram will be taken in
                     * % columns. Otherwise, hist will collapse the values across the
                     * % entire set of surrogates and return xc which is 1 x lagslength
                     * % instead of lagslength x (numsurrogates+1)
                     * xcmat = concatenate(xcmat,nan);
                     */
                    mlfAssign(
                      &xcmat,
                      mlfConcatenate(mclVv(xcmat, "xcmat"), _mxarray10_, NULL));
                /*
                 * end
                 */
                }
                /*
                 * % do histogram of the values in xcreps
                 * xc = histcie(xcmat,lagbins2,'DropLast');        
                 */
                mlfAssign(
                  &xc,
                  mlfNHistcie(
                    1,
                    NULL,
                    mclVv(xcmat, "xcmat"),
                    mclVv(lagbins2, "lagbins2"),
                    _mxarray11_,
                    NULL));
                /*
                 * 
                 * % now compute the shift predictor
                 * % subtract rep number for the second cell so that the 1st rep of
                 * % cell 1 is compared to 2nd rep of cell 2
                 * sperepnums2 = mod(erepnums2 - 2,reps) + 1;
                 */
                mlfAssign(
                  &sperepnums2,
                  mclPlus(
                    mlfMod(
                      mclMinus(mclVv(erepnums2, "erepnums2"), _mxarray1_),
                      mclVv(reps, "reps")),
                    _mxarray0_));
                /*
                 * % subtract the rep number for the central window for the second cell
                 * sprepnums2 = mod(repnums2 - 2,reps) + 1;
                 */
                mlfAssign(
                  &sprepnums2,
                  mclPlus(
                    mlfMod(
                      mclMinus(mclVv(repnums2, "repnums2"), _mxarray1_),
                      mclVv(reps, "reps")),
                    _mxarray0_));
                /*
                 * % find the reps that have spikes in the central window of at least 1 cell
                 * % can't use union since it will include reps with spikes in the central
                 * % window for 1 cell and no spikes in the extended window for the other
                 * [spintrepnums1,spsria1,spsrib2] = intersect(erepnums1,sprepnums2);
                 */
                mlfAssign(
                  &spintrepnums1,
                  mlfNIntersect(
                    3,
                    &spsria1,
                    &spsrib2,
                    mclVv(erepnums1, "erepnums1"),
                    mclVv(sprepnums2, "sprepnums2"),
                    NULL));
                /*
                 * [spintrepnums2,spsrib1,spsria2] = intersect(sperepnums2,repnums1);
                 */
                mlfAssign(
                  &spintrepnums2,
                  mlfNIntersect(
                    3,
                    &spsrib1,
                    &spsria2,
                    mclVv(sperepnums2, "sperepnums2"),
                    mclVv(repnums1, "repnums1"),
                    NULL));
                /*
                 * sprepsboth = unique([spintrepnums1; spintrepnums2]);
                 */
                mlfAssign(
                  &sprepsboth,
                  mlfNUnique(
                    1,
                    NULL,
                    NULL,
                    mlfVertcat(
                      mclVv(spintrepnums1, "spintrepnums1"),
                      mclVv(spintrepnums2, "spintrepnums2"),
                      NULL),
                    NULL));
                /*
                 * % save rep numbers so we can compute synchrony only on these reps
                 * sprepnvec = vecr(sprepsboth);
                 */
                mlfAssign(&sprepnvec, mlfVecr(mclVv(sprepsboth, "sprepsboth")));
                /*
                 * sprepnl = length(sprepnvec);
                 */
                mlfAssign(
                  &sprepnl,
                  mlfScalar(mclLengthInt(mclVv(sprepnvec, "sprepnvec"))));
                /*
                 * 
                 * if(sprepnl>0)
                 */
                if (mclGtBool(mclVv(sprepnl, "sprepnl"), _mxarray4_)) {
                    /*
                     * xcreps = cell(sprepnl,1);
                     */
                    mlfAssign(
                      &xcreps,
                      mlfCell(mclVv(sprepnl, "sprepnl"), _mxarray0_, NULL));
                    /*
                     * % set the limits for the indices in trialspikebin
                     * startlimits(:,[EXT2COL CEN2COL]) = circshift(startlimits(:,[EXT2COL CEN2COL]),-1);
                     */
                    mclArrayAssign2(
                      &startlimits,
                      mlfCircshift(
                        mclArrayRef2(
                          mclVv(startlimits, "startlimits"),
                          mlfCreateColonIndex(),
                          mlfHorzcat(
                            mclVv(EXT2COL, "EXT2COL"),
                            mclVv(CEN2COL, "CEN2COL"),
                            NULL)),
                        _mxarray9_),
                      mlfCreateColonIndex(),
                      mlfHorzcat(
                        mclVv(EXT2COL, "EXT2COL"),
                        mclVv(CEN2COL, "CEN2COL"),
                        NULL));
                    /*
                     * endlimits(:,[EXT2COL CEN2COL]) = circshift(endlimits(:,[EXT2COL CEN2COL]),-1);
                     */
                    mclArrayAssign2(
                      &endlimits,
                      mlfCircshift(
                        mclArrayRef2(
                          mclVv(endlimits, "endlimits"),
                          mlfCreateColonIndex(),
                          mlfHorzcat(
                            mclVv(EXT2COL, "EXT2COL"),
                            mclVv(CEN2COL, "CEN2COL"),
                            NULL)),
                        _mxarray9_),
                      mlfCreateColonIndex(),
                      mlfHorzcat(
                        mclVv(EXT2COL, "EXT2COL"),
                        mclVv(CEN2COL, "CEN2COL"),
                        NULL));
                    /*
                     * 
                     * for repi = 1:sprepnl
                     */
                    {
                        int v_ = mclForIntStart(1);
                        int e_ = mclForIntEnd(mclVv(sprepnl, "sprepnl"));
                        if (v_ > e_) {
                            mlfAssign(&repi, _mxarray6_);
                        } else {
                            /*
                             * % get rep number
                             * repn = sprepsboth(repi);
                             * % get spikes times for the central window
                             * c1lindices = startlimits(repn,CEN1COL):endlimits(repn,CEN1COL);
                             * c2lindices = startlimits(repn,CEN2COL):endlimits(repn,CEN2COL);
                             * sptimes1 = ctrialspikebin(c1lindices);
                             * sptimes2 = ctrialspikebin(c2lindices);
                             * nsc1 = length(c1lindices);
                             * nsc2 = length(c2lindices);
                             * 
                             * % get the spike times for the extended window
                             * ec1lindices = startlimits(repn,EXT1COL):endlimits(repn,EXT1COL);
                             * ec2lindices = startlimits(repn,EXT2COL):endlimits(repn,EXT2COL);
                             * espt1 = trialspikebin(ec1lindices);
                             * espt2 = trialspikebin(ec2lindices);
                             * 
                             * % get the spike times for the extended window after removing the
                             * % spikes in the central window. This makes computing the xcorr 
                             * % easier since we don't have to subtract out the double-counting 
                             * % of the central spikes
                             * excsptimes1 = setdiff(espt1,sptimes1);
                             * excsptimes2 = setdiff(espt2,sptimes2);
                             * excnsc1 = length(excsptimes1);
                             * excnsc2 = length(excsptimes2);
                             * 
                             * % check spike times
                             * % 				concat([full([trialspikebin(ec1lindices) trialspiken(ec1lindices) repspike(ec1lindices)])], ...
                             * % 						[full([trialspikebin(ec2lindices) trialspiken(ec2lindices) repspike(ec2lindices)])], ...
                             * % 						[full([ctrialspikebin(c1lindices) crepspike(c1lindices)])], ...
                             * % 						[full([ctrialspikebin(c2lindices) crepspike(c2lindices)])],'Columnwise')
                             * 
                             * % create matrix to compute cross-correlation between central 
                             * % windows of cells 1 and 2
                             * mat1 = [sptimes1(:) repmat(-1,nsc1,1)];
                             * mat2 = [ones(1,nsc2); sptimes2(:)'];
                             * % compute xcorr for the data
                             * xcvals = mat1 * mat2;
                             * 
                             * % compute xcorr between central window of cell 1 and extended
                             * % window of cell 2
                             * emat2 = [ones(1,excnsc2); excsptimes2(:)'];
                             * excvals1 = mat1 * emat2;
                             * 
                             * % compute xcorr between central window of cell 2 and extended
                             * % window of cell 1
                             * emat1 = [excsptimes1(:) repmat(-1,excnsc1,1)];
                             * excvals2 = emat1 * mat2;
                             * xcreps{repi} = [xcvals(:); excvals1(:); excvals2(:)];
                             * end % for repi = 1:repnl
                             */
                            for (; ; ) {
                                mlfAssign(
                                  &repn,
                                  mclIntArrayRef1(
                                    mclVv(sprepsboth, "sprepsboth"), v_));
                                mlfAssign(
                                  &c1lindices,
                                  mlfColon(
                                    mclArrayRef2(
                                      mclVv(startlimits, "startlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(CEN1COL, "CEN1COL")),
                                    mclArrayRef2(
                                      mclVv(endlimits, "endlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(CEN1COL, "CEN1COL")),
                                    NULL));
                                mlfAssign(
                                  &c2lindices,
                                  mlfColon(
                                    mclArrayRef2(
                                      mclVv(startlimits, "startlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(CEN2COL, "CEN2COL")),
                                    mclArrayRef2(
                                      mclVv(endlimits, "endlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(CEN2COL, "CEN2COL")),
                                    NULL));
                                mlfAssign(
                                  &sptimes1,
                                  mclArrayRef1(
                                    mclVv(ctrialspikebin, "ctrialspikebin"),
                                    mclVv(c1lindices, "c1lindices")));
                                mlfAssign(
                                  &sptimes2,
                                  mclArrayRef1(
                                    mclVv(ctrialspikebin, "ctrialspikebin"),
                                    mclVv(c2lindices, "c2lindices")));
                                mlfAssign(
                                  &nsc1,
                                  mlfScalar(
                                    mclLengthInt(
                                      mclVv(c1lindices, "c1lindices"))));
                                mlfAssign(
                                  &nsc2,
                                  mlfScalar(
                                    mclLengthInt(
                                      mclVv(c2lindices, "c2lindices"))));
                                mlfAssign(
                                  &ec1lindices,
                                  mlfColon(
                                    mclArrayRef2(
                                      mclVv(startlimits, "startlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(EXT1COL, "EXT1COL")),
                                    mclArrayRef2(
                                      mclVv(endlimits, "endlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(EXT1COL, "EXT1COL")),
                                    NULL));
                                mlfAssign(
                                  &ec2lindices,
                                  mlfColon(
                                    mclArrayRef2(
                                      mclVv(startlimits, "startlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(EXT2COL, "EXT2COL")),
                                    mclArrayRef2(
                                      mclVv(endlimits, "endlimits"),
                                      mclVv(repn, "repn"),
                                      mclVv(EXT2COL, "EXT2COL")),
                                    NULL));
                                mlfAssign(
                                  &espt1,
                                  mclArrayRef1(
                                    mclVv(trialspikebin, "trialspikebin"),
                                    mclVv(ec1lindices, "ec1lindices")));
                                mlfAssign(
                                  &espt2,
                                  mclArrayRef1(
                                    mclVv(trialspikebin, "trialspikebin"),
                                    mclVv(ec2lindices, "ec2lindices")));
                                mlfAssign(
                                  &excsptimes1,
                                  mlfNSetdiff(
                                    1,
                                    NULL,
                                    mclVv(espt1, "espt1"),
                                    mclVv(sptimes1, "sptimes1"),
                                    NULL));
                                mlfAssign(
                                  &excsptimes2,
                                  mlfNSetdiff(
                                    1,
                                    NULL,
                                    mclVv(espt2, "espt2"),
                                    mclVv(sptimes2, "sptimes2"),
                                    NULL));
                                mlfAssign(
                                  &excnsc1,
                                  mlfScalar(
                                    mclLengthInt(
                                      mclVv(excsptimes1, "excsptimes1"))));
                                mlfAssign(
                                  &excnsc2,
                                  mlfScalar(
                                    mclLengthInt(
                                      mclVv(excsptimes2, "excsptimes2"))));
                                mlfAssign(
                                  &mat1,
                                  mlfHorzcat(
                                    mclArrayRef1(
                                      mclVv(sptimes1, "sptimes1"),
                                      mlfCreateColonIndex()),
                                    mlfRepmat(
                                      _mxarray9_,
                                      mclVv(nsc1, "nsc1"),
                                      _mxarray0_),
                                    NULL));
                                mlfAssign(
                                  &mat2,
                                  mlfVertcat(
                                    mlfOnes(
                                      _mxarray0_, mclVv(nsc2, "nsc2"), NULL),
                                    mlfCtranspose(
                                      mclArrayRef1(
                                        mclVv(sptimes2, "sptimes2"),
                                        mlfCreateColonIndex())),
                                    NULL));
                                mlfAssign(
                                  &xcvals,
                                  mclMtimes(
                                    mclVv(mat1, "mat1"), mclVv(mat2, "mat2")));
                                mlfAssign(
                                  &emat2,
                                  mlfVertcat(
                                    mlfOnes(
                                      _mxarray0_,
                                      mclVv(excnsc2, "excnsc2"),
                                      NULL),
                                    mlfCtranspose(
                                      mclArrayRef1(
                                        mclVv(excsptimes2, "excsptimes2"),
                                        mlfCreateColonIndex())),
                                    NULL));
                                mlfAssign(
                                  &excvals1,
                                  mclMtimes(
                                    mclVv(mat1, "mat1"),
                                    mclVv(emat2, "emat2")));
                                mlfAssign(
                                  &emat1,
                                  mlfHorzcat(
                                    mclArrayRef1(
                                      mclVv(excsptimes1, "excsptimes1"),
                                      mlfCreateColonIndex()),
                                    mlfRepmat(
                                      _mxarray9_,
                                      mclVv(excnsc1, "excnsc1"),
                                      _mxarray0_),
                                    NULL));
                                mlfAssign(
                                  &excvals2,
                                  mclMtimes(
                                    mclVv(emat1, "emat1"),
                                    mclVv(mat2, "mat2")));
                                mlfIndexAssign(
                                  &xcreps,
                                  "{?}",
                                  mlfScalar(v_),
                                  mlfVertcat(
                                    mclArrayRef1(
                                      mclVv(xcvals, "xcvals"),
                                      mlfCreateColonIndex()),
                                    mclArrayRef1(
                                      mclVv(excvals1, "excvals1"),
                                      mlfCreateColonIndex()),
                                    mclArrayRef1(
                                      mclVv(excvals2, "excvals2"),
                                      mlfCreateColonIndex()),
                                    NULL));
                                if (v_ == e_) {
                                    break;
                                }
                                ++v_;
                            }
                            mlfAssign(&repi, mlfScalar(v_));
                        }
                    }
                    /*
                     * % do histogram of the values in xcreps
                     * xc(:,SHIFTCOL) = vecc(histcie(full(cell2mat(xcreps)),lagbins2,'DropLast'));
                     */
                    mclArrayAssign2(
                      &xc,
                      mlfVecc(
                        mlfNHistcie(
                          1,
                          NULL,
                          mlfFull(mlfCell2mat(mclVv(xcreps, "xcreps"))),
                          mclVv(lagbins2, "lagbins2"),
                          _mxarray11_,
                          NULL)),
                      mlfCreateColonIndex(),
                      mclVv(SHIFTCOL, "SHIFTCOL"));
                /*
                 * end % if(~isempty(trialspiken) && (repnl>0))
                 */
                }
                /*
                 * 
                 * % now shuffle the spike times for the remaining reps that are
                 * % not in repnvec
                 * % create cell array to store cell1limits and cell2limits and its
                 * % corresponding c1l1 and c2l1
                 * clcell = {cell1limits,cell2limits};
                 */
                mlfAssign(
                  &clcell,
                  mlfCellhcat(
                    mclVv(cell1limits, "cell1limits"),
                    mclVv(cell2limits, "cell2limits"),
                    NULL));
                /*
                 * cl1cell = {c1l1,c2l1};
                 */
                mlfAssign(
                  &cl1cell,
                  mlfCellhcat(mclVv(c1l1, "c1l1"), mclVv(c2l1, "c2l1"), NULL));
                /*
                 * eclcell = {ecell1limits,ecell2limits};
                 */
                mlfAssign(
                  &eclcell,
                  mlfCellhcat(
                    mclVv(ecell1limits, "ecell1limits"),
                    mclVv(ecell2limits, "ecell2limits"),
                    NULL));
                /*
                 * ecl1cell = {ec1l1,ec2l1};
                 */
                mlfAssign(
                  &ecl1cell,
                  mlfCellhcat(
                    mclVv(ec1l1, "ec1l1"), mclVv(ec2l1, "ec2l1"), NULL));
                /*
                 * repnumscell = {repnums1,repnums2};
                 */
                mlfAssign(
                  &repnumscell,
                  mlfCellhcat(
                    mclVv(repnums1, "repnums1"),
                    mclVv(repnums2, "repnums2"),
                    NULL));
                /*
                 * erepnumscell = {erepnums1,erepnums2};
                 */
                mlfAssign(
                  &erepnumscell,
                  mlfCellhcat(
                    mclVv(erepnums1, "erepnums1"),
                    mclVv(erepnums2, "erepnums2"),
                    NULL));
                /*
                 * for celli = cellvec
                 */
                {
                    mclForLoopIterator viter__3;
                    for (mclForStart(
                           &viter__3, mclVv(cellvec, "cellvec"), NULL, NULL);
                         mclForNext(&viter__3, &celli);
                         ) {
                        mlfAssign(&_T0_, mclVv(celli, "celli"));
                        /*
                         * % determine the repetitions that are not covered by repnvec,
                         * % i.e. reps that do not have spikes in both cells
                         * [repsnopair{celli},srinopair] = setdiff(repnumscell{celli},repsboth);
                         */
                        mclFeval(
                          mlfIndexVarargout(
                            &repsnopair, "{?}", _T0_, &srinopair, "", NULL),
                          mlxSetdiff,
                          mlfIndexRef(
                            mclVv(repnumscell, "repnumscell"),
                            "{?}",
                            mclVv(celli, "celli")),
                          mclVv(repsboth, "repsboth"),
                          NULL);
                        /*
                         * % sri1nopair = srinopair + 1;
                         * % find the indices for the extended windows
                         * [discard,esrinopair] = intersect(erepnumscell{celli},repsnopair{celli});
                         */
                        mclFeval(
                          mlfVarargout(&discard, &esrinopair, NULL),
                          mlxIntersect,
                          mlfIndexRef(
                            mclVv(erepnumscell, "erepnumscell"),
                            "{?}",
                            mclVv(celli, "celli")),
                          mlfIndexRef(
                            mclVv(repsnopair, "repsnopair"),
                            "{?}",
                            mclVv(celli, "celli")),
                          NULL);
                        /*
                         * % esri1nopair = esrinopair + 1;
                         * a = repsnopair{celli};
                         */
                        mlfAssign(
                          &a,
                          mlfIndexRef(
                            mclVv(repsnopair, "repsnopair"),
                            "{?}",
                            mclVv(celli, "celli")));
                        /*
                         * for repi = 1:length(repsnopair{celli})
                         */
                        {
                            int v_ = mclForIntStart(1);
                            int e_
                              = mclForIntEnd(
                                  mclFeval(
                                    mclValueVarargout(),
                                    mlxLength,
                                    mlfIndexRef(
                                      mclVv(repsnopair, "repsnopair"),
                                      "{?}",
                                      mclVv(celli, "celli")),
                                    NULL));
                            if (v_ > e_) {
                                mlfAssign(&repi, _mxarray6_);
                            } else {
                                /*
                                 * repi2 = repnl+repi;
                                 * % get the spike times for the central window
                                 * srepi = srinopair(repi);
                                 * clindices = cl1cell{celli}(srepi) : clcell{celli}(srepi);
                                 * sptime = ctrialspikebin(clindices);
                                 * sptrains{celli}{repi2}(:,DATACOL) = sptime;
                                 * nspikes = length(sptime);
                                 * % get the spike times for the extended window
                                 * esrepi = esrinopair(repi);
                                 * eclindices = ecl1cell{celli}(esrepi) : eclcell{celli}(esrepi);
                                 * esptime = trialspikebin(eclindices);
                                 * enspikes = length(esptime);
                                 * % get the spike times for the extended window without the
                                 * % spikes from the central window
                                 * excsptime = setdiff(esptime,sptime);
                                 * 
                                 * % check the spike times
                                 * % concat(full([sptime crepspike(clindices)]),full([esptime repspike(eclindices)]),'Columnwise')
                                 * 
                                 * % only do while loop if there is more than 1 spike
                                 * % and the intervals between the spikes is such that 
                                 * % there will be overlaps
                                 * bspikes = enspikes > 1;
                                 * bisi = ~isempty(find(diff(esptime)<shufflevecl));
                                 * % check if there are isi's smaller than 1 in the extended
                                 * % window
                                 * nref = length(find(diff(excsptime)<1));
                                 * if( bspikes && bisi)
                                 * nfc = 1;
                                 * else
                                 * nfc = 0;
                                 * end
                                 * 
                                 * % loop over number of surrogates
                                 * for surri = surrvec
                                 * surri1 = surrvec1(surri);
                                 * % check if we need to check for refractory period 
                                 * % violations
                                 * if(nfc)
                                 * % make sure we go into the while loop
                                 * nisi0 = 1;
                                 * % check for refractory period violations
                                 * while(nisi0)
                                 * % get random numbers for cell i
                                 * randnums = floor(rand(nspikes,1) ...
                                 * * shufflevecl) - shuffle;
                                 * % add random numbers to spike times for cell i and
                                 * % sort the spike times
                                 * spt1 = sort(sptime+randnums);
                                 * % take the diff and find 0's
                                 * % no need to worry about negative numbers since
                                 * % the first number is always going to be larger
                                 * % or equal to nref
                                 * nisi0 = length(find(diff(spt1)==0)) - nref;
                                 * end % while(nisi0)
                                 * sptrains{celli}{repi2}(:,surri1) = spt1;
                                 * else % if(nfc)
                                 * % don't have to worry about overlap between the 
                                 * % shuffled spikes so just add the shuffle to the
                                 * % spike times
                                 * sptrains{celli}{repi2}(:,surri1) = sort(sptime + ...
                                 * floor(rand(nspikes,1) ...
                                 * * shufflevecl) - shuffle);
                                 * end % if(nfc)
                                 * end % for surri = surrvec
                                 * end % for repi = 1:length(repsnopair{celli})
                                 */
                                for (; ; ) {
                                    mlfAssign(
                                      &repi2,
                                      mclPlus(
                                        mclVv(repnl, "repnl"), mlfScalar(v_)));
                                    mlfAssign(
                                      &srepi,
                                      mclIntArrayRef1(
                                        mclVv(srinopair, "srinopair"), v_));
                                    mlfAssign(
                                      &clindices,
                                      mclFeval(
                                        mclValueVarargout(),
                                        mlxColon,
                                        mlfIndexRef(
                                          mclVv(cl1cell, "cl1cell"),
                                          "{?}(?)",
                                          mclVv(celli, "celli"),
                                          mclVv(srepi, "srepi")),
                                        mlfIndexRef(
                                          mclVv(clcell, "clcell"),
                                          "{?}(?)",
                                          mclVv(celli, "celli"),
                                          mclVv(srepi, "srepi")),
                                        NULL));
                                    mlfAssign(
                                      &sptime,
                                      mclArrayRef1(
                                        mclVv(ctrialspikebin, "ctrialspikebin"),
                                        mclVv(clindices, "clindices")));
                                    mlfIndexAssign(
                                      &sptrains,
                                      "{?}{?}(?,?)",
                                      mclVv(celli, "celli"),
                                      mclVv(repi2, "repi2"),
                                      mlfCreateColonIndex(),
                                      mclVv(DATACOL, "DATACOL"),
                                      mclVv(sptime, "sptime"));
                                    mlfAssign(
                                      &nspikes,
                                      mlfScalar(
                                        mclLengthInt(mclVv(sptime, "sptime"))));
                                    mlfAssign(
                                      &esrepi,
                                      mclIntArrayRef1(
                                        mclVv(esrinopair, "esrinopair"), v_));
                                    mlfAssign(
                                      &eclindices,
                                      mclFeval(
                                        mclValueVarargout(),
                                        mlxColon,
                                        mlfIndexRef(
                                          mclVv(ecl1cell, "ecl1cell"),
                                          "{?}(?)",
                                          mclVv(celli, "celli"),
                                          mclVv(esrepi, "esrepi")),
                                        mlfIndexRef(
                                          mclVv(eclcell, "eclcell"),
                                          "{?}(?)",
                                          mclVv(celli, "celli"),
                                          mclVv(esrepi, "esrepi")),
                                        NULL));
                                    mlfAssign(
                                      &esptime,
                                      mclArrayRef1(
                                        mclVv(trialspikebin, "trialspikebin"),
                                        mclVv(eclindices, "eclindices")));
                                    mlfAssign(
                                      &enspikes,
                                      mlfScalar(
                                        mclLengthInt(
                                          mclVv(esptime, "esptime"))));
                                    mlfAssign(
                                      &excsptime,
                                      mlfNSetdiff(
                                        1,
                                        NULL,
                                        mclVv(esptime, "esptime"),
                                        mclVv(sptime, "sptime"),
                                        NULL));
                                    mlfAssign(
                                      &bspikes,
                                      mclGt(
                                        mclVv(enspikes, "enspikes"),
                                        _mxarray0_));
                                    mlfAssign(
                                      &bisi,
                                      mclNot(
                                        mlfIsempty(
                                          mlfFind(
                                            NULL,
                                            NULL,
                                            mclLt(
                                              mlfDiff(
                                                mclVv(esptime, "esptime"),
                                                NULL,
                                                NULL),
                                              mclVv(
                                                shufflevecl,
                                                "shufflevecl"))))));
                                    mlfAssign(
                                      &nref,
                                      mlfScalar(
                                        mclLengthInt(
                                          mlfFind(
                                            NULL,
                                            NULL,
                                            mclLt(
                                              mlfDiff(
                                                mclVv(excsptime, "excsptime"),
                                                NULL,
                                                NULL),
                                              _mxarray0_)))));
                                    if (mclScalarToBool(
                                          mclVv(bspikes, "bspikes"))
                                        && mclScalarToBool(
                                             mclVv(bisi, "bisi"))) {
                                        mlfAssign(&nfc, _mxarray0_);
                                    } else {
                                        mlfAssign(&nfc, _mxarray4_);
                                    }
                                    {
                                        mclForLoopIterator viter__4;
                                        for (mclForStart(
                                               &viter__4,
                                               mclVv(surrvec, "surrvec"),
                                               NULL,
                                               NULL);
                                             mclForNext(&viter__4, &surri);
                                             ) {
                                            mlfAssign(
                                              &surri1,
                                              mclArrayRef1(
                                                mclVv(surrvec1, "surrvec1"),
                                                mclVv(surri, "surri")));
                                            if (mlfTobool(mclVv(nfc, "nfc"))) {
                                                mlfAssign(&nisi0, _mxarray0_);
                                                while (mlfTobool(
                                                         mclVv(
                                                           nisi0, "nisi0"))) {
                                                    mlfAssign(
                                                      &randnums,
                                                      mclMinus(
                                                        mlfFloor(
                                                          mclMtimes(
                                                            mlfNRand(
                                                              1,
                                                              mclVv(
                                                                nspikes,
                                                                "nspikes"),
                                                              _mxarray0_,
                                                              NULL),
                                                            mclVv(
                                                              shufflevecl,
                                                              "shufflevecl"))),
                                                        mclVv(
                                                          shuffle, "shuffle")));
                                                    mlfAssign(
                                                      &spt1,
                                                      mlfSort(
                                                        NULL,
                                                        mclPlus(
                                                          mclVv(
                                                            sptime, "sptime"),
                                                          mclVv(
                                                            randnums,
                                                            "randnums")),
                                                        NULL));
                                                    mlfAssign(
                                                      &nisi0,
                                                      mclMinus(
                                                        mlfScalar(
                                                          mclLengthInt(
                                                            mlfFind(
                                                              NULL,
                                                              NULL,
                                                              mclEq(
                                                                mlfDiff(
                                                                  mclVv(
                                                                    spt1,
                                                                    "spt1"),
                                                                  NULL,
                                                                  NULL),
                                                                _mxarray4_)))),
                                                        mclVv(nref, "nref")));
                                                }
                                                mlfIndexAssign(
                                                  &sptrains,
                                                  "{?}{?}(?,?)",
                                                  mclVv(celli, "celli"),
                                                  mclVv(repi2, "repi2"),
                                                  mlfCreateColonIndex(),
                                                  mclVv(surri1, "surri1"),
                                                  mclVv(spt1, "spt1"));
                                            } else {
                                                mlfIndexAssign(
                                                  &sptrains,
                                                  "{?}{?}(?,?)",
                                                  mclVv(celli, "celli"),
                                                  mclVv(repi2, "repi2"),
                                                  mlfCreateColonIndex(),
                                                  mclVv(surri1, "surri1"),
                                                  mlfSort(
                                                    NULL,
                                                    mclMinus(
                                                      mclPlus(
                                                        mclVv(sptime, "sptime"),
                                                        mlfFloor(
                                                          mclMtimes(
                                                            mlfNRand(
                                                              1,
                                                              mclVv(
                                                                nspikes,
                                                                "nspikes"),
                                                              _mxarray0_,
                                                              NULL),
                                                            mclVv(
                                                              shufflevecl,
                                                              "shufflevecl")))),
                                                      mclVv(
                                                        shuffle, "shuffle")),
                                                    NULL));
                                            }
                                        }
                                        mclDestroyForLoopIterator(viter__4);
                                    }
                                    if (v_ == e_) {
                                        break;
                                    }
                                    ++v_;
                                }
                                mlfAssign(&repi, mlfScalar(v_));
                            }
                        }
                    /*
                     * end % for celli = [1 2]
                     */
                    }
                    mclDestroyForLoopIterator(viter__3);
                }
                /*
                 * 
                 * repswspike1 = [repsboth; repsnopair{1}(:)];
                 */
                mlfAssign(
                  &repswspike1,
                  mlfVertcat(
                    mclVv(repsboth, "repsboth"),
                    mlfIndexRef(
                      mclVv(repsnopair, "repsnopair"),
                      "{?}(?)",
                      _mxarray0_,
                      mlfCreateColonIndex()),
                    NULL));
                /*
                 * repswspike2 = [repsboth; repsnopair{2}(:)];
                 */
                mlfAssign(
                  &repswspike2,
                  mlfVertcat(
                    mclVv(repsboth, "repsboth"),
                    mlfIndexRef(
                      mclVv(repsnopair, "repsnopair"),
                      "{?}(?)",
                      _mxarray1_,
                      mlfCreateColonIndex()),
                    NULL));
            /*
             * end % if(~isempty(trialspiken) && ~isempty(repnvec))
             */
            }
            /*
             * % convert xc to sparse matrix to save space
             * xc = sparse(xc);
             */
            mlfAssign(
              &xc, mlfSparse(mclVv(xc, "xc"), NULL, NULL, NULL, NULL, NULL));
            /*
             * % save data to file
             * save([surrprefix num2str(windown,'%04d')],'xc','sptrains', ...
             * 'repswspike1','repswspike2');
             */
            mlfSave(
              mlfHorzcat(
                mclVv(surrprefix, "surrprefix"),
                mlfNum2str(mclVv(windown, "windown"), _mxarray13_),
                NULL),
              "w",
              "xc",
              xc,
              "sptrains",
              sptrains,
              "repswspike1",
              repswspike1,
              "repswspike2",
              repswspike2,
              NULL);
        /*
         * end % for windown = startwindow:endwindow
         */
        }
        mclDestroyForLoopIterator(viter__);
    }
    mxDestroyArray(DATACOL);
    mxDestroyArray(SHIFTCOL);
    mxDestroyArray(EXT1COL);
    mxDestroyArray(EXT2COL);
    mxDestroyArray(CEN1COL);
    mxDestroyArray(CEN2COL);
    mxDestroyArray(ans);
    mxDestroyArray(cellpairdatafile);
    mxDestroyArray(startwindow);
    mxDestroyArray(endwindow);
    mxDestroyArray(numwindows);
    mxDestroyArray(windown);
    mxDestroyArray(startbins2);
    mxDestroyArray(fmin);
    mxDestroyArray(endbins2);
    mxDestroyArray(fmax);
    mxDestroyArray(binidx);
    mxDestroyArray(trialspiken);
    mxDestroyArray(repspike);
    mxDestroyArray(emptytsn);
    mxDestroyArray(binidxsize);
    mxDestroyArray(trialspikebin);
    mxDestroyArray(eurepspike);
    mxDestroyArray(eurepspikea);
    mxDestroyArray(eurepspikeb);
    mxDestroyArray(reps);
    mxDestroyArray(ecell2i);
    mxDestroyArray(ecell1i);
    mxDestroyArray(erepnums1);
    mxDestroyArray(erepnums2);
    mxDestroyArray(startbins);
    mxDestroyArray(fmin2);
    mxDestroyArray(endbins);
    mxDestroyArray(fmax2);
    mxDestroyArray(ctsbi);
    mxDestroyArray(ctrialspikebin);
    mxDestroyArray(crepspike);
    mxDestroyArray(urepspike);
    mxDestroyArray(urepspikea);
    mxDestroyArray(urepspikeb);
    mxDestroyArray(cell2i);
    mxDestroyArray(cell1i);
    mxDestroyArray(repnums1);
    mxDestroyArray(repnums2);
    mxDestroyArray(intrepnums1);
    mxDestroyArray(sria1);
    mxDestroyArray(srib2);
    mxDestroyArray(intrepnums2);
    mxDestroyArray(srib1);
    mxDestroyArray(sria2);
    mxDestroyArray(repsboth);
    mxDestroyArray(repnvec);
    mxDestroyArray(repnl);
    mxDestroyArray(xctemp);
    mxDestroyArray(xc);
    mxDestroyArray(sptrains);
    mxDestroyArray(repswspike1);
    mxDestroyArray(repswspike2);
    mxDestroyArray(xcreps);
    mxDestroyArray(lrepn1);
    mxDestroyArray(lrepn2);
    mxDestroyArray(spiketrains);
    mxDestroyArray(elrepn1);
    mxDestroyArray(elrepn2);
    mxDestroyArray(esptrains);
    mxDestroyArray(ecell1limits);
    mxDestroyArray(ecell2limits);
    mxDestroyArray(ec1l1);
    mxDestroyArray(ec2l1);
    mxDestroyArray(cell1limits);
    mxDestroyArray(cell2limits);
    mxDestroyArray(blr1);
    mxDestroyArray(c1l1);
    mxDestroyArray(c2l1);
    mxDestroyArray(tmpstartlimits);
    mxDestroyArray(startlimits);
    mxDestroyArray(tmpendlimits);
    mxDestroyArray(endlimits);
    mxDestroyArray(repi);
    mxDestroyArray(repn);
    mxDestroyArray(c1lindices);
    mxDestroyArray(c2lindices);
    mxDestroyArray(sptimes);
    mxDestroyArray(nsc1);
    mxDestroyArray(nsc2);
    mxDestroyArray(nspikescells);
    mxDestroyArray(ec1lindices);
    mxDestroyArray(ec2lindices);
    mxDestroyArray(espt1);
    mxDestroyArray(espt2);
    mxDestroyArray(esptimes);
    mxDestroyArray(ensc1);
    mxDestroyArray(ensc2);
    mxDestroyArray(enspikescells);
    mxDestroyArray(excspt1);
    mxDestroyArray(excspt2);
    mxDestroyArray(excsptimes);
    mxDestroyArray(excnsc1);
    mxDestroyArray(excnsc2);
    mxDestroyArray(excnspikescells);
    mxDestroyArray(mat1);
    mxDestroyArray(mat2);
    mxDestroyArray(xcvals);
    mxDestroyArray(nxcvals);
    mxDestroyArray(emat2);
    mxDestroyArray(excvals1);
    mxDestroyArray(nexcvals1);
    mxDestroyArray(emat1);
    mxDestroyArray(excvals2);
    mxDestroyArray(nexcvals2);
    mxDestroyArray(numSurrogates2);
    mxDestroyArray(xc1);
    mxDestroyArray(cellvec);
    mxDestroyArray(cellvecr);
    mxDestroyArray(noskip1);
    mxDestroyArray(numSurrogates1);
    mxDestroyArray(noskip2);
    mxDestroyArray(needrefcheck);
    mxDestroyArray(nfc);
    mxDestroyArray(celli);
    mxDestroyArray(bspikes);
    mxDestroyArray(shufflevecl);
    mxDestroyArray(bisi);
    mxDestroyArray(nref);
    mxDestroyArray(surri);
    mxDestroyArray(surrvec);
    mxDestroyArray(surrvec1);
    mxDestroyArray(surri1);
    mxDestroyArray(surrvec2);
    mxDestroyArray(surri2);
    mxDestroyArray(nisi0);
    mxDestroyArray(shuffle);
    mxDestroyArray(randnums);
    mxDestroyArray(spt1);
    mxDestroyArray(xcmat);
    mxDestroyArray(lagbins2);
    mxDestroyArray(sperepnums2);
    mxDestroyArray(sprepnums2);
    mxDestroyArray(spintrepnums1);
    mxDestroyArray(spsria1);
    mxDestroyArray(spsrib2);
    mxDestroyArray(spintrepnums2);
    mxDestroyArray(spsrib1);
    mxDestroyArray(spsria2);
    mxDestroyArray(sprepsboth);
    mxDestroyArray(sprepnvec);
    mxDestroyArray(sprepnl);
    mxDestroyArray(sptimes1);
    mxDestroyArray(sptimes2);
    mxDestroyArray(excsptimes1);
    mxDestroyArray(excsptimes2);
    mxDestroyArray(clcell);
    mxDestroyArray(cl1cell);
    mxDestroyArray(eclcell);
    mxDestroyArray(ecl1cell);
    mxDestroyArray(repnumscell);
    mxDestroyArray(erepnumscell);
    mxDestroyArray(repsnopair);
    mxDestroyArray(_T0_);
    mxDestroyArray(srinopair);
    mxDestroyArray(discard);
    mxDestroyArray(esrinopair);
    mxDestroyArray(a);
    mxDestroyArray(repi2);
    mxDestroyArray(srepi);
    mxDestroyArray(clindices);
    mxDestroyArray(sptime);
    mxDestroyArray(nspikes);
    mxDestroyArray(esrepi);
    mxDestroyArray(eclindices);
    mxDestroyArray(esptime);
    mxDestroyArray(enspikes);
    mxDestroyArray(excsptime);
    mxDestroyArray(surrprefix);
    mxDestroyArray(varargin);
    mxDestroyArray(sdatafile);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
}
