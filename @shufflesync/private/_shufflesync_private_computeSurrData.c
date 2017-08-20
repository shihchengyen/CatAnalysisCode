/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 31 15:20:38 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "computeSurrData" 
 */
#include "_shufflesync_private_computeSurrData.h"
#include "concat.h"
#include "concatenate.h"
#include "convmtx.h"
#include "histcie.h"
#include "jbtest.h"
#include "kstest2.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "lillietest.h"
#include "nptDir.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;
static mxArray * _mxarray2_;
static mxArray * _mxarray3_;

static mxChar _array5_[5] = { '*', '.', 'm', 'a', 't' };
static mxArray * _mxarray4_;

static mxChar _array7_[37] = { 'I', 'n', 'c', 'o', 'm', 'p', 'l', 'e', 't', 'e',
                               ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f',
                               ' ', 's', 'u', 'r', 'r', 'o', 'g', 'a', 't', 'e',
                               ' ', 'f', 'i', 'l', 'e', 's', '!' };
static mxArray * _mxarray6_;
static double _ieee_nan_;
static mxArray * _mxarray8_;
static mxArray * _mxarray9_;

static double _array11_[2] = { 1.0, 2.0 };
static mxArray * _mxarray10_;

static mxChar _array13_[3] = { 'o', 'f', 'f' };
static mxArray * _mxarray12_;

static mxChar _array15_[19] = { 'M', 'A', 'T', 'L', 'A', 'B', ':',
                                'd', 'i', 'v', 'i', 'd', 'e', 'B',
                                'y', 'Z', 'e', 'r', 'o' };
static mxArray * _mxarray14_;
static mxArray * _mxarray16_;

static mxChar _array18_[4] = { '%', '0', '4', 'd' };
static mxArray * _mxarray17_;

static mxChar _array20_[10] = { 'C', 'o', 'l', 'u', 'm',
                                'n', 'w', 'i', 's', 'e' };
static mxArray * _mxarray19_;
static mxArray * _mxarray21_;

static mxChar _array23_[8] = { 'D', 'r', 'o', 'p', 'L', 'a', 's', 't' };
static mxArray * _mxarray22_;

static mxChar _array25_[8] = { 'D', 'a', 't', 'a', 'C', 'o', 'l', 's' };
static mxArray * _mxarray24_;

static double _array27_[2] = { 1.0, -1.0 };
static mxArray * _mxarray26_;
static mxArray * _mxarray28_;
static mxArray * _mxarray29_;

static mxChar _array31_[2] = { 'o', 'n' };
static mxArray * _mxarray30_;

void InitializeModule__shufflesync_private_computeSurrData(void) {
    _mxarray0_ = mclInitializeDouble(1.0);
    _mxarray1_ = mclInitializeDouble(2.0);
    _mxarray2_ = mclInitializeDouble(3.0);
    _mxarray3_ = mclInitializeDouble(4.0);
    _mxarray4_ = mclInitializeString(5, _array5_);
    _mxarray6_ = mclInitializeString(37, _array7_);
    _ieee_nan_ = mclGetNaN();
    _mxarray8_ = mclInitializeDouble(_ieee_nan_);
    _mxarray9_ = mclInitializeDouble(6.0);
    _mxarray10_ = mclInitializeDoubleVector(1, 2, _array11_);
    _mxarray12_ = mclInitializeString(3, _array13_);
    _mxarray14_ = mclInitializeString(19, _array15_);
    _mxarray16_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray17_ = mclInitializeString(4, _array18_);
    _mxarray19_ = mclInitializeString(10, _array20_);
    _mxarray21_ = mclInitializeDouble(.5);
    _mxarray22_ = mclInitializeString(8, _array23_);
    _mxarray24_ = mclInitializeString(8, _array25_);
    _mxarray26_ = mclInitializeDoubleVector(1, 2, _array27_);
    _mxarray28_ = mclInitializeDouble(-1.0);
    _mxarray29_ = mclInitializeDouble(0.0);
    _mxarray30_ = mclInitializeString(2, _array31_);
}

void TerminateModule__shufflesync_private_computeSurrData(void) {
    mxDestroyArray(_mxarray30_);
    mxDestroyArray(_mxarray29_);
    mxDestroyArray(_mxarray28_);
    mxDestroyArray(_mxarray26_);
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray22_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * M_shufflesync_private_computeSurrData(mxArray * * surrmean,
                                                       mxArray * * surrstd,
                                                       mxArray * * surrhist,
                                                       mxArray * * zscore,
                                                       int nargout_,
                                                       mxArray * sdatafile);

_mexLocalFunctionTable _local_function_table__shufflesync_private_computeSurrData
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlf_shufflesync_private_computeSurrData" contains the normal
 * interface for the "@shufflesync/private/computeSurrData" M-function from
 * file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/@shufflesync/private/computeSurrD
 * ata.m" (lines 1-350). This function processes any input arguments and passes
 * them to the implementation version of the function, appearing above.
 */
mxArray * mlf_shufflesync_private_computeSurrData(mxArray * * surrmean,
                                                  mxArray * * surrstd,
                                                  mxArray * * surrhist,
                                                  mxArray * * zscore,
                                                  mxArray * sdatafile) {
    int nargout = 1;
    mxArray * data = NULL;
    mxArray * surrmean__ = NULL;
    mxArray * surrstd__ = NULL;
    mxArray * surrhist__ = NULL;
    mxArray * zscore__ = NULL;
    mlfEnterNewContext(4, 1, surrmean, surrstd, surrhist, zscore, sdatafile);
    if (surrmean != NULL) {
        ++nargout;
    }
    if (surrstd != NULL) {
        ++nargout;
    }
    if (surrhist != NULL) {
        ++nargout;
    }
    if (zscore != NULL) {
        ++nargout;
    }
    data
      = M_shufflesync_private_computeSurrData(
          &surrmean__, &surrstd__, &surrhist__, &zscore__, nargout, sdatafile);
    mlfRestorePreviousContext(
      4, 1, surrmean, surrstd, surrhist, zscore, sdatafile);
    if (surrmean != NULL) {
        mclCopyOutputArg(surrmean, surrmean__);
    } else {
        mxDestroyArray(surrmean__);
    }
    if (surrstd != NULL) {
        mclCopyOutputArg(surrstd, surrstd__);
    } else {
        mxDestroyArray(surrstd__);
    }
    if (surrhist != NULL) {
        mclCopyOutputArg(surrhist, surrhist__);
    } else {
        mxDestroyArray(surrhist__);
    }
    if (zscore != NULL) {
        mclCopyOutputArg(zscore, zscore__);
    } else {
        mxDestroyArray(zscore__);
    }
    return mlfReturnValue(data);
}

/*
 * The function "mlx_shufflesync_private_computeSurrData" contains the feval
 * interface for the "@shufflesync/private/computeSurrData" M-function from
 * file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/@shufflesync/private/computeSurrD
 * ata.m" (lines 1-350). The feval function calls the implementation version of
 * @shufflesync/private/computeSurrData through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
void mlx_shufflesync_private_computeSurrData(int nlhs,
                                             mxArray * plhs[],
                                             int nrhs,
                                             mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[5];
    int i;
    if (nlhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: @shufflesync/private/computeSurrData Line: "
            "1 Column: 1 The function \"@shufflesync/private/computeSurrData\""
            " was called with more than the declared number of outputs (5)."),
          NULL);
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: @shufflesync/private/computeSurrData Line: "
            "1 Column: 1 The function \"@shufflesync/private/computeSurrData\""
            " was called with more than the declared number of inputs (1)."),
          NULL);
    }
    for (i = 0; i < 5; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0]
      = M_shufflesync_private_computeSurrData(
          &mplhs[1], &mplhs[2], &mplhs[3], &mplhs[4], nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 5 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 5; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "M_shufflesync_private_computeSurrData" is the implementation
 * version of the "@shufflesync/private/computeSurrData" M-function from file
 * "/Users/syen/Documents/ShihCheng/Matlab/Cat/@shufflesync/private/computeSurrD
 * ata.m" (lines 1-350). It contains the actual compiled code for that
 * M-function. It is a static function and must only be called from one of the
 * interface functions, appearing below.
 */
/*
 * function [data,surrmean,surrstd,surrhist,zscore] = computeSurrData(sdatafile)
 */
static mxArray * M_shufflesync_private_computeSurrData(mxArray * * surrmean,
                                                       mxArray * * surrstd,
                                                       mxArray * * surrhist,
                                                       mxArray * * zscore,
                                                       int nargout_,
                                                       mxArray * sdatafile) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(
          &_local_function_table__shufflesync_private_computeSurrData);
    mxArray * data = NULL;
    mxArray * resultfile = NULL;
    mxArray * m2 = NULL;
    mxArray * m1 = NULL;
    mxArray * xc = NULL;
    mxArray * repi = NULL;
    mxArray * stdiff = NULL;
    mxArray * binidx2 = NULL;
    mxArray * binidx = NULL;
    mxArray * spmarkmat2 = NULL;
    mxArray * spmarkmat1 = NULL;
    mxArray * markmat2 = NULL;
    mxArray * markmat1 = NULL;
    mxArray * mat2 = NULL;
    mxArray * mat1 = NULL;
    mxArray * nspikes = NULL;
    mxArray * binidxsize = NULL;
    mxArray * sprepvec2 = NULL;
    mxArray * repvec2 = NULL;
    mxArray * sprepvec = NULL;
    mxArray * repvec = NULL;
    mxArray * reps = NULL;
    mxArray * dzscore = NULL;
    mxArray * si = NULL;
    mxArray * bins = NULL;
    mxArray * SurrHistStep = NULL;
    mxArray * SurrHistMin = NULL;
    mxArray * maxn = NULL;
    mxArray * surrvec1 = NULL;
    mxArray * loopcounts = NULL;
    mxArray * srtimes = NULL;
    mxArray * surrtimes = NULL;
    mxArray * sptimes = NULL;
    mxArray * sptrains = NULL;
    mxArray * cellvec = NULL;
    mxArray * cidx = NULL;
    mxArray * loopbins = NULL;
    mxArray * sspstd = NULL;
    mxArray * sspm = NULL;
    mxArray * sspsthmat = NULL;
    mxArray * dspsth2 = NULL;
    mxArray * dspsth1 = NULL;
    mxArray * spsth = NULL;
    mxArray * psth2 = NULL;
    mxArray * psth1 = NULL;
    mxArray * psth = NULL;
    mxArray * psthhistcbins = NULL;
    mxArray * spsthbins = NULL;
    mxArray * endbins = NULL;
    mxArray * startbins = NULL;
    mxArray * smat = NULL;
    mxArray * l1spt2 = NULL;
    mxArray * l1spt1 = NULL;
    mxArray * l1 = NULL;
    mxArray * l1sp2 = NULL;
    mxArray * l1sp1 = NULL;
    mxArray * l0spt2 = NULL;
    mxArray * l0spt1 = NULL;
    mxArray * l0 = NULL;
    mxArray * l0sp2 = NULL;
    mxArray * l0sp1 = NULL;
    mxArray * spt = NULL;
    mxArray * dummy = NULL;
    mxArray * testdata = NULL;
    mxArray * lagi = NULL;
    mxArray * sstd = NULL;
    mxArray * smean = NULL;
    mxArray * xcsdata = NULL;
    mxArray * xcdata = NULL;
    mxArray * nsyncspikes = NULL;
    mxArray * l = NULL;
    mxArray * sidx = NULL;
    mxArray * dspsthz = NULL;
    mxArray * defaultspt = NULL;
    mxArray * smoothmtx = NULL;
    mxArray * surridx2 = NULL;
    mxArray * surridx1 = NULL;
    mxArray * didx = NULL;
    mxArray * cwinidx2 = NULL;
    mxArray * cwinidx = NULL;
    mxArray * ps2ind2 = NULL;
    mxArray * ps2ind = NULL;
    mxArray * ps1ind2 = NULL;
    mxArray * ps1ind = NULL;
    mxArray * lpsthbins = NULL;
    mxArray * lcentralwin = NULL;
    mxArray * halfwinsize = NULL;
    mxArray * ws1 = NULL;
    mxArray * psthwin = NULL;
    mxArray * ns2 = NULL;
    mxArray * ns1 = NULL;
    mxArray * psths = NULL;
    mxArray * psthbins = NULL;
    mxArray * pbins = NULL;
    mxArray * dvec = NULL;
    mxArray * windowlength = NULL;
    mxArray * bintimes = NULL;
    mxArray * extrabins = NULL;
    mxArray * bsize = NULL;
    mxArray * bt1l = NULL;
    mxArray * bintimes1 = NULL;
    mxArray * subbins = NULL;
    mxArray * sadd = NULL;
    mxArray * shuffle = NULL;
    mxArray * limat = NULL;
    mxArray * jbmat = NULL;
    mxArray * lagsvec = NULL;
    mxArray * hpval = NULL;
    mxArray * dxcstd = NULL;
    mxArray * dxcmean = NULL;
    mxArray * dxcz = NULL;
    mxArray * dxcorr = NULL;
    mxArray * srind2 = NULL;
    mxArray * srind = NULL;
    mxArray * numsurrvec = NULL;
    mxArray * centerbins = NULL;
    mxArray * NumCentralBins = NULL;
    mxArray * llcenter = NULL;
    mxArray * smstmp = NULL;
    mxArray * sumvec = NULL;
    mxArray * lagslength = NULL;
    mxArray * sdata = NULL;
    mxArray * numSurrogates = NULL;
    mxArray * shiftp = NULL;
    mxArray * numwindows = NULL;
    mxArray * surrnum = NULL;
    mxArray * surrlist = NULL;
    mxArray * surrprefix = NULL;
    mxArray * cellpairdatafile = NULL;
    mxArray * ans = NULL;
    mxArray * KSNSPIKES = NULL;
    mxArray * SURR_COL = NULL;
    mxArray * SHIFT_COL = NULL;
    mxArray * DATA_COL = NULL;
    mclCopyArray(&sdatafile);
    /*
     * 
     * % set constants
     * DATA_COL = 1;
     */
    mlfAssign(&DATA_COL, _mxarray0_);
    /*
     * SHIFT_COL = 2;
     */
    mlfAssign(&SHIFT_COL, _mxarray1_);
    /*
     * SURR_COL = 3;
     */
    mlfAssign(&SURR_COL, _mxarray2_);
    /*
     * % the kstest2 is not acurate if n1*n2/(n1+n2) < 4. Since
     * % the number of surrogate spikes is always going to be 1000
     * % times the number of real spikes, this simplifies to
     * % 1000n^2/1001n < 4, which further simplies to 
     * % n < 4001/1000, which is pretty close to 4
     * KSNSPIKES = 4;
     */
    mlfAssign(&KSNSPIKES, _mxarray3_);
    /*
     * 
     * % load the surrogate data file
     * load(sdatafile);
     */
    {
        mxArray * name_ = mclInitialize(mclVa(sdatafile, "sdatafile"));
        mclLoadConditional(
          name_,
          "DATA_COL",
          &DATA_COL,
          "KSNSPIKES",
          &KSNSPIKES,
          "NumCentralBins",
          &NumCentralBins,
          "SHIFT_COL",
          &SHIFT_COL,
          "SURR_COL",
          &SURR_COL,
          "SurrHistMin",
          &SurrHistMin,
          "SurrHistStep",
          &SurrHistStep,
          "ans",
          &ans,
          "binidx",
          &binidx,
          "binidx2",
          &binidx2,
          "binidxsize",
          &binidxsize,
          "bins",
          &bins,
          "bintimes",
          &bintimes,
          "bintimes1",
          &bintimes1,
          "bsize",
          &bsize,
          NULL);
        {
            mxArray * name_0 = mclInitialize(name_);
            mclLoadConditional(
              name_0,
              "bt1l",
              &bt1l,
              "cellpairdatafile",
              &cellpairdatafile,
              "cellvec",
              &cellvec,
              "centerbins",
              &centerbins,
              "cidx",
              &cidx,
              "cwinidx",
              &cwinidx,
              "cwinidx2",
              &cwinidx2,
              "data",
              &data,
              "defaultspt",
              &defaultspt,
              "didx",
              &didx,
              "dspsth1",
              &dspsth1,
              "dspsth2",
              &dspsth2,
              "dspsthz",
              &dspsthz,
              "dummy",
              &dummy,
              "dvec",
              &dvec,
              NULL);
            {
                mxArray * name_1 = mclInitialize(name_0);
                mclLoadConditional(
                  name_1,
                  "dxcmean",
                  &dxcmean,
                  "dxcorr",
                  &dxcorr,
                  "dxcstd",
                  &dxcstd,
                  "dxcz",
                  &dxcz,
                  "dzscore",
                  &dzscore,
                  "endbins",
                  &endbins,
                  "extrabins",
                  &extrabins,
                  "halfwinsize",
                  &halfwinsize,
                  "hpval",
                  &hpval,
                  "jbmat",
                  &jbmat,
                  "l",
                  &l,
                  "l0",
                  &l0,
                  "l0sp1",
                  &l0sp1,
                  "l0sp2",
                  &l0sp2,
                  "l0spt1",
                  &l0spt1,
                  NULL);
                {
                    mxArray * name_2 = mclInitialize(name_1);
                    mclLoadConditional(
                      name_2,
                      "l0spt2",
                      &l0spt2,
                      "l1",
                      &l1,
                      "l1sp1",
                      &l1sp1,
                      "l1sp2",
                      &l1sp2,
                      "l1spt1",
                      &l1spt1,
                      "l1spt2",
                      &l1spt2,
                      "lagi",
                      &lagi,
                      "lagslength",
                      &lagslength,
                      "lagsvec",
                      &lagsvec,
                      "lcentralwin",
                      &lcentralwin,
                      "limat",
                      &limat,
                      "llcenter",
                      &llcenter,
                      "loopbins",
                      &loopbins,
                      "loopcounts",
                      &loopcounts,
                      "lpsthbins",
                      &lpsthbins,
                      NULL);
                    {
                        mxArray * name_3 = mclInitialize(name_2);
                        mclLoadConditional(
                          name_3,
                          "m1",
                          &m1,
                          "m2",
                          &m2,
                          "markmat1",
                          &markmat1,
                          "markmat2",
                          &markmat2,
                          "mat1",
                          &mat1,
                          "mat2",
                          &mat2,
                          "maxn",
                          &maxn,
                          "ns1",
                          &ns1,
                          "ns2",
                          &ns2,
                          "nspikes",
                          &nspikes,
                          "nsyncspikes",
                          &nsyncspikes,
                          "numSurrogates",
                          &numSurrogates,
                          "numsurrvec",
                          &numsurrvec,
                          "numwindows",
                          &numwindows,
                          "pbins",
                          &pbins,
                          NULL);
                        {
                            mxArray * name_4 = mclInitialize(name_3);
                            mclLoadConditional(
                              name_4,
                              "ps1ind",
                              &ps1ind,
                              "ps1ind2",
                              &ps1ind2,
                              "ps2ind",
                              &ps2ind,
                              "ps2ind2",
                              &ps2ind2,
                              "psth",
                              &psth,
                              "psth1",
                              &psth1,
                              "psth2",
                              &psth2,
                              "psthbins",
                              &psthbins,
                              "psthhistcbins",
                              &psthhistcbins,
                              "psths",
                              &psths,
                              "psthwin",
                              &psthwin,
                              "repi",
                              &repi,
                              "reps",
                              &reps,
                              "repvec",
                              &repvec,
                              "repvec2",
                              &repvec2,
                              NULL);
                            {
                                mxArray * name_5 = mclInitialize(name_4);
                                mclLoadConditional(
                                  name_5,
                                  "resultfile",
                                  &resultfile,
                                  "sadd",
                                  &sadd,
                                  "sdata",
                                  &sdata,
                                  "sdatafile",
                                  &sdatafile,
                                  "shiftp",
                                  &shiftp,
                                  "shuffle",
                                  &shuffle,
                                  "si",
                                  &si,
                                  "sidx",
                                  &sidx,
                                  "smat",
                                  &smat,
                                  "smean",
                                  &smean,
                                  "smoothmtx",
                                  &smoothmtx,
                                  "smstmp",
                                  &smstmp,
                                  "spmarkmat1",
                                  &spmarkmat1,
                                  "spmarkmat2",
                                  &spmarkmat2,
                                  "sprepvec",
                                  &sprepvec,
                                  NULL);
                                {
                                    mxArray * name_6 = mclInitialize(name_5);
                                    mclLoadConditional(
                                      name_6,
                                      "sprepvec2",
                                      &sprepvec2,
                                      "spsth",
                                      &spsth,
                                      "spsthbins",
                                      &spsthbins,
                                      "spt",
                                      &spt,
                                      "sptimes",
                                      &sptimes,
                                      "sptrains",
                                      &sptrains,
                                      "srind",
                                      &srind,
                                      "srind2",
                                      &srind2,
                                      "srtimes",
                                      &srtimes,
                                      "sspm",
                                      &sspm,
                                      "sspstd",
                                      &sspstd,
                                      "sspsthmat",
                                      &sspsthmat,
                                      "sstd",
                                      &sstd,
                                      "startbins",
                                      &startbins,
                                      "stdiff",
                                      &stdiff,
                                      NULL);
                                    {
                                        mxArray * name_7
                                          = mclInitialize(name_6);
                                        mclLoadConditional(
                                          name_7,
                                          "subbins",
                                          &subbins,
                                          "sumvec",
                                          &sumvec,
                                          "surrhist",
                                          surrhist,
                                          "surridx1",
                                          &surridx1,
                                          "surridx2",
                                          &surridx2,
                                          "surrlist",
                                          &surrlist,
                                          "surrmean",
                                          surrmean,
                                          "surrnum",
                                          &surrnum,
                                          "surrprefix",
                                          &surrprefix,
                                          "surrstd",
                                          surrstd,
                                          "surrtimes",
                                          &surrtimes,
                                          "surrvec1",
                                          &surrvec1,
                                          "testdata",
                                          &testdata,
                                          "windowlength",
                                          &windowlength,
                                          "ws1",
                                          &ws1,
                                          NULL);
                                        mclLoadConditional(
                                          name_7,
                                          "xc",
                                          &xc,
                                          "xcdata",
                                          &xcdata,
                                          "xcsdata",
                                          &xcsdata,
                                          "zscore",
                                          zscore,
                                          NULL);
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
     * % load the data file
     * load(cellpairdatafile);
     */
    {
        mxArray * name_
          = mclInitialize(mclVv(cellpairdatafile, "cellpairdatafile"));
        mclLoadConditional(
          name_,
          "DATA_COL",
          &DATA_COL,
          "KSNSPIKES",
          &KSNSPIKES,
          "NumCentralBins",
          &NumCentralBins,
          "SHIFT_COL",
          &SHIFT_COL,
          "SURR_COL",
          &SURR_COL,
          "SurrHistMin",
          &SurrHistMin,
          "SurrHistStep",
          &SurrHistStep,
          "ans",
          &ans,
          "binidx",
          &binidx,
          "binidx2",
          &binidx2,
          "binidxsize",
          &binidxsize,
          "bins",
          &bins,
          "bintimes",
          &bintimes,
          "bintimes1",
          &bintimes1,
          "bsize",
          &bsize,
          NULL);
        {
            mxArray * name_8 = mclInitialize(name_);
            mclLoadConditional(
              name_8,
              "bt1l",
              &bt1l,
              "cellpairdatafile",
              &cellpairdatafile,
              "cellvec",
              &cellvec,
              "centerbins",
              &centerbins,
              "cidx",
              &cidx,
              "cwinidx",
              &cwinidx,
              "cwinidx2",
              &cwinidx2,
              "data",
              &data,
              "defaultspt",
              &defaultspt,
              "didx",
              &didx,
              "dspsth1",
              &dspsth1,
              "dspsth2",
              &dspsth2,
              "dspsthz",
              &dspsthz,
              "dummy",
              &dummy,
              "dvec",
              &dvec,
              NULL);
            {
                mxArray * name_9 = mclInitialize(name_8);
                mclLoadConditional(
                  name_9,
                  "dxcmean",
                  &dxcmean,
                  "dxcorr",
                  &dxcorr,
                  "dxcstd",
                  &dxcstd,
                  "dxcz",
                  &dxcz,
                  "dzscore",
                  &dzscore,
                  "endbins",
                  &endbins,
                  "extrabins",
                  &extrabins,
                  "halfwinsize",
                  &halfwinsize,
                  "hpval",
                  &hpval,
                  "jbmat",
                  &jbmat,
                  "l",
                  &l,
                  "l0",
                  &l0,
                  "l0sp1",
                  &l0sp1,
                  "l0sp2",
                  &l0sp2,
                  "l0spt1",
                  &l0spt1,
                  NULL);
                {
                    mxArray * name_00 = mclInitialize(name_9);
                    mclLoadConditional(
                      name_00,
                      "l0spt2",
                      &l0spt2,
                      "l1",
                      &l1,
                      "l1sp1",
                      &l1sp1,
                      "l1sp2",
                      &l1sp2,
                      "l1spt1",
                      &l1spt1,
                      "l1spt2",
                      &l1spt2,
                      "lagi",
                      &lagi,
                      "lagslength",
                      &lagslength,
                      "lagsvec",
                      &lagsvec,
                      "lcentralwin",
                      &lcentralwin,
                      "limat",
                      &limat,
                      "llcenter",
                      &llcenter,
                      "loopbins",
                      &loopbins,
                      "loopcounts",
                      &loopcounts,
                      "lpsthbins",
                      &lpsthbins,
                      NULL);
                    {
                        mxArray * name_01 = mclInitialize(name_00);
                        mclLoadConditional(
                          name_01,
                          "m1",
                          &m1,
                          "m2",
                          &m2,
                          "markmat1",
                          &markmat1,
                          "markmat2",
                          &markmat2,
                          "mat1",
                          &mat1,
                          "mat2",
                          &mat2,
                          "maxn",
                          &maxn,
                          "ns1",
                          &ns1,
                          "ns2",
                          &ns2,
                          "nspikes",
                          &nspikes,
                          "nsyncspikes",
                          &nsyncspikes,
                          "numSurrogates",
                          &numSurrogates,
                          "numsurrvec",
                          &numsurrvec,
                          "numwindows",
                          &numwindows,
                          "pbins",
                          &pbins,
                          NULL);
                        {
                            mxArray * name_02 = mclInitialize(name_01);
                            mclLoadConditional(
                              name_02,
                              "ps1ind",
                              &ps1ind,
                              "ps1ind2",
                              &ps1ind2,
                              "ps2ind",
                              &ps2ind,
                              "ps2ind2",
                              &ps2ind2,
                              "psth",
                              &psth,
                              "psth1",
                              &psth1,
                              "psth2",
                              &psth2,
                              "psthbins",
                              &psthbins,
                              "psthhistcbins",
                              &psthhistcbins,
                              "psths",
                              &psths,
                              "psthwin",
                              &psthwin,
                              "repi",
                              &repi,
                              "reps",
                              &reps,
                              "repvec",
                              &repvec,
                              "repvec2",
                              &repvec2,
                              NULL);
                            {
                                mxArray * name_03 = mclInitialize(name_02);
                                mclLoadConditional(
                                  name_03,
                                  "resultfile",
                                  &resultfile,
                                  "sadd",
                                  &sadd,
                                  "sdata",
                                  &sdata,
                                  "sdatafile",
                                  &sdatafile,
                                  "shiftp",
                                  &shiftp,
                                  "shuffle",
                                  &shuffle,
                                  "si",
                                  &si,
                                  "sidx",
                                  &sidx,
                                  "smat",
                                  &smat,
                                  "smean",
                                  &smean,
                                  "smoothmtx",
                                  &smoothmtx,
                                  "smstmp",
                                  &smstmp,
                                  "spmarkmat1",
                                  &spmarkmat1,
                                  "spmarkmat2",
                                  &spmarkmat2,
                                  "sprepvec",
                                  &sprepvec,
                                  NULL);
                                {
                                    mxArray * name_04 = mclInitialize(name_03);
                                    mclLoadConditional(
                                      name_04,
                                      "sprepvec2",
                                      &sprepvec2,
                                      "spsth",
                                      &spsth,
                                      "spsthbins",
                                      &spsthbins,
                                      "spt",
                                      &spt,
                                      "sptimes",
                                      &sptimes,
                                      "sptrains",
                                      &sptrains,
                                      "srind",
                                      &srind,
                                      "srind2",
                                      &srind2,
                                      "srtimes",
                                      &srtimes,
                                      "sspm",
                                      &sspm,
                                      "sspstd",
                                      &sspstd,
                                      "sspsthmat",
                                      &sspsthmat,
                                      "sstd",
                                      &sstd,
                                      "startbins",
                                      &startbins,
                                      "stdiff",
                                      &stdiff,
                                      NULL);
                                    {
                                        mxArray * name_05
                                          = mclInitialize(name_04);
                                        mclLoadConditional(
                                          name_05,
                                          "subbins",
                                          &subbins,
                                          "sumvec",
                                          &sumvec,
                                          "surrhist",
                                          surrhist,
                                          "surridx1",
                                          &surridx1,
                                          "surridx2",
                                          &surridx2,
                                          "surrlist",
                                          &surrlist,
                                          "surrmean",
                                          surrmean,
                                          "surrnum",
                                          &surrnum,
                                          "surrprefix",
                                          &surrprefix,
                                          "surrstd",
                                          surrstd,
                                          "surrtimes",
                                          &surrtimes,
                                          "surrvec1",
                                          &surrvec1,
                                          "testdata",
                                          &testdata,
                                          "windowlength",
                                          &windowlength,
                                          "ws1",
                                          &ws1,
                                          NULL);
                                        mclLoadConditional(
                                          name_05,
                                          "xc",
                                          &xc,
                                          "xcdata",
                                          &xcdata,
                                          "xcsdata",
                                          &xcsdata,
                                          "zscore",
                                          zscore,
                                          NULL);
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
                mxDestroyArray(name_9);
            }
            mxDestroyArray(name_8);
        }
        mxDestroyArray(name_);
    }
    /*
     * 
     * % get list of surrogate files
     * surrlist = nptDir([surrprefix '*.mat']);
     */
    mlfAssign(
      &surrlist,
      mlfNptDir(
        mlfHorzcat(mclVv(surrprefix, "surrprefix"), _mxarray4_, NULL), NULL));
    /*
     * % check to make sure we have the complete set of files
     * surrnum = length(surrlist);
     */
    mlfAssign(&surrnum, mlfScalar(mclLengthInt(mclVv(surrlist, "surrlist"))));
    /*
     * if(surrnum<numwindows)
     */
    if (mclLtBool(mclVv(surrnum, "surrnum"), mclVv(numwindows, "numwindows"))) {
        /*
         * error('Incomplete number of surrogate files!');
         */
        mlfError(_mxarray6_, NULL);
    /*
     * else
     */
    } else {
        /*
         * % allocate memory for output arguments
         * data = zeros(numwindows,1);
         */
        mlfAssign(
          &data, mlfZeros(mclVv(numwindows, "numwindows"), _mxarray0_, NULL));
        /*
         * shiftp = data;
         */
        mlfAssign(&shiftp, mclVv(data, "data"));
        /*
         * sdata = zeros(numSurrogates,numwindows);
         */
        mlfAssign(
          &sdata,
          mlfZeros(
            mclVv(numSurrogates, "numSurrogates"),
            mclVv(numwindows, "numwindows"),
            NULL));
        /*
         * % create vector for summing appropriate bins in xcorr vector
         * % create vector using windowlength
         * sumvec = zeros(1,lagslength);
         */
        mlfAssign(
          &sumvec, mlfZeros(_mxarray0_, mclVv(lagslength, "lagslength"), NULL));
        /*
         * % allocate memory for storing mean and std of surrogate xcorr
         * smstmp = sumvec';
         */
        mlfAssign(&smstmp, mlfCtranspose(mclVv(sumvec, "sumvec")));
        /*
         * % find the center of lagslength
         * llcenter = ceil(lagslength/2);
         */
        mlfAssign(
          &llcenter,
          mlfCeil(mclMrdivide(mclVv(lagslength, "lagslength"), _mxarray1_)));
        /*
         * centerbins = (NumCentralBins-1)/2;
         */
        mlfAssign(
          &centerbins,
          mclMrdivide(
            mclMinus(mclVv(NumCentralBins, "NumCentralBins"), _mxarray0_),
            _mxarray1_));
        /*
         * % create ones in places that correspond to NumCentralBins
         * sumvec((llcenter-centerbins):(llcenter+centerbins)) = 1;
         */
        mclArrayAssign1(
          &sumvec,
          _mxarray0_,
          mlfColon(
            mclMinus(
              mclVv(llcenter, "llcenter"), mclVv(centerbins, "centerbins")),
            mclPlus(
              mclVv(llcenter, "llcenter"), mclVv(centerbins, "centerbins")),
            NULL));
        /*
         * % generate indices for the surrogates so we don't have to keep
         * % generating it inside the for loop
         * numsurrvec = 1:numSurrogates;
         */
        mlfAssign(
          &numsurrvec,
          mlfColon(_mxarray0_, mclVv(numSurrogates, "numSurrogates"), NULL));
        /*
         * srind = numsurrvec + 2;
         */
        mlfAssign(&srind, mclPlus(mclVv(numsurrvec, "numsurrvec"), _mxarray1_));
        /*
         * srind2 = srind - 1;
         */
        mlfAssign(&srind2, mclMinus(mclVv(srind, "srind"), _mxarray0_));
        /*
         * % create memory for saving xcorr data
         * dxcorr = zeros(lagslength,numwindows);
         */
        mlfAssign(
          &dxcorr,
          mlfZeros(
            mclVv(lagslength, "lagslength"),
            mclVv(numwindows, "numwindows"),
            NULL));
        /*
         * dxcz = repmat(nan,lagslength,numwindows);
         */
        mlfAssign(
          &dxcz,
          mlfRepmat(
            _mxarray8_,
            mclVv(lagslength, "lagslength"),
            mclVv(numwindows, "numwindows")));
        /*
         * dxcmean = dxcz;
         */
        mlfAssign(&dxcmean, mclVv(dxcz, "dxcz"));
        /*
         * dxcstd = dxcz;
         */
        mlfAssign(&dxcstd, mclVv(dxcz, "dxcz"));
        /*
         * % create matrix to save result of testing psth distributions
         * hpval = repmat(nan,2,numwindows);
         */
        mlfAssign(
          &hpval,
          mlfRepmat(_mxarray8_, _mxarray1_, mclVv(numwindows, "numwindows")));
        /*
         * % generate vector of lagslength
         * lagsvec = 1:lagslength;
         */
        mlfAssign(
          &lagsvec,
          mlfColon(_mxarray0_, mclVv(lagslength, "lagslength"), NULL));
        /*
         * % generate matrix for result of jbtest
         * jbmat = repmat(nan,lagslength,numwindows);
         */
        mlfAssign(
          &jbmat,
          mlfRepmat(
            _mxarray8_,
            mclVv(lagslength, "lagslength"),
            mclVv(numwindows, "numwindows")));
        /*
         * % generate matrix for result of lillietest
         * limat = jbmat;
         */
        mlfAssign(&limat, mclVv(jbmat, "jbmat"));
        /*
         * 
         * % get vector that picks out just the data without the shuffle edges
         * sadd = shuffle + 1;
         */
        mlfAssign(&sadd, mclPlus(mclVv(shuffle, "shuffle"), _mxarray0_));
        /*
         * % get the time of the subbins
         * bintimes1 = subbins;
         */
        mlfAssign(&bintimes1, mclVv(subbins, "subbins"));
        /*
         * % get length of bintimes1
         * bt1l = length(bintimes1);
         */
        mlfAssign(
          &bt1l, mlfScalar(mclLengthInt(mclVv(bintimes1, "bintimes1"))));
        /*
         * % get binsize so we can add the necessary bins for the shuffles that 
         * % go off the edge
         * bsize = bintimes1(2)-bintimes1(1);
         */
        mlfAssign(
          &bsize,
          mclMinus(
            mclIntArrayRef1(mclVv(bintimes1, "bintimes1"), 2),
            mclIntArrayRef1(mclVv(bintimes1, "bintimes1"), 1)));
        /*
         * % get extra bins
         * extrabins = (1:shuffle)'*bsize;
         */
        mlfAssign(
          &extrabins,
          mlf_times_transpose(
            mlfColon(_mxarray0_, mclVv(shuffle, "shuffle"), NULL),
            mclVv(bsize, "bsize"),
            _mxarray9_));
        /*
         * bintimes = [-flipud(extrabins); bintimes1; bintimes1(bt1l) + extrabins];
         */
        mlfAssign(
          &bintimes,
          mlfVertcat(
            mclUminus(mlfFlipud(mclVv(extrabins, "extrabins"))),
            mclVv(bintimes1, "bintimes1"),
            mclPlus(
              mclArrayRef1(mclVv(bintimes1, "bintimes1"), mclVv(bt1l, "bt1l")),
              mclVv(extrabins, "extrabins")),
            NULL));
        /*
         * dvec = sadd:(windowlength+shuffle);
         */
        mlfAssign(
          &dvec,
          mlfColon(
            mclVv(sadd, "sadd"),
            mclPlus(
              mclVv(windowlength, "windowlength"), mclVv(shuffle, "shuffle")),
            NULL));
        /*
         * pbins = windowlength+(2*shuffle);
         */
        mlfAssign(
          &pbins,
          mclPlus(
            mclVv(windowlength, "windowlength"),
            mclMtimes(_mxarray1_, mclVv(shuffle, "shuffle"))));
        /*
         * psthbins = repmat(nan,pbins,numwindows);
         */
        mlfAssign(
          &psthbins,
          mlfRepmat(
            _mxarray8_,
            mclVv(pbins, "pbins"),
            mclVv(numwindows, "numwindows")));
        /*
         * % set second index first so we don't have to change memory size
         * psths{2} = psthbins;
         */
        mlfIndexAssign(&psths, "{?}", _mxarray1_, mclVv(psthbins, "psthbins"));
        /*
         * psths{1} = psthbins;
         */
        mlfIndexAssign(&psths, "{?}", _mxarray0_, mclVv(psthbins, "psthbins"));
        /*
         * 
         * ns1 = numSurrogates + 1;
         */
        mlfAssign(
          &ns1, mclPlus(mclVv(numSurrogates, "numSurrogates"), _mxarray0_));
        /*
         * ns2 = numSurrogates + 2;
         */
        mlfAssign(
          &ns2, mclPlus(mclVv(numSurrogates, "numSurrogates"), _mxarray1_));
        /*
         * ws1 = psthwin - 1;
         */
        mlfAssign(&ws1, mclMinus(mclVv(psthwin, "psthwin"), _mxarray0_));
        /*
         * % get half window size
         * halfwinsize = ws1 / 2;
         */
        mlfAssign(&halfwinsize, mclMrdivide(mclVv(ws1, "ws1"), _mxarray1_));
        /*
         * % get length of central window
         * lcentralwin = windowlength;
         */
        mlfAssign(&lcentralwin, mclVv(windowlength, "windowlength"));
        /*
         * % compute length of psthbins and save for use later
         * lpsthbins = lcentralwin + ws1;
         */
        mlfAssign(
          &lpsthbins,
          mclPlus(mclVv(lcentralwin, "lcentralwin"), mclVv(ws1, "ws1")));
        /*
         * ps1ind = [numsurrvec ns1];
         */
        mlfAssign(
          &ps1ind,
          mlfHorzcat(mclVv(numsurrvec, "numsurrvec"), mclVv(ns1, "ns1"), NULL));
        /*
         * ps1ind2 = (2 * ns1) + [1 2];
         */
        mlfAssign(
          &ps1ind2,
          mclPlus(mclMtimes(_mxarray1_, mclVv(ns1, "ns1")), _mxarray10_));
        /*
         * ps2ind = ps1ind + ns1;
         */
        mlfAssign(&ps2ind, mclPlus(mclVv(ps1ind, "ps1ind"), mclVv(ns1, "ns1")));
        /*
         * ps2ind2 = ps1ind2 + 2;
         */
        mlfAssign(&ps2ind2, mclPlus(mclVv(ps1ind2, "ps1ind2"), _mxarray1_));
        /*
         * % get indices representing window
         * cwinidx = 1:lcentralwin;
         */
        mlfAssign(
          &cwinidx,
          mlfColon(_mxarray0_, mclVv(lcentralwin, "lcentralwin"), NULL));
        /*
         * cwinidx2 = cwinidx+lcentralwin;
         */
        mlfAssign(
          &cwinidx2,
          mclPlus(
            mclVv(cwinidx, "cwinidx"), mclVv(lcentralwin, "lcentralwin")));
        /*
         * % compute the relevant indices
         * didx = cwinidx + ws1;
         */
        mlfAssign(&didx, mclPlus(mclVv(cwinidx, "cwinidx"), mclVv(ws1, "ws1")));
        /*
         * surridx1 = 2:ns1;
         */
        mlfAssign(&surridx1, mlfColon(_mxarray1_, mclVv(ns1, "ns1"), NULL));
        /*
         * surridx2 = surridx1 + ns1;
         */
        mlfAssign(
          &surridx2, mclPlus(mclVv(surridx1, "surridx1"), mclVv(ns1, "ns1")));
        /*
         * smoothmtx = convmtx(ones(psthwin,1),lpsthbins);
         */
        mlfAssign(
          &smoothmtx,
          mlfConvmtx(
            mlfOnes(mclVv(psthwin, "psthwin"), _mxarray0_, NULL),
            mclVv(lpsthbins, "lpsthbins")));
        /*
         * defaultspt = repmat(nan,1,ns1);
         */
        mlfAssign(
          &defaultspt, mlfRepmat(_mxarray8_, _mxarray0_, mclVv(ns1, "ns1")));
        /*
         * dspsthz = zeros(numwindows,2*lcentralwin);
         */
        mlfAssign(
          &dspsthz,
          mlfZeros(
            mclVv(numwindows, "numwindows"),
            mclMtimes(_mxarray1_, mclVv(lcentralwin, "lcentralwin")),
            NULL));
        /*
         * 
         * % pdf for truncated normal distribution
         * % taken from example from Mathworks:
         * % http://www.mathworks.com/products/statistics/demos.html?file=/products/demos/shipping/stats/customdist2demo.html
         * % but since we have values of 0, we normalize by the cummulative
         * % pdf for the value of -1
         * % pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma) ./ (1-normcdf(-1,mu,sigma));
         * 
         * % turn off warning
         * warning off MATLAB:divideByZero
         */
        mclPrintAns(&ans, mlfNWarning(0, NULL, _mxarray12_, _mxarray14_, NULL));
        /*
         * for sidx = 1:numwindows
         */
        {
            int v_ = mclForIntStart(1);
            int e_ = mclForIntEnd(mclVv(numwindows, "numwindows"));
            if (v_ > e_) {
                mlfAssign(&sidx, _mxarray16_);
            } else {
                /*
                 * l = load([surrprefix num2str(sidx,'%04d')]);
                 * % if there were no spikes or if there were no reps in which there
                 * % were spikes in both cells, skip this analysis since there would
                 * % no xcorr values to look at.
                 * if(~isempty(l.sptrains))
                 * % *** Compute Synchronous Spikes ***
                 * % sum counts within maxlags for data and surrogates by creating 
                 * % appropriate vector
                 * nsyncspikes = sumvec * l.xc;
                 * % separate the surrogates from the data
                 * data(sidx) = nsyncspikes(DATA_COL);
                 * shiftp(sidx) = nsyncspikes(SHIFT_COL);
                 * sdata(:,sidx) = nsyncspikes(srind)';
                 * 
                 * % *** Compute XCORR Z-Scores and test for Normality ***
                 * % save the xcorr data so we don't have to compute it in the plot 
                 * % function
                 * xcdata = l.xc(:,DATA_COL);
                 * % get surrogate data
                 * xcsdata = l.xc(:,srind)';
                 * % initialize smean and sstd to zeros
                 * smean = smstmp;
                 * sstd = smstmp;
                 * % compute the jbtest to test for presence of normal distributions
                 * for lagi = lagsvec
                 * testdata = full(xcsdata(:,lagi));
                 * [dummy,jbmat(lagi,sidx)] = jbtest(testdata);
                 * [dummy,limat(lagi,sidx)] = lillietest(testdata);
                 * %                 if(sum(testdata)>0)
                 * %                     start = [mean(testdata) std(testdata)];
                 * %                     try
                 * %                         [phat,pci] = mle(testdata,'pdf',pdf_truncnorm,'start',start,'lower',[-Inf 0]);
                 * %                         % fprintf('%d mean: %f std: %f\n',lagi,phat(1),phat(2));
                 * %                         smean(lagi) = phat(1);
                 * %                         sstd(lagi) = phat(2);
                 * %                     catch
                 * %                         fprintf('%d %s\n',lagi,lasterr);
                 * %                         % phat = nan;
                 * %                     end
                 * %                 %else
                 * %                 %    fprintf('%d Skipped\n',lagi);
                 * %                 end
                 * end
                 * % compute the mean and std of the surrogates
                 * smean = mean(xcsdata)';
                 * sstd = std(xcsdata)';
                 * % convert data to z-score using the mean and std of the surrogates
                 * dxcz(:,sidx) = (xcdata-smean)./sstd;
                 * dxcorr(:,sidx) = xcdata;
                 * dxcmean(:,sidx) = smean;
                 * dxcstd(:,sidx) = sstd;	 
                 * 
                 * % *** Test smoothed PSTH for differences between data and surrogates ***
                 * % load frame before and frame after, which is 2 windows before 
                 * % and after since we are stepping in half frame increments
                 * % convert cell arrays to matrices
                 * % store in cell array since there is a for-loop later which
                 * % will make use of that data
                 * if(~isempty(l.sptrains{2}))
                 * spt{2} = full(cell2mat(l.sptrains{2}));
                 * else
                 * spt{2} = defaultspt;
                 * end
                 * if(~isempty(l.sptrains{1}))
                 * spt{1} = full(cell2mat(l.sptrains{1}));
                 * else
                 * spt{1} = defaultspt;
                 * end
                 * l0sp1 = [];
                 * l0sp2 = [];
                 * if(sidx>2)
                 * l0 = load([surrprefix num2str(sidx-2,'%04d')]);
                 * if(~isempty(l0.sptrains))
                 * if(~isempty(l0.sptrains{1}))
                 * l0spt1 = full(cell2mat(l0.sptrains{1}));
                 * l0sp1 = l0spt1(:,1);
                 * end
                 * if(~isempty(l0.sptrains{2}))
                 * l0spt2 = full(cell2mat(l0.sptrains{2}));
                 * l0sp2 = l0spt2(:,1);
                 * end
                 * end
                 * end
                 * l1sp1 = [];
                 * l1sp2 = [];
                 * if(sidx<(numwindows-1))
                 * l1 = load([surrprefix num2str(sidx+2,'%04d')]);
                 * if(~isempty(l1.sptrains))
                 * if(~isempty(l1.sptrains{1}))
                 * l1spt1 = full(cell2mat(l1.sptrains{1}));
                 * l1sp1 = l1spt1(:,1);
                 * end
                 * if(~isempty(l1.sptrains{2}))
                 * l1spt2 = full(cell2mat(l1.sptrains{2}));
                 * l1sp2 = l1spt2(:,1);		
                 * end
                 * end
                 * end
                 * % concatenate data together in order to perform histogram just
                 * % once
                 * smat = concat(spt{1},spt{2},l0sp1,l1sp1, ...
                 * l0sp2,l1sp2,'Columnwise');
                 * % shift by 0.5 at the beginning and the end so we can use
                 * % histcie instead of hist as hist includes all values outside
                 * % the range in the first and last bins
                 * spsthbins = (startbins(sidx)-halfwinsize):(endbins(sidx)+halfwinsize);
                 * psthhistcbins = [spsthbins-0.5 spsthbins(lpsthbins)+0.5];
                 * psth = histcie(smat,psthhistcbins,'DropLast','DataCols');
                 * % add data from frame before and frame after before computing
                 * % smoothed PSTH
                 * psth1 = psth(:,ps1ind) + repmat(sum(psth(:,ps1ind2),2),1,ns1);
                 * psth2 = psth(:,ps2ind) + repmat(sum(psth(:,ps2ind2),2),1,ns1);
                 * % now compute the sliding window smoothed psth
                 * spsth = smoothmtx * [psth1 psth2];
                 * % grab the data
                 * dspsth1 = spsth(didx,1);
                 * dspsth2 = spsth(didx,ns2);
                 * % grab the surrogates
                 * sspsthmat = [spsth(didx,surridx1)' spsth(didx,surridx2)'];
                 * % compute mean and std of the surrogates
                 * sspm = mean(sspsthmat);
                 * sspstd = std(sspsthmat);
                 * % compute z-score for the data
                 * dspsthz(sidx,:) = ([dspsth1' dspsth2'] - sspm) ./ sspstd;
                 * 
                 * loopbins = ( (startbins(sidx)-shuffle):(endbins(sidx)+shuffle) )';
                 * for cidx = cellvec
                 * % get the spike times of the first cell
                 * % call the full function to make sure any sparse matrices 
                 * % gets converted to full matrices otherwise kstest2 will
                 * % not work. Sparse matrices are possible now since binidx 
                 * % is saved as a sparse matrix in shufflesync and the code
                 * % in shufflesyncsurr will index into binidx and save it 
                 * % directly to l.sptrains.
                 * sptrains = spt{cidx};
                 * if(~isempty(sptrains))
                 * % get the spike times of the data
                 * sptimes = sptrains(:,DATA_COL);
                 * % the kstest2 is not acurate if n1*n2/(n1+n2) < 4. Since
                 * % the number of surrogate spikes is always going to be 1000
                 * % times the number of real spikes, this simplifies to
                 * % 1000n^2/1001n < 4, which further simplies to 
                 * % n < 4001/1000, which is pretty close to 4
                 * if(length(sptimes)>KSNSPIKES)
                 * % get the spike times of the surrogates
                 * surrtimes = sptrains(:,srind2);
                 * srtimes = surrtimes(:);
                 * % do test for difference
                 * [dummy,hpval(cidx,sidx)] = kstest2(sptimes,srtimes);
                 * end
                 * 
                 * % check to see if we need to pad sptrains so the histogram is
                 * % taken in the right dimension
                 * if(size(sptrains,1)==1)
                 * sptrains = concatenate(sptrains,nan);
                 * end
                 * loopcounts = hist(sptrains,loopbins);
                 * psths{cidx}(:,sidx) = mean(loopcounts(:,surrvec1),2);
                 * end
                 * end % for cidx = cellvec
                 * psthbins(:,sidx) = loopbins;
                 * end % if(~isempty(l.sptrains))
                 * end % for sidx = 1:numwindows
                 */
                for (; ; ) {
                    mlfAssign(
                      &l,
                      mlfLoadStruct(
                        mlfHorzcat(
                          mclVv(surrprefix, "surrprefix"),
                          mlfNum2str(mlfScalar(v_), _mxarray17_),
                          NULL),
                        NULL));
                    if (mclNotBool(
                          mclFeval(
                            mclValueVarargout(),
                            mlxIsempty,
                            mlfIndexRef(mclVv(l, "l"), ".sptrains"),
                            NULL))) {
                        mlfAssign(
                          &nsyncspikes,
                          mclFeval(
                            mclValueVarargout(),
                            mlxMtimes,
                            mclVv(sumvec, "sumvec"),
                            mlfIndexRef(mclVv(l, "l"), ".xc"),
                            NULL));
                        mclIntArrayAssign1(
                          &data,
                          mclArrayRef1(
                            mclVv(nsyncspikes, "nsyncspikes"),
                            mclVv(DATA_COL, "DATA_COL")),
                          v_);
                        mclIntArrayAssign1(
                          &shiftp,
                          mclArrayRef1(
                            mclVv(nsyncspikes, "nsyncspikes"),
                            mclVv(SHIFT_COL, "SHIFT_COL")),
                          v_);
                        mclArrayAssign2(
                          &sdata,
                          mlfCtranspose(
                            mclArrayRef1(
                              mclVv(nsyncspikes, "nsyncspikes"),
                              mclVv(srind, "srind"))),
                          mlfCreateColonIndex(),
                          mlfScalar(v_));
                        mlfAssign(
                          &xcdata,
                          mlfIndexRef(
                            mclVv(l, "l"),
                            ".xc(?,?)",
                            mlfCreateColonIndex(),
                            mclVv(DATA_COL, "DATA_COL")));
                        mlfAssign(
                          &xcsdata,
                          mclFeval(
                            mclValueVarargout(),
                            mlxCtranspose,
                            mlfIndexRef(
                              mclVv(l, "l"),
                              ".xc(?,?)",
                              mlfCreateColonIndex(),
                              mclVv(srind, "srind")),
                            NULL));
                        mlfAssign(&smean, mclVv(smstmp, "smstmp"));
                        mlfAssign(&sstd, mclVv(smstmp, "smstmp"));
                        {
                            mclForLoopIterator viter__;
                            for (mclForStart(
                                   &viter__,
                                   mclVv(lagsvec, "lagsvec"),
                                   NULL,
                                   NULL);
                                 mclForNext(&viter__, &lagi);
                                 ) {
                                mlfAssign(
                                  &testdata,
                                  mlfFull(
                                    mclArrayRef2(
                                      mclVv(xcsdata, "xcsdata"),
                                      mlfCreateColonIndex(),
                                      mclVv(lagi, "lagi"))));
                                mclFeval(
                                  mlfIndexVarargout(
                                    &dummy, "",
                                    &jbmat,
                                    "(?,?)",
                                    mclVv(lagi, "lagi"),
                                    mlfScalar(v_),
                                    NULL),
                                  mlxJbtest,
                                  mclVv(testdata, "testdata"),
                                  NULL);
                                mclFeval(
                                  mlfIndexVarargout(
                                    &dummy, "",
                                    &limat,
                                    "(?,?)",
                                    mclVv(lagi, "lagi"),
                                    mlfScalar(v_),
                                    NULL),
                                  mlxLillietest,
                                  mclVv(testdata, "testdata"),
                                  NULL);
                            }
                            mclDestroyForLoopIterator(viter__);
                        }
                        mlfAssign(
                          &smean,
                          mlfCtranspose(
                            mlfMean(mclVv(xcsdata, "xcsdata"), NULL)));
                        mlfAssign(
                          &sstd,
                          mlfCtranspose(
                            mlfStd(mclVv(xcsdata, "xcsdata"), NULL, NULL)));
                        mclArrayAssign2(
                          &dxcz,
                          mclRdivide(
                            mclMinus(
                              mclVv(xcdata, "xcdata"), mclVv(smean, "smean")),
                            mclVv(sstd, "sstd")),
                          mlfCreateColonIndex(),
                          mlfScalar(v_));
                        mclArrayAssign2(
                          &dxcorr,
                          mclVv(xcdata, "xcdata"),
                          mlfCreateColonIndex(),
                          mlfScalar(v_));
                        mclArrayAssign2(
                          &dxcmean,
                          mclVv(smean, "smean"),
                          mlfCreateColonIndex(),
                          mlfScalar(v_));
                        mclArrayAssign2(
                          &dxcstd,
                          mclVv(sstd, "sstd"),
                          mlfCreateColonIndex(),
                          mlfScalar(v_));
                        if (mclNotBool(
                              mclFeval(
                                mclValueVarargout(),
                                mlxIsempty,
                                mlfIndexRef(
                                  mclVv(l, "l"), ".sptrains{?}", _mxarray1_),
                                NULL))) {
                            mlfIndexAssign(
                              &spt,
                              "{?}",
                              _mxarray1_,
                              mlfFull(
                                mclFeval(
                                  mclValueVarargout(),
                                  mlxCell2mat,
                                  mlfIndexRef(
                                    mclVv(l, "l"), ".sptrains{?}", _mxarray1_),
                                  NULL)));
                        } else {
                            mlfIndexAssign(
                              &spt,
                              "{?}",
                              _mxarray1_,
                              mclVv(defaultspt, "defaultspt"));
                        }
                        if (mclNotBool(
                              mclFeval(
                                mclValueVarargout(),
                                mlxIsempty,
                                mlfIndexRef(
                                  mclVv(l, "l"), ".sptrains{?}", _mxarray0_),
                                NULL))) {
                            mlfIndexAssign(
                              &spt,
                              "{?}",
                              _mxarray0_,
                              mlfFull(
                                mclFeval(
                                  mclValueVarargout(),
                                  mlxCell2mat,
                                  mlfIndexRef(
                                    mclVv(l, "l"), ".sptrains{?}", _mxarray0_),
                                  NULL)));
                        } else {
                            mlfIndexAssign(
                              &spt,
                              "{?}",
                              _mxarray0_,
                              mclVv(defaultspt, "defaultspt"));
                        }
                        mlfAssign(&l0sp1, _mxarray16_);
                        mlfAssign(&l0sp2, _mxarray16_);
                        if (mclGtBool(mlfScalar(v_), _mxarray1_)) {
                            mlfAssign(
                              &l0,
                              mlfLoadStruct(
                                mlfHorzcat(
                                  mclVv(surrprefix, "surrprefix"),
                                  mlfNum2str(mlfScalar(v_ - 2), _mxarray17_),
                                  NULL),
                                NULL));
                            if (mclNotBool(
                                  mclFeval(
                                    mclValueVarargout(),
                                    mlxIsempty,
                                    mlfIndexRef(mclVv(l0, "l0"), ".sptrains"),
                                    NULL))) {
                                if (mclNotBool(
                                      mclFeval(
                                        mclValueVarargout(),
                                        mlxIsempty,
                                        mlfIndexRef(
                                          mclVv(l0, "l0"),
                                          ".sptrains{?}",
                                          _mxarray0_),
                                        NULL))) {
                                    mlfAssign(
                                      &l0spt1,
                                      mlfFull(
                                        mclFeval(
                                          mclValueVarargout(),
                                          mlxCell2mat,
                                          mlfIndexRef(
                                            mclVv(l0, "l0"),
                                            ".sptrains{?}",
                                            _mxarray0_),
                                          NULL)));
                                    mlfAssign(
                                      &l0sp1,
                                      mclArrayRef2(
                                        mclVv(l0spt1, "l0spt1"),
                                        mlfCreateColonIndex(),
                                        _mxarray0_));
                                }
                                if (mclNotBool(
                                      mclFeval(
                                        mclValueVarargout(),
                                        mlxIsempty,
                                        mlfIndexRef(
                                          mclVv(l0, "l0"),
                                          ".sptrains{?}",
                                          _mxarray1_),
                                        NULL))) {
                                    mlfAssign(
                                      &l0spt2,
                                      mlfFull(
                                        mclFeval(
                                          mclValueVarargout(),
                                          mlxCell2mat,
                                          mlfIndexRef(
                                            mclVv(l0, "l0"),
                                            ".sptrains{?}",
                                            _mxarray1_),
                                          NULL)));
                                    mlfAssign(
                                      &l0sp2,
                                      mclArrayRef2(
                                        mclVv(l0spt2, "l0spt2"),
                                        mlfCreateColonIndex(),
                                        _mxarray0_));
                                }
                            }
                        }
                        mlfAssign(&l1sp1, _mxarray16_);
                        mlfAssign(&l1sp2, _mxarray16_);
                        if (mclLtBool(
                              mlfScalar(v_),
                              mclMinus(
                                mclVv(numwindows, "numwindows"), _mxarray0_))) {
                            mlfAssign(
                              &l1,
                              mlfLoadStruct(
                                mlfHorzcat(
                                  mclVv(surrprefix, "surrprefix"),
                                  mlfNum2str(mlfScalar(v_ + 2), _mxarray17_),
                                  NULL),
                                NULL));
                            if (mclNotBool(
                                  mclFeval(
                                    mclValueVarargout(),
                                    mlxIsempty,
                                    mlfIndexRef(mclVv(l1, "l1"), ".sptrains"),
                                    NULL))) {
                                if (mclNotBool(
                                      mclFeval(
                                        mclValueVarargout(),
                                        mlxIsempty,
                                        mlfIndexRef(
                                          mclVv(l1, "l1"),
                                          ".sptrains{?}",
                                          _mxarray0_),
                                        NULL))) {
                                    mlfAssign(
                                      &l1spt1,
                                      mlfFull(
                                        mclFeval(
                                          mclValueVarargout(),
                                          mlxCell2mat,
                                          mlfIndexRef(
                                            mclVv(l1, "l1"),
                                            ".sptrains{?}",
                                            _mxarray0_),
                                          NULL)));
                                    mlfAssign(
                                      &l1sp1,
                                      mclArrayRef2(
                                        mclVv(l1spt1, "l1spt1"),
                                        mlfCreateColonIndex(),
                                        _mxarray0_));
                                }
                                if (mclNotBool(
                                      mclFeval(
                                        mclValueVarargout(),
                                        mlxIsempty,
                                        mlfIndexRef(
                                          mclVv(l1, "l1"),
                                          ".sptrains{?}",
                                          _mxarray1_),
                                        NULL))) {
                                    mlfAssign(
                                      &l1spt2,
                                      mlfFull(
                                        mclFeval(
                                          mclValueVarargout(),
                                          mlxCell2mat,
                                          mlfIndexRef(
                                            mclVv(l1, "l1"),
                                            ".sptrains{?}",
                                            _mxarray1_),
                                          NULL)));
                                    mlfAssign(
                                      &l1sp2,
                                      mclArrayRef2(
                                        mclVv(l1spt2, "l1spt2"),
                                        mlfCreateColonIndex(),
                                        _mxarray0_));
                                }
                            }
                        }
                        mlfAssign(
                          &smat,
                          mclFeval(
                            mclValueVarargout(),
                            mlxConcat,
                            mlfIndexRef(mclVv(spt, "spt"), "{?}", _mxarray0_),
                            mlfIndexRef(mclVv(spt, "spt"), "{?}", _mxarray1_),
                            mclVv(l0sp1, "l0sp1"),
                            mclVv(l1sp1, "l1sp1"),
                            mclVv(l0sp2, "l0sp2"),
                            mclVv(l1sp2, "l1sp2"),
                            _mxarray19_,
                            NULL));
                        mlfAssign(
                          &spsthbins,
                          mlfColon(
                            mclMinus(
                              mclIntArrayRef1(
                                mclVv(startbins, "startbins"), v_),
                              mclVv(halfwinsize, "halfwinsize")),
                            mclPlus(
                              mclIntArrayRef1(mclVv(endbins, "endbins"), v_),
                              mclVv(halfwinsize, "halfwinsize")),
                            NULL));
                        mlfAssign(
                          &psthhistcbins,
                          mlfHorzcat(
                            mclMinus(
                              mclVv(spsthbins, "spsthbins"), _mxarray21_),
                            mclPlus(
                              mclArrayRef1(
                                mclVv(spsthbins, "spsthbins"),
                                mclVv(lpsthbins, "lpsthbins")),
                              _mxarray21_),
                            NULL));
                        mlfAssign(
                          &psth,
                          mlfNHistcie(
                            1,
                            NULL,
                            mclVv(smat, "smat"),
                            mclVv(psthhistcbins, "psthhistcbins"),
                            _mxarray22_,
                            _mxarray24_,
                            NULL));
                        mlfAssign(
                          &psth1,
                          mclPlus(
                            mclArrayRef2(
                              mclVv(psth, "psth"),
                              mlfCreateColonIndex(),
                              mclVv(ps1ind, "ps1ind")),
                            mlfRepmat(
                              mlfSum(
                                mclArrayRef2(
                                  mclVv(psth, "psth"),
                                  mlfCreateColonIndex(),
                                  mclVv(ps1ind2, "ps1ind2")),
                                _mxarray1_),
                              _mxarray0_,
                              mclVv(ns1, "ns1"))));
                        mlfAssign(
                          &psth2,
                          mclPlus(
                            mclArrayRef2(
                              mclVv(psth, "psth"),
                              mlfCreateColonIndex(),
                              mclVv(ps2ind, "ps2ind")),
                            mlfRepmat(
                              mlfSum(
                                mclArrayRef2(
                                  mclVv(psth, "psth"),
                                  mlfCreateColonIndex(),
                                  mclVv(ps2ind2, "ps2ind2")),
                                _mxarray1_),
                              _mxarray0_,
                              mclVv(ns1, "ns1"))));
                        mlfAssign(
                          &spsth,
                          mclMtimes(
                            mclVv(smoothmtx, "smoothmtx"),
                            mlfHorzcat(
                              mclVv(psth1, "psth1"),
                              mclVv(psth2, "psth2"),
                              NULL)));
                        mlfAssign(
                          &dspsth1,
                          mclArrayRef2(
                            mclVv(spsth, "spsth"),
                            mclVv(didx, "didx"),
                            _mxarray0_));
                        mlfAssign(
                          &dspsth2,
                          mclArrayRef2(
                            mclVv(spsth, "spsth"),
                            mclVv(didx, "didx"),
                            mclVv(ns2, "ns2")));
                        mlfAssign(
                          &sspsthmat,
                          mlfHorzcat(
                            mlfCtranspose(
                              mclArrayRef2(
                                mclVv(spsth, "spsth"),
                                mclVv(didx, "didx"),
                                mclVv(surridx1, "surridx1"))),
                            mlfCtranspose(
                              mclArrayRef2(
                                mclVv(spsth, "spsth"),
                                mclVv(didx, "didx"),
                                mclVv(surridx2, "surridx2"))),
                            NULL));
                        mlfAssign(
                          &sspm, mlfMean(mclVv(sspsthmat, "sspsthmat"), NULL));
                        mlfAssign(
                          &sspstd,
                          mlfStd(mclVv(sspsthmat, "sspsthmat"), NULL, NULL));
                        mclArrayAssign2(
                          &dspsthz,
                          mclRdivide(
                            mclMinus(
                              mlfHorzcat(
                                mlfCtranspose(mclVv(dspsth1, "dspsth1")),
                                mlfCtranspose(mclVv(dspsth2, "dspsth2")),
                                NULL),
                              mclVv(sspm, "sspm")),
                            mclVv(sspstd, "sspstd")),
                          mlfScalar(v_),
                          mlfCreateColonIndex());
                        mlfAssign(
                          &loopbins,
                          mlfCtranspose(
                            mlfColon(
                              mclMinus(
                                mclIntArrayRef1(
                                  mclVv(startbins, "startbins"), v_),
                                mclVv(shuffle, "shuffle")),
                              mclPlus(
                                mclIntArrayRef1(mclVv(endbins, "endbins"), v_),
                                mclVv(shuffle, "shuffle")),
                              NULL)));
                        {
                            mclForLoopIterator viter__;
                            for (mclForStart(
                                   &viter__,
                                   mclVv(cellvec, "cellvec"),
                                   NULL,
                                   NULL);
                                 mclForNext(&viter__, &cidx);
                                 ) {
                                mlfAssign(
                                  &sptrains,
                                  mlfIndexRef(
                                    mclVv(spt, "spt"),
                                    "{?}",
                                    mclVv(cidx, "cidx")));
                                if (mclNotBool(
                                      mlfIsempty(
                                        mclVv(sptrains, "sptrains")))) {
                                    mlfAssign(
                                      &sptimes,
                                      mclArrayRef2(
                                        mclVv(sptrains, "sptrains"),
                                        mlfCreateColonIndex(),
                                        mclVv(DATA_COL, "DATA_COL")));
                                    if (mclGtBool(
                                          mlfScalar(
                                            mclLengthInt(
                                              mclVv(sptimes, "sptimes"))),
                                          mclVv(KSNSPIKES, "KSNSPIKES"))) {
                                        mlfAssign(
                                          &surrtimes,
                                          mclArrayRef2(
                                            mclVv(sptrains, "sptrains"),
                                            mlfCreateColonIndex(),
                                            mclVv(srind2, "srind2")));
                                        mlfAssign(
                                          &srtimes,
                                          mclArrayRef1(
                                            mclVv(surrtimes, "surrtimes"),
                                            mlfCreateColonIndex()));
                                        mclFeval(
                                          mlfIndexVarargout(
                                            &dummy, "",
                                            &hpval,
                                            "(?,?)",
                                            mclVv(cidx, "cidx"),
                                            mlfScalar(v_),
                                            NULL),
                                          mlxKstest2,
                                          mclVv(sptimes, "sptimes"),
                                          mclVv(srtimes, "srtimes"),
                                          NULL);
                                    }
                                    if (mclEqBool(
                                          mlfSize(
                                            mclValueVarargout(),
                                            mclVv(sptrains, "sptrains"),
                                            _mxarray0_),
                                          _mxarray0_)) {
                                        mlfAssign(
                                          &sptrains,
                                          mlfConcatenate(
                                            mclVv(sptrains, "sptrains"),
                                            _mxarray8_,
                                            NULL));
                                    }
                                    mlfAssign(
                                      &loopcounts,
                                      mlfNHist(
                                        1,
                                        NULL,
                                        mclVv(sptrains, "sptrains"),
                                        mclVv(loopbins, "loopbins")));
                                    mlfIndexAssign(
                                      &psths,
                                      "{?}(?,?)",
                                      mclVv(cidx, "cidx"),
                                      mlfCreateColonIndex(),
                                      mlfScalar(v_),
                                      mlfMean(
                                        mclArrayRef2(
                                          mclVv(loopcounts, "loopcounts"),
                                          mlfCreateColonIndex(),
                                          mclVv(surrvec1, "surrvec1")),
                                        _mxarray1_));
                                }
                            }
                            mclDestroyForLoopIterator(viter__);
                        }
                        mclArrayAssign2(
                          &psthbins,
                          mclVv(loopbins, "loopbins"),
                          mlfCreateColonIndex(),
                          mlfScalar(v_));
                    }
                    if (v_ == e_) {
                        break;
                    }
                    ++v_;
                }
                mlfAssign(&sidx, mlfScalar(v_));
            }
        }
        /*
         * 
         * % compute histogram of surrogates
         * surrmean = mean(sdata)';
         */
        mlfAssign(
          surrmean, mlfCtranspose(mlfMean(mclVv(sdata, "sdata"), NULL)));
        /*
         * surrstd = std(sdata)';
         */
        mlfAssign(
          surrstd, mlfCtranspose(mlfStd(mclVv(sdata, "sdata"), NULL, NULL)));
        /*
         * % set bins so it covers the max of the data as well as of the surrogates
         * maxn = max([sdata(:); data(:); shiftp(:)]);
         */
        mlfAssign(
          &maxn,
          mlfMax(
            NULL,
            mlfVertcat(
              mclArrayRef1(mclVv(sdata, "sdata"), mlfCreateColonIndex()),
              mclArrayRef1(mclVv(data, "data"), mlfCreateColonIndex()),
              mclArrayRef1(mclVv(shiftp, "shiftp"), mlfCreateColonIndex()),
              NULL),
            NULL,
            NULL));
        /*
         * bins = SurrHistMin:SurrHistStep:(ceil(maxn/SurrHistStep) ...
         */
        mlfAssign(
          &bins,
          mlfColon(
            mclVv(SurrHistMin, "SurrHistMin"),
            mclVv(SurrHistStep, "SurrHistStep"),
            mclMtimes(
              mlfCeil(
                mclMrdivide(
                  mclVv(maxn, "maxn"), mclVv(SurrHistStep, "SurrHistStep"))),
              mclVv(SurrHistStep, "SurrHistStep"))));
        /*
         * * SurrHistStep);
         * % do hist since values are integers
         * surrhist = hist(sdata,bins);
         */
        mlfAssign(
          surrhist,
          mlfNHist(1, NULL, mclVv(sdata, "sdata"), mclVv(bins, "bins")));
        /*
         * % find bins with std not equal to 0
         * si = find(surrstd);
         */
        mlfAssign(&si, mlfFind(NULL, NULL, mclVv(*surrstd, "surrstd")));
        /*
         * % create z-score vector consisting of nans
         * dzscore = repmat(nan,numwindows,1);
         */
        mlfAssign(
          &dzscore,
          mlfRepmat(_mxarray8_, mclVv(numwindows, "numwindows"), _mxarray0_));
        /*
         * % fill bins that don't have std of 0
         * dzscore(si) = (data(si) - surrmean(si)) ./ surrstd(si);
         */
        mclArrayAssign1(
          &dzscore,
          mclRdivide(
            mclMinus(
              mclArrayRef1(mclVv(data, "data"), mclVv(si, "si")),
              mclArrayRef1(mclVv(*surrmean, "surrmean"), mclVv(si, "si"))),
            mclArrayRef1(mclVv(*surrstd, "surrstd"), mclVv(si, "si"))),
          mclVv(si, "si"));
        /*
         * % convert matrices with no nan's to sparse matrices to save space
         * data = sparse(data);
         */
        mlfAssign(
          &data, mlfSparse(mclVv(data, "data"), NULL, NULL, NULL, NULL, NULL));
        /*
         * shiftp = sparse(shiftp);
         */
        mlfAssign(
          &shiftp,
          mlfSparse(mclVv(shiftp, "shiftp"), NULL, NULL, NULL, NULL, NULL));
        /*
         * surrmean = sparse(surrmean);
         */
        mlfAssign(
          surrmean,
          mlfSparse(
            mclVv(*surrmean, "surrmean"), NULL, NULL, NULL, NULL, NULL));
        /*
         * surrstd = sparse(surrstd);
         */
        mlfAssign(
          surrstd,
          mlfSparse(mclVv(*surrstd, "surrstd"), NULL, NULL, NULL, NULL, NULL));
        /*
         * surrhist = sparse(surrhist);
         */
        mlfAssign(
          surrhist,
          mlfSparse(
            mclVv(*surrhist, "surrhist"), NULL, NULL, NULL, NULL, NULL));
        /*
         * dxcorr = sparse(dxcorr);
         */
        mlfAssign(
          &dxcorr,
          mlfSparse(mclVv(dxcorr, "dxcorr"), NULL, NULL, NULL, NULL, NULL));
        /*
         * 
         * % calculate the synchronous spikes for the raster display
         * repvec = 1:reps;
         */
        mlfAssign(&repvec, mlfColon(_mxarray0_, mclVv(reps, "reps"), NULL));
        /*
         * sprepvec = circshift(repvec,[1 -1]);
         */
        mlfAssign(
          &sprepvec, mlfCircshift(mclVv(repvec, "repvec"), _mxarray26_));
        /*
         * repvec2 = repvec + reps;
         */
        mlfAssign(
          &repvec2, mclPlus(mclVv(repvec, "repvec"), mclVv(reps, "reps")));
        /*
         * sprepvec2 = sprepvec + reps;
         */
        mlfAssign(
          &sprepvec2,
          mclPlus(mclVv(sprepvec, "sprepvec"), mclVv(reps, "reps")));
        /*
         * nspikes = binidxsize(1);
         */
        mlfAssign(
          &nspikes, mclIntArrayRef1(mclVv(binidxsize, "binidxsize"), 1));
        /*
         * mat1 = repmat(-1,nspikes,2);
         */
        mlfAssign(
          &mat1, mlfRepmat(_mxarray28_, mclVv(nspikes, "nspikes"), _mxarray1_));
        /*
         * mat2 = ones(2,nspikes);
         */
        mlfAssign(&mat2, mlfOnes(_mxarray1_, mclVv(nspikes, "nspikes"), NULL));
        /*
         * markmat1 = zeros(nspikes,reps);
         */
        mlfAssign(
          &markmat1,
          mlfZeros(mclVv(nspikes, "nspikes"), mclVv(reps, "reps"), NULL));
        /*
         * markmat2 = markmat1;
         */
        mlfAssign(&markmat2, mclVv(markmat1, "markmat1"));
        /*
         * spmarkmat1 = markmat1;
         */
        mlfAssign(&spmarkmat1, mclVv(markmat1, "markmat1"));
        /*
         * spmarkmat2 = markmat1;
         */
        mlfAssign(&spmarkmat2, mclVv(markmat1, "markmat1"));
        /*
         * % replace 0's in binidx with nan's so we can distinguish
         * % 0's due to synchronous spikes
         * binidx2 = binidx;
         */
        mlfAssign(&binidx2, mclVv(binidx, "binidx"));
        /*
         * binidx2(binidx2==0) = nan;
         */
        mclArrayAssign1(
          &binidx2, _mxarray8_, mclEq(mclVv(binidx2, "binidx2"), _mxarray29_));
        /*
         * % find the spikes that are within NumCentralBins. Add 1 so that we
         * % can just use the < operator instead of the <= operator 
         * stdiff = (NumCentralBins+1)/2;
         */
        mlfAssign(
          &stdiff,
          mclMrdivide(
            mclPlus(mclVv(NumCentralBins, "NumCentralBins"), _mxarray0_),
            _mxarray1_));
        /*
         * for repi = repvec
         */
        {
            mclForLoopIterator viter__;
            for (mclForStart(&viter__, mclVv(repvec, "repvec"), NULL, NULL);
                 mclForNext(&viter__, &repi);
                 ) {
                /*
                 * % get spikes from each cell
                 * mat1(:,1) = binidx2(:,repi);
                 */
                mclArrayAssign2(
                  &mat1,
                  mclArrayRef2(
                    mclVv(binidx2, "binidx2"),
                    mlfCreateColonIndex(),
                    mclVv(repi, "repi")),
                  mlfCreateColonIndex(),
                  _mxarray0_);
                /*
                 * mat2(2,:) = binidx2(:,repvec2(repi))';
                 */
                mclArrayAssign2(
                  &mat2,
                  mlfCtranspose(
                    mclArrayRef2(
                      mclVv(binidx2, "binidx2"),
                      mlfCreateColonIndex(),
                      mclArrayRef1(
                        mclVv(repvec2, "repvec2"), mclVv(repi, "repi")))),
                  _mxarray1_,
                  mlfCreateColonIndex());
                /*
                 * % compute xcorr and take the absolute value so we can more easily
                 * % search for differences that are within NumCentralBins
                 * xc = abs(mat1 * mat2);
                 */
                mlfAssign(
                  &xc,
                  mlfAbs(mclMtimes(mclVv(mat1, "mat1"), mclVv(mat2, "mat2"))));
                /*
                 * % find the spikes that are smaller than stdiff
                 * [m1,m2] = find(xc<stdiff);
                 */
                mlfAssign(
                  &m1,
                  mlfFind(
                    &m2,
                    NULL,
                    mclLt(mclVv(xc, "xc"), mclVv(stdiff, "stdiff"))));
                /*
                 * markmat1(m1,repi) = 1;
                 */
                mclArrayAssign2(
                  &markmat1, _mxarray0_, mclVv(m1, "m1"), mclVv(repi, "repi"));
                /*
                 * markmat2(m2,repi) = 1;
                 */
                mclArrayAssign2(
                  &markmat2, _mxarray0_, mclVv(m2, "m2"), mclVv(repi, "repi"));
                /*
                 * % get shift predictor spikes
                 * mat2(2,:) = binidx2(:,sprepvec2(repi))';
                 */
                mclArrayAssign2(
                  &mat2,
                  mlfCtranspose(
                    mclArrayRef2(
                      mclVv(binidx2, "binidx2"),
                      mlfCreateColonIndex(),
                      mclArrayRef1(
                        mclVv(sprepvec2, "sprepvec2"), mclVv(repi, "repi")))),
                  _mxarray1_,
                  mlfCreateColonIndex());
                /*
                 * xc = abs(mat1 * mat2);
                 */
                mlfAssign(
                  &xc,
                  mlfAbs(mclMtimes(mclVv(mat1, "mat1"), mclVv(mat2, "mat2"))));
                /*
                 * % find the spikes that are smaller than stdiff
                 * [m1,m2] = find(xc<stdiff);
                 */
                mlfAssign(
                  &m1,
                  mlfFind(
                    &m2,
                    NULL,
                    mclLt(mclVv(xc, "xc"), mclVv(stdiff, "stdiff"))));
                /*
                 * spmarkmat1(m1,repi) = 1;
                 */
                mclArrayAssign2(
                  &spmarkmat1,
                  _mxarray0_,
                  mclVv(m1, "m1"),
                  mclVv(repi, "repi"));
                /*
                 * spmarkmat2(m2,sprepvec(repi)) = 1;
                 */
                mclArrayAssign2(
                  &spmarkmat2,
                  _mxarray0_,
                  mclVv(m2, "m2"),
                  mclArrayRef1(
                    mclVv(sprepvec, "sprepvec"), mclVv(repi, "repi")));
            /*
             * end
             */
            }
            mclDestroyForLoopIterator(viter__);
        }
        /*
         * % convert to sparse matrices to save space
         * markmat1 = sparse(markmat1);
         */
        mlfAssign(
          &markmat1,
          mlfSparse(mclVv(markmat1, "markmat1"), NULL, NULL, NULL, NULL, NULL));
        /*
         * spmarkmat1 = sparse(spmarkmat1);
         */
        mlfAssign(
          &spmarkmat1,
          mlfSparse(
            mclVv(spmarkmat1, "spmarkmat1"), NULL, NULL, NULL, NULL, NULL));
        /*
         * markmat2 = sparse(markmat2);
         */
        mlfAssign(
          &markmat2,
          mlfSparse(mclVv(markmat2, "markmat2"), NULL, NULL, NULL, NULL, NULL));
        /*
         * spmarkmat2 = sparse(spmarkmat2);
         */
        mlfAssign(
          &spmarkmat2,
          mlfSparse(
            mclVv(spmarkmat2, "spmarkmat2"), NULL, NULL, NULL, NULL, NULL));
        /*
         * warning on MATLAB:divideByZero
         */
        mclPrintAns(&ans, mlfNWarning(0, NULL, _mxarray30_, _mxarray14_, NULL));
        /*
         * 
         * % save data
         * save(resultfile,'data','shiftp','surrmean','surrstd','surrhist', ...
         * 'dzscore','dxcorr','dxcmean','dxcstd','dxcz','hpval', ...
         * 'markmat1','spmarkmat1','markmat2','spmarkmat2','psthbins', ...
         * 'psths','jbmat','limat','dspsthz');
         */
        {
            mxArray * name_ = mclInitialize(mclVv(resultfile, "resultfile"));
            mlfSave(
              name_,
              "w",
              "data",
              data,
              "shiftp",
              shiftp,
              "surrmean",
              *surrmean,
              "surrstd",
              *surrstd,
              "surrhist",
              *surrhist,
              "dzscore",
              dzscore,
              "dxcorr",
              dxcorr,
              "dxcmean",
              dxcmean,
              "dxcstd",
              dxcstd,
              "dxcz",
              dxcz,
              "hpval",
              hpval,
              "markmat1",
              markmat1,
              "spmarkmat1",
              spmarkmat1,
              "markmat2",
              markmat2,
              "spmarkmat2",
              spmarkmat2,
              NULL);
            mlfSave(
              name_,
              "u",
              "psthbins",
              psthbins,
              "psths",
              psths,
              "jbmat",
              jbmat,
              "limat",
              limat,
              "dspsthz",
              dspsthz,
              NULL);
            mxDestroyArray(name_);
        }
    /*
     * end % if(surrnum<numwindows)
     */
    }
    mclValidateOutput(
      data, 1, nargout_, "data", "@shufflesync/private/computeSurrData");
    mclValidateOutput(
      *surrmean,
      2,
      nargout_,
      "surrmean",
      "@shufflesync/private/computeSurrData");
    mclValidateOutput(
      *surrstd, 3, nargout_, "surrstd", "@shufflesync/private/computeSurrData");
    mclValidateOutput(
      *surrhist,
      4,
      nargout_,
      "surrhist",
      "@shufflesync/private/computeSurrData");
    mclValidateOutput(
      *zscore, 5, nargout_, "zscore", "@shufflesync/private/computeSurrData");
    mxDestroyArray(DATA_COL);
    mxDestroyArray(SHIFT_COL);
    mxDestroyArray(SURR_COL);
    mxDestroyArray(KSNSPIKES);
    mxDestroyArray(ans);
    mxDestroyArray(cellpairdatafile);
    mxDestroyArray(surrprefix);
    mxDestroyArray(surrlist);
    mxDestroyArray(surrnum);
    mxDestroyArray(numwindows);
    mxDestroyArray(shiftp);
    mxDestroyArray(numSurrogates);
    mxDestroyArray(sdata);
    mxDestroyArray(lagslength);
    mxDestroyArray(sumvec);
    mxDestroyArray(smstmp);
    mxDestroyArray(llcenter);
    mxDestroyArray(NumCentralBins);
    mxDestroyArray(centerbins);
    mxDestroyArray(numsurrvec);
    mxDestroyArray(srind);
    mxDestroyArray(srind2);
    mxDestroyArray(dxcorr);
    mxDestroyArray(dxcz);
    mxDestroyArray(dxcmean);
    mxDestroyArray(dxcstd);
    mxDestroyArray(hpval);
    mxDestroyArray(lagsvec);
    mxDestroyArray(jbmat);
    mxDestroyArray(limat);
    mxDestroyArray(shuffle);
    mxDestroyArray(sadd);
    mxDestroyArray(subbins);
    mxDestroyArray(bintimes1);
    mxDestroyArray(bt1l);
    mxDestroyArray(bsize);
    mxDestroyArray(extrabins);
    mxDestroyArray(bintimes);
    mxDestroyArray(windowlength);
    mxDestroyArray(dvec);
    mxDestroyArray(pbins);
    mxDestroyArray(psthbins);
    mxDestroyArray(psths);
    mxDestroyArray(ns1);
    mxDestroyArray(ns2);
    mxDestroyArray(psthwin);
    mxDestroyArray(ws1);
    mxDestroyArray(halfwinsize);
    mxDestroyArray(lcentralwin);
    mxDestroyArray(lpsthbins);
    mxDestroyArray(ps1ind);
    mxDestroyArray(ps1ind2);
    mxDestroyArray(ps2ind);
    mxDestroyArray(ps2ind2);
    mxDestroyArray(cwinidx);
    mxDestroyArray(cwinidx2);
    mxDestroyArray(didx);
    mxDestroyArray(surridx1);
    mxDestroyArray(surridx2);
    mxDestroyArray(smoothmtx);
    mxDestroyArray(defaultspt);
    mxDestroyArray(dspsthz);
    mxDestroyArray(sidx);
    mxDestroyArray(l);
    mxDestroyArray(nsyncspikes);
    mxDestroyArray(xcdata);
    mxDestroyArray(xcsdata);
    mxDestroyArray(smean);
    mxDestroyArray(sstd);
    mxDestroyArray(lagi);
    mxDestroyArray(testdata);
    mxDestroyArray(dummy);
    mxDestroyArray(spt);
    mxDestroyArray(l0sp1);
    mxDestroyArray(l0sp2);
    mxDestroyArray(l0);
    mxDestroyArray(l0spt1);
    mxDestroyArray(l0spt2);
    mxDestroyArray(l1sp1);
    mxDestroyArray(l1sp2);
    mxDestroyArray(l1);
    mxDestroyArray(l1spt1);
    mxDestroyArray(l1spt2);
    mxDestroyArray(smat);
    mxDestroyArray(startbins);
    mxDestroyArray(endbins);
    mxDestroyArray(spsthbins);
    mxDestroyArray(psthhistcbins);
    mxDestroyArray(psth);
    mxDestroyArray(psth1);
    mxDestroyArray(psth2);
    mxDestroyArray(spsth);
    mxDestroyArray(dspsth1);
    mxDestroyArray(dspsth2);
    mxDestroyArray(sspsthmat);
    mxDestroyArray(sspm);
    mxDestroyArray(sspstd);
    mxDestroyArray(loopbins);
    mxDestroyArray(cidx);
    mxDestroyArray(cellvec);
    mxDestroyArray(sptrains);
    mxDestroyArray(sptimes);
    mxDestroyArray(surrtimes);
    mxDestroyArray(srtimes);
    mxDestroyArray(loopcounts);
    mxDestroyArray(surrvec1);
    mxDestroyArray(maxn);
    mxDestroyArray(SurrHistMin);
    mxDestroyArray(SurrHistStep);
    mxDestroyArray(bins);
    mxDestroyArray(si);
    mxDestroyArray(dzscore);
    mxDestroyArray(reps);
    mxDestroyArray(repvec);
    mxDestroyArray(sprepvec);
    mxDestroyArray(repvec2);
    mxDestroyArray(sprepvec2);
    mxDestroyArray(binidxsize);
    mxDestroyArray(nspikes);
    mxDestroyArray(mat1);
    mxDestroyArray(mat2);
    mxDestroyArray(markmat1);
    mxDestroyArray(markmat2);
    mxDestroyArray(spmarkmat1);
    mxDestroyArray(spmarkmat2);
    mxDestroyArray(binidx);
    mxDestroyArray(binidx2);
    mxDestroyArray(stdiff);
    mxDestroyArray(repi);
    mxDestroyArray(xc);
    mxDestroyArray(m1);
    mxDestroyArray(m2);
    mxDestroyArray(resultfile);
    mxDestroyArray(sdatafile);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return data;
}
