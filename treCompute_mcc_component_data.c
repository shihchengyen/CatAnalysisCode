/*
 * MATLAB Compiler: 4.3 (R14SP3)
 * Date: Wed Jan 17 16:44:31 2007
 * Arguments: "-B" "macro_default" "-m" "-W" "main" "-T" "link:exe"
 * "treCompute" 
 */

#include "mclmcr.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_treCompute_session_key[] = {
        '5', '1', '0', '4', 'B', '3', '4', '3', '5', '3', '0', '3', 'E', 'D',
        '2', '5', 'A', 'E', '0', '9', '4', '5', 'E', 'A', '3', '3', '5', '2',
        '4', '7', '7', '0', '2', 'A', '7', '4', '4', '7', '9', '6', 'E', '3',
        '9', '4', '1', '6', '2', '1', '6', '8', 'F', 'F', 'D', '6', '3', '4',
        '9', '8', '9', 'A', '0', '7', '8', 'B', '4', '7', '5', '1', '6', '6',
        'D', '8', '3', '9', '6', '1', '2', '9', '3', '1', 'E', 'B', '9', '8',
        '6', '1', '5', '9', '3', 'B', 'C', 'D', 'B', 'D', '8', '1', '8', 'C',
        '4', '1', '2', 'B', '1', '1', '9', '0', 'B', '0', '3', 'A', 'F', '7',
        'D', '2', '3', '7', '7', '6', '3', '5', '6', 'A', '9', 'F', '1', 'C',
        '4', 'F', 'F', 'B', 'F', '1', '6', 'B', '8', 'F', '9', '4', '0', '5',
        '5', 'A', 'C', 'D', '4', '6', '7', '8', 'C', 'E', '8', '0', '1', '1',
        '3', '2', 'B', '1', '4', '0', '2', '1', '4', '8', 'B', '7', '7', '1',
        '7', 'E', 'A', 'D', '2', 'D', 'A', '8', 'E', '6', '8', '0', '6', '6',
        '9', 'E', 'C', 'A', '3', '6', 'B', '7', '7', '4', '9', '1', 'B', 'D',
        '8', '5', '8', 'D', '8', '2', '8', '8', '2', '7', 'D', '8', '2', 'E',
        '0', '0', 'F', 'B', '3', '8', 'B', '5', '8', '0', 'E', '0', '5', 'A',
        'E', 'F', '4', '6', '0', '1', 'D', 'E', 'D', '6', '6', '5', 'F', 'A',
        'E', '3', '2', '2', '0', 'A', 'C', 'A', '2', '5', '9', '0', '9', '1',
        '3', '0', '3', '7', '\0'};

const unsigned char __MCC_treCompute_public_key[] = {
        '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9',
        '2', 'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1',
        '0', '1', '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B',
        '0', '0', '3', '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1',
        '0', '0', 'C', '4', '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3',
        'A', '5', '2', '0', '6', '5', '8', 'F', '6', 'F', '8', 'E', '0', '1',
        '3', '8', 'C', '4', '3', '1', '5', 'B', '4', '3', '1', '5', '2', '7',
        '7', 'E', 'D', '3', 'F', '7', 'D', 'A', 'E', '5', '3', '0', '9', '9',
        'D', 'B', '0', '8', 'E', 'E', '5', '8', '9', 'F', '8', '0', '4', 'D',
        '4', 'B', '9', '8', '1', '3', '2', '6', 'A', '5', '2', 'C', 'C', 'E',
        '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4', 'D', '0', '8', '5',
        'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2', 'E', 'D', 'E',
        '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6', '3', '7',
        '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E', '6',
        '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
        '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1',
        'B', 'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9',
        '9', '0', '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0',
        'B', '6', '1', 'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B',
        '5', '8', 'F', 'C', '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6',
        'E', 'B', '7', 'E', 'C', 'D', '3', '1', '7', '8', 'B', '5', '6', 'A',
        'B', '0', 'F', 'A', '0', '6', 'D', 'D', '6', '4', '9', '6', '7', 'C',
        'B', '1', '4', '9', 'E', '5', '0', '2', '0', '1', '1', '1', '\0'};

static const char * MCC_treCompute_matlabpath_data[] = 
    { "treCompute/", "toolbox/compiler/deploy/",
      "Users/syen/Documents/ShihCheng/Matlab/npt/scripts/miscellaneous/",
      "Users/syen/Documents/ShihCheng/Matlab/npt/scripts/view/",
      "Users/syen/Documents/ShihCheng/Matlab/MatlabCentral/",
      "Users/syen/Documents/ShihCheng/Matlab/fvt/view/",
      "Users/syen/Documents/ShihCheng/Matlab/tftb/",
      "$TOOLBOXMATLABDIR/general/", "$TOOLBOXMATLABDIR/ops/",
      "$TOOLBOXMATLABDIR/lang/", "$TOOLBOXMATLABDIR/elmat/",
      "$TOOLBOXMATLABDIR/elfun/", "$TOOLBOXMATLABDIR/specfun/",
      "$TOOLBOXMATLABDIR/matfun/", "$TOOLBOXMATLABDIR/datafun/",
      "$TOOLBOXMATLABDIR/polyfun/", "$TOOLBOXMATLABDIR/funfun/",
      "$TOOLBOXMATLABDIR/sparfun/", "$TOOLBOXMATLABDIR/scribe/",
      "$TOOLBOXMATLABDIR/graph2d/", "$TOOLBOXMATLABDIR/graph3d/",
      "$TOOLBOXMATLABDIR/specgraph/", "$TOOLBOXMATLABDIR/graphics/",
      "$TOOLBOXMATLABDIR/uitools/", "$TOOLBOXMATLABDIR/strfun/",
      "$TOOLBOXMATLABDIR/imagesci/", "$TOOLBOXMATLABDIR/iofun/",
      "$TOOLBOXMATLABDIR/audiovideo/", "$TOOLBOXMATLABDIR/timefun/",
      "$TOOLBOXMATLABDIR/datatypes/", "$TOOLBOXMATLABDIR/verctrl/",
      "$TOOLBOXMATLABDIR/codetools/", "$TOOLBOXMATLABDIR/helptools/",
      "$TOOLBOXMATLABDIR/demos/", "$TOOLBOXMATLABDIR/timeseries/",
      "$TOOLBOXMATLABDIR/hds/", "toolbox/local/", "toolbox/compiler/",
      "toolbox/optim/", "toolbox/signal/signal/", "toolbox/nnet/nnet/",
      "toolbox/nnet/nnutils/" };

static const char * MCC_treCompute_classpath_data[] = 
    { "" };

static const char * MCC_treCompute_libpath_data[] = 
    { "" };

static const char * MCC_treCompute_app_opts_data[] = 
    { "" };

static const char * MCC_treCompute_run_opts_data[] = 
    { "" };

static const char * MCC_treCompute_warning_state_data[] = 
    { "" };


mclComponentData __MCC_treCompute_component_data = { 

    /* Public key data */
    __MCC_treCompute_public_key,

    /* Component name */
    "treCompute",

    /* Component Root */
    "",

    /* Application key data */
    __MCC_treCompute_session_key,

    /* Component's MATLAB Path */
    MCC_treCompute_matlabpath_data,

    /* Number of directories in the MATLAB Path */
    42,

    /* Component's Java class path */
    MCC_treCompute_classpath_data,
    /* Number of directories in the Java class path */
    0,

    /* Component's load library path (for extra shared libraries) */
    MCC_treCompute_libpath_data,
    /* Number of directories in the load library path */
    0,

    /* MCR instance-specific runtime options */
    MCC_treCompute_app_opts_data,
    /* Number of MCR instance-specific runtime options */
    0,

    /* MCR global runtime options */
    MCC_treCompute_run_opts_data,
    /* Number of MCR global runtime options */
    0,
    
    /* Component preferences directory */
    "treCompute_C00EB8A89CC934B1D959DD322FA58D9D",

    /* MCR warning status data */
    MCC_treCompute_warning_state_data,
    /* Number of MCR warning status modifiers */
    0,

    /* Path to component - evaluated at runtime */
    NULL

};

#ifdef __cplusplus
}
#endif


