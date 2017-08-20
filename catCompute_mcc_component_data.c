/*
 * MATLAB Compiler: 4.3 (R14SP3)
 * Date: Wed Jan 17 16:46:15 2007
 * Arguments: "-B" "macro_default" "-m" "-W" "main" "-T" "link:exe"
 * "catCompute" 
 */

#include "mclmcr.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_catCompute_session_key[] = {
        'A', '6', 'D', 'C', '0', '6', '5', 'A', 'B', '2', 'D', '6', '8', '5',
        'B', '0', '6', '9', '6', '6', 'D', 'E', '8', '0', '5', 'D', '3', '7',
        '3', 'F', '2', 'C', '6', '4', '1', '6', 'E', 'F', '8', 'D', '8', 'B',
        '4', '6', 'F', '0', '6', '0', '4', '4', '8', '6', '9', 'D', 'C', 'E',
        '5', 'A', '1', 'C', 'D', 'F', 'B', '5', '3', 'C', '4', 'F', '2', '4',
        '1', 'D', 'F', 'B', '3', '4', '7', 'C', '9', 'D', 'F', '5', 'F', '9',
        '2', '6', '6', 'E', '1', '8', '3', '5', '9', '8', '8', '1', '0', '2',
        'C', '4', 'C', '2', '1', 'D', 'B', '3', 'B', '7', '6', '3', '1', '7',
        '5', '0', '7', '1', 'C', 'C', 'C', 'F', '3', '5', '8', 'C', '4', '3',
        'D', '7', 'A', 'A', 'B', 'E', 'F', '8', 'A', '5', '8', 'C', '7', '6',
        '9', '6', '3', '4', 'D', 'B', 'E', '4', '2', 'B', '2', '4', '9', 'C',
        'E', 'D', 'A', '3', '4', '5', '1', '3', '9', '9', '5', '9', '1', 'A',
        '0', '9', '4', '0', '9', 'A', '5', '2', 'D', '6', '6', 'D', '0', 'D',
        '5', '9', '4', 'A', 'F', '4', '5', '5', '1', '2', '3', '1', '7', '9',
        '3', '2', '4', 'C', 'D', 'D', '9', '9', '9', 'B', 'B', '3', '5', '7',
        'A', '4', 'B', '3', '9', '9', 'B', '6', 'D', '6', '8', '4', '6', 'E',
        'B', '4', 'E', '8', '5', '8', 'A', '8', '3', 'A', '8', '4', '7', '3',
        'B', 'A', '6', '6', 'B', '4', '7', '4', 'B', 'A', '1', '1', '2', '6',
        '1', '3', '9', 'F', '\0'};

const unsigned char __MCC_catCompute_public_key[] = {
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

static const char * MCC_catCompute_matlabpath_data[] = 
    { "catCompute/", "toolbox/compiler/deploy/",
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

static const char * MCC_catCompute_classpath_data[] = 
    { "" };

static const char * MCC_catCompute_libpath_data[] = 
    { "" };

static const char * MCC_catCompute_app_opts_data[] = 
    { "" };

static const char * MCC_catCompute_run_opts_data[] = 
    { "" };

static const char * MCC_catCompute_warning_state_data[] = 
    { "" };


mclComponentData __MCC_catCompute_component_data = { 

    /* Public key data */
    __MCC_catCompute_public_key,

    /* Component name */
    "catCompute",

    /* Component Root */
    "",

    /* Application key data */
    __MCC_catCompute_session_key,

    /* Component's MATLAB Path */
    MCC_catCompute_matlabpath_data,

    /* Number of directories in the MATLAB Path */
    42,

    /* Component's Java class path */
    MCC_catCompute_classpath_data,
    /* Number of directories in the Java class path */
    0,

    /* Component's load library path (for extra shared libraries) */
    MCC_catCompute_libpath_data,
    /* Number of directories in the load library path */
    0,

    /* MCR instance-specific runtime options */
    MCC_catCompute_app_opts_data,
    /* Number of MCR instance-specific runtime options */
    0,

    /* MCR global runtime options */
    MCC_catCompute_run_opts_data,
    /* Number of MCR global runtime options */
    0,
    
    /* Component preferences directory */
    "catCompute_1D6DB158874DE8A9E7A1BE0128736A05",

    /* MCR warning status data */
    MCC_catCompute_warning_state_data,
    /* Number of MCR warning status modifiers */
    0,

    /* Path to component - evaluated at runtime */
    NULL

};

#ifdef __cplusplus
}
#endif


