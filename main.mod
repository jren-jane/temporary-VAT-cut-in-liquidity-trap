// *******************************************************
//                  Main Dynare Module
// *******************************************************

// I modified this code which was originally written by 
// Drautzburg and Uhlig (2015)

// *******************************************************
//				LOAD PARAMETERS WHICH MAY CHANGE
//				
//				ALSO SETS OPTIONS FOR SIMULATION
// *******************************************************s

    @#include "main_config.mod"

// *******************************************************
//				DECLARING VARIABLES AND PARAMETERS
// *******************************************************

    @#include "var_par_declaration.mod"

// *******************************************************
//				PARAMETER INPUT FOR SW PARAMETERS
// *******************************************************

    @#include "parameters_setting.mod"

// parameters derived from steady state restrictions
    @#include "steadystate_restrictions.mod"

// *******************************************************
//				MODEL
// *******************************************************

    @#include "modelequations_ZLB.mod"

// *******************************************************
//				SHOCKS
// *******************************************************

    @#include "shocks.mod"

// *******************************************************
//				INITIAL CONDITIONS FOR OBS EQNS
// *******************************************************

    @#include "initial_conditions.mod"

// *******************************************************
//				THE PROPER SIMULATION
// *******************************************************

    simul(periods=@{SIMULATIONPERIOD},lmmcp);
    rplot pinf;
    rplot y;
    rplot pinf_mc;
    rplot pinf_mi;
    rplot r;
    rplot c;
    rplot inve;
    rplot k;
    rplot rk;
    rplot w;
    rplot lab;