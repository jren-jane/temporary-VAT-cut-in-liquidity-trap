// *******************************************************************
//
// VARIABLE AND PARAMETER DEFINITION
//
// *******************************************************************

// I modified this code which was originally written by 
// Drautzburg and Uhlig (2015)

// *******************************************************************
//
// MODEL ENDOGENOUS VARIABLES
//
// *******************************************************************

var 

    // observation equations
    labobs robs pinfobs dy dc dinve dw rspreadobs 

    // observation equations added for the open economy               
    xobs dx dm                                     
                                                                   
    // technolgy shock, government bond wedge shock, exogenous government spending, investment technology shock, monetary policy shock, private bond wedge
    a b g qs ms rspread 

    // markup shocks
    ewma epinfma epinfma_mc epinfma_mi epinfma_x

// sticky price-wage economy variables                                                                    

    // capital utilization, return on capital, capital, private capital, consumption, investment, output
    zcap rk k pk c inve y 

    // labor, wage, federal funds rate, private capital, bond holdings, interest rate by Taylor rule, real exchange rate, foreign asset
    lab w r kp bonds rtaylor x as               

    // marginal cost
    mc mc_mc mc_mi mc_x

    // price and wage markup 
    spinf spinf_mc spinf_mi spinf_x sw

    // inflation
    pinf pinf_c pinf_mc pinf_mi pinf_x pinf_bar

    // relative prices
    gamma_mcd gamma_mid gamma_cd gamma_id gamma_xstar gamma_f 

    // change in terms of trade
    dS

    // UIP
    phi_tilde

    // foreign economy variables
    pinfstar ystar rstar

    // Rule of Thumb vs. RA endogenous variables
    cRoT cRA

    // Aggregate transfers, Agent-specific transfers                                                                
    sRoT sRA                                                                

// flexible price-wage economy variables

    // exogenous given processes such as government spending are the same

    // utilization rate, return on capital, capital, private capital, consumption, investment, output
    zcapf rkf kf pkf cf invef yf                                     

    // labor, wage, federal funds rate, shadow price of private capital, bonds, real exchange rate, foreign asset
    labf wf rrf kpf bondsf xf asf

    // markups
    mcf mc_mcf mc_mif mc_xf

    // markup shocks
    spinff spinf_mcf spinf_mif spinf_xf swf

    // inflation
    pinff pinf_cf pinf_mcf pinf_mif pinf_xf 

    // relative prices
    gamma_mcdf gamma_midf gamma_cdf gamma_idf gamma_xstarf gamma_ff

    // change in terms of trade
    dSf

   // foreign economy variables
    pinfstarf ystarf rrstarf

    // Rule of Thumb vs. RA endogenous variables
    cRoTf cRAf

    // Aggregate transfers, Agent-specific transfers                                                              
    sRoTf sRAf;                                                             

// *******************************************************************
//
// MODEL EXOGENOUS VARIABLES
//
// *******************************************************************

varexo ea eb eg eqs em erspread                                             // genuine iid shocks
    epinf ew epinf_mc epinf_mi epinf_x                                      // markup shocks
    epinf_bar                                                               // varying inflation target
    epinfstar eystar erstar                                                 // foreign VAR shocks
    ephitilde                                                               // foreign bond premium shock
	zerointerest zerointerest_gradual                                       // impulses to policy shocks
    dtau_v dtau_vmc dtau_vmi dtau_vx dtau_vf;                               // predetermined VAT shock


// *******************************************************************
//
// MODEL PARAMETERS
//
// *******************************************************************


parameters 

// Kimball curvature for wage and price
    curvw curvp 

// labor share in production, government spending share in output
    calfa cgy cg

// capital utilization, discount rate, investment adjustment cost, depreciation rate, inverse of the intertemporal elasticity of substitution
    czcap cbeta csadjcost cdelta csigma 

// habit persistence, parameter in the production function that equals price markup, 
    chabb cfc csigl clandaw 

// price markup, detrended discount rate, steady states for the interest rate, the inflation rate, return on capital, wage
    clandap cbetabar cr cpie crk cw 

// markup parameters, steady states
    clambda_mc clambda_mi   // pass on values to ceta_mc and ceta_mi

// markup parameters, MA processes
    cmaw cmap cmap_mc cmap_mi cmap_x

// growth rate, which differs from the relative prices
    cgamma

// monetary policy parameters
    crpi crdy cry crr crdpi crx

// estimated parameters
    ctrend consteb conster constelab constepinf constebeta constex

// steady states with various ratios
    cikbar cik clk cky ciy ccy crkky cwly cwlc cwhlc cystary ccmc cimi ccmim ccdc cidi

// Calvo parameter
    cprobw cprobp cprob_mc cprob_mi cprob_x cprob_e

// inflation indexation
    cindw cindp cind_mc cind_mi cind_x

// steady states of elasticity of substitution
    ceta_c ceta_i ceta_f
    ceta_mc ceta_mi

// steady states of relative prices
    cgamma_id cgamma_cd cgamma_dc cgamma_mcc cgamma_cmc cgamma_imi cgamma_mii

// share of consumption and investment
    comega_c comega_i

// UIP, foreign bond premium
   cphi_a

// tax rates
    ctau_k ctau_l ctau_v ctau_vmc ctau_vmi ctau_vx

// bonds to GDP, revenue to GDP
    cby crevy

// persistence parameters
    crhoa crhob crhog crhoqs crhoms   

// markup persistence
    crhopinf crhow crhopinf_mc crhopinf_mi crhopinf_x

// Phillips curve parameter, varying inflation target shock parameter
    crhopi

// bond premium
    crhophi_tilde

// Fraction Rule of Thumb consumers, steady state subsidies to GDP, consumption to GDP of different agents, transfer to GDP of different agents
    cphi csy ccRoTy ccRAy

// new shock processes
    consterspread crhospread

// Parameters from exogenous VAR
    ForLag111 ,ForLag112 ,ForLag113 ,ForLag121 ,ForLag122 ,ForLag123 ,ForLag131 ,ForLag132 ,ForLag133 ,
    ForLag211 ,ForLag212 ,ForLag213 ,ForLag221 ,ForLag222 ,ForLag223 ,ForLag231 ,ForLag232 ,ForLag233 ,
    ForLag311 ,ForLag312 ,ForLag313 ,ForLag321 ,ForLag322 ,ForLag323 ,ForLag331 ,ForLag332 ,ForLag333 ,
    ForLag411 ,ForLag412 ,ForLag413 ,ForLag421 ,ForLag422 ,ForLag423 ,ForLag431 ,ForLag432 ,ForLag433 ,
    ForShock11 ,ForShock12 ,ForShock13 ,ForShock21 ,ForShock22 ,ForShock23 ,ForShock31 ,ForShock32 ,ForShock33;


