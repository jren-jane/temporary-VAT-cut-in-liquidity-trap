
// *******************************************************************
//
// fixed/calibrated parameters: SW / CCWT parameters
//
// *******************************************************************

// tax rates

    ctau_l=0.28;
    ctau_k=0.36;
    ctau_v=0.15;

    cphi=0.25; // share of RoT consumers
    consteb=0; // demeand data
    cby = 2.52; // government bonds-GDP ratio

    cdelta=0.0145; //depreciation rate
    clandaw=150/100 ; // SS markup labor market
    cg=0.193; //exogenous spending GDP-ratio
    curvp=10; //curvature Kimball aggregator goods market
    curvw=10; //curvature Kimball aggregator labor market


// *******************************************************************
//
// fixed/calibrated parameters: ALLV parameters
//
// *******************************************************************

  // foreign VAT
    ctau_vmc = 0.1;
    ctau_vmi = 0.1;
    ctau_vx = 0.1;

// investment and consumption share
    comega_i=0.55; //imported investment share
    comega_c=0.31; //imported consumption share

    ceta_c=5.00;  //substitution elasticity
    crhopi=0.975; //inflation target persistence in the Phillips curves

// *******************************************************************
//
// estimated parameters
//
// *******************************************************************

// measurement equations

  ctrend=0.8000; //quarterly trend growth rate to GDP
  constebeta=0.1968;
  constepinf=0.5617; //quarterly SS inflation rate
  constelab=2.0593;
  constex=0.5000;
  consterspread=0.3874;

// various DU and SW parameters

  calfa=0.0100; //labor share in production
  csigma=0.2500;//intertemporal elasticity of substitution
  cfc=1.2412; // parameter in the production function
  cgy=0.4462; // parameter in the government spending process
  csadjcost= 5.3476; //investment adjustment cost
  chabb=    0.4767;  // habit persistence 
  cprobw=   0.6829;  //calvo parameter labor market
  csigl=    0.2500; //labor supply elasticity
  cprobp=   0.5000; //calvo parameter goods market
  cindw=    0.2804; //indexation labor market
  cindp=    0.2276; //indexation goods market
  czcap=    0.4723;//capital utilization

// monetary policy

  crpi=     2.0443; //Taylor rule reaction to inflation
  crr=      0.8103;//Taylor rule interest rate smoothing
  cry=      0.0882;//Taylor rule long run reaction to output gap
  crdy=     0.2247;//Taylor rule short run reaction to output gap

//shock parameters

  crhoa=    0.9870;
  crhob=    0.9976;
  crhog=    0.8634;
  crhoqs=   0.1414;
  crhoms=   0.2801; 
  crhospread=0.8809;
  crhopinf= 0.8246; //price markup persistence
  crhow=    0.3338; //wage markup persistence
  cmap =    0.7500; //MA process
  cmaw  =   0.4188; //MA process


//substitution elasticity

  ceta_i=0.0403; //substitution elasticity investment
  ceta_f=3.0000; //substitution elasticity foreign

//Calvo parameter

  cprob_mc=0.4937;
  cprob_mi=0.5013;
  cprob_x=0.4537;
  cprob_e=0.7870;

//indexation parameter

  cind_mc=0.2287;
  cind_mi=0.5023;
  cind_x=0.2971;

//markups

  clambda_mc=1.7955; // steady states, determineds ceta_mc and ceta_mi, and then cgamma_mc and cgamma_mi
  clambda_mi=0.0584;

  crhopinf_mc=0.8749; // for the MA processes
  crhopinf_mi=0.5077;
  crhopinf_x=0.8688;

  cmap_mc=0.8497; // for the MA processes
  cmap_mi=0.5024; 
  cmap_x=0.7439;  

// monetary policy parameters

  crr  =0.7870;
  crpi =1.0300;
  crdpi=0.1310;
  crx  =0.2011;
  cry  =0.0743;
  crdy =0.0910;

// foreign bond premium

  cphi_a=0.0385;
  crhophi_tilde=0.8996;


// *******************************************************
//	         PARAMETER INPUT FOR FOREIGN VAR
// *******************************************************

  ForLag111 = 0.3768043;
  ForLag211 = 0.0736594;
  ForLag311 = 0.1843169;
  ForLag411 = 0.251153;

  ForLag112 = 0.113258;
  ForLag212 = 0.0181732;
  ForLag312 = 0.0366196;
  ForLag412 = 0.0001788 ;

  ForLag113 = -.1090663;
  ForLag213 = 0.2314199;
  ForLag313 = -0.2251971;
  ForLag413 = 0.0657894;

  ForLag121 = -0.4791982 ;
  ForLag221 = 0.1209666 ;
  ForLag321 = 0.2539076 ;
  ForLag421 = 0.3217595;

  ForLag122 = 0.2549007;
  ForLag222 = 0.2787732;
  ForLag322 = 0.0510823;
  ForLag422 = -0.0295298;

  ForLag123 = 0.8270144;
  ForLag223 = -0.8922233;
  ForLag323 = 0.3042302;
  ForLag423 = -0.1867868;

  ForLag131 = -0.0138056;
  ForLag231 = 0.0226833;
  ForLag331 = 0.0251806;
  ForLag431 = -0.0168382;

  ForLag132 = 0.0127795;
  ForLag232 = -0.0019983;
  ForLag332 = 0.024528;
  ForLag432 = -0.0009625;

  ForLag133 = 1.378601;
  ForLag233 = -0.2744204;
  ForLag333 = -0.1252314;
  ForLag433 = -0.0228878;

  ForShock11=0.1899;
  ForShock12=0;
  ForShock13=0;
  ForShock21=0;
  ForShock22=0.5912;
  ForShock23=0;
  ForShock31=-0.0035;
  ForShock32=0.0351;
  ForShock33=0.0841;

