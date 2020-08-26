// *******************************************************
//                  Main Dynare Module
// *******************************************************


    @#define betamean="0.30"
    @#define betastd="0.2"

// this corresponds to beta(1.5, 4.5) in standard notation



// *******************************************************************
//
// MODEL ENDOGENOUS VARIABLES
//
// *******************************************************************

var 
    
    // observation equations
    labobs robs pinfobs dy dc dinve dw 

//rspreadobs        
    
    // observation equations added for the open economy               
    //xobs dx dm                                                              

    // technolgy shock, government bond wedge shock, exogenous government spending, investment technology shock, monetary policy shock, private bond spread
    a b g qs ms rspread 

    // markup shocks
    ewma epinfma epinfma_mc epinfma_mi epinfma_x

// sticky price-wage economy variables                                                                    

    // capital utilization, return on capital, capital, private capital, consumption, investment, output, transfers
    zcap rk k pk c inve y s

    // labor, wage, federal funds rate, private capital, bond holdings, interest rate by Taylor rule, real exchange rate, foreign asset
    lab w r kp bonds x as             
// rtaylor  

    // marginal cost
    mc mc_mc mc_mi mc_x

    // price and wage markup 
    spinf spinf_mc spinf_mi spinf_x sw

    // inflation
    pinf pinf_c pinf_mc pinf_mi pinf_x 
// pinf_bar

    // relative prices
    gamma_mcd gamma_mid gamma_cd gamma_id gamma_xstar gamma_f 

    // change in terms of trade
    dS

    // UIP
    phi_tilde

    // foreign economy variables
    pinfstar ystar rstar
                                                        


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

    // transfers                                                             
    sf;                                                             


// *******************************************************************
//
// MODEL EXOGENOUS VARIABLES
//
// *******************************************************************

varexo ea eb eg eqs em erspread                                             // genuine iid shocks
    epinf ew epinf_mc epinf_mi epinf_x                                      // markup shocks
//    epinf_bar                                                               // varying inflation target
    ephitilde;                                                              // foreign bond premium shock

 
// *******************************************************************
//
// MODEL PARAMETERS
//
// *******************************************************************

parameters 

// Kimball curvature for wage and price
    curvw curvp 

// labor share in production, parameter in the government spending process, government spending share in output
    calfa cgy cg

// capital utilization, investment adjustment cost, depreciation rate, inverse of the intertemporal elasticity of substitution
    czcap csadjcost cdelta csigma 

// habit persistence, parameter in the production function that equals price markup, labor supply elasticity, wage markup
    chabb cfc csigl clandaw 

// markup parameters, steady states
    clambda_mc clambda_mi                                                   // pass on values to ceta_mc and ceta_mi

// markup parameters, MA processes
    cmaw cmap cmap_mc cmap_mi cmap_x

// monetary policy parameters
    crpi crdy cry crr crdpi crx

// estimated parameters
    ctrend constelab constepinf constebeta 
// constex

// Calvo parameter
    cprobw cprobp cprob_mc cprob_mi cprob_x cprob_e

// inflation indexation
    cindw cindp cind_mc cind_mi cind_x

// steady states of elasticity of substitution
    ceta_c ceta_i ceta_f 

// share of consumption and investment
    comega_c comega_i

// UIP, foreign bond premium
   cphi_a

// tax rates
    ctau_k ctau_l ctau_v ctau_vmc ctau_vmi ctau_vx

// bonds to GDP
    cby 

// persistence parameters
    crhoa crhob crhog crhoqs crhoms   

// markup persistence
    crhopinf crhow crhopinf_mc crhopinf_mi crhopinf_x

// Phillips curve parameter, varying inflation target shock parameter
    crhopi

// bond premium
    crhophi_tilde

// new shock processes
    consterspread crhospread

// Parameters from exogenous VAR
    ForLag111 ,ForLag112 ,ForLag113 ,ForLag121 ,ForLag122 ,ForLag123 ,ForLag131 ,ForLag132 ,ForLag133 ,
    ForLag211 ,ForLag212 ,ForLag213 ,ForLag221 ,ForLag222 ,ForLag223 ,ForLag231 ,ForLag232 ,ForLag233 ,
    ForLag311 ,ForLag312 ,ForLag313 ,ForLag321 ,ForLag322 ,ForLag323 ,ForLag331 ,ForLag332 ,ForLag333 ,
    ForLag411 ,ForLag412 ,ForLag413 ,ForLag421 ,ForLag422 ,ForLag423 ,ForLag431 ,ForLag432 ,ForLag433 ,
    ForShock11 ,ForShock12 ,ForShock13 ,ForShock21 ,ForShock22 ,ForShock23 ,ForShock31 ,ForShock32 ,ForShock33;


// *******************************************************
//				PARAMETER INPUT FOR SW PARAMETERS
// *******************************************************
   
    cby = 2.52;

    ctau_l=0.2;  // UK base rate
    ctau_k=0.28; // UK base rate
    ctau_v=0.175;// UK rate before changes

    crhospread=0.9;

    consterspread=0;
 //   constex=0.5;

// *******************************************************************
//
// SW / CCWT parameters
//
// *******************************************************************

// fixed parameters

// fixed parameters
  cdelta  =0.025; //depreciation rate
  clandaw=1.500 ; // SS markup labor market
  cg=0.193; //exogenous spending GDP-ratio
  curvp=10; //curvature Kimball aggregator goods market
  curvw=10; //curvature Kimball aggregator labor market

// estimated parameters initialisation

  ctrend=0.8000; //quarterly trend growth rate to GDP
  constebeta=0.1657;
  constepinf=0.5000; //quarterly SS inflation rate
  constelab=0;
// constelab=0.5509;

  calfa=0.2419; //labor share in production
  csigma=0.8381;//intertemporal elasticity of substitution
  cfc=1.6797; // parameter in the production function
  cgy=0.5187; // parameter in the government spending process
  csadjcost= 5.7606; //investment adjustment cost
  chabb=    0.7133;  // habit persistence 
  cprobw=   0.7061;  //calvo parameter labor market
  csigl=    1.8383; //labor supply elasticity
  cprobp=   0.6523; //calvo parameter goods market
  cindw=    0.5845; //indexation labor market
  cindp=    0.2432; //indexation goods market
  czcap=    0.5462;//capital utilization

//  crpi=     2.0443; //Taylor rule reaction to inflation
//  crr=      0.8103;//Taylor rule interest rate smoothing
//  cry=      0.0882;//Taylor rule long run reaction to output gap
//  crdy=     0.2247;//Taylor rule short run reaction to output gap

//shock parameters

  crhoa=    0.9577;
  crhob=    0.2194;
  crhog=    0.9767;
  crhoqs=   0.7113;
  crhoms=   0.1479; 
  crhopinf= 0.8895; //price markup persistence
  crhow=    0.9688; //wage markup persistence
  cmap =    0.7010;
  cmaw  =   0.8503;

// *******************************************************************
//
// ALLV parameters
//
// *******************************************************************

// fixed parameters
  comega_i=0.55; //imported investment share
  comega_c=0.31; //imported consumption share
  crhopi=0.975; //inflation target persistence in the Phillips curves
  ceta_c=5.00;  //substitution elasticity


// elasticity of substitution
  ceta_i=1.696; //substitution elasticity investment
  ceta_f=1.486; //substitution elasticity foreign

//Calvo parameter

  cprob_mc=0.444;
  cprob_mi=0.721;
  cprob_x=0.612;
  cprob_e=0.787;

//indexation parameter

  cind_mc=0.220;
  cind_mi=0.231;
  cind_x=0.185;

//markups

  clambda_mc=1.633; // steady states, determines ceta_mc and ceta_mi, and then cgamma_mc and cgamma_mi
  clambda_mi=1.275;

  cmap_mc=0.7; // for the MA processes
  cmap_mi=0.7; 
  cmap_x=0.7;  

  crhopinf_mc=0.970; // for the MA processes
  crhopinf_mi=0.963;
  crhopinf_x=0.886;

  cphi_a=0.252;

  crhophi_tilde=0.955;

// monetary policy parameters
  crr  =0.881;
  crpi =1.73;
  crdpi=0.310;
  crx  =-0.009;
  cry  =0.104;
  crdy =0.128;

// foreign VAT
  ctau_vmc = 0.1;
  ctau_vmi = 0.1;
  ctau_vx = 0.1;

// *******************************************************
//	         PARAMETER INPUT FOR FOREIGN VAR
// *******************************************************


// The foreign economy is described by a VAR with 4 lags

  ForLag111 = 0.3768043;
  ForLag211 = 0.0736594;
  ForLag311 = 0.1843169;
  ForLag411 = 0.251153;

  ForLag112 = 0.113258;
  ForLag212 = 0.0181732;
  ForLag312 = 0.0366196;
  ForLag412 = 0.0001788;

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


// *******************************************************
//				MODEL
// *******************************************************

// *************************************************
// Model equations for extended Smets-Wouters model
// *************************************************

model; 

// parameter dependencies from steady state conditions

// *******************************************************************
//
// ALLV parameters
//
// *******************************************************************

       #ceta_mc=clambda_mc/(clambda_mc-1); //Substitution elasticity among imported consumption goods
       #ceta_mi=clambda_mi/(clambda_mi-1);//Substitution elasticity among imported investment goods
       #cgamma_id=((1-comega_i)+comega_i*(ceta_mi/(ceta_mi-1))^(1-ceta_i))^(1/(1-ceta_i)); //A.11
       #cgamma_cd=((1-comega_c)+comega_c*(ceta_mc/(ceta_mc-1))^(1-ceta_c))^(1/(1-ceta_c)); //A.9
       #cgamma_dc=1/cgamma_cd;
       #cgamma_cmc=((1-comega_c)*((ceta_mc-1)/ceta_mc)^(1-ceta_c)+comega_c)^(1/(1-ceta_c)); //A.10
       #cgamma_mcc=1/cgamma_cmc;
       #cgamma_imi=((1-comega_i)*((ceta_mi-1)/ceta_mi)^(1-ceta_i)+comega_i)^(1/(1-ceta_i)); //A.12
       #cgamma_mii=1/cgamma_imi;
       #ccmc = comega_c*cgamma_cmc^ceta_c; //A.26
       #cimi = comega_i*cgamma_imi^ceta_i; //A.26
       #ccdc = (1-comega_c)*cgamma_cd^ceta_c; //A.25
       #cidi = (1-comega_i)*cgamma_id^ceta_i; //A.25
 
// *******************************************************************
//
// SW / CCWT parameters
//
// *******************************************************************

       #clandap=cfc; // B.18, price markup lambda
       #cpie=1+constepinf/100; 
       #cgamma=1+ctrend/100 ; 
       #cbeta=1/(1+constebeta/100);
       #cbetabar=cbeta*cgamma^(-csigma); // normalization definition
       #cr=cpie/(cbeta*cgamma^(-csigma)); 
       #crk=(cbeta^(-1)*(cgamma^csigma) -(1-cdelta) - cdelta*ctau_k) /(1-ctau_k); // B.33c
       #cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa)); // B.11
       #cikbar=(1-(1-cdelta)/cgamma); // B.30
       #cik=(1-(1-cdelta)/cgamma)*cgamma; // B.30
       #clk=((1-calfa)/calfa)*(crk/cw); // B.9
       #cky=cfc*(clk)^(calfa-1); // B.8
       #ciy=cik*cky;
       #ccy=(1-cg- ((1-comega_i)*cgamma_id^ceta_i + comega_i*(cgamma_imi)^ceta_i) *cik*cky) / ((1-comega_c)*cgamma_cd^ceta_c + comega_c*cgamma_cmc^ceta_c); // B.62 modified resource constraint
       #crkky=crk*cky;
       #cwlc=(1-calfa)/calfa*crk*cky/ccy;
       #cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
       #cwly=1-crk*cky; // labor and capital for intermediate firms, B.16, B.18
       #conster=(cr-1)*100;
       #csy=ctau_l*cwlc*ccy+ctau_k*(crk-cdelta)*cky+(ctau_v/(1+ctau_v))*ccy+((cpie*cgamma/cr)-1)*cby/(cpie*cgamma)-cg; // government budget
       #cystary = comega_c*cgamma_cmc^ceta_c*ccy + comega_i*cgamma_imi^ceta_i*ciy; // from ALLV
       #ccmim = (ccmc*ccy) / (cimi*ciy); // from ALLV


 // *************************************************
// FLEXIBLE PRICE-WAGE ECONOMY
// *************************************************


        // instead of monetary policy
          pinff = 0;

          mcf = 0;
          mc_mcf = 0;
          mc_mif = 0;
          mc_xf = 0;

        // B.65. Marginal cost is zero in the flexible price economy
	      mcf =  calfa*rkf+(1-calfa)*(wf) - 1*a  ;

        // B.74e. HH's FOC for capacity utilization
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;

        // B.64. Capital-labor ratio
	      rkf =  (wf)+labf-kf ;

        // B.69. Capital services
	      kf =  kpf(-1)+zcapf ;

        // B.74d. HH's FOC for private investment
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*(pkf-gamma_idf) ) +qs ;

        // B.74bc. HH's FOC for private capital. The shadow price for private capital
	      pkf = -rrf-b-rspread+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*0*b +(crk*(1-ctau_k)/(crk*(1-ctau_k)+cdelta*ctau_k+(1-cdelta)))*rkf(1) 
			+  ((1-cdelta)/(crk*(1-ctau_k)+cdelta*ctau_k+(1-cdelta)))*pkf(1) ;

        // B.73 and B.74. RA HH's FOC for consumption and bond holdings
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) 
			+((csigma-1)*(1-ctau_l)*(cwlc/clandaw)/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) 
			- (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+b -gamma_cdf) + 0*b ;
            //the zero premultiplying b refers to the scaling of the shock only. 

        // B.63. Production function
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );

        // B.71. Wage-setting
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) + 1*swf;

        // B.68. Law of motion for capital
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs*(1+cbetabar*cgamma) ;

        // imported consumption goods
          mc_mcf = -mc_xf-gamma_xstarf-gamma_mcdf;

        // imported investment goods
          mc_mif = -mc_xf-gamma_xstarf-gamma_midf;      

        // relative prices
          gamma_mcdf=gamma_mcdf(-1)+pinf_mcf-pinff;
          gamma_midf=gamma_midf(-1)+pinf_mif-pinff;
          gamma_xstarf=gamma_xstarf(-1)+pinf_xf-pinfstarf;
          mc_xf=mc_xf(-1)+pinff-pinf_xf-dSf;
          gamma_ff=mc_xf+gamma_xstarf;
          gamma_cdf=comega_c*(cgamma_mcc)^(1-ceta_c)*gamma_mcdf;
          gamma_idf=comega_i*(cgamma_mii)^(1-ceta_i)*gamma_midf;


        // A B.10 UIP condition
          dSf(+1)-(rrf+b-rrstarf)-cphi_a*asf+phi_tilde=0;

        // A B.11 aggregate resource constraint
          (1-comega_c)*cgamma_cd^ceta_c*ccy*(cf+ceta_c*gamma_cdf)
            +(1-comega_i)*cgamma_id^ceta_i*ciy*(invef+ceta_i*gamma_idf)
            +g+cystary*(ystarf-ceta_f*gamma_xstarf)
            =yf-1*crkky*zcapf;

        // A B.17 net foreign assets
           cgamma*cpie*asf=cgamma*cpie*cr*(-cystary*mc_xf-ceta_f*cystary*gamma_xstarf+cystary*ystarf+(ccmc*ccy+cimi*ciy)*gamma_ff
               -(ccmc*ccy*(-ceta_c*(1-comega_c)*(cgamma_cd)^(-(1-ceta_c))*gamma_mcdf +cf)
               +cimi*ciy*(-ceta_i*(1-comega_i)*(cgamma_id)^(-(1-ceta_i))*gamma_midf+invef)))
               +cr*asf(-1);

        // CPI_Inflation
          pinf_cf=((1-comega_c)*(cgamma_dc)^(1-ceta_c))*pinff+((comega_c)*(cgamma_mcc)^(1-ceta_c))*pinf_mcf;

        // real exchange rate
          xf=-comega_c*(cgamma_cmc)^(-(1-ceta_c))*gamma_mcdf-gamma_xstarf-mc_xf;

        // B.76. Government budget constraint
 //         (1/cgamma)*(cgamma*g+(cby/cpie)*(bondsf(-1)-0)
 //             -cgamma*ctau_k*cky*(crk*rkf+(crk-cdelta)* kpf(-1)) - ctau_l*cwlc*ccy*(wf+labf)*cgamma 
 //             -(ctau_v/(1+ctau_v))*ccy*cf*cgamma )+ csy*sf = 0;

       ctau_l*cwlc*ccy*(wf+labf)+ (cby/(cr))*(bondsf-rrf-b)
        =(1/cgamma)*(cgamma*(g+csy*sf)+(cby/cpie)*(bondsf(-1)-0)-(ctau_v/(1+ctau_v))*ccy*cf*cgamma
                    -cgamma*ctau_k*cky*(crk*rkf+(crk-cdelta)* kpf(-1)) );



// *************************************************
// STICKY PRICE-WAGE ECONOMY
// *************************************************


        // B.65. Marginal cost
	      mc =  calfa*rk+(1-calfa)*(w) - 1*a ;

        // B.74e. HH's FOC for capital utilization
	      zcap =  (1/(czcap/(1-czcap)))* rk ;

        // B.64. Capital-labor ratio
	      rk =  w+lab-k ;

        // B.69. Capital services
	      k =  kp(-1)+zcap ;

        // B.74d. HH's FOC for investment
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*(pk-gamma_id) ) +qs ;

        // B.74bc. HH's FOC for capital. Shadow price for private capital
	      pk = -r+pinf(1)-b-rspread+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*0*b +(crk*(1-ctau_k)/(crk*(1-ctau_k)+cdelta*ctau_k+(1-cdelta)))*rk(1) 
			+  ((1-cdelta)/(crk*(1-ctau_k)+cdelta*ctau_k+(1-cdelta)))*pk(1) ;

        // RA HH's FOC for consumption
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +
			((csigma-1)*(1-ctau_l)*(cwlc/clandaw)/(csigma*(1+chabb/cgamma)))*(lab-lab(+1))
			- (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) +b -gamma_cd) +0*b ;
            //the zero premultiplying b refers to the scaling of the shock only. 
            // notice that the SIGN of b is flipped (consinstently) compared to the original SW code: now it matches the write-up.

        // UIP condition
          dS(+1)-(b+r-rstar)-cphi_a*as+phi_tilde=0;

        // aggregate resource constraint
          (1-comega_c)*cgamma_cd^ceta_c*ccy*(c+ceta_c*gamma_cd)
            +(1-comega_i)*cgamma_id^ceta_i*ciy*(inve+ceta_i*gamma_id)
            +g+cystary*(ystar-ceta_f*gamma_xstar)
            =y-1*crkky*zcap;

        // net foreign assets
           cgamma*cpie*as=cgamma*cpie*cr*(-cystary*mc_x-ceta_f*cystary*gamma_xstar+cystary*ystar+(ccmc*ccy+cimi*ciy)*gamma_f
               -(ccmc*ccy*(-ceta_c*(1-comega_c)*(cgamma_cd)^(-(1-ceta_c))*gamma_mcd +c)
               +cimi*ciy*(-ceta_i*(1-comega_i)*(cgamma_id)^(-(1-ceta_i))*gamma_mid+inve)))
               +cr*as(-1);

        // B.63. Production function
	      y = cfc*( calfa*k+(1-calfa)*lab +a );

        // domestic Phillips curve
	     pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
           +((1-ctau_v)*(1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ;
   
         
        // Phillips curve for imported consumption goods
//	     pinf_mc-pinf_bar =  (1/(1+cbetabar*cgamma*cind_mc)) * ( cbetabar*cgamma*(pinf_mc(1)-crhopi*pinf_bar) +cind_mc*(pinf_mc(-1)-pinf_bar) -cbetabar*cgamma*cind_mc*(1-crhopi)*pinf_bar
//           +((1-cprob_mc)*(1-cbetabar*cgamma*cprob_mc)/cprob_mc)/((cfc-1)*curvp+1)*(mc_mc)  )  + spinf_mc;
          
	     pinf_mc =  (1/(1+cbetabar*cgamma*cind_mc)) * ( cbetabar*cgamma*pinf_mc(1) +cind_mc*pinf_mc(-1) 
           +((1-ctau_vmc)*(1-cprob_mc)*(1-cbetabar*cgamma*cprob_mc)/cprob_mc)/((cfc-1)*curvp+1)*(mc_mc)  )  + spinf_mc ;

        // marginal cost for imported consumption goods
         mc_mc=-mc_x-gamma_xstar-gamma_mcd;

        // Phillips curve for imported investment goods
//	     pinf_mi-pinf_bar =  (1/(1+cbetabar*cgamma*cind_mi)) * ( cbetabar*cgamma*(pinf_mi(1)-crhopi*pinf_bar) +cind_mi*(pinf_mi(-1)-pinf_bar) -cbetabar*cgamma*cind_mi*(1-crhopi)*pinf_bar
//           +((1-cprob_mi)*(1-cbetabar*cgamma*cprob_mi)/cprob_mi)/((cfc-1)*curvp+1)*(mc_mi)  )  + spinf_mi;

	     pinf_mi =  (1/(1+cbetabar*cgamma*cind_mi)) * ( cbetabar*cgamma*pinf_mi(1) +cind_mi*pinf_mi(-1) 
           +((1-ctau_vmi)*(1-cprob_mi)*(1-cbetabar*cgamma*cprob_mi)/cprob_mi)/((cfc-1)*curvp+1)*(mc_mi)  )  + spinf_mi ;

        // marginal cost for imported investment goods
         mc_mi=-mc_x-gamma_xstar-gamma_mid;

        // Phillips curve for exporting firms
//	     pinf_x-pinf_bar =  (1/(1+cbetabar*cgamma*cind_x)) * ( cbetabar*cgamma*(pinf_x(1)-crhopi*pinf_bar) +cind_x*(pinf_x(-1)-pinf_bar) -cbetabar*cgamma*cind_x*(1-crhopi)*pinf_bar
//           +((1-cprob_x)*(1-cbetabar*cgamma*cprob_x)/cprob_x)/((cfc-1)*curvp+1)*(mc_x)  )  + spinf_x;

	     pinf_x =  (1/(1+cbetabar*cgamma*cind_x)) * ( cbetabar*cgamma*pinf_x(1) +cind_x*pinf_x(-1) 
           +((1-ctau_vx)*(1-cprob_x)*(1-cbetabar*cgamma*cprob_x)/cprob_x)/((cfc-1)*curvp+1)*(mc_x)  )  + spinf_x ;

        // B.84 Taylor rule

//           r = rtaylor;

//            rtaylor = (1-crr)*pinf_bar
//                    +crpi*(1-crr)*(pinf_c(-1)-pinf_bar)
//                    +cry*(1-crr)*(y(-1)-yf(-1))
//                    +crx*(1-crr)*(x(-1)-xf(-1))
//                    +crdpi*(pinf_c-pinf_c(-1))
//                    +crdy*(y-yf-y(-1)+yf(-1))
//                    +crr*rtaylor(-1)
//                    +ms ;

         // B.84. Monetary policy--Taylor rule
          r=(crpi*(1-crr)*pinf+cry*(1-crr)*(y-yf)+crdy*(y-yf-y(-1)+yf(-1))+crr*r(-1))+crx*(1-crr)*(x(-1)-xf(-1))+ms ;

        // CPI_Inflation
          pinf_c=((1-comega_c)*(cgamma_dc)^(1-ceta_c))*pinf+((comega_c)*(cgamma_mcc)^(1-ceta_c))*pinf_mc;

        // real exchange rate
          x=-comega_c*(cgamma_cmc)^(-(1-ceta_c))*gamma_mcd-gamma_xstar-mc_x;

        // B.72 -- Union prices according to RA, but all agents get union wage and profit
        // wage determination
//		w =  (1/(1+cbetabar*cgamma))*w(-1)							
//               +((cbetabar*cgamma)/(1+cbetabar*cgamma))*w(1)			
//               +(cindw/(1+cbetabar*cgamma))*(pinf_c(-1)-pinf_bar)				
//               -(cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*(pinf_c-crhopi*pinf_bar)		
//               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*(pinf(1)-crhopi*pinf_bar)				
//               -1/(1+cbetabar*cgamma)*(pinf-pinf_bar)
//               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)			
//			*(1/((clandaw-1)*curvw+1))*									
//               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) 	
//			-w) 										                	
//               + 1*sw ;				

               // B.72 -- Union prices according to RA, but all agents get union wage and profit
		w =  (1/(1+cbetabar*cgamma))*w(-1)							//ok
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)					//ok
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)						//ok
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf			//ok
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)				//ok
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)			//ok
			*(1/((clandaw-1)*curvw+1))*									//what's this? markup. clandaw=1+lambda_w
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) 	//ok
			-w) 										                	//ok
               + 1*sw ;										            	//seems to be rescaled						             

        // relative prices
        gamma_mcd=gamma_mcd(-1)+pinf_mc-pinf;
        gamma_mid=gamma_mid(-1)+pinf_mi-pinf;
        gamma_xstar=gamma_xstar(-1)+pinf_x-pinfstar;
        mc_x=mc_x(-1)+pinf-pinf_x-dS;
        gamma_f=mc_x+gamma_xstar;
        gamma_cd=comega_c*(cgamma_mcc)^(1-ceta_c)*gamma_mcd;
        gamma_id=comega_i*(cgamma_mii)^(1-ceta_i)*gamma_mid;

        // B.76. Government budget constraint
 //       (1/cgamma)*(cgamma*g+(cby/cpie)*(bonds(-1)-pinf)
 //       -cgamma*ctau_k*cky*(crk*rk+(crk-cdelta)* kp(-1)) - ctau_l*cwlc*ccy*(w+lab)*cgamma
 //       -(ctau_v/(1+ctau_v))*ccy*c*cgamma)+ csy*s = 0;

      ctau_l*cwlc*ccy*(w+lab)+ (cby/(cr))*(bonds-r-b)=(1/cgamma)*(cgamma*(g+csy*s)+(cby/cpie)*(bonds(-1)-pinf)-(ctau_v/(1+ctau_v))*ccy*c*cgamma
                    -cgamma*ctau_k*cky*(crk*rk+(crk-cdelta)* kp(-1))  );

         // B.68. Law of motion for capital
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs*(1+cbetabar*cgamma) ;


// *************************************************
// FOREIGN VAR
// *************************************************

//Foreign VAR
   pinfstar  = ForLag111*pinfstar(-1)  + ForLag112*ystar(-1) + ForLag113*rstar(-1)
             +ForLag211*pinfstar(-2)  + ForLag212*ystar(-2) + ForLag213*rstar(-2)
             +ForLag311*pinfstar(-3)  + ForLag312*ystar(-3) + ForLag313*rstar(-3)
      		+ForLag411*pinfstar(-4)  + ForLag412*ystar(-4) + ForLag413*rstar(-4) ;


   ystar = ForLag121*pinfstar(-1)  + ForLag122*ystar(-1) + ForLag123*rstar(-1)
             +ForLag221*pinfstar(-2)  + ForLag222*ystar(-2) + ForLag223*rstar(-2)
             +ForLag321*pinfstar(-3)  + ForLag322*ystar(-3) + ForLag323*rstar(-3)
      		+ForLag421*pinfstar(-4)  + ForLag422*ystar(-4) + ForLag423*rstar(-4) ;



   rstar  = ForLag131*pinfstar(-1)  + ForLag132*ystar(-1) + ForLag133*rstar(-1)
             +ForLag231*pinfstar(-2)  + ForLag232*ystar(-2) + ForLag233*rstar(-2)
             +ForLag331*pinfstar(-3)  + ForLag332*ystar(-3) + ForLag333*rstar(-3)
      		+ForLag431*pinfstar(-4)  + ForLag432*ystar(-4) + ForLag433*rstar(-4) ;



//Foreign VAR for flexible economy (added)
   pinfstarf  = 0;

   ystarf = ForLag121*pinfstarf(-1)  + ForLag122*ystarf(-1) + ForLag123*rrstarf(-1)
             +ForLag221*pinfstarf(-2)  + ForLag222*ystarf(-2) + ForLag223*rrstarf(-2)
             +ForLag321*pinfstarf(-3)  + ForLag322*ystarf(-3) + ForLag323*rrstarf(-3)
      		+ForLag421*pinfstarf(-4)  + ForLag422*ystarf(-4) + ForLag423*rrstarf(-4) ;

   rrstarf  = ForLag131*pinfstarf(-1)  + ForLag132*ystarf(-1) + ForLag133*rrstarf(-1)
             +ForLag231*pinfstarf(-2)  + ForLag232*ystarf(-2) + ForLag233*rrstarf(-2)
             +ForLag331*pinfstarf(-3)  + ForLag332*ystarf(-3) + ForLag333*rrstarf(-3)
      		+ForLag431*pinfstarf(-4)  + ForLag432*ystarf(-4) + ForLag433*rrstarf(-4) ;



// *************************************************
// EXOGENOUS PROCESSES - STOCHASTIC
// *************************************************

	      a = crhoa*a(-1)  + (ea); // B.85a. technology
	      b = crhob*b(-1) + eb; // B.85i. government bond wedge
	      g = crhog*(g(-1)) + eg + cgy*(ea); // B.85cd. autonomous government spending
                                             // the last term describes net exports
	      qs = crhoqs*qs(-1) + eqs; // B.85k. investment technology/value of private investment
	      ms = crhoms*ms(-1) + em; // B.85b. monetary policy

        // The price and wage markups follow ARMA processes. 
        // The inclusion of the MA term is designed to capture the high-frequency fluctuations in inflation

	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1); 
	          epinfma=epinf;

 	      spinf_mc = crhopinf_mc*spinf_mc(-1) + epinfma_mc - cmap_mc*epinfma_mc(-1); 
	          epinfma_mc=epinf_mc;

	      spinf_mi = crhopinf_mi*spinf_mi(-1) + epinfma_mi - cmap_mi*epinfma_mi(-1); 
	          epinfma_mi=epinf_mi;

	      spinf_x = crhopinf_x*spinf_x(-1) + epinfma_x - cmap_x*epinfma_x(-1); 
	          epinfma_x=epinf_x; 

	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=ew; 

        // Markups are zero in the flexible price economy
          spinff = 0;
          spinf_mcf = 0;
          spinf_mif = 0;
          spinf_xf = 0;
          swf = 0;
       
//          pinf_bar=crhopi*pinf_bar(-1)+epinf_bar; // the moving inflation target

          rspread=crhospread*1*rspread(-1)+erspread; // B.85j. private bond spread

          phi_tilde=crhophi_tilde*phi_tilde(-1)+ephitilde; // foreign bond premium

        // Set debt change to zero. Interpretation: spending offset by changes in lump-sum taxes.
          bondsf = 0;
          bonds = 0;

// *************************************************
// MEASUREMENT EQUATIONS
// *************************************************

        dy=y-y(-1)+ctrend; // B.89a. output, scaled by 100
        dc=c-c(-1) + 1/((1/ccdc)+(1/ccmc))*ceta_c*(1/cgamma_cmc-1/cgamma_cd)*(gamma_mcd-gamma_mcd(-1)) + ctrend; // consumption growth
        dinve=inve-inve(-1) + 1/((1/cidi)+(1/cimi))*ceta_i*(1/cgamma_imi-1/cgamma_id)*(gamma_mid-gamma_mid(-1)) + ctrend; // investment growth
        dw=w-w(-1)+ctrend; // B.89e. wage, scaled by 100
        pinfobs = 1*(pinf) + constepinf; // B.89f. inflation
        robs =    1*(r) + conster; // B.89g. bank rate
 //       rspreadobs = rspread+consterspread; // B.89h. private bond wedge
        labobs = lab + constelab; // B.89i. labor

      //added for the open economy
//        xobs = x - constex; // real exchange rate
//        dx = (ystar-ceta_f*gamma_xstar) - (ystar(-1)-ceta_f*gamma_xstar(-1)) + ctrend; // real exports
//        dm = ccmim / (ccmim + 1) * ((c - c(-1)) - ceta_c*(1-comega_c)*(cgamma_cd)^(-(1-ceta_c)) * (gamma_mcd-gamma_mcd(-1)))
//            + 1 / (ccmim + 1) * ((inve - inve(-1)) - ceta_i*(1-comega_i)*(cgamma_id)^(-(1-ceta_i)) * (gamma_mid-gamma_mid(-1)))
//            + ctrend; // real imports

end; 

// *******************************************************
//				INITIAL CONDITIONS FOR OBS EQNS
// *******************************************************

initval;

//      robs = 0.8;
//      pinfobs = 0.5;
//      dy = 0.8;
//      dc = 0.8;
//      dinve = 0.8;
//      dw = 0.8;

dy=ctrend;
dc=ctrend;
dinve=ctrend;
dw=ctrend;
pinfobs = constepinf;
robs = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100;
labobs = constelab;

    // The rest do not need to be specified because the data is demeaned

    // added for the open economy
//      xobs = -0.5;
//      dx = 0.8;
//      dm = 0.8;
      
    
   end;

   check;
   
   steady; 



// *******************************************************
//				SHOCKS
// *******************************************************

    shocks;

   
        var ea;
        stderr 0.5730;
        var eb;
        stderr 2.1670;
        var eg;
        stderr 0.7724;
        var eqs;
        stderr 1.3669;
        var em;
        stderr 0.1273;
        var epinf;
        stderr 0.4315;
        var ew;
        stderr 0.5299;
        var erspread;
        stderr 0.2688;

        // added for the open economy
        //consumption import markup
        var epinf_mc;
        stderr 2.882;
        //investment import markup
        var epinf_mi;
        stderr 0.354;
        //export markup
        var epinf_x;
        stderr 1.124;
        //inflation target shock
//        var epinf_bar;
//        stderr 0.053;
        //risk premium shock
        var ephitilde;
        stderr 0.183;


    end;

stoch_simul(irf=20,nograph);

close all;

        estimated_params;
        // PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
        // PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
            stderr ea,0.4565,0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr eb,0.9871,0.025,5,INV_GAMMA_PDF,0.1,2;
            stderr eg,0.4058,0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr eqs,1.1966,0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr em,0.2184,0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr epinf,0.2096,0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr ew,0.2558,0.01,3,INV_GAMMA_PDF,0.1,2;
//            stderr erspread,0.0775,0.025,5,INV_GAMMA_PDF,0.1,2; //0.195 is the entire sample mean

        // shocks added for the open economy
            stderr epinf_mc, 2.882, 0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr epinf_mi, 0.354, 0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr epinf_x, 1.124, 0.01,3,INV_GAMMA_PDF,0.1,2;
//            stderr epinf_bar, 0.053, 0.01,3,INV_GAMMA_PDF,0.1,2;
            stderr ephitilde, 0.183, 0.01,3,INV_GAMMA_PDF,0.1,2;
 
        // persistence parameters
            crhoa,0.9594 ,.01,.9999,BETA_PDF,0.5,0.20;
            crhob,0.7134,.01,.9999,BETA_PDF,0.5,0.20;
            crhog,0.9828,.01,.9999,BETA_PDF,0.5,0.20;
            crhoqs,0.5531,.01,.9999,BETA_PDF,0.5,0.20;
            crhoms,0.2172,.01,.9999,BETA_PDF,0.5,0.20;
            crhopinf,0.9909,.01,.9999,BETA_PDF,0.5,0.20;
            crhow,0.9755,.001,.9999,BETA_PDF,0.5,0.20;
//            crhospread,.9220,.01,.9999,BETA_PDF,0.5,0.20; // .9894 is the sample AR(1) coefficient 
            crhophi_tilde,0.955,.01,.9999,BETA_PDF,0.5,0.20;
            crhopinf_mc,0.970,0.01,.9999,BETA_PDF,0.5,0.2;
            crhopinf_mi,0.963,0.01,.9999,BETA_PDF,0.5,0.2;
            crhopinf_x,0.886,0.01,.9999,BETA_PDF,0.5,0.2;
       
        // markups
            clambda_mc,1.633,0.01,10,INV_GAMMA_PDF,0.1,2;
            clambda_mi,1.275,0.01,10,INV_GAMMA_PDF,0.1,2;

        // markup parameters, MA processes
            cmap,0.9095,0.01,.9999,BETA_PDF,0.5,0.2;
            cmaw,0.9318,0.01,.9999,BETA_PDF,0.5,0.2;
            cmap_mc,0.9,0.01,.9999,BETA_PDF,0.5,0.2;
            cmap_mi,0.9,0.01,.9999,BETA_PDF,0.5,0.2;
            cmap_x,0.9,0.01,.9999,BETA_PDF,0.5,0.2;

        // substitution elasticity
            ceta_i,1.696,0.01,5,INV_GAMMA_PDF,0.1,2;
            ceta_f,1.486,0.01,5,INV_GAMMA_PDF,0.1,2;

        // various DU and SW parameters
            csadjcost,5.0068,2,15,NORMAL_PDF,4,1.5;
            cgy,0.2678,0.01,2.0,NORMAL_PDF,0.5,0.25;
            calfa,0.2419,0.01,1.0,NORMAL_PDF,0.3,0.05;
            csigl,1.8216,0.25,10,NORMAL_PDF,2,0.75;
            czcap,0.5169,0.01,1,BETA_PDF,0.5,0.15;
            cfc,1.5,1.0,3,NORMAL_PDF,1.25,0.125; 
            csigma,1.2013,0.01,3,NORMAL_PDF,1.50,0.375;
            chabb,0.8019 ,0.001,0.9999,BETA_PDF,0.7,0.1;
        
        // Calvo parameters
            cprobp,0.5430,0.3,0.93,BETA_PDF,0.5,0.10;
            cprobw,0.7844,0.3,0.95,BETA_PDF,0.5,0.1;
            cprob_mc, 0.463,0.3,0.95,BETA_PDF,0.5,0.1;
            cprob_mi, 0.740,0.3,0.95,BETA_PDF,0.5,0.1;
            cprob_x, 0.639,0.3,0.95,BETA_PDF,0.5,0.1;

        // inflation indexation
            cindw,0.5430,0.01,0.99,BETA_PDF,0.5,0.15;
            cindp,0.2430,0.01,0.99,BETA_PDF,0.5,0.15;
            cind_mc,0.220,0.01,0.99,BETA_PDF,0.5,0.15;
            cind_mi,0.231,0.01,0.99,BETA_PDF,0.5,0.15;
            cind_x,0.185,0.01,0.99,BETA_PDF,0.5,0.15;

        // monetary policy parameters
            crpi,1.8675,1.03,3,NORMAL_PDF,1.5,0.25;
            crr,0.8948,0.5,0.975,BETA_PDF,0.75,0.10;
            cry, 0.0674,0.001,0.5,NORMAL_PDF,0.125,0.05;
            crdy, 0.1801,0.001,0.5,NORMAL_PDF,0.125,0.05;
//            crdpi, 0.310,0.001,0.5,NORMAL_PDF,0.125,0.05;
            crx, -0.009,-0.01,0.5,NORMAL_PDF,0.125,0.05;

        // measurement equations   
            constepinf,0.6795,0.1,2.0,GAMMA_PDF,0.625,0.1;//20;
            constebeta, 0.1011,0.01,2.0,GAMMA_PDF,0.25,0.1;//0.20;
            constelab,-0.4086,-10.0,10.0,NORMAL_PDF,0.0,2.0;
 //           constex,0.5,-10.0,10.0,NORMAL_PDF,0.0,2.0;     
            ctrend,0.4940,0.1,0.8,NORMAL_PDF,0.4,0.10;
 //           consterspread,0.4551,0.1,2.0,GAMMA_PDF,0.5,0.1;

        // foreign bond risk premium
            cphi_a,0.252,0.01,3,INV_GAMMA_PDF,0.1,2;

        end;

    write_latex_prior_table;

    varobs dy dc dinve labobs pinfobs dw robs;
//rspreadobs;
//xobs dx dm;


    // 1995 Q1 - 2016 Q4
    estimation(order=1, optim=('MaxIter',200), datafile=data, mode_file=Posterior_Generation_mode, mode_compute=5, first_obs=1,presample=4,lik_init=2,prefilter=0,mh_replic=20000,mh_nblocks=1,mh_jscale=0.275,mh_drop=0.2 , mode_check , forecast=12, nodiagnostic, tex) labobs robs pinfobs dy dc dinve dw ewma epinfma spinf sw zcapf rkf kf pkf cf invef yf labf wf rrf kpf bondsf mc zcap rk k pk c inve y lab pinf w r a b g qs ms kp bonds;
    //nobs=84,

    stoch_simul(irf=20,nograph);
    shock_decomp = 1;
    shock_decomposition y;

    collect_latex_files;

