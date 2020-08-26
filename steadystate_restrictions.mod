

      // steady states of the original model
      // which appear as parameters in the log-linearized model

// *******************************************************************
//
// ALLV parameters
//
// *******************************************************************

       ceta_mc=clambda_mc/(clambda_mc-1); //Substitution elasticity among imported consumption goods
       ceta_mi=clambda_mi/(clambda_mi-1);//Substitution elasticity among imported investment goods
       cgamma_id=((1-comega_i)+comega_i*(ceta_mi/(ceta_mi-1))^(1-ceta_i))^(1/(1-ceta_i)); //A.11
       cgamma_cd=((1-comega_c)+comega_c*(ceta_mc/(ceta_mc-1))^(1-ceta_c))^(1/(1-ceta_c)); //A.9
       cgamma_dc=1/cgamma_cd;
       cgamma_cmc=((1-comega_c)*((ceta_mc-1)/ceta_mc)^(1-ceta_c)+comega_c)^(1/(1-ceta_c)); //A.10
       cgamma_mcc=1/cgamma_cmc;
       cgamma_imi=((1-comega_i)*((ceta_mi-1)/ceta_mi)^(1-ceta_i)+comega_i)^(1/(1-ceta_i)); //A.12
       cgamma_mii=1/cgamma_imi;


       ccmc = comega_c*cgamma_cmc^ceta_c;
       cimi = comega_i*cgamma_imi^ceta_i;
       ccdc = (1-comega_c)*cgamma_cd^ceta_c; //A.25
       cidi = (1-comega_i)*cgamma_id^ceta_i; //A.25

// *******************************************************************
//
// SW / CCWT parameters
//
// *******************************************************************

       clandap=cfc; // B.18, price markup lambda
       cpie=1+constepinf/100; 
       cgamma=1+ctrend/100 ; 
       cbeta=1/(1+constebeta/100);

       cbetabar=cbeta*cgamma^(-csigma); // normalization definition
       cr=cpie/(cbeta*cgamma^(-csigma)); 

       crk=(cbeta^(-1)*(cgamma^csigma) -(1-cdelta) - cdelta*ctau_k) * cgamma_id/(1-ctau_k); // B.33c

       cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa)); // B.11

       cikbar=(1-(1-cdelta)/cgamma); // B.30
       cik=(1-(1-cdelta)/cgamma)*cgamma; // B.30
       clk=((1-calfa)/calfa)*(crk/cw); // B.9
       cky=cfc*(clk)^(calfa-1); // B.8
       ciy=cik*cky;
       ccy=1-cg-cik*cky; // B.62
       crkky=crk*cky;

       cwlc=(1-calfa)/calfa*crk*cky/ccy;
       cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;

       cwly=1-crk*cky; // labor and capital for intermediate firms, B.16, B.18

       crevy=(ctau_l*cwlc*ccy+ctau_k*(crk-cdelta)*cky);

       conster=(cr-1)*100;

       csy=ctau_l*cwlc*ccy+ctau_k*(crk-cdelta)*cky+(ctau_v/(1+ctau_v))*ccy+((cpie*cgamma/cr)-1)*cby/(cpie*cgamma)-cg; // government budget
       ccRoTy=csy+(1-ctau_l)*cwlc*ccy; // B.36
       ccRAy=(ccy-cphi*ccRoTy)/(1-cphi); // B.60


       cystary = comega_c*cgamma_cmc^ceta_c*ccy + comega_i*cgamma_imi^ceta_i*ciy; // from ALLV

       ccmim = (ccmc*ccy) / (cimi*ciy);



 




