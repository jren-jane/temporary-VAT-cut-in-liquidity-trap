
// *******************************************************
//				MODEL
// *******************************************************


// *************************************************
// Model equations for extended Smets-Wouters model
// *************************************************

model; 


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
	      cRAf = (chabb/cgamma)/(1+chabb/cgamma)*cRAf(-1) + (1/(1+chabb/cgamma))*cRAf(+1) 
			+((csigma-1)*(1-ctau_l)*(cwlc/clandaw)*(ccy/ccRAy)/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) 
			- (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+b -gamma_cdf) + 0*b ;
            //the zero premultiplying b refers to the scaling of the shock only. 

        // B.86. Total consumption is a population and consumption share weighted average
          cf=(ccRAy/ccy)*(1-cphi)*cRAf+(ccRoTy/ccy)*cphi*cRoTf;

        // B.63. Production function
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );

        // B.71. Wage-setting
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) + 1*swf;

        // B.75. Rule of Thumb agents consume what they get:
          cRoTf =  cwlc*(ccy/ccRoTy)*((1-ctau_l)*(wf+labf))+ (csy/ccRoTy)*sRoTf+(1/ccRoTy)*(yf/clandap);

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


        // UIP condition
          dSf(+1)-(rrf-rrstarf)-cphi_a*asf+phi_tilde=0;

        // aggregate resource constraint
          (1-comega_c)*cgamma_cd^ceta_c*ccy*(cf+ceta_c*gamma_cdf)
            +(1-comega_i)*cgamma_id^ceta_i*ciy*(invef+ceta_i*gamma_idf)
            +g+cystary*(ystarf-ceta_f*gamma_xstarf)
            =yf-1*crkky*zcapf;

        // net foreign assets
           cgamma*cpie*asf=cgamma*cpie*cr*(-cystary*mc_xf-ceta_f*cystary*gamma_xstarf+cystary*ystarf+(ccmc*ccy+cimi*ciy)*gamma_ff
               -(ccmc*ccy*(-ceta_c*(1-comega_c)*(cgamma_cd)^(-(1-ceta_c))*gamma_mcdf +cf)
               +cimi*ciy*(-ceta_i*(1-comega_i)*(cgamma_id)^(-(1-ceta_i))*gamma_midf+invef)))
               +cr*asf(-1);

        // CPI_Inflation
          pinf_cf=((1-comega_c)*(cgamma_dc)^(1-ceta_c))*pinff+((comega_c)*(cgamma_mcc)^(1-ceta_c))*pinf_mcf;

        // real exchange rate
          xf=-comega_c*(cgamma_cmc)^(-(1-ceta_c))*gamma_mcdf-gamma_xstarf-mc_xf;

        // B.76. Government budget constraint
        (1/cgamma)*(cgamma*g+(cby/cpie)*(bondsf(-1)-0)
        -cgamma*ctau_k*cky*(crk*rkf+(crk-cdelta)* kpf(-1)) - ctau_l*cwlc*ccy*(wf+labf)*cgamma 
        -(ctau_v/(1+ctau_v))*ccy*cf*cgamma - ccy*cgamma*dtau_vf/((1+ctau_v)^2))
        + csy*((1-cphi)*sRAf+cphi*sRoTf) = 0;


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
	      cRA = (chabb/cgamma)/(1+chabb/cgamma)*cRA(-1) + (1/(1+chabb/cgamma))*cRA(+1) +
			((csigma-1)*(1-ctau_l)*(cwlc/clandaw)*(ccy/ccRAy)/(csigma*(1+chabb/cgamma)))*(lab-lab(+1))
			- (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) +b -gamma_cd) +0*b ;
            //the zero premultiplying b refers to the scaling of the shock only. 
            // notice that the SIGN of b is flipped (consinstently) compared to the original SW code: now it matches the write-up.

        // UIP condition
          dS(+1)-(b+r-rstar)-cphi_a*as+phi_tilde=0;

        // B.86. Total consumption is a population and consumption share weighted average
          c=(ccRAy/ccy)*(1-cphi)*cRA+(ccRoTy/ccy)*cphi*cRoT;

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
//	     pinf-pinf_bar =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*(pinf(1)-crhopi*pinf_bar) +cindp*(pinf(-1)-pinf_bar) -cbetabar*cgamma*cindp*(1-crhopi)*pinf_bar
//           +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc))    + spinf 
//           +((1-cprobp)/cprobp) * (ctau_v/(1+ctau_v)) * (1/(1+cbetabar*cgamma*cindp)) * (dtau_v - cprobp*cgamma*cbetabar*dtau_v(+1));

	     pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
           +((1-ctau_v)*(1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf
           +((1-cprobp)/cprobp) * ((1-cbetabar*cgamma*cindp)/(1+cbetabar*cgamma*cindp)) * (ctau_v/(1+ctau_v))* dtau_v;

        // Phillips curve for imported consumption goods
//	     pinf_mc-pinf_bar =  (1/(1+cbetabar*cgamma*cind_mc)) * ( cbetabar*cgamma*(pinf_mc(1)-crhopi*pinf_bar) +cind_mc*(pinf_mc(-1)-pinf_bar) -cbetabar*cgamma*cind_mc*(1-crhopi)*pinf_bar
//           +((1-cprob_mc)*(1-cbetabar*cgamma*cprob_mc)/cprob_mc)/((cfc-1)*curvp+1)*(mc_mc)  )  + spinf_mc
//           +((1-cprob_mc)/cprob_mc) * (ctau_vmc/(1+ctau_vmc)) * (1/(1+cbetabar*cgamma*cind_mc)) * (dtau_vmc - cprob_mc*cgamma*cbetabar*dtau_vmc(+1));

	     pinf_mc =  (1/(1+cbetabar*cgamma*cind_mc)) * ( cbetabar*cgamma*pinf_mc(1) +cind_mc*pinf_mc(-1) 
           +((1-ctau_vmc)*(1-cprob_mc)*(1-cbetabar*cgamma*cprob_mc)/cprob_mc)/((cfc-1)*curvp+1)*(mc_mc)  )  + spinf_mc
           +((1-cprob_mc)/cprob_mc) * ((1-cbetabar*cgamma*cind_mc)/(1+cbetabar*cgamma*cind_mc)) * (ctau_vmc/(1+ctau_vmc))* dtau_vmc;

        
        // imported consumption goods
         mc_mc=-mc_x-gamma_xstar-gamma_mcd;

        // Phillips curve for imported investment goods
//	     pinf_mi-pinf_bar =  (1/(1+cbetabar*cgamma*cind_mi)) * ( cbetabar*cgamma*(pinf_mi(1)-crhopi*pinf_bar) +cind_mi*(pinf_mi(-1)-pinf_bar) -cbetabar*cgamma*cind_mi*(1-crhopi)*pinf_bar
//           +((1-cprob_mi)*(1-cbetabar*cgamma*cprob_mi)/cprob_mi)/((cfc-1)*curvp+1)*(mc_mi)  )  + spinf_mi
//           +((1-cprob_mi)/cprob_mi) * (ctau_vmi/(1+ctau_vmi)) * (1/(1+cbetabar*cgamma*cind_mi)) * (dtau_vmi - cprob_mi*cgamma*cbetabar*dtau_vmi(+1));

	     pinf_mi =  (1/(1+cbetabar*cgamma*cind_mi)) * ( cbetabar*cgamma*pinf_mi(1) +cind_mi*pinf_mi(-1) 
           +((1-ctau_vmi)*(1-cprob_mi)*(1-cbetabar*cgamma*cprob_mi)/cprob_mi)/((cfc-1)*curvp+1)*(mc_mi)  )  + spinf_mi
           +((1-cprob_mi)/cprob_mi) * ((1-cbetabar*cgamma*cind_mi)/(1+cbetabar*cgamma*cind_mi)) * (ctau_vmi/(1+ctau_vmi))* dtau_vmi;

        // imported investment goods
         mc_mi=-mc_x-gamma_xstar-gamma_mid;

        // Phillips curve for exporting firms
//	     pinf_x-pinf_bar =  (1/(1+cbetabar*cgamma*cind_x)) * ( cbetabar*cgamma*(pinf_x(1)-crhopi*pinf_bar) +cind_x*(pinf_x(-1)-pinf_bar) -cbetabar*cgamma*cind_x*(1-crhopi)*pinf_bar
//           +((1-cprob_x)*(1-cbetabar*cgamma*cprob_x)/cprob_x)/((cfc-1)*curvp+1)*(mc_x)  )  + spinf_x
//           +((1-cprob_x)/cprob_x) * (ctau_vx/(1+ctau_vx)) * (1/(1+cbetabar*cgamma*cind_x)) * (dtau_vx - cprob_x*cgamma*cbetabar*dtau_vx(+1));

	     pinf_x =  (1/(1+cbetabar*cgamma*cind_x)) * ( cbetabar*cgamma*pinf_x(1) +cind_x*pinf_x(-1) 
           +((1-ctau_vx)*(1-cprob_x)*(1-cbetabar*cgamma*cprob_x)/cprob_x)/((cfc-1)*curvp+1)*(mc_x)  )  + spinf_x
           +((1-cprob_x)/cprob_x) * ((1-cbetabar*cgamma*cind_x)/(1+cbetabar*cgamma*cind_x)) * (ctau_vx/(1+ctau_vx))* dtau_vx;


        // B.84 and an adjustment for the zero interest provision

        // Consider the different interest rates and implement the different ZLB options
        // FEDFUNDS Option 3 is the max{ . } implementation
        // FEDFUNDS Option 2 is a modifified version of CCWT where the transition to standard interest rates is gradual
        // FEDFUNDS Option 1 is CCWT and also the default. Need to set zerointerest to zero to get a standard Taylor rule

@#if ZERO_FEDFUNDS_OPT!=3
     
        r=(1-zerointerest)*rtaylor+zerointerest*r(-1);

            rtaylor = (1-zerointerest_gradual)*
                    ((1-crr)*pinf_bar
                    +crpi*(1-crr)*(pinf_c(-1)-pinf_bar)
                    +cry*(1-crr)*(y(-1)-yf(-1))
                    +crx*(1-crr)*(x(-1)-xf(-1))
                    +crdpi*(pinf_c-pinf_c(-1))
                    +crdy*(y-yf-y(-1)+yf(-1))
                    +crr*rtaylor(-1))
                    +zerointerest_gradual*rtaylor(-1)
                    +ms ;
                   

@#else
// max operator in the DU version
        //r=max(rtaylor, (1-cr+0.000625)*@{CR_SCALE});

// I am using a toolkit

          r = rtaylor;
          [mcp='r > -1.944781619515523']


            rtaylor = (1-crr)*pinf_bar
                    +crpi*(1-crr)*(pinf_c(-1)-pinf_bar)
                    +cry*(1-crr)*(y(-1)-yf(-1))
                    +crx*(1-crr)*(x(-1)-xf(-1))
                    +crdpi*(pinf_c-pinf_c(-1))
                    +crdy*(y-yf-y(-1)+yf(-1))
                    +crr*rtaylor(-1)
                    +ms ;

@#endif

        // CPI_Inflation
          pinf_c=((1-comega_c)*(cgamma_dc)^(1-ceta_c))*pinf+((comega_c)*(cgamma_mcc)^(1-ceta_c))*pinf_mc;

        // real exchange rate
          x=-comega_c*(cgamma_cmc)^(-(1-ceta_c))*gamma_mcd-gamma_xstar-mc_x;
          
        // B.75. Rule of Thumb agents consume what they get:
          cRoT = cwlc*(ccy/ccRoTy)*((1-ctau_l)*(w+lab))+(csy/ccRoTy)*sRoT+(1/ccRoTy)*(y/clandap-mc);

        // B.72 wage determination
		w =  (1/(1+cbetabar*cgamma))*w(-1)							
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)			
               +(cindw/(1+cbetabar*cgamma))*(pinf_c(-1)-pinf_bar)				
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*(pinf_c-crhopi*pinf_bar)		
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*(pinf(1)-crhopi*pinf_bar)				
               -1/(1+cbetabar*cgamma)*(pinf-pinf_bar)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)			
			*(1/((clandaw-1)*curvw+1))*									
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1)-w) 										                	
               + 1*sw ;										             

        // relative prices
        gamma_mcd=gamma_mcd(-1)+pinf_mc-pinf;
        gamma_mid=gamma_mid(-1)+pinf_mi-pinf;
        gamma_xstar=gamma_xstar(-1)+pinf_x-pinfstar;
        mc_x=mc_x(-1)+pinf-pinf_x-dS;
        gamma_f=mc_x+gamma_xstar;
        gamma_cd=comega_c*(cgamma_mcc)^(1-ceta_c)*gamma_mcd;
        gamma_id=comega_i*(cgamma_mii)^(1-ceta_i)*gamma_mid;

        // B.76. Government budget constraint
        (1/cgamma)*(cgamma*g+(cby/cpie)*(bonds(-1)-pinf)
        -cgamma*ctau_k*cky*(crk*rk+(crk-cdelta)* kp(-1)) - ctau_l*cwlc*ccy*(w+lab)*cgamma
        -(ctau_v/(1+ctau_v))*ccy*c*cgamma - ccy*cgamma*dtau_v/((1+ctau_v)^2))
        + csy*((1-cphi)*sRA+cphi*sRoT) = 0;

         // B.68. Law of motion for capital
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs*(1+cbetabar*cgamma) ;


// *************************************************
// FOREIGN VAR
// *************************************************

//Foreign VAR
   pinfstar  = ForLag111*pinfstar(-1)  + ForLag112*ystar(-1) + ForLag113*rstar(-1)
             +ForLag211*pinfstar(-2)  + ForLag212*ystar(-2) + ForLag213*rstar(-2)
             +ForLag311*pinfstar(-3)  + ForLag312*ystar(-3) + ForLag313*rstar(-3)
      		+ForLag411*pinfstar(-4)  + ForLag412*ystar(-4) + ForLag413*rstar(-4) 
            + ForShock11*epinfstar + ForShock12*eystar + ForShock13*erstar ;

   ystar = ForLag121*pinfstar(-1)  + ForLag122*ystar(-1) + ForLag123*rstar(-1)
             +ForLag221*pinfstar(-2)  + ForLag222*ystar(-2) + ForLag223*rstar(-2)
             +ForLag321*pinfstar(-3)  + ForLag322*ystar(-3) + ForLag323*rstar(-3)
      		+ForLag421*pinfstar(-4)  + ForLag422*ystar(-4) + ForLag423*rstar(-4) 
            + ForShock21*epinfstar + ForShock22*eystar + ForShock23*erstar ;

   rstar  = ForLag131*pinfstar(-1)  + ForLag132*ystar(-1) + ForLag133*rstar(-1)
             +ForLag231*pinfstar(-2)  + ForLag232*ystar(-2) + ForLag233*rstar(-2)
             +ForLag331*pinfstar(-3)  + ForLag332*ystar(-3) + ForLag333*rstar(-3)
      		+ForLag431*pinfstar(-4)  + ForLag432*ystar(-4) + ForLag433*rstar(-4) 
            + ForShock31*epinfstar + ForShock32*eystar + ForShock33*erstar ;


//Foreign VAR for flexible economy (added)
   pinfstarf  = 0;

   ystarf = ForLag121*pinfstarf(-1)  + ForLag122*ystarf(-1) + ForLag123*rrstarf(-1)
             +ForLag221*pinfstarf(-2)  + ForLag222*ystarf(-2) + ForLag223*rrstarf(-2)
             +ForLag321*pinfstarf(-3)  + ForLag322*ystarf(-3) + ForLag323*rrstarf(-3)
      		+ForLag421*pinfstarf(-4)  + ForLag422*ystarf(-4) + ForLag423*rrstarf(-4) 
            + 0*epinfstar + ForShock22*eystar + ForShock23*erstar ;

   rrstarf  = ForLag131*pinfstarf(-1)  + ForLag132*ystarf(-1) + ForLag133*rrstarf(-1)
             +ForLag231*pinfstarf(-2)  + ForLag232*ystarf(-2) + ForLag233*rrstarf(-2)
             +ForLag331*pinfstarf(-3)  + ForLag332*ystarf(-3) + ForLag333*rrstarf(-3)
      		+ForLag431*pinfstarf(-4)  + ForLag432*ystarf(-4) + ForLag433*rrstarf(-4) 
            + 0*epinfstar + ForShock32*eystar + ForShock33*erstar ;


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

        // Markup shocks in the flexible price-wage economy equals zero
          spinff = 0;
          spinf_mcf = 0;
          spinf_mif = 0;
          spinf_xf = 0;
          swf = 0;

          pinf_bar=crhopi*pinf_bar(-1)+epinf_bar; // moving inflation target

          rspread=crhospread*1*rspread(-1)+erspread; // B.85j. private bond spread

          phi_tilde=crhophi_tilde*phi_tilde(-1)+ephitilde; // foreign bond premium

          sRoTf=@{cpsi_sf}*sRAf; // allowing for a free ratio between transfers to two types of HH's
          sRoT=@{cpsi_s}*sRA;

        // Set debt change to zero. Interpretation: spending offset by changes in lump-sum taxes.
        // Bonds only exist to determine the interest rate
          bondsf = 0;
          bonds = 0;

// *************************************************
// MEASUREMENT EQUATIONS
// *************************************************

        dy=y-y(-1)+ctrend; // B.89a. output, scaled by 100
        dc=c-c(-1)+ctrend; // B.89b. consumption, scaled by 100
        dinve=inve-inve(-1)+ctrend; // B.89c. private investment, scaled by 100
        dw=w-w(-1)+ctrend; // B.89e. wage, scaled by 100
        pinfobs = 1*(pinf) + constepinf; // B.89f. inflation
        robs =    1*(r) + conster; // B.89g. bank rate
        rspreadobs = rspread+consterspread; // B.89h. private bond wedge
        labobs = lab + constelab; // B.89i. labor
        
        xobs = x - constex; // real exchange rate
        dx = (ystar-ceta_f*gamma_xstar) - (ystar(-1)-ceta_f*gamma_xstar(-1)) + ctrend; // real exports
        dm = ccmim / (ccmim + 1) * ((c - c(-1)) - ceta_c*(1-comega_c)*(cgamma_cd)^(-(1-ceta_c)) * (gamma_mcd-gamma_mcd(-1)))
            + 1 / (ccmim + 1) * ((inve - inve(-1)) - ceta_i*(1-comega_i)*(cgamma_id)^(-(1-ceta_i)) * (gamma_mid-gamma_mid(-1)))
            + ctrend; // real imports

end; 


