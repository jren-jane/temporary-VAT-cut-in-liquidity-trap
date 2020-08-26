// I modified this code which was originally written by 
// Drautzburg and Uhlig (2015)

shocks;

// Gov. bond premium
       var eb;  
       periods 1;
       values @{BONDSHOCK};

// Corporate gov. bond premium
       var erspread;  
       periods 1;
       values @{SPREADSHOCK};

// VAT cut
       var dtau_v;
       periods @{CUT_START}:@{CUT_END};
       values -0.05; // The shocks are scaled with the steady state rates

// Modelling the ZLB by force
@#if ZERO_FEDFUNDS_OPT==1
  var zerointerest;
  periods 1:@{NO_ZERO_QTRS};
  values 1;
@#endif

@#if ZERO_FEDFUNDS_OPT==2
  var zerointerest_gradual;
  periods 1:@{NO_ZERO_QTRS};
  values 1;
@#endif

end;