DATA_SECTION

  init_int nobs
  init_vector nexposed(1,15)
  init_vector TBL(1,15)
  init_vector Kill(1,15)

PARAMETER_SECTION

  objective_function_value f
  init_bounded_number c(0,1)
  init_bounded_number d(0,50)
  init_bounded_number g(-1,25)
  vector prob(1,nobs)    // per capita mort prob
  sdreport_number r_c
  sdreport_number r_d
  sdreport_number r_g
PROCEDURE_SECTION

  dvariable fpen=0.0;         // penalty variable
  // power-Ricker
  prob = c*pow(elem_prod(TBL/d,exp(1-TBL/d)),g);
  // penalties: constrain 0.001 <= prob <= 0.999
  prob = posfun(prob,0.001,fpen);
  f += 1000*fpen;
  prob = 1-posfun(1-prob,0.001,fpen);
  f += 1000*fpen;
  // binomial negative log-likelihood
  f -= sum( log_comb(nexposed,Kill)+
            elem_prod(Kill,log(prob))+
            elem_prod(nexposed-Kill,log(1-prob)));
  r_c=c;
  r_d=d;
  r_g=g;
