data {
  int N;
  array[N] int y;
  array[N] int n;
  // for one an two-stage indicators 
  // these are the main effects of interest
  int P1;
  matrix[N, P1] X1;
  // other terms bar interactions and ext proph related
  int P2;
  matrix[N, P2] X2;
  // extra fields for extended proph
  int P3;
  matrix[N, P3] X3;
  // fields for interactions
  int P4;
  matrix[N, P4] X4;
  // auxiliary terms
  // observed probability of one-stage
  real pr_one;
  // number preferring 1/2 stage under rev
  int N1;
  int N2;
  array[N1] int ix1;
  array[N2] int ix2;
  // priors
  int prior_only;
  real sd_a0;
  real sd_b1;
  real sd_b2;
  real sd_b3;
  real sd_b4;
}
transformed data {
}
parameters{
  real a0;
  // surgical effects
  vector[P1] b1;
  // other terms
  vector[P2] b2;
  // ext prophy
  vector[P3] b3;
  // interactions
  vector[P4] b4;
}
transformed parameters{
  vector[N] eta;
  for(i in 1:N){
    // X1[i,2] should be an indicator of whether two-stage has been done
    if(X1[i,2] == 0){
      eta[i] = a0 + X1[i,]*b1 + X2[i,]*b2 +             X4[i, ]*b4 ;    
    }
    if(X1[i,2] == 1){
      eta[i] = a0 + X1[i,]*b1 + X2[i,]*b2 + X3[i,]*b3 + X4[i, ]*b4 ;     
    }
  }
}
model{
  target += normal_lpdf(a0 | 0, sd_a0);
  target += normal_lpdf(b1 | 0, sd_b1);
  target += normal_lpdf(b2 | 0, sd_b2);
  target += normal_lpdf(b3 | 0, sd_b3);
  target += normal_lpdf(b4 | 0, sd_b4);
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }
}
generated quantities{
  vector[N] p;    
  real marg_p0;                                               
  real marg_p1;       
  real rd;                                                    
                                                              
  vector[N] cond_p0;                                          
  vector[N] cond_p1;                                               
                                                              
  // to get to pate rather than sate                          
  // Bayesian bootstrap (weights)                             
  vector[N] w = dirichlet_rng(to_vector(n)); 
                                                              
  p = inv_logit(eta);
  
  // drop out the terms relating to sa == 1 and sa == 2
  // drop out last two terms that are only applicable for sa == 2
  // dair - no X4 since sa = 0 hence zero contribution.
  cond_p0 = inv_logit(a0 +         X2*b2); 
  // one
  cond_p1[ix1] = inv_logit(a0 + b1[1] + X2[ix1,]*b2 +         X4[ix1,]*b4); 
  // two
  cond_p1[ix2] = inv_logit(a0 + b1[2] + X2[ix2,]*b2 + X3[ix2,]*b3 + X4[ix2,]*b4);                          
                                                              
  // taking average over bayesian bootstrap weights           
  marg_p0 = w' * cond_p0;                                  
  marg_p1 = w' * cond_p1;    

  rd = marg_p1 - marg_p0;
}
