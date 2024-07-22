data {
  int N;
  array[N] int y;
  matrix[N, 4] X;
  
  
  int prior_only;
  real mu_a0;
  real mu_b;
  
  real sd_a0;
  real sd_b;
  real r_se;
}

parameters{
  real a0;
  vector[4] b;
  real<lower = 0> s_e;
  // real<lower = 0> s_sa;
  // vector[J] z;
}
transformed parameters{
  vector[N] mu;
  mu = a0 + X * b;
}
model{
  target += normal_lpdf(a0 | mu_a0, sd_a0);
  target += normal_lpdf(b | mu_b, sd_b);
  target += exponential_lpdf(s_e | r_se);
  if(!prior_only){
    target += normal_lpdf(y | mu, s_e);
  }
}
generated quantities{
  // real marg_p0;
  // real marg_p1;
  // real rd;
  // 
  // vector[N] cond_p0;    
  // vector[N] cond_p1;    
  // 
  // // to get to pate rather than sate
  // // Bayesian bootstrap (weights)
  // vector[N] wgts = dirichlet_rng(rep_vector(1, N));    
  //   
  // cond_p0 = a0 + 0    + X[, 2:4] * b[2:4];
  // cond_p1 = a0 + X * b;
  //   
  // // taking average over bayesian bootstrap weights    
  // marg_p0 = wgts' * cond_p0;    
  // marg_p1 = wgts' * cond_p1;    
  // rd = marg_p1 - marg_p0;    
  
}
