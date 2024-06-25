data {
  int N;
  array[N] int y;
  int P;
  matrix[N, P] X;
 
  int prior_only;
  real sd_a0;
  real sd_b;
}

parameters{
  real a0;
  vector[P] b;
}
transformed parameters{
  vector[N] eta;
  eta = a0 + X * b;
}
model{
  target += normal_lpdf(a0 | 0, sd_a0);
  target += normal_lpdf(b | 0, sd_b);
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta);
  }
}
generated quantities{
  real marg_p0;
  real marg_p1;
  real rd;

  vector[N] cond_p0;
  vector[N] cond_p1;

  // to get to pate rather than sate
  // Bayesian bootstrap (weights)
  vector[N] wgts = dirichlet_rng(rep_vector(1, N));

  // first term in X relates to intervention of interest
  cond_p0 = inv_logit(a0 + 0    + X[, 3:P] * b[3:P]);
  cond_p1 = inv_logit(a0 + X * b);

  // taking average over bayesian bootstrap weights
  marg_p0 = wgts' * cond_p0;
  marg_p1 = wgts' * cond_p1;
  rd = marg_p1 - marg_p0;
  
}
