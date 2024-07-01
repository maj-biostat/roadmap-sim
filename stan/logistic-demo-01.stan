data {
  int<lower=0> N;
  array[N] int y;
  array[N] int n;
  array[N] int trt;
  array[N] int sex;
}
parameters {
  real b0;
  real b_trt;
  real b_sex;
}
transformed parameters{
  vector[N] eta;
  for(i in 1:N){
    eta[i] = b0 + b_trt * trt[i] + b_sex * sex[i];
  }
}
model {
  target += normal_lpdf(b0 | 0, 1.5);
  target += normal_lpdf(b_trt  | 0, 1);
  target += normal_lpdf(b_sex  | 0, 1);
  // long form can be abbreviated.
  target += binomial_logit_lpmf(y | n, eta) ;
}

