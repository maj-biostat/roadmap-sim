// partial pooling
// estimates a treatment effect for each silo assuming exchangeability
data {
  int N;
  array[N] int y;
  array[N] int n;
  vector[N] x;
  array[N] int s;
}
parameters {
  real a;
  vector[3] b_z;
  real<lower=0> sd_b;
}
transformed parameters{
  vector[3] b;
  b = b_z * sd_b;
}
model{
  target += normal_lpdf(a | 0, 1.5);
  target += normal_lpdf(b_z | 0, 1);
  target += exponential_lpdf(sd_b | 100);
  for(i in 1:N){
    target += binomial_logit_lpmf(y[i] | n[i], a + x[i] * b[s[i]]); 
  } 
}
generated quantities{
  matrix[3,2] p;
  real b_avg;
  
  b_avg = a + normal_rng(0, sd_b);
  
  
  for(i in 1:3){
    p[i,1] = inv_logit(a);
    p[i,2] = inv_logit(a + b[i]);
  }
}
