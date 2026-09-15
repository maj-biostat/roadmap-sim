// group all your data up first
// instructions:
data{ 
  
  int N;
  
  array[N] int y;
  array[N] int n;
  
  int K_reg;
  int K_d4;
  
  array[N] int reg;
  array[N] int d4;
  
  // priors
  vector[2] pri_b_0;
  vector[3] pri_b_d4;
  
  int prior_only;
}
transformed data{
  
}
parameters{
  real b_0;
  real mu_reg;
  vector[K_reg-1] z_reg;
  real<lower=0> sig_reg;
  vector[K_d4-1] b_d4_raw;
}
transformed parameters{
  vector[K_reg] b_reg;
  vector[K_d4] b_d4;
  
  b_d4[1] = 0.0;
  
  b_reg[1] = 0.0;
  b_reg[2:K_reg] = mu_reg + z_reg * sig_reg;
  b_d4[2:K_d4] = b_d4_raw;
  
} 
model{
  target += logistic_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  
  target += normal_lpdf(mu_reg | 0, 3);
  target += normal_lpdf(z_reg | 0, 1);
  target += exponential_lpdf(sig_reg | 1);
  
  target += student_t_lpdf(b_d4_raw | pri_b_d4[1], pri_b_d4[2], pri_b_d4[3]);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, b_0 + b_reg[reg] + b_d4[d4]);  
  }

}
generated quantities{
  
}
