data {
  int N;
  array[N] int y;
  array[N] int n;
  matrix[N, 1] X;
  int prior_only;
  real sd_a0;
  real sd_b;
}
parameters{
  real a0;
  vector[1] b;
}
transformed parameters{
  vector[N] eta;
  eta = a0 + X * b;
}
model{
  target += normal_lpdf(a0 | 0, sd_a0);
  target += normal_lpdf(b | 0, sd_b);  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);
  }
}
generated quantities{   
  real marg_p0;                                               
  real marg_p1;    
  real rd;                                                    

  marg_p0 = inv_logit(a0); 
  marg_p1 = inv_logit(a0 + b[1]);                      

  rd = marg_p1 - marg_p0;
}
