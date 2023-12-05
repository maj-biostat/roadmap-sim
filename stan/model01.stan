data {
  int N;
  array[N] int y;
  array[N] int n;
  // covariates excludes intercept
  int P;
  matrix[N,P] X;
}
parameters {
  real b0;
  vector[P] b;  
}
model{
  target += normal_lpdf(b0 | 0, 1.5);
  target += normal_lpdf(b | 0, 1);
  target += binomial_logit_lpmf(y | n, b0 + X * b); 
}
generated quantities{
  vector[N] p_y;
  p_y = inv_logit(b0 + X * b);
}
