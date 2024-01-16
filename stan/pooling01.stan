// pooled
// estimates a single treatment effect averaged across all silos
data {
  int N;
  array[N] int y;
  array[N] int n;
  vector[N] x;
}
parameters {
  real a;
  real b;
}
model{
  target += normal_lpdf(a | 0, 1.5);
  target += normal_lpdf(b | 0, 1);
  for(i in 1:N){
    target += binomial_logit_lpmf(y[i] | n[i], a + x[i] * b); 
  }
}
generated quantities{
  vector[2] p;
  p[1] = inv_logit(a);
  p[2] = inv_logit(a + b);
}
