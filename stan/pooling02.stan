// no pooling
// estimates a treatment effect for each silo
data {
  int N;
  array[N] int y;
  array[N] int n;
  vector[N] x;
  array[N] int s;
}
parameters {
  real a;
  vector[3] b;
}
model{
  target += normal_lpdf(a | 0, 1.5);
  target += normal_lpdf(b | 0, 1);
  for(i in 1:N){
    target += binomial_logit_lpmf(y[i] | n[i], a + x[i] * b[s[i]]); 
  }
  
}
generated quantities{
  matrix[3,2] p;
  
  for(i in 1:3){
    p[i,1] = inv_logit(a);
    p[i,2] = inv_logit(a + b[i]);
  }
}
