// Model for:
// independent estimates of treatment arm by subgroup means
data {
  int N;
  array[N] int y;
  int P;
  matrix[N,P] X;
  real s;
  int prior_only;
}
transformed data {
}
parameters{
  vector[P] b;
}
transformed parameters{ 
  vector[N] eta = X*b;
}
model{
  target += normal_lpdf(b | 0, s);
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta);  
  }
}
generated quantities{
}
