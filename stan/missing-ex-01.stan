// Model for:
data {
  int N_obs;
  array[N_obs] int y;
  array[N_obs] int n;
  int P;
  matrix[N_obs,P] X_obs;
  int prior_only;
}
transformed data {
}
parameters{
  vector[P] b;
}
transformed parameters{ 
  vector[N_obs] eta = X_obs*b;
}
model{
  target += normal_lpdf(b | 0, 3);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }
}
generated quantities{
  
}
