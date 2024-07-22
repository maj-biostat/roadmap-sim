// Model for:
data {
  int N;
  array[N] int y;
  int J;
  array[N] int j;
  real pri_r;
  int prior_only;
}
transformed data {
}
parameters{
  real mu;
  real<lower=0> s_y;
  vector[J] z_eta;
}
transformed parameters{ 
  vector[J] eta = mu + z_eta * s_y;
}
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(z_eta | 0, 1);
  target += exponential_lpdf(s_y | pri_r);
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta[j]);  
  }
}
generated quantities{
}
