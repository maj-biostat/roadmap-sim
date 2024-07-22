// Model for:
data {
  int N;
  array[N] int y;
  int J;
  array[N] int j;
  real pri_s_norm;
  int prior_only;
}
transformed data {
}
parameters{
  real mu;
  vector[J] z_eta;
}
transformed parameters{ 
  vector[J] eta = mu + z_eta;
}
model{
  // grand mean
  target += logistic_lpdf(mu | 0, 1);
  // offset
  target += normal_lpdf(z_eta | 0, pri_s_norm);
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta[j]);  
  }
}
generated quantities{
}
