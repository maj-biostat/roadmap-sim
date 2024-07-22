// Model for:
data {
  int N;
  array[N] int y;
  int J;
  array[N] int j;
  
  // alpha controls the mix between L1 and L2 penalties. 
  // when alpha = 1, we have pure lasso regression (double exp)
  // when alpha = 0, we have pure ridge regression (normal with 1/lambda scale).
  
  // ridge is just an informative normal prior with scale 1/(2*lambda)
  // which tends to shrink values across the board
  // lasso is using a double exponential prior that will tend to shrink small
  // values towards zero and leave large values as they are
  
  real<lower=0, upper=1> pri_alpha;  // mixing parameter between L1 and L2
  real<lower=0> pri_lambda;  // overall regularization strength
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
  target += logistic_lpdf(mu | 0, 1);
  // implicit uniform prior on z_eta transformed by the following
  // obviously this is just updating the log-density.
  // the implied probability density is proportional to
  // exp(-lambda * (alpha * sum(abs(beta)) + (1 - alpha) * 0.5 * dot_product(beta, beta)))
  target += -pri_lambda * (pri_alpha * sum(abs(z_eta)) +
    (1 - pri_alpha) * 0.5 * dot_product(z_eta, z_eta));
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta[j]);  
  }
}
generated quantities{
}
