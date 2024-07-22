// Model for:
// independent estimates of treatment arm by subgroup means
data {
  int N;
  array[N] int y;
  // arm index
  array[N] int x1trt;
  // dimension of design matrix
  int ncX1des;
  int nrX1des;
  matrix[nrX1des, ncX1des] X1des;
  vector[ncX1des] sx1;
  int prior_only;
}
transformed data {
  // build full design matrices
  matrix[N, ncX1des] X1 = X1des[x1trt];
}
parameters{
  real mu;
  vector[ncX1des] bx1;
}
transformed parameters{ 
  vector[N] eta = mu + X1*bx1;
}
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(bx1 | 0, sx1);
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta);  
  }
}
generated quantities{
}
