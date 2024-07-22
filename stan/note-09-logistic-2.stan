// Model for:
// independent estimates of treatment arm by subgroup means
data {
  int N;
  array[N] int y;
  array[N] int n;
  // arm index
  array[N, 3] int trt;
  // factor 1 design matrix
  int ncXades;
  int nrXades;
  matrix[nrXades, ncXades] Xades;
  vector[ncXades] sa;
  // factor 2 design matrix
  int ncXbdes;
  int nrXbdes;
  matrix[nrXbdes, ncXbdes] Xbdes;
  vector[ncXbdes] sb;
  // factor 1/2 interaction design matrix
  int ncXabdes;
  int nrXabdes;
  matrix[nrXabdes, ncXabdes] Xabdes;
  vector[ncXabdes] sab;
  
  int prior_only;
}
transformed data {
  // build full design matrices
  matrix[N, ncXades] Xa = Xades[trt[,1]];
  matrix[N, ncXbdes] Xb = Xbdes[trt[,2]];
  // indexes the relevant col in the design matrix
  matrix[N, ncXabdes] Xab = Xabdes[trt[,3]];
}
parameters{
  real mu;
  vector[ncXades] ba;  
  vector[ncXbdes] bb;
  vector[ncXabdes] bab;
}
transformed parameters{ 
  vector[N] eta = mu + Xa*ba + Xb*bb + Xab*bab;
}
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(ba | 0, sa);
  target += normal_lpdf(bb | 0, sb);
  target += normal_lpdf(bab | 0, sab);
  
  if(!prior_only){
    // target += bernoulli_logit_lpmf(y | eta);  
    target += binomial_logit_lpmf(y | n, eta);  
  }
}
generated quantities{
}
