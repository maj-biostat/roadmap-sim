// make it work, make it right, make it quicker
data {
  int N;
  array[N] int y;
  // intervention arms (one less that total num arms)
  int J;
  array[N] int j;
  // subgroups of interest creating effect heterogeneity
  int K;
  array[N] int k;
  int prior_only;
}
transformed data {
}
parameters{
  
  // overall trt effect across all intervention types
  real mu;
  // b/w trt type variation
  real<lower=0> s;
  vector[J] z;
  
  // within trt variation due to group membership 
  // (shared across  all trt types)
  real<lower=0> s_j;
  matrix[J, K] z_j;
  
}
transformed parameters{
  
  vector[N] eta;
  
  // no reference group, just model variation on response over all arms
  // assuming exchangeability
  
  // trt x groups
  matrix[J, K] mu_j_k;
  // average effect by trt type
  vector[J] mu_j;
  // within each joint type we have some mean effect of revision
  mu_j = mu + z * s;
  // within trt type offsets
  for(i in 1:J){
    mu_j_k[i, ] = to_row_vector(mu_j[i] + z_j[i, ] .* s_j);
  }
  
  for(i in 1:N){
    eta[i] = mu_j_k[j[i], k[i]];
  }
  
}
model{
  
  target += std_normal_lpdf(mu);
  
  target += student_t_lpdf(s | 3, 0, 1) - 
     1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(z);
  
  target += student_t_lpdf(s_j | 3, 0, 1) - 
     1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(to_vector(z_j));
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta);  
  }
  
}
generated quantities{

}
