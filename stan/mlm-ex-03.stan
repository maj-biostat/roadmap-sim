// Model for:
// independent estimates of treatment arm means
// shared within group variance for deflections from trt arm means 
// due to subgroup membership
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
  real mu;
  // between treatment arm variation 
  real<lower=0> s;
  vector[J] z_j;
  // parameters to construct offsets from trt arm means
  real<lower=0> s_j;
  matrix[J,K] z_j_k;
}
transformed parameters{
  vector[J] mu_j;
  matrix[J, K] mu_j_k;
  vector[N] eta;
  
  mu_j = mu + z_j .* s; 
  for(i in 1:J){
    mu_j_k[i, ] = to_row_vector(mu_j[i] + z_j_k[i, ] .* s_j);  
  }
  for(i in 1:N){
    eta[i] = mu_j_k[j[i], k[i]];
  }
}
model{
  // overall mean
  target += logistic_lpdf(mu | 0, 1);
  // b/w trt arm variation
  target += std_normal_lpdf(to_vector(z_j));   
  target += student_t_lpdf(s | 3, 0, 2) - 
     1 * student_t_lccdf(0 | 3, 0, 2);
  // w/in trt arm variation   
  target += std_normal_lpdf(to_vector(z_j_k));      
  target += student_t_lpdf(s_j | 3, 0, 2) - 
     1 * student_t_lccdf(0 | 3, 0, 2);
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta);  
  }
}
generated quantities{
}
