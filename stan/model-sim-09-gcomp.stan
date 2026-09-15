// group all your data up first
// instructions:
data{ 
  
  int N;
  
  array[N] int y;
  
  int K_trt;
  int K_cov1;
  int K_cov2;
  
  array[N] int trt;
  array[N] int cov1;
  array[N] int cov2;
  
  // priors
  vector[2] pri_b_0;
  vector[2] pri_b_trt;
  vector[2] pri_b_cov1;
  vector[2] pri_b_cov2;
}
transformed data{
  
}
parameters{
  real b_0;
  vector[K_trt-1] b_trt_raw;
  vector[K_cov1-1] b_cov1_raw;
  vector[K_cov2-1] b_cov2_raw;
}
transformed parameters{
  vector[K_trt] b_trt;
  vector[K_cov1] b_cov1;
  vector[K_cov2] b_cov2;
  
  b_trt[1] = 0.0;
  b_cov1[1] = 0.0;
  b_cov2[1] = 0.0;
  
  b_trt[2:K_trt] = b_trt_raw;
  b_cov1[2:K_cov1] = b_cov1_raw;
  b_cov2[2:K_cov2] = b_cov2_raw;
  
} 
model{
  target += logistic_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  target += normal_lpdf(b_trt_raw | pri_b_trt[1], pri_b_trt[2]);
  target += normal_lpdf(b_cov1_raw | pri_b_cov1[1], pri_b_cov1[2]);
  target += normal_lpdf(b_cov2_raw | pri_b_cov2[1], pri_b_cov2[2]);
  
  target += bernoulli_logit_lpmf(y | b_0 + b_trt[trt] + b_cov1[cov1] + b_cov2[cov2]);  

}
generated quantities{
  
  vector[N] w = dirichlet_rng(rep_vector(1.0, N));
  
  vector[N] p_1 = inv_logit(b_0 + b_trt[1] + b_cov1[cov1] + b_cov2[cov2]);
  vector[N] p_2 = inv_logit(b_0 + b_trt[2] + b_cov1[cov1] + b_cov2[cov2]);
  
  real theta_1 = w' * p_1   ;
  real theta_2 = w' * p_2   ;
  
  real mu_1 = mean(p_1)   ;
  real mu_2 = mean(p_2)  ;
  
  
}
