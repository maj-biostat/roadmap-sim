data {
  int<lower=1> N;
  array[N] int y;
  array[N] int n;

  array[N] int domA;
  array[N] int s;
  array[N] int j;
  
  // for generated quantities block, much simpler to set this up in R
  // rather than stan
  array[3] int N_s;

  array[N_s[1]] int ix_s1;
  array[N_s[2]] int ix_s2;
  array[N_s[3]] int ix_s3;

  array[2] int N_j;
  array[N_j[1]] int ix_j1;
  array[N_j[2]] int ix_j2;
  
}
transformed data{
  
  array[N_s[1]] int n_s1 = n[ix_s1];
  array[N_s[2]] int n_s2 = n[ix_s2];
  array[N_s[3]] int n_s3 = n[ix_s3];
  
  array[N_j[1]] int n_j1 = n[ix_j1];
  array[N_j[2]] int n_j2 = n[ix_j2];  
}
parameters {
  real alpha;
  
  vector[2] bA_raw;
  vector[2] bs_raw;
  real bj_raw;

  real<lower=0.0> s_domA_sg;
  matrix[3,3] z_As;
  matrix[3,2] z_Aj;
}
transformed parameters{
  vector[3] bA;
  vector[3] bs;
  vector[2] bj;
  matrix[3,3] gAs;
  matrix[3,2] gAj;
  
  bA[1] = 0.0;
  bA[2:3] = bA_raw;
  bs[1] = 0.0;
  bs[2:3] = bs_raw;
  bj[1] = 0.0;
  bj[2] = bj_raw;
  
  gAs = s_domA_sg * z_As;
  gAj = s_domA_sg * z_Aj;
  
  vector[N] eta;
  for(i in 1:N){
    eta[i] = alpha + bA[domA[i]] + 
      bs[s[i]] + gAs[domA[i],s[i]] + 
      bj[j[i]] + gAj[domA[i],j[i]];
  }
}
model {
  
  target += logistic_lpdf(alpha | 0, 1);
  target += std_normal_lpdf(bA_raw);
  target += std_normal_lpdf(bs_raw);
  target += std_normal_lpdf(bj_raw);
  target += std_normal_lpdf(to_vector(z_As));
  target += std_normal_lpdf(to_vector(z_Aj));
  // single variance across all subgroups
  target += exponential_lpdf(s_domA_sg | 1);
  // likelihood
  target += binomial_logit_lpmf(y | n, eta);
  
}
generated quantities{
  
  // matrix[3,3] p_s;
  // matrix[3,2] p_j;
  // 
  vector[N] p1;
  vector[N] p2;
  vector[N] p3;

  for(i in 1:N){
    // treatment arm 1
    p1[i] = inv_logit(alpha + bA[1] + bs[s[i]] + gAs[1,s[i]] + bj[j[i]] + gAj[1,j[i]]);
    // treatment arm 2
    p2[i] = inv_logit(alpha + bA[2] + bs[s[i]] + gAs[2,s[i]] + bj[j[i]] + gAj[2,j[i]]);
    // treatment arm 3
    p3[i] = inv_logit(alpha + bA[3] + bs[s[i]] + gAs[3,s[i]] + bj[j[i]] + gAj[3,j[i]]);
  }
  
  vector[N_s[1]] w_s1 = dirichlet_rng(to_vector(n_s1));
  vector[N_s[2]] w_s2 = dirichlet_rng(to_vector(n_s2));
  vector[N_s[3]] w_s3 = dirichlet_rng(to_vector(n_s3));

  // mean risk silo1 treatment arms 1, 2, 3
  vector[3] mu_p_dA_s1 = [ w_s1' * p1[ix_s1], w_s1' * p2[ix_s1], w_s1' * p3[ix_s1] ]';
  // mean risk silo2 treatment arms 1, 2, 3
  vector[3] mu_p_dA_s2 = [ w_s2' * p1[ix_s2], w_s2' * p2[ix_s2], w_s2' * p3[ix_s2] ]';
  // mean risk silo3 treatment arms 1, 2, 3
  vector[3] mu_p_dA_s3 = [ w_s3' * p1[ix_s3], w_s3' * p2[ix_s3], w_s3' * p3[ix_s3] ]';
  
  // mean risk differences for silo1 (rev(1) vs dair, rev(2) vs dair)
  vector[2] mu_rd_s1 = [ mu_p_dA_s1[2] - mu_p_dA_s1[1], mu_p_dA_s1[3] - mu_p_dA_s1[1] ]';
  // mean risk differences for silo2 (rev(1) vs dair, rev(2) vs dair)
  vector[2] mu_rd_s2 = [ mu_p_dA_s2[2] - mu_p_dA_s2[1], mu_p_dA_s2[3] - mu_p_dA_s2[1] ]';
  // mean risk differences for silo3 (rev(1) vs dair, rev(2) vs dair)
  vector[2] mu_rd_s3 = [ mu_p_dA_s3[2] - mu_p_dA_s3[1], mu_p_dA_s3[3] - mu_p_dA_s3[1] ]';

  vector[N_j[1]] w_j1 = dirichlet_rng(to_vector(n_j1));
  vector[N_j[2]] w_j2 = dirichlet_rng(to_vector(n_j2));

  // mean risk joint1 treatment arms 1, 2, 3
  vector[3] mu_p_dA_j1 = [ w_j1' * p1[ix_j1], w_j1' * p2[ix_j1], w_j1' * p3[ix_j1] ]';
  // mean risk joint2 treatment arms 1, 2, 3
  vector[3] mu_p_dA_j2 = [ w_j2' * p1[ix_j2], w_j2' * p2[ix_j2], w_j2' * p3[ix_j2] ]'   ;

  // risk differences of interest (rev(1) vs dair, rev(2) vs dair)
  vector[2] mu_rd_j1 = [ mu_p_dA_j1[2] - mu_p_dA_j1[1], mu_p_dA_j1[3] - mu_p_dA_j1[1] ]';
  vector[2] mu_rd_j2 = [ mu_p_dA_j2[2] - mu_p_dA_j2[1], mu_p_dA_j2[3] - mu_p_dA_j2[1] ]';
  
  
}





