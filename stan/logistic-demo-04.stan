// make it work, make it right, make it quicker
data {
  int N;
  array[N] int y;
  // idx into silo
  array[N] int silo;
  // idx into joint
  array[N] int jnt;
  
  // surg type pref (dair, one, two)
  array[N] int u_d1;
  // idx into domain 1 trts
  array[N] int i_d1;
}
transformed data {
}
parameters{
  real b0;
  
  // silo
  matrix[3, 2] b_s;
  // baseline adj for surg pref (one-stage = 2, two-stage = 3)
  vector[2] b_u_d1_raw;
  
  // domain 1
  
  // non-rand effects (dair/hip, one/hip, ..., two/knee)
  vector[6] b_d1_non_rand;
  
  // overall mean for revision effect across one-stage (hip), two-stage (hip)
  // one-stage (knee), two-stage (knee) groups
  real mu_d1;
  real<lower=0> s_mu_d1;
  vector[2] z_mu_d1;
  // within joint variation
  real<lower=0> s_b_d1;
  vector[4] z_b_d1;
  
}
transformed parameters{
  
  vector[N] eta;
  vector[3] b_u_d1;
  
  // domain 1 effects, strat by joint
  matrix[6, 2] b_d1;
  // average effect (offset) of one-stage irrespective of joint
  vector[2] mu_d1_j;
  
  // baseline pref adj 
  b_u_d1[1] = 0.0;
  b_u_d1[2:3] = b_u_d1_raw;
  
  // within each joint type we have some mean effect of revision
  mu_d1_j = mu_d1 + z_mu_d1 * s_mu_d1;
  
  // non-rand effects for those recv dair, one, two in hip and knee
  b_d1[1, ] = to_row_vector(b_d1_non_rand[1:2]);
  b_d1[2, ] = to_row_vector(b_d1_non_rand[3:4]);
  b_d1[3, ] = to_row_vector(b_d1_non_rand[5:6]);
  // reference group for hip and knee
  b_d1[4, ] = rep_row_vector(0.0, 2);
  // around the joint means, we have the effects of one-stage and two-stage rev
  b_d1[5, ] = to_row_vector(mu_d1_j + z_b_d1[1:2] * s_b_d1);
  b_d1[6, ] = to_row_vector(mu_d1_j + z_b_d1[3:4] * s_b_d1);

  for(i in 1:N){
    eta[i] = b0 + 
      b_s[silo[i], jnt[i]] + b_u_d1[u_d1[i]] +
      b_d1[i_d1[i], jnt[i]];
  }
  
}
model{
  
  target += logistic_lpdf(b0 | 0, 1);
  
  // silo/baseline adj
  target += std_normal_lpdf(to_vector(b_s));
  target += std_normal_lpdf(b_u_d1_raw);
  
  // domain 1 effect components
  target += std_normal_lpdf(b_d1_non_rand);
  
  target += std_normal_lpdf(mu_d1);
  target += student_t_lpdf(s_mu_d1 | 3, 0, 2) - 
    1 * student_t_lccdf(0 | 3, 0, 2);
  // target += exponential_lpdf(s_mu_d1 | 1);
  target += std_normal_lpdf(z_mu_d1);
  
  target += student_t_lpdf(s_b_d1 | 3, 0, 2) - 
    1 * student_t_lccdf(0 | 3, 0, 2);
  // target += exponential_lpdf(s_b_d1 | 1);
  target += std_normal_lpdf(z_b_d1);
  
  target += bernoulli_logit_lpmf(y | eta);
  
}
generated quantities{
  
}
