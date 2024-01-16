data {
  // early
  int N_e;
  array[N_e] int e_su;
  array[N_e] int e_y;
  array[N_e] int e_n;
  vector[N_e] e_ec;
  vector[N_e] e_ecp; // 1 - e_ec
  array[N_e] int e_c;
  
  int N_l;
  array[N_l] int l_su;
  array[N_l] int l_y;
  array[N_l] int l_n;
  vector[N_l] l_ec;
  vector[N_l] l_ecp; // 1 - l_ec
  array[N_l] int l_c;
  vector[N_l] l_ea;
  vector[N_l] l_eap;
  array[N_l] int l_a;
  vector[N_l] l_eb;
  vector[N_l] l_eb1p;
  vector[N_l] l_eb2p;
  array[N_l] int l_b;
  
  int N_c;
  array[N_c] int c_su;
  array[N_c] int c_y;
  array[N_c] int c_n;
  vector[N_c] c_ec;
  vector[N_c] c_ecp; // 1 - c_ec
  array[N_c] int c_c;
  vector[N_c] c_ea;
  vector[N_c] c_eap;
  array[N_c] int c_a;
  vector[N_c] c_eb;
  vector[N_c] c_eb1p;
  vector[N_c] c_eb2p;
  array[N_c] int c_b;
  
}
transformed data {
}
parameters {
  vector[6] alpha;
  real gamma_b;
  real gamma_c;
  real b_a_l_raw;
  real b_b1_l_raw;
  real b_b2_l_raw;
  real b_a_c_raw;
  real b_b1_c_raw;
  real b_b2_c_raw;
  real b_c_raw;
}
transformed parameters{
  // Include non-randomised items in parameter vector but set to zero 
  // see b_c[3] below. Gives a simpler way to build up model without index
  // overrun or conditionals for building linear predictor.
  
  vector[2] b_a_l; // dair vs revision (late silo)
  vector[3] b_b1_l; // note length - dair vs revision (late silo, recvd one-stage) 
  vector[3] b_b2_l; // note length - dair vs revision (late silo, recvd two-stage)
  
  vector[2] b_a_c; // dair vs revision (chronic silo)
  vector[2] b_b1_c; // dair vs revision (chronic silo, recvd one-stage)
  vector[2] b_b2_c; // dair vs revision (chronic silo, recvd two-stage)
  
  vector[3] b_c; // note length
  
  b_a_l[1] = 0.0;
  b_a_l[2] = b_a_l_raw;
  
  b_b1_l[1] = 0.0;
  b_b1_l[2] = b_b1_l_raw;
  b_b1_l[3] = 0.0;
  
  b_b2_l[1] = 0.0;
  b_b2_l[2] = b_b2_l_raw;
  b_b2_l[3] = 0.0;
  
  b_a_c[1] = 0.0;
  b_a_c[2] = b_a_c_raw;
  
  b_b1_c[1] = 0.0;
  b_b1_c[2] = b_b1_c_raw;
  
  b_b2_c[1] = 0.0;
  b_b2_c[2] = b_b2_c_raw;
  
  b_c[1] = 0.0;
  b_c[2] = b_c_raw;
  // handles other, but this can be ignored in post-processing
  b_c[3] = 0.0;
}
model{
  target += normal_lpdf(alpha | 0, 1.5);
  
  target += std_normal_lpdf(gamma_b);
  target += std_normal_lpdf(gamma_c);
  // all silos
  target += std_normal_lpdf(b_c_raw);
  // late silo
  target += std_normal_lpdf(b_a_l_raw);
  target += std_normal_lpdf(b_b1_l_raw);
  target += std_normal_lpdf(b_b2_l_raw);
  // chronic silo
  target += std_normal_lpdf(b_a_c_raw);
  target += std_normal_lpdf(b_b1_c_raw);
  target += std_normal_lpdf(b_b2_c_raw);
  
  
  // early
  target += binomial_logit_lpmf(e_y | e_n, alpha[e_su] + e_ec * gamma_c + e_ecp .* b_c[e_c]) ; 
  // late
  target += binomial_logit_lpmf(l_y | l_n, alpha[l_su] + l_eb * gamma_b + l_ec * gamma_c + 
                                  l_eap .* b_a_l[l_a] +
                                  l_eb1p .* b_b1_l[l_b] + l_eb2p .* b_b1_l[l_b] +
                                  l_ecp .* b_c[l_c]) ; 
  // chronic
  target += binomial_logit_lpmf(c_y | c_n, alpha[c_su] + c_eb * gamma_b + c_ec * gamma_c + 
                                  c_eap .* b_a_c[c_a] +
                                  c_eb1p .* b_b1_c[c_b] + c_eb2p .* b_b1_c[c_b] +
                                  c_ecp .* b_c[c_c]) ; 
  
}
generated quantities{
  // vector[6] p_y;
  // p_y = inv_logit(alpha);
}
