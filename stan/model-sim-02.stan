data {
  // early
  int N_e;
  array[N_e] int e_su;
  array[N_e] int e_y;
  array[N_e] int e_n;
  vector[N_e] e_ec; // membership
  vector[N_e] e_ecp; // non-membership, ie 1 - e_ec
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
  vector[N_l] l_eb1;
  vector[N_l] l_eb2;
  vector[N_l] l_ebp;
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
  vector[N_c] c_eb1;
  vector[N_c] c_eb2;
  vector[N_c] c_ebp;
  array[N_c] int c_b;
  
  // priors
  real pri_sig_b_c;
  real pri_sig_a_l;
  real pri_sig_b1_l;
  real pri_sig_b2_l;
  real pri_sig_a_c;
  real pri_sig_b1_c;
  real pri_sig_b2_c;
  
}
transformed data {
  int N = N_e + N_l + N_c;
}
parameters {
  vector[6] alpha;
  // real gamma_b;
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
  
  // target += std_normal_lpdf(gamma_b);
  target += std_normal_lpdf(gamma_c);
  // all silos
  target += normal_lpdf(b_c_raw | 0, pri_sig_b_c);
  // late silo
  target += normal_lpdf(b_a_l_raw | 0, pri_sig_a_l);
  target += normal_lpdf(b_b1_l_raw | 0, pri_sig_b1_l);
  target += normal_lpdf(b_b2_l_raw | 0, pri_sig_b2_l);
  // chronic silo
  target += normal_lpdf(b_a_c_raw | 0, pri_sig_a_c);
  target += normal_lpdf(b_b1_c_raw | 0, pri_sig_b1_c);
  target += normal_lpdf(b_b2_c_raw | 0, pri_sig_b2_c);
  
  // likelihood chunks pertaining to each silo
  target += binomial_logit_lpmf(e_y | e_n, alpha[e_su] +
                                      e_ecp * gamma_c + 
                                      e_ec .* b_c[e_c]) ;      
                                    
  target += binomial_logit_lpmf(l_y | l_n, alpha[l_su] + 
                                    l_ecp * gamma_c +
                                    l_ea .* b_a_l[l_a] +
                                    l_eb1 .* b_b1_l[l_b] +  
                                    l_eb2 .* b_b2_l[l_b] +
                                    l_ec .* b_c[l_c]
                                    ) ;    
                                    
  target += binomial_logit_lpmf(c_y | c_n, alpha[c_su] +
                                      c_ecp * gamma_c +
                                      c_ea .* b_a_c[c_a] +
                                      c_eb1 .* b_b1_c[c_b] + 
                                      c_eb2 .* b_b2_c[c_b] +
                                      c_ec .* b_c[c_c]) ; 
}
generated quantities{
  vector[N_e] eta_e ;
  vector[N_l] eta_l ;
  vector[N_c] eta_c ;
  vector[N] eta;
  
  eta_e = alpha[e_su] + e_ecp * gamma_c + e_ec .* b_c[e_c];
  eta_l = alpha[l_su] + l_ecp * gamma_c +
    l_ea .* b_a_l[l_a] + 
    l_eb1 .* b_b1_l[l_b] + l_eb2 .* b_b2_l[l_b] +
    l_ec .* b_c[l_c];
  eta_c = alpha[c_su] + c_ecp * gamma_c +
    c_ea .* b_a_c[c_a] +
    c_eb1 .* b_b1_c[c_b] + c_eb2 .* b_b2_c[c_b] +
    c_ec .* b_c[c_c];
    
  eta = append_row(eta_e, append_row(eta_l, eta_c));

}
