

data {
  int N;
  array[N] int y;
  array[N] int n;
  // variation due to silo/joint
  vector[N] l1;
  vector[N] l2;
  // reveal
  vector[N] er;
  vector[N] ed;
  vector[N] ef;
  // surgery
  vector[N] r;
  // duration
  vector[N] d;
  // vector[N] rp; // revision was recvd
  vector[N] srp0; // dair recvd
  vector[N] srp1; // one-stage revision recvd
  vector[N] srp2; // two-stage revision recvd
  // choice
  vector[N] f;
  
  vector[2] pri_m_sd;
  vector[9] pri_b_sd;
  int prior_only;
}
transformed data{
  vector[N] erx;
  vector[N] erx_srp1;
  vector[N] erx_srp2;
  vector[N] er_r;
  vector[N] er_r_srp1;
  vector[N] er_r_srp2;
  vector[N] ed_d_srp1;
  vector[N] ed_d_srp2;
  vector[N] ef_f;
  
  erx = (1-er) ;
  erx_srp1 = (1-er) .* srp1;
  erx_srp2 = (1-er) .* srp2;
  
  er_r = er .* r;
  
  er_r_srp1 = er_r .* srp1;
  er_r_srp2 = er_r .* srp2;

// you can remove rp srp1 indicates one-stage happened, srp2 indicates two stage
// the rp becomes redundant.
  ed_d_srp1 = ed .* d .* srp1;
  ed_d_srp2 = ed .* d .* srp2;

  ef_f = ef .* f;
}
parameters {
  real a0;
  vector[2] m;
  vector[9] b;
}
transformed parameters{
  vector[N] eta;
 
  eta = a0 + 
      m[1]*l1 + m[2]*l2 +  
      (b[1]*erx + b[2]*erx_srp1 + b[3]*erx_srp2)  +
      (b[4]*er_r_srp1 + b[5]*er_r_srp2) +
      // b[6] * (1 - ed) + 
      (b[6]*ed_d_srp1 + b[7]*ed_d_srp2)  + 
      b[8] * (1 - ef) + (b[9]*ef_f);
      
}
model{
  target += logistic_lpdf(a0 | 0, 1);
  target += normal_lpdf(m | 0, pri_m_sd);
  target += normal_lpdf(b | 0, pri_b_sd);
  
  // likelihood chunks pertaining to each silo
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta) ;      
  }
  
  
}
generated quantities{

  // vector[N] eta_r_0;
  // vector[N] eta_r_1;
  // vector[N] eta_d_0;
  // vector[N] eta_d_1;
  // vector[N] eta_f_0;
  // vector[N] eta_f_1;
  // 
  // {
  // 
  //   // predictions on log-odds scale setting all participants to the level
  //   // of interest, e.g. assume all have r = 0
  //   eta_r_0   =  a0 +
  //     m[1]*l1   + m[2]*l2 +
  //     (b[1] * erx_srp0 + b[2] * erx_srp1 + b[3] * erx_srp2) +
  //     // b[6] * (1 - ed  ) + 
  //     (b[6]*ed_d_srp1   + b[7]*ed_d_srp2)   +
  //     b[8] * (1 - ef  ) + (b[9]*ef_f)  ;
  //   // r = 1
  //   // (only those revealed to surgery domain contribute due to the er term)
  //   eta_r_1   =  a0 +
  //     m[1]*l1   + m[2]*l2 +
  //     // note the use of er here and not er_r, i.e. the r is set to 1 for
  //     // everyone
  //     (b[1] * erx_srp0 + b[2] * erx_srp1 + b[3] * erx_srp2) +
  //     (b[4]*er .* srp1 + b[5]*er .* srp2) +
  //     // b[6] * (1 - ed  ) + 
  //     (b[6]*ed_d_srp1   + b[7]*ed_d_srp2)   +
  //     b[8] * (1 - ef  ) + (b[9]*ef_f);
  // 
  //   // duration, d = 0
  //   eta_d_0   =  a0 +
  //     m[1]*l1 + m[2]*l2 +
  //     (b[1] * erx_srp0 + b[2] * erx_srp1 + b[3] * erx_srp2) +
  //     (b[4]*er_r_srp1 + b[5]*er_r_srp2) +
  //     // b[6] * (1 - ed) +
  //     b[8] * (1 - ef) + (b[9]*ef_f);
  //   // d = 1
  //   eta_d_1   =  a0 +
  //     m[1]*l1 + m[2]*l2 +
  //     (b[1] * erx_srp0 + b[2] * erx_srp1 + b[3] * erx_srp2) +
  //     (b[4]*er_r_srp1 + b[5]*er_r_srp2) +
  //     // note the use of ed and rp here and not ed_rp_d_srp1, i.e. the d is
  //     // set to 1 for everyone
  //     // b[6] * (1 - ed) + 
  //     (b[6]*ed .* srp1 + b[7]*ed .* srp2)  +
  //     b[8] * (1 - ef) + (b[9]*ef_f);
  // 
  //   // choice, f = 0
  //   eta_f_0   =  a0 +
  //     m[1]*l1 + m[2]*l2 +
  //     (b[1] * erx_srp0 + b[2] * erx_srp1 + b[3] * erx_srp2) +
  //     (b[4]*er_r_srp1 + b[5]*er_r_srp2) +
  //     // b[6] * (1 - ed) + 
  //     (b[6]*ed_d_srp1 + b[7]*ed_d_srp2)  +
  //     b[8] * (1 - ef)  ;
  //   // f = 1
  //   eta_f_1   =  a0 +
  //     m[1]*l1 + m[2]*l2 +
  //     (b[1] * erx_srp0 + b[2] * erx_srp1 + b[3] * erx_srp2) +
  //     (b[4]*er_r_srp1 + b[5]*er_r_srp2) +
  //     // b[6] * (1 - ed) + 
  //     (b[6]*ed_d_srp1 + b[7]*ed_d_srp2)  +
  //     // note the use of ef here and not ef_f
  //     b[8] * (1 - ef) + (b[9]*ef);
  // 
  // }

}
