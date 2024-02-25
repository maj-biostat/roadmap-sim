

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
  vector[N] rp; // revision was recvd
  vector[N] srp2; // two-stage revision was recvd
  // choice
  vector[N] f;
  
  vector[2] pri_m_sd;
  vector[8] pri_b_sd;
  int prior_only;
}
transformed data{
  vector[N] er_r;
  vector[N] er_r_srp2;
  vector[N] ed_rp_d;
  vector[N] ed_rp_d_srp2;
  vector[N] ef_f;
  
  er_r = er .* r;
  er_r_srp2 = er_r .* srp2;
  
  ed_rp_d = ed .* rp .* d;
  ed_rp_d_srp2 = ed_rp_d .* srp2;
  
  ef_f = ef .* f;
}
parameters {
  real a0;
  vector[2] m;
  vector[8] b;
}
transformed parameters{
  vector[N] eta;
 
  eta = a0 + 
      m[1]*l1 + m[2]*l2 +  
      b[1] * (1 - er) + (b[2]*er_r + b[3]*er_r_srp2) +
      b[4] * (1 - ed) + (b[5]*ed_rp_d + b[6]*ed_rp_d_srp2)  + 
      b[7] * (1 - ef) + (b[8]*ef_f);
      
}
model{
  target += normal_lpdf(a0 | 0, 1.5);
  target += normal_lpdf(m | 0, pri_m_sd);
  target += normal_lpdf(b | 0, pri_b_sd);
  
  // likelihood chunks pertaining to each silo
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta) ;      
  }
  
  
}
generated quantities{

  vector[N] eta_r_0;
  vector[N] eta_r_1;
  vector[N] eta_d_0;
  vector[N] eta_d_1;
  vector[N] eta_f_0;
  vector[N] eta_f_1;

  {
 
    // r = 0
    eta_r_0   =  a0 +
      m[1]*l1   + m[2]*l2 +
      b[1] * (1 - er  ) + 
      b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   +
      b[7] * (1 - ef  ) + (b[8]*ef_f)  ;
      
    eta_r_1   =  a0 +
      m[1]*l1   + m[2]*l2 +
      // note the use of er here and not er_r
      b[1] * (1 - er  ) + (b[2]*er + b[3]*er .* srp2) +
      b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   +
      b[7] * (1 - ef  ) + (b[8]*ef_f);

    // d
    eta_d_0   =  a0 + 
      m[1]*l1 + m[2]*l2 +  
      b[1] * (1 - er) + (b[2]*er_r + b[3]*er_r_srp2) +
      b[4] * (1 - ed) + 
      b[7] * (1 - ef) + (b[8]*ef_f);

    eta_d_1   =  a0 + 
      m[1]*l1 + m[2]*l2 +  
      b[1] * (1 - er) + (b[2]*er_r + b[3]*er_r_srp2) +
      // note the use of ed and rp here and not ed_rp_d
      b[4] * (1 - ed) + (b[5]*ed .* rp + b[6]*ed .* rp .* srp2)  + 
      b[7] * (1 - ef) + (b[8]*ef_f);

    // f
    eta_f_0   =  a0 +
      m[1]*l1 + m[2]*l2 +
      b[1] * (1 - er) + (b[2]*er_r + b[3]*er_r_srp2) +
      b[4] * (1 - ed) + (b[5]*ed_rp_d + b[6]*ed_rp_d_srp2)  +
      b[7] * (1 - ef)  ;

    eta_f_1   =  a0 +
      m[1]*l1 + m[2]*l2 +
      b[1] * (1 - er) + (b[2]*er_r + b[3]*er_r_srp2) +
      b[4] * (1 - ed) + (b[5]*ed_rp_d + b[6]*ed_rp_d_srp2)  +
      // note the use of ef here and not ef_f
      b[7] * (1 - ef) + (b[8]*ef);
      
    // b_r = sum( (eta_r1 - eta_r0) .* to_vector(n)  ) / sum(n);
    
  }
  
  
  
  

}
