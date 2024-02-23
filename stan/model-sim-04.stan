

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
  target += normal_lpdf(m | 0, 1);
  target += normal_lpdf(b | 0, 1);
  
  // likelihood chunks pertaining to each silo
  target += binomial_logit_lpmf(y | n, eta) ;    
  
}
generated quantities{
  
  // real b_r;
  // real b_d1;
  // real b_d2;
  // real b_f;
  vector[N] eta_r_0;
  vector[N] eta_r_1;
  vector[N] eta_d1_0;
  vector[N] eta_d1_1;
  vector[N] eta_d2_0;
  vector[N] eta_d2_1;
  vector[N] eta_f_0;
  vector[N] eta_f_1;

  {
    vector[sum(n)] delta;
    int is = 1;
    int ie = n[1];
    
    // vector[N] eta_r1;
    // vector[N] eta_d0_1;
    // vector[N] eta_d1_1;
    // vector[N] eta_d0_2;
    // vector[N] eta_d1_2;
    // vector[N] eta_f0;
    // vector[N] eta_f1;
       
    // r = 0
    eta_r0   =  a0 +
      m[1]*l1   + m[2]*l2 +
      b[1] * (1 - er  ) + 
      b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   +
      b[7] * (1 - ef  ) + (b[8]*ef_f)  ;
      
    eta_r1   =  a0 +
      m[1]*l1   + m[2]*l2 +
      b[1] * (1 - er  ) + (b[2]*er + b[3]*er .* srp2) +
      b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   +
      b[7] * (1 - ef  ) + (b[8]*ef_f);
      
      
    
    // for(i in 1:N){
    //   // print(i, ", ", is, ", ", ie, ", ", eta_r1[i] - eta_r0[i], ", ", to_int(n[i]));
    //   delta[is:ie] = rep_vector(eta_r1[i] - eta_r0[i], to_int(n[i]));
    //   is = ie + 1;
    //   if(i < N){
    //     ie = sum(n[1:(i+1)]);  
    //   }
    // }
 
    // r = 1
    // vector[N] er_r_1_srp2 = er .* srp2;
    // eta_r1   =  a0 + 
    //   m[1]*l1   + m[1]*l2   + m[3]*j   + 
    //   m[4]*l1_j   + m[5]*l2_j   +
    //   // we use er instead of er_r b'coz all that were randomised into the 
    //   // surg domain are being set to have received rev
    //   // note that we use er .* srp2 instead of er_r_srp2
    //   b[1] * (1 - er  ) + (b[2]*er + b[3]*er .* srp2) +
    //   b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   + 
    //   b[7] * (1 - ef  ) + (b[8]*ef_f)  ;

    // d = 0, srp = 0 (long)
    // eta_d0_1   = a0 + 
    //   m[1]*l1   + m[1]*l2   + m[3]*j   + 
    //   m[4]*l1_j   + m[5]*l2_j   + 
    //   b[1] * (1 - er  ) + (b[2]*er_r)   +
    //   b[4] * (1 - ed  ) + 
    //   b[7] * (1 - ef  ) + (b[8]*ef_f)  ;
    //   
    // // d = 1, srp = 0 (short)
    // eta_d1_1   = a0 + 
    //   m[1]*l1   + m[1]*l2   + m[3]*j   + 
    //   m[4]*l1_j   + m[5]*l2_j   + 
    //   b[1] * (1 - er  ) + (b[2]*er_r)   +
    //   b[4] * (1 - ed  ) + (b[5]*ed_rp_d)    + 
    //   b[7] * (1 - ef  ) + (b[8]*ef_f)  ;
    // 
    // // d = 0, srp = 1 (long)
    // eta_d0_2   = a0 + 
    //   m[1]*l1   + m[1]*l2   + m[3]*j   + 
    //   m[4]*l1_j   + m[5]*l2_j   + 
    //   b[1] * (1 - er  ) + (b[2]*er_r + b[3]*er_r_srp2) +
    //   b[4] * (1 - ed  ) + 
    //   b[7] * (1 - ef  ) + (b[8]*ef_f)  ;
    //   
    // // d = 1, srp = 1 (short)
    // eta_d1_2   = a0 + 
    //   m[1]*l1   + m[1]*l2   + m[3]*j   + 
    //   m[4]*l1_j   + m[5]*l2_j   + 
    //   b[1] * (1 - er  ) + (b[2]*er_r + b[3]*er_r_srp2) +
    //   b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   + 
    //   b[7] * (1 - ef  ) + (b[8]*ef_f)  ;
    // 
    // // f = 0
    // eta_f0   = a0 + 
    //   m[1]*l1   + m[1]*l2   + m[3]*j   + 
    //   m[4]*l1_j   + m[5]*l2_j   + 
    //   b[1] * (1 - er  ) + (b[2]*er_r + b[3]*er_r_srp2) +
    //   b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   + 
    //   b[7] * (1 - ef  ) ;
    //   
    // // f = 1
    // eta_f1   = a0 + 
    //   m[1]*l1   + m[1]*l2   + m[3]*j   + 
    //   m[4]*l1_j   + m[5]*l2_j   + 
    //   b[1] * (1 - er  ) + (b[2]*er_r + b[3]*er_r_srp2) +
    //   b[4] * (1 - ed  ) + (b[5]*ed_rp_d   + b[6]*ed_rp_d_srp2)   + 
    //   b[7] * (1 - ef  ) + (b[8]*ef_f)  ;
    
    // b_r = mean(eta_r1 - eta_r0);
    // b_d1 = mean(eta_d1_1 - eta_d0_1);
    // b_d2 = mean(eta_d1_2 - eta_d0_2);
    // b_f = mean(eta_f1 - eta_f0);
    
    // b_r = sum( (eta_r1 - eta_r0) .* to_vector(n)  ) / sum(n);
    
  }
  
  
  
  

}
