// group all your data up first
// instructions:
data{ 
  int N;
  
  array[N] int y;
  array[N] int n;
  
  int K_silo;
  int K_jnt;
  int K_pref;
  int K_d1;
  int K_d2;
  int K_d3;
  int K_d4;
  
  array[N] int silo;
  array[N] int jnt;
  array[N] int pref;
  array[N] int d1;
  array[N] int d2;
  array[N] int d3;
  array[N] int d4;
  
  // int N_ix_d1_1;
  // int N_ix_d1_2;
  // int N_ix_d1_3;
  // 
  // array[N_ix_d1_1] int ix_d1_1;
  // array[N_ix_d1_2] int ix_d1_2;
  // array[N_ix_d1_3] int ix_d1_3;
  
  // g-comp setup  
  
  // int nrd1p;
  // array[nrd1p] int d1_s;
  // array[nrd1p] int d1_j;
  // array[nrd1p] int d1_p;
  // array[nrd1p] int d1_d2;
  // array[nrd1p] int d1_d3;
  // array[nrd1p] int d1_d4;
  // array[nrd1p] int nd1p;
  // 
  // int nrd2p;
  // array[nrd2p] int d2_s;
  // array[nrd2p] int d2_j;
  // array[nrd2p] int d2_p;
  // array[nrd2p] int d2_d1;
  // array[nrd2p] int nd2p;
  
  int prior_only;
}
transformed data{
}
parameters{
  real mu;
  vector[K_silo-1] bs_raw;
  vector[K_jnt-1] bj_raw;
  vector[K_pref-1] bp_raw;
  vector[K_d1-1] bd1_raw;
  vector[K_d2-1] bd2_raw;
  vector[K_d3-1] bd3_raw;
  vector[K_d4-1] bd4_raw;
}
transformed parameters{
  vector[K_silo] bs;
  vector[K_jnt] bj;
  vector[K_pref] bp;
  vector[K_d1] bd1;
  vector[K_d2] bd2;
  vector[K_d3] bd3;
  vector[K_d4] bd4;
  
  vector[N] eta;
  
  bs[1] = 0.0;
  bj[1] = 0.0;
  bp[1] = 0.0;
  bd1[1] = 0.0;
  bd2[1] = 0.0;
  bd3[1] = 0.0;
  bd4[1] = 0.0;
  
  bs[2:K_silo] = bs_raw;
  bj[2:K_jnt] = bj_raw;
  bp[2:K_pref] = bp_raw;
  bd1[2:K_d1] = bd1_raw;
  bd2[2:K_d2] = bd2_raw;
  bd3[2:K_d3] = bd3_raw;
  bd4[2:K_d4] = bd4_raw;

  eta = mu + bs[silo] + bj[jnt] + bp[pref] + 
    bd1[d1] + bd2[d2] + bd3[d3] + bd4[d4];  
    
  // eta[ix_d1_2] = mu + bs[silo[ix_d1_2]] + bj[jnt[ix_d1_2]] + bp[pref[ix_d1_2]] + 
  //   bd1[d1[ix_d1_2]] + bd2[d2[ix_d1_2]]; 
  //   
  // eta[ix_d1_3] = mu + bs[silo[ix_d1_3]] + bj[jnt[ix_d1_3]] + bp[pref[ix_d1_3]] + 
  //   bd1[d1[ix_d1_3]] + bd2[d2[ix_d1_3]];  

} 
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(bs_raw | 0, 1);
  target += normal_lpdf(bj_raw | 0, 1);
  target += normal_lpdf(bp_raw | 0, 1);
  target += normal_lpdf(bd1_raw | 0, 1);
  target += normal_lpdf(bd2_raw | 0, 1);
  target += normal_lpdf(bd3_raw | 0, 1);
  target += normal_lpdf(bd4_raw | 0, 1);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  vector[K_silo-1] bd1_delta;
  real bd2_delta;
  real bd3_delta;
  real bd4_delta;

  bd1_delta[1] = bd1[2] - bd1[1];
  bd1_delta[2] = bd1[3] - bd1[1];
  bd2_delta = bd2[3] - bd2[2];
  bd3_delta = bd3[3] - bd3[2];
  bd4_delta = bd4[3] - bd4[2];

  // vector[nrd1p] wgtsd1 = dirichlet_rng(to_vector(nd1p));

  // bd2 is at reference level for all groups to attempt to get a fair comparison
  // bd3 cannot be included for any  again to try to get a fair comparison
  // vector[nrd1p] eta_d1_1 = mu + bs[d1_s] + bj[d1_j] + bp[d1_p] + bd1[1] + bd4[d1_d4] ;
  // vector[nrd1p] eta_d1_2 = mu + bs[d1_s] + bj[d1_j] + bp[d1_p] + bd1[2] + bd4[d1_d4] ;
  // vector[nrd1p] eta_d1_3 = mu + bs[d1_s] + bj[d1_j] + bp[d1_p] + bd1[3] + bd4[d1_d4] ;
  // 
  // real bd1_1 = wgtsd1' * eta_d1_1   ;
  // real bd1_2 = wgtsd1' * eta_d1_2   ;
  // real bd1_3 = wgtsd1' * eta_d1_3   ;
  // 
  // vector[K_silo-1] bd1_gamma;
  // 
  // bd1_gamma[1] = bd1_2 - bd1_1;
  // bd1_gamma[2] = bd1_3 - bd1_1;
  
}
