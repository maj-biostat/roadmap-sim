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
  
  array[N] int silo;
  array[N] int jnt;
  array[N] int pref;
  array[N] int d1;
  array[N] int d2;
  
  int N_ix_d1_1;
  int N_ix_d1_2;
  int N_ix_d1_3;
  
  array[N_ix_d1_1] int ix_d1_1;
  array[N_ix_d1_2] int ix_d1_2;
  array[N_ix_d1_3] int ix_d1_3;
  
  // g-comp setup  
  
  int nrd1p;
  array[nrd1p] int d1_s;
  array[nrd1p] int d1_j;
  array[nrd1p] int d1_p;
  array[nrd1p] int d1_d2;
  array[nrd1p] int nd1p;
  
  int nrd2p;
  array[nrd2p] int d2_s;
  array[nrd2p] int d2_j;
  array[nrd2p] int d2_p;
  array[nrd2p] int d2_d1;
  array[nrd2p] int nd2p;
  
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
}
transformed parameters{
  vector[K_silo] bs;
  vector[K_jnt] bj;
  vector[K_pref] bp;
  vector[K_d1] bd1;
  vector[K_d2] bd2;
  
  vector[N] eta;
  
  bs[1] = 0.0;
  bj[1] = 0.0;
  bp[1] = 0.0;
  bd1[1] = 0.0;
  bd2[1] = 0.0;
  
  bs[2:K_silo] = bs_raw;
  bj[2:K_jnt] = bj_raw;
  bp[2:K_pref] = bp_raw;
  bd1[2:K_d1] = bd1_raw;
  bd2[2:K_d2] = bd2_raw;

  // Approach this way due to identifiability issue.
  // If d1 = 1 then d2 = 1
  // If d1 = 2 then d2 = 2 or d2 = 3.
  // i.e. there is nothing to identify the parameter d1=1 from d2=1.
  // Need to have "not randomised to A/B duration" be conditional on d1 != 1.
  // For d1=1, d2 only equal to 1 therefore for d1 = 1, bd2 drops out of the 
  // linear predictor. bd2 estimated only for the case where d1 = 2 (onestage) 
  // or d1 = 3 (twostage).
  eta[ix_d1_1] = mu + bs[silo[ix_d1_1]] + bj[jnt[ix_d1_1]] + bp[pref[ix_d1_1]] + 
    bd1[d1[ix_d1_1]] ;  
  eta[ix_d1_2] = mu + bs[silo[ix_d1_2]] + bj[jnt[ix_d1_2]] + bp[pref[ix_d1_2]] + 
    bd1[d1[ix_d1_2]] + bd2[d2[ix_d1_2]];  
  eta[ix_d1_3] = mu + bs[silo[ix_d1_3]] + bj[jnt[ix_d1_3]] + bp[pref[ix_d1_3]] + 
    bd1[d1[ix_d1_3]] + bd2[d2[ix_d1_3]];  

} 
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(bs_raw | 0, 1);
  target += normal_lpdf(bj_raw | 0, 1);
  target += normal_lpdf(bp_raw | 0, 1);
  target += normal_lpdf(bd1_raw | 0, 1);
  target += normal_lpdf(bd2_raw | 0, 1);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  vector[N] wgts1 = dirichlet_rng(to_vector(n));
  vector[nrd1p] wgtsd1 = dirichlet_rng(to_vector(nd1p));
  vector[nrd2p] wgtsd2 = dirichlet_rng(to_vector(nd2p));
  
  // no bd2 for d1 = 1
  vector[nrd1p] eta_d1_1 = mu + bs[d1_s] + bj[d1_j] + bp[d1_p] + bd1[1] ; 
  // d2 is logically constrained from d1
  vector[nrd1p] eta_d1_2 = mu + bs[d1_s] + bj[d1_j] + bp[d1_p] + bd1[2] + bd2[1];
  vector[nrd1p] eta_d1_3 = mu + bs[d1_s] + bj[d1_j] + bp[d1_p] + bd1[3] + bd2[1];
  
  vector[nrd2p] eta_d2_2 = mu + bs[d2_s] + bj[d2_j] + bp[d2_p] + bd1[d2_d1] + bd2[2];
  vector[nrd2p] eta_d2_3 = mu + bs[d2_s] + bj[d2_j] + bp[d2_p] + bd1[d2_d1] + bd2[3];
  
  real mu_pop = wgts1' * eta;
  
  real bd1_1 = wgtsd1' * eta_d1_1 - mu_pop  ;
  real bd1_2 = wgtsd1' * eta_d1_2 - mu_pop  ;
  real bd1_3 = wgtsd1' * eta_d1_3 - mu_pop  ;
  
  real bd2_2 = wgtsd2' * eta_d2_2 - mu_pop  ;
  real bd2_3 = wgtsd2' * eta_d2_3 - mu_pop  ;
  
}
