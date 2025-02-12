// group all your data up first
// instructions:
data{ 
  
  // full model setup
  int N;
  
  array[N] int y;
  array[N] int n;
  
  // total indexes per covariate
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
  
  // subset to specific parts of the sample in order to focus on 
  // the comparisons of interest, e.g. surgical is based on the 
  // late-acute group only and the effect of interest is obtained 
  // via g-computation
  int N_d1;
  array[N_d1] int d1_s;
  array[N_d1] int d1_j;
  array[N_d1] int d1_p;
  // assignments to each treatment level within the domain 1 cohort
  array[N_d1] int d1_d1;
  array[N_d1] int d1_d2;
  array[N_d1] int d1_d3;
  array[N_d1] int d1_d4;
  array[N_d1] int n_d1;
  int N_d1_p1;
  int N_d1_p2;
  array[N_d1_p1] int ix_d1_p1;
  array[N_d1_p2] int ix_d1_p2;
  array[N_d1_p1] int n_d1_p1;
  array[N_d1_p2] int n_d1_p2;
  
  // d2 is evaluated for those receiving one-stage revision
  int N_d2;
  array[N_d2] int d2_s;
  array[N_d2] int d2_j;
  array[N_d2] int d2_p;
  array[N_d2] int d2_d1;
  array[N_d2] int d2_d2;
  array[N_d2] int d2_d3;
  array[N_d2] int d2_d4;
  array[N_d2] int n_d2;
  
  // d3 is evaluated for those receiving two-stage revision
  int N_d3;
  array[N_d3] int d3_s;
  array[N_d3] int d3_j;
  array[N_d3] int d3_p;
  array[N_d3] int d3_d1;
  array[N_d3] int d3_d2;
  array[N_d3] int d3_d3;
  array[N_d3] int d3_d4;
  array[N_d3] int n_d3;
  
  // d4 is evaluated for those receiving two-stage revision
  int N_d4;
  array[N_d4] int d4_s;
  array[N_d4] int d4_j;
  array[N_d4] int d4_p;
  array[N_d4] int d4_d1;
  array[N_d4] int d4_d2;
  array[N_d4] int d4_d3;
  array[N_d4] int d4_d4;
  array[N_d4] int n_d4;
  
  // priors
  vector[2] pri_mu;
  vector[2] pri_b_silo;
  vector[2] pri_b_jnt;
  vector[2] pri_b_prf;
  vector[2] pri_b_trt;
  
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

} 
model{
  target += logistic_lpdf(mu | pri_mu[1], pri_mu[2]);
  target += normal_lpdf(bs_raw | pri_b_silo[1], pri_b_silo[2]);
  target += normal_lpdf(bj_raw | pri_b_jnt[1], pri_b_jnt[2]);
  target += normal_lpdf(bp_raw | pri_b_prf[1], pri_b_prf[2]);
  target += normal_lpdf(bd1_raw | pri_b_trt[1], pri_b_trt[2]);
  target += normal_lpdf(bd2_raw | pri_b_trt[1], pri_b_trt[2]);
  target += normal_lpdf(bd3_raw | pri_b_trt[1], pri_b_trt[2]);
  target += normal_lpdf(bd4_raw | pri_b_trt[1], pri_b_trt[2]);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  
  array[N] real y_pred;
  real p_hat_pred;

  for (i in 1:N) {
    y_pred[i] = bernoulli_logit_rng(eta[i]);
  }
  p_hat_pred = mean(y_pred);

  // Surgery domain (D1) comparisons of interest are revision relative to dair
  // restricted to the late acute group. Here I just compute the direct 
  // one-stage relative to dair 
  // two-stage relative to dair
  // comparisons are those are the ones that can be validated directly relative to
  // the data generation process.
  vector[N_d1] wgtsd1 = dirichlet_rng(to_vector(n_d1));  
  vector[N_d1_p1] wgtsd1_p1 = dirichlet_rng(to_vector(n_d1_p1));  
  vector[N_d1_p2] wgtsd1_p2 = dirichlet_rng(to_vector(n_d1_p2));  
  // First level of bd2 selected since this is the only comparison that 
  // applies to all. First level is reference, i.e. equals zero.
  // Similarly, first level of bd3 selected since this is the only comparison 
  // that applies to all. Again, first level is reference, i.e. equals zero.  
  vector[N_d1] eta_d1_1 = mu + bs[d1_s] + bj[d1_j] + bp[d1_p] + bd1[1] + bd2[1] + bd3[1] + bd4[d1_d4];
  // Assignment to revision will either end up being one or two-stage.
  // those that had the preference for one get one and those that had pref for
  // two get two. Thus we use different weights (since these are subsets of 
  // our late acute cohort).
  // I have set bp explicitly here but it would be ok to just use the data passed
  // in as all records should have been selected based on the required preference.
  vector[N_d1_p1] eta_d1_2 = mu + bs[d1_s[ix_d1_p1]] + bj[d1_j[ix_d1_p1]] + bp[1] + 
    bd1[2] + bd2[1] + bd3[1] + bd4[d1_d4[ix_d1_p1]];
  vector[N_d1_p2] eta_d1_3 = mu + bs[d1_s[ix_d1_p2]] + bj[d1_j[ix_d1_p2]] + bp[2] + 
    bd1[3] + bd2[1] + bd3[1] + bd4[d1_d4[ix_d1_p2]]; 
  real nu_d1_1 = wgtsd1' * eta_d1_1   ;           
  real nu_d1_2 = wgtsd1_p1' * eta_d1_2   ;           
  real nu_d1_3 = wgtsd1_p2' * eta_d1_3   ;  
  
  // effect of interest is the weight sum of the revision effects (wgts being 
  // the observed proportion of preferring each type) relative to dair
  // multiplication by 1.0 is required to cast to float.
  real prf_1 = (N_d1_p1 * 1.0/N_d1);
  real prf_2 = (N_d1_p2 * 1.0/N_d1);
  real nu_d1_23 = ((N_d1_p1 * 1.0/N_d1) * nu_d1_2) + ((N_d1_p2 * 1.0/N_d1) * nu_d1_3) ; 
  real lor_d1 = nu_d1_23 - nu_d1_1;
  // on the risk scale we would have
  real p_d1_1 = inv_logit(nu_d1_1);
  real p_d1_23 = inv_logit(nu_d1_23);
  real rd_d1 = p_d1_23 - p_d1_1;

  
  // AB duration domain (D2) comparisons of interest are 6 wks relative to 12 wks
  // for those that have one-stage revision
  vector[N_d2] wgtsd2 = dirichlet_rng(to_vector(n_d2)); 
  // In the data passed to stan, d1 needs to be set to 2 here for all units 
  // (corresponding to one-stage), similarly d3 (ext-proph) needs to be 1 here 
  // for all units (corresponding to non-rand intervention). This is so that we
  // are comparing apples with apples.
  vector[N_d2] eta_d2_2 = mu + bs[d2_s] + bj[d2_j] + bp[d2_p] + bd1[d2_d1] + bd2[2] + bd3[d2_d3] + bd4[d2_d4];
  vector[N_d2] eta_d2_3 = mu + bs[d2_s] + bj[d2_j] + bp[d2_p] + bd1[d2_d1] + bd2[3] + bd3[d2_d3] + bd4[d2_d4];    
  real nu_d2_2 = wgtsd2' * eta_d2_2   ;           
  real nu_d2_3 = wgtsd2' * eta_d2_3   ; 
  // only one comparison of interest         
  real lor_d2 = nu_d2_3 - nu_d2_2;  
  real p_d2_2 = inv_logit(nu_d2_2);
  real p_d2_3 = inv_logit(nu_d2_3);
  real rd_d2 = p_d2_3 - p_d2_2;
  
  // Ext prophylaxis domain (D3) comparisons of interest are 12 wks relative to none
  // for those that have two-stage revision
  vector[N_d3] wgtsd3 = dirichlet_rng(to_vector(n_d3)); 
  // In the data passed to stan, d1 needs to be set to 3 here for all units 
  // (corresponding to two-stage), similarly d3 needs to be 1 here for all units 
  // (corresponding to non-rand). This is so that we are comparing apples with apples.
  vector[N_d3] eta_d3_2 = mu + bs[d3_s] + bj[d3_j] + bp[d3_p] + bd1[d3_d1] + bd2[d3_d2] + bd3[2] + bd4[d3_d4];
  vector[N_d3] eta_d3_3 = mu + bs[d3_s] + bj[d3_j] + bp[d3_p] + bd1[d3_d1] + bd2[d3_d2] + bd3[3] + bd4[d3_d4];  
  real nu_d3_2 = wgtsd3' * eta_d3_2   ;           
  real nu_d3_3 = wgtsd3' * eta_d3_3   ; 
  // only one comparison of interest     
  real lor_d3 = nu_d3_3 - nu_d3_2;  
  real p_d3_2 = inv_logit(nu_d3_2);
  real p_d3_3 = inv_logit(nu_d3_3);
  real rd_d3 = p_d3_3 - p_d3_2;
   
  // Ab choice domain (D4) comparisons of interest are rif to none
  vector[N_d4] wgtsd4 = dirichlet_rng(to_vector(n_d4)); 
  vector[N_d4] eta_d4_2 = mu + bs[d4_s] + bj[d4_j] + bp[d4_p] + bd1[d4_d1] + bd2[d4_d2] + bd3[d4_d3] + bd4[2];
  vector[N_d4] eta_d4_3 = mu + bs[d4_s] + bj[d4_j] + bp[d4_p] + bd1[d4_d1] + bd2[d4_d2] + bd3[d4_d3] + bd4[3];    
  real nu_d4_2 = wgtsd4' * eta_d4_2   ;           
  real nu_d4_3 = wgtsd4' * eta_d4_3   ; 
  // only one comparison of interest   
  real lor_d4 = nu_d4_3 - nu_d4_2;  
  real p_d4_2 = inv_logit(nu_d4_2);
  real p_d4_3 = inv_logit(nu_d4_3);
  real rd_d4 = p_d4_3 - p_d4_2;
  
  
}
