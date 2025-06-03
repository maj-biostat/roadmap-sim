// group all your data up first
// instructions:
data{ 
  
  int N;
  
  array[N] int y;
  array[N] int n;
  
  // Parameter indexes
  array[N] int s;
  array[N] int pref;
  // silo x intv - converted via K_d1 and silo to index b/w 1:9, see d1_ix
  array[N] int d1;
  // intv + non-rand 
  array[N] int d2;
  array[N] int d3;
  array[N] int d4;
  
  // Total indexes per covariate, e.g. number of silos, d1 interventions etc
  int K_s;
  int K_p;
  int K_d1;
  int K_d2;
  int K_d3;
  int K_d4;
  
  // G-comp related subsets
  
  // subset to specific parts of the sample in order to focus on 
  // the comparisons of interest, e.g. surgical is based on the 
  // late-acute group only and the effect of interest is obtained 
  // via g-computation
  int N_d1;
  array[N_d1] int d1_s;
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
  real prop_p1; // proportion with preference towards one-stage
  real prop_p2; // proportion with preference towards two-stage
  
  // d2 is evaluated for those receiving one-stage revision
  int N_d2;
  array[N_d2] int d2_s;
  array[N_d2] int d2_p;
  array[N_d2] int d2_d1;
  array[N_d2] int d2_d2;
  array[N_d2] int d2_d3;
  array[N_d2] int d2_d4;
  array[N_d2] int n_d2;
  
  // d3 is evaluated for those receiving two-stage revision
  int N_d3;
  array[N_d3] int d3_s;
  array[N_d3] int d3_p;
  array[N_d3] int d3_d1;
  array[N_d3] int d3_d2;
  array[N_d3] int d3_d3;
  array[N_d3] int d3_d4;
  array[N_d3] int n_d3;
  
  // d4 is evaluated for those receiving two-stage revision
  int N_d4;
  array[N_d4] int d4_s;
  array[N_d4] int d4_p;
  array[N_d4] int d4_d1;
  array[N_d4] int d4_d2;
  array[N_d4] int d4_d3;
  array[N_d4] int d4_d4;
  array[N_d4] int n_d4;
  
  // priors
  vector[2] pri_mu;
  vector[2] pri_bs;
  vector[2] pri_bp;
  vector[2] pri_b1;
  vector[2] pri_b2;
  vector[2] pri_b3;
  vector[2] pri_b4;
  
  int prior_only;
}
transformed data{
  array[N] int d1_ix;
  
  array[N_d2] int d2_d1_ix;
  array[N_d3] int d3_d1_ix;
  array[N_d4] int d4_d1_ix;
  
  
  for(i in 1:N){
    d1_ix[i] = d1[i] + (K_d1 * (s[i] - 1));  
  } 
  // for g-comp to pick up the correct surgical domain parameter
  for(i in 1:N_d2){
    // d2_d1 should be one-stage (2) for everything since we are conditioning on 
    // one-stage but this will convert to a silo specific index for one-stage
    d2_d1_ix[i] = d2_d1[i] + (K_d1 * (d2_s[i] - 1));  
  }
  for(i in 1:N_d3){
    // As above but for two-stage. 
    // d3_d1 should be two-stage (2) for everything since we are conditioning on 
    // two-stage but this will convert to a silo specific index for two-stage
    d3_d1_ix[i] = d3_d1[i] + (K_d1 * (d3_s[i] - 1));  
  }
  for(i in 1:N_d4){
    d4_d1_ix[i] = d4_d1[i] + (K_d1 * (d4_s[i] - 1));  
  }
}
parameters{
  real mu;
  vector[K_s-1] bs_raw;
  vector[K_p-1] bp_raw;
  vector[(K_d1 * K_s) - 1] bd1_raw;
  vector[K_d2-1] bd2_raw;
  vector[K_d3-1] bd3_raw;
  vector[K_d4-1] bd4_raw;
}
transformed parameters{
  vector[K_s] bs;
  vector[K_p] bp;
  // The surgical domain needs a parameter to account for the fact that the
  // non randomised comparisons of dair/rev(1)/rev(2) are correctly reflected
  // in the linear predictor for the duration domains.
  // For example, suppose that there is no effect of revision in the early silo
  // (for whatever reason) but in the late acute group (our randomised comparison
  // for the surgical domain) there is an effect of revision. 
  // We account for the non-randomised entry into the surgical domain via the
  // silo parameters as these are identical to an indicator for non-randomised
  // surgical intervention. 
  // Here we try to account for silo specific surgical domain effects (even if
  // simply due to the non-randomised nature of the comparison for early and 
  // chronic) and their influence on the duration domains.
  // Below I declare the silo by domain effects for surgical intervention.
  vector[K_d1 * K_s] bd1;
  vector[K_d2] bd2;
  vector[K_d3] bd3;
  vector[K_d4] bd4;
  
  vector[N] eta;
  
  bs[1] = 0.0;
  bp[1] = 0.0;
  bd1[1] = 0.0;
  bd2[1] = 0.0;
  bd3[1] = 0.0;
  bd4[1] = 0.0;
  
  bs[2:K_s] = bs_raw;
  bp[2:K_p] = bp_raw;
  bd1[2:(K_d1 * K_s)] = bd1_raw;
  bd2[2:K_d2] = bd2_raw;
  bd3[2:K_d3] = bd3_raw;
  bd4[2:K_d4] = bd4_raw;
  
  for(i in 1:N){
    // dair
    if(d1[i] == 1){
      // exclude both bd2 and bd3    
      eta[i] = mu + bs[s[i]] + bp[pref[i]] + bd1[d1_ix[i]] + bd4[d4[i]]; 
    } else if (d1[i] == 2){
      // rev(1)
      eta[i] = mu + bs[s[i]] + bp[pref[i]] + bd1[d1_ix[i]] + bd2[d2[i]] + bd4[d4[i]];  
    } else {
      // rev(2)
      eta[i] = mu + bs[s[i]] + bp[pref[i]] + bd1[d1_ix[i]] + bd3[d3[i]] + bd4[d4[i]];  
    }
  }
  

} 
model{
  target += logistic_lpdf(mu | pri_mu[1], pri_mu[2]);
  target += normal_lpdf(bs_raw | pri_bs[1], pri_bs[2]);
  target += normal_lpdf(bp_raw | pri_bp[1], pri_bp[2]);
  target += normal_lpdf(bd1_raw | pri_b1[1], pri_b1[2]);
  target += normal_lpdf(bd2_raw | pri_b2[1], pri_b2[2]);
  target += normal_lpdf(bd3_raw | pri_b3[1], pri_b3[2]);
  target += normal_lpdf(bd4_raw | pri_b4[1], pri_b4[2]);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  
  
  // Surgery domain (D1) comparisons of interest are revision relative to dair
  // restricted to the late acute group. 
  vector[N_d1] wgtsd1 = dirichlet_rng(to_vector(n_d1));
  vector[N_d1_p1] wgtsd1_p1 = dirichlet_rng(to_vector(n_d1_p1));
  vector[N_d1_p2] wgtsd1_p2 = dirichlet_rng(to_vector(n_d1_p2));
  
  // For dair and rev(2), d2 is undefined and omitted from model.
  // For dair and rev(1), d3 is undefined and omitted from model.
  
  // The comparison is dair (with whatever is recvd for d2 backbone abx) and
  // by defn no d3 (ext proph) vs rev.
  // Rev is decomposed into rev(1) and rev(2).
  // We take rev(1) as including d2 (non-randomised) and d3 is undef
  // We take rev(2) as including d3 (non-randomised) and d2 is undef

  // An alternative might be to consider bd2 and bd3 when given as usual.
  // However, this will induce effects in d1 even if d1 shows no effect but
  // either of the duration domains has an effect.
  
  vector[N_d1] mu_d1_1 = inv_logit(mu + bs[2] + bd1[4] + bp[d1_p] + bd4[d1_d4]);
  
  // Assignment to revision will either end up being one or two-stage.
  // those that had the preference for one get one and those that had pref for
  // two get two. Thus we use different weights (since these are subsets of
  // our late acute cohort).
  // I have set bp explicitly here but it would be ok to just use the data passed
  // in as all records should have been selected based on the required preference.
  vector[N_d1_p1] mu_d1_2 = inv_logit(mu + bs[2] + bp[1] + bd1[5] + bd2[1] + bd4[d1_d4[ix_d1_p1]]);
  vector[N_d1_p2] mu_d1_3 = inv_logit(mu + bs[2] + bp[2] + bd1[6] + bd3[1] + bd4[d1_d4[ix_d1_p2]]);
  real p_d1_1 = wgtsd1' * mu_d1_1   ;
  real p_d1_2 = wgtsd1_p1' * mu_d1_2   ;
  real p_d1_3 = wgtsd1_p2' * mu_d1_3   ;

  // our interpretation of revision:
  real p_d1_23 = (prop_p1 * p_d1_2) + (prop_p2 * p_d1_3) ;
  real rd_d1 = p_d1_23 - p_d1_1;
  real lor_d1 = log((p_d1_23 * (1-p_d1_1)) / ( (1 - p_d1_23) * p_d1_1));

  // AB duration domain (D2) comparisons of interest are 6 wks relative to 12 wks
  // for those that have one-stage revision
  vector[N_d2] wgtsd2 = dirichlet_rng(to_vector(n_d2));
  // In the data passed to stan, d1 needs to be set to one-stage revision
  // for all units. However, this is now a silo specific view to account for the
  // possibility of differential surgical effects (for whatever reason).
  // So, we index the bd1 using d2_d1_ix to pick up either the non-rand or rand
  // comparison from the surg domain.
  // Additionally, d3 (ext-proph) is undefined for all units here and so is 
  // omitted.
  // If we do not do these conditioning steps, then we would not be making
  // logical comparisons based on the design constraints, e.g. if you have
  // one stage revision, then it is logically impossible to receive ext proph
  // based on the design rules.
  vector[N_d2] mu_d2_2 = inv_logit(mu + bs[d2_s] + bp[d2_p] + bd1[d2_d1_ix] + bd2[2] + bd4[d2_d4] );
  vector[N_d2] mu_d2_3 = inv_logit(mu + bs[d2_s] + bp[d2_p] + bd1[d2_d1_ix] + bd2[3] + bd4[d2_d4] );
  real p_d2_2 = wgtsd2' * mu_d2_2 ;
  real p_d2_3 = wgtsd2' * mu_d2_3 ;

  real rd_d2 = p_d2_3 - p_d2_2;
  real lor_d2 = log((p_d2_3 * (1-p_d2_2)) / ( (1 - p_d2_3) * p_d2_2));

  // Ext prophylaxis domain (D3) comparisons of interest are 12 wks relative to none
  // for those that have two-stage revision
  vector[N_d3] wgtsd3 = dirichlet_rng(to_vector(n_d3));
  // In the data passed to stan, d1 needs to be set to two-stage revision
  // for all units but now this is a silo specific view so we index the bd1
  // parameter to pick up either the rand or non-rand comparison from the surg
  // domain.
  // Similarly, d2 is now underfined so omitted.
  vector[N_d3] mu_d3_2 = inv_logit(mu + bs[d3_s] + bp[d3_p] + bd1[d3_d1_ix] + bd3[2] + bd4[d3_d4]);
  vector[N_d3] mu_d3_3 = inv_logit(mu + bs[d3_s] + bp[d3_p] + bd1[d3_d1_ix] + bd3[3] + bd4[d3_d4]);
  real p_d3_2 = wgtsd3' * mu_d3_2 ;
  real p_d3_3 = wgtsd3' * mu_d3_3 ;

  real rd_d3 = p_d3_3 - p_d3_2;
  real lor_d3 = log((p_d3_3 * (1-p_d3_2)) / ( (1 - p_d3_3) * p_d3_2));

  // Ab choice domain (D4) comparisons of interest are rif to none
  vector[N_d4] wgtsd4 = dirichlet_rng(to_vector(n_d4));
  vector[N_d4] mu_d4_2;
  vector[N_d4] mu_d4_3;

  for(i in 1:N_d4){
    
    if(d4_d1[i] == 1){
      // exclude both bd2 and bd3
      mu_d4_2[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd4[2]);
      mu_d4_3[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd4[3]);
      
    } else if (d4_d1[i] == 2){
      // exclude the bd3 term
      mu_d4_2[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd2[d4_d2[i]] + bd4[2]);
      mu_d4_3[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd2[d4_d2[i]] + bd4[3]);
      
    } else {  // d4_d1 is 3
      // exclude the bd2 term
      mu_d4_2[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd3[d4_d3[i]] + bd4[2]);
      mu_d4_3[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd3[d4_d3[i]] + bd4[3]);
      
    }
    
  }
  real p_d4_2 = wgtsd4' * mu_d4_2   ;
  real p_d4_3 = wgtsd4' * mu_d4_3   ;

  real rd_d4 = p_d4_3 - p_d4_2;
  real lor_d4 = log((p_d4_3 * (1-p_d4_2)) / ( (1 - p_d4_3) * p_d4_2));
  
  
}
