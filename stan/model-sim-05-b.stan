

// group all your data up first
// instructions:
data{ 
  // full model setup
  int N;
  
  array[N] int y;
  array[N] int n;
  
  // total indexes per covariate
  int K_silo;
  int K_d1;
  
  array[N] int silo;
  array[N] int d1;
  
  // priors
  vector[2] pri_mu;
  vector[2] pri_b_silo;
  
  int prior_only;
}
transformed data{
}
parameters{
  // real mu;
  // vector[K_silo-1] bs_raw;
  matrix[K_d1, K_silo] bd1_std; 
  
  vector[K_d1] bd1_nu;                    // location of beta[ , j]
  vector<lower=0>[K_d1] bd1_tau;              // scale of beta[ , j]
  cholesky_factor_corr[K_d1] bd1_LOmega;     // Cholesky of correlation of beta[ , j]
}
transformed parameters{
  // vector[K_silo] bs;
  vector[N] eta;
  
  // for K_silo = 3 rep_matrix(bd1_nu, K_silo) will produce:
  // b1 b1 b1
  // b2 b2 b2
  // b3 b3 b3
  // ie domain parameter rows by silo cols
  // diag_pre_multiply returns the product of the diagonal matrix formed from 
  // the vector v and the matrix m, i.e., diag_matrix(v) * m.
  matrix[K_d1, K_silo] bd1 = rep_matrix(bd1_nu, K_silo)
                      + diag_pre_multiply(bd1_tau, bd1_LOmega) * bd1_std;
                      
  // matrix[K_d1, K_silo] bd1 = diag_pre_multiply(bd1_tau, bd1_LOmega) * bd1_std;                   

  // bs[1] = 0.0;
  // 
  // bs[2:K_silo] = bs_raw;

  for(i in 1:N){
    eta[i] = bd1[d1[i], silo[i]] ;   
  }
  

} 
model{
  // target += logistic_lpdf(mu | pri_mu[1], pri_mu[2]);
  // target += normal_lpdf(bs_raw | pri_b_silo[1], pri_b_silo[2]);
  
  target += normal_lpdf(bd1_nu | 0, 1);
  
  target += normal_lpdf(to_vector(bd1_std) | 0, 1);
  target += exponential_lpdf(bd1_tau | 1);
  target += lkj_corr_cholesky_lpdf(bd1_LOmega | 1);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  
}
