// group all your data up first
// instructions:
data{ 
  int N;
  
  array[N] int y;
  array[N] int n;
  
  array[N] int silo;
  array[N] int jnt;
  array[N] int pref;
  array[N] int d1;
  array[N] int d2;
  // array[N] int g1;

  // design matrices
  // silo
  int nrXs;
  int ncXs;
  matrix[nrXs, ncXs] Xsdes;
  vector[ncXs] ss;
  // joint
  int nrXj;
  int ncXj;
  matrix[nrXj, ncXj] Xjdes;
  vector[ncXj] sj;
  // prev rev
  int nrXp;
  int ncXp;
  matrix[nrXp, ncXp] Xpdes;
  vector[ncXp] sp;
  // domain 1 
  int nrXd1;
  int ncXd1;
  matrix[nrXd1, ncXd1] Xd1des;
  vector[ncXd1] sd1;
  // silo by domain 1 interaction - 
  // target silo specific effects for surgery not an average effect
  // int ncXg1;
  // int nrXg1;
  // matrix[nrXg1, ncXg1] Xg1des;
  // vector[ncXg1] sg1;
  
  // domain 2 
  int nrXd2;
  int ncXd2;
  matrix[nrXd2, ncXd2] Xd2des;
  vector[ncXd2] sd2;
  
  int nrp;
  int ncp;
  matrix[nrp, ncp] Xp1;
  matrix[nrp, ncp] Xp2;
  array[nrp] int nd2;
  
  int prior_only;
}
transformed data{
  
  // build full design matrices
  matrix[N, ncXs] Xs = Xsdes[silo,];
  matrix[N, ncXj] Xj = Xjdes[jnt,];
  matrix[N, ncXj] Xp = Xpdes[pref,];
  matrix[N, ncXd1] Xd1 = Xd1des[d1,];
  // indexes the relevant col in the design matrix
  // matrix[N, ncXg1] Xg1 = Xg1des[g1,];
  matrix[N, ncXd2] Xd2 = Xd2des[d2,];

}
parameters{
  real mu;
  vector[ncXs] bs;
  vector[ncXj] bj;
  vector[ncXp] bp;
  vector[ncXd1] bd1;
  vector[ncXd2] bd2;
  // vector[ncXg1] bg1;
}
transformed parameters{
  vector[N] eta = mu + Xs*bs + Xj*bj + Xp*bp + Xd1*bd1 + Xd2*bd2  ;
  vector[ncXs + ncXj + ncXp + ncXd1 + ncXd2] b;
  b[1:ncXs] = bs;
  b[(ncXs + 1):(ncXs+ncXj)] = bj;
  b[(ncXs+ncXj + 1):(ncXs+ncXj+ncXp)] = bp;
  b[(ncXs+ncXj+ncXp + 1):(ncXs+ncXj+ncXp+ncXd1)] = bd1;
  b[(ncXs+ncXj+ncXp+ncXd1 + 1):(ncXs+ncXj+ncXp+ncXd1+ncXd2)] = bd2;
  
} 
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(bs | 0, ss);
  target += normal_lpdf(bj | 0, sj);
  target += normal_lpdf(bp | 0, sp);
  target += normal_lpdf(bd1 | 0, sd1);
  target += normal_lpdf(bd2 | 0, sd2);
  
  // target += normal_lpdf(bg1 | 0, sg1);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  vector[N] wgts1 = dirichlet_rng(to_vector(n));  
  vector[nrp] wgts2 = dirichlet_rng(to_vector(nd2));     
  vector[nrp] eta_d2 = mu + Xp1 * b;
  vector[nrp] eta_d3 = mu + Xp2 * b;
  
  real mu_pop = wgts1' * eta;
  real bd2_2 = wgts2' * eta_d2 - mu_pop  ;
  real bd2_3 = wgts2' * eta_d3 - mu_pop  ;  

}
