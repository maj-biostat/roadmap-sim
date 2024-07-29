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
  // domain 2 
  int nrXd2;
  int ncXd2;
  matrix[nrXd2, ncXd2] Xd2des;
  vector[ncXd2] sd2;
  
  // g-comp setup  
  int nrd1p;
  int ncd1p;
  matrix[nrd1p, ncd1p] Xd1p1;
  matrix[nrd1p, ncd1p] Xd1p2;
  array[nrd1p] int nd1p;
  // 
  int nrd2p;
  int ncd2p;
  matrix[nrd2p, ncd2p] Xd2p1;
  matrix[nrd2p, ncd2p] Xd2p2;
  array[nrd2p] int nd2p;
  
  int prior_only;
}
transformed data{
  
  // build full design matrices
  matrix[N, ncXs] Xs = Xsdes[silo,];
  matrix[N, ncXj] Xj = Xjdes[jnt,];
  matrix[N, ncXj] Xp = Xpdes[pref,];
  matrix[N, ncXd1] Xd1 = Xd1des[d1,];
  matrix[N, ncXd2] Xd2 = Xd2des[d2,];

}
parameters{
  real mu;
  vector[ncXs] bs;
  vector[ncXj] bj;
  vector[ncXp] bp;
  vector[ncXd1] bd1;
  vector[ncXd2] bd2;
}
transformed parameters{
  vector[N] eta = mu + Xs*bs + Xj*bj + Xp*bp + Xd1*bd1 + Xd2*bd2 ;
  // vector[N] eta = mu + Xs*bs + Xj*bj + Xp*bp + Xd1*bd1 + Xd2*bd2  ;
  // vector[ncXs + ncXj + ncXp + ncXd1 + ncXd2] b;
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
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  vector[N] wgts1 = dirichlet_rng(to_vector(n));
  vector[nrd1p] wgtsd1 = dirichlet_rng(to_vector(nd1p));
  vector[nrd2p] wgtsd2 = dirichlet_rng(to_vector(nd2p));
  // 
  vector[nrd1p] eta_d1_1 = mu + Xd1p1 * b;
  vector[nrd1p] eta_d1_2 = mu + Xd1p2 * b;
  vector[nrd1p] eta_d1_3 = mu + Xd1p2 * b;
  // 
  vector[nrd2p] eta_d2_2 = mu + Xd2p1 * b;
  vector[nrd2p] eta_d2_3 = mu + Xd2p2 * b;
  // 
  real mu_pop = wgts1' * eta;
  // 
  real bd1_1 = wgtsd1' * eta_d1_1 - mu_pop  ;
  real bd1_2 = wgtsd1' * eta_d1_2 - mu_pop  ;
  real bd1_3 = wgtsd1' * eta_d1_3 - mu_pop  ;
  // 
  real bd2_2 = wgtsd2' * eta_d2_2 - mu_pop  ;
  real bd2_3 = wgtsd2' * eta_d2_3 - mu_pop  ;

}
