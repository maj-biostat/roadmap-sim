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
  array[N] int g1;

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
  int ncXg1;
  int nrXg1;
  matrix[nrXg1, ncXg1] Xg1des;
  vector[ncXg1] sg1;
  
  int prior_only;
}
transformed data{
  
  // build full design matrices
  matrix[N, ncXs] Xs = Xsdes[silo,];
  matrix[N, ncXj] Xj = Xjdes[jnt,];
  matrix[N, ncXj] Xp = Xpdes[pref,];
  matrix[N, ncXd1] Xd1 = Xd1des[d1,];

  // indexes the relevant col in the design matrix
  matrix[N, ncXg1] Xg1 = Xg1des[g1,];

}
parameters{
  real mu;
  vector[ncXs] bs;
  vector[ncXj] bj;
  vector[ncXp] bp;
  vector[ncXd1] bd1;
  vector[ncXg1] bg1;
}
transformed parameters{
  vector[N] eta = mu + Xs*bs + Xj*bj + Xp*bp + Xd1*bd1 + Xg1*bg1  ;
} 
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(bs | 0, ss);
  target += normal_lpdf(bj | 0, sj);
  target += normal_lpdf(bp | 0, sp);
  target += normal_lpdf(bd1 | 0, sd1);
  target += normal_lpdf(bg1 | 0, sg1);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
}
