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
  array[N] int g1;
  array[N] int g2;

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
  
  // domain 2 
  int nrXd2;
  int ncXd2;
  matrix[nrXd2, ncXd2] Xd2des;
  vector[ncXd2] sd2;
  
  // silo by domain 2 interaction - 
  // target surg specific effects for d2 not an average effect
  int ncXg2;
  int nrXg2;
  matrix[nrXg2, ncXg2] Xg2des;
  vector[ncXg2] sg2;
  
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
  matrix[N, ncXd2] Xd2 = Xd2des[d2,];
  matrix[N, ncXg2] Xg2 = Xg2des[g2,];

}
parameters{
  real mu;
  vector[ncXs] bs;
  vector[ncXj] bj;
  vector[ncXp] bp;
  vector[ncXd1] bd1;
  vector[ncXd2] bd2;
  vector[ncXg1] bg1;
  vector[ncXg2] bg2;
}
transformed parameters{
  vector[N] eta = mu + Xs*bs + Xj*bj + Xp*bp + Xd1*bd1 + Xd2*bd2 + Xg1*bg1 + Xg2*bg2 ;
} 
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(bs | 0, ss);
  target += normal_lpdf(bj | 0, sj);
  target += normal_lpdf(bp | 0, sp);
  target += normal_lpdf(bd1 | 0, sd1);
  target += normal_lpdf(bd2 | 0, sd2);
  
  target += normal_lpdf(bg1 | 0, sg1);
  target += normal_lpdf(bg2 | 0, sg2);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
}
