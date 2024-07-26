// group all your data up first
// instructions:
data{ 
  int N;
  array[N] int y;
  array[N] int n;
  
  // int S;
  array[N] int silo;
  array[N] int jnt;
  array[N] int pref;
  // arm index domain 1
  array[N] int d1;
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
  // surgical preference
  int nrXp;
  int ncXp;
  matrix[nrXp, ncXp] Xpdes;
  vector[ncXp] sp;
  // domain 1 surgery
  int nrXd1;
  int ncXd1;
  matrix[nrXd1, ncXd1] Xd1des;
  vector[ncXd1] sd1;
  
  
  int prior_only;
}
transformed data{
  
  // build full design matrices
  matrix[N, ncXs] Xs = Xsdes[silo];
  matrix[N, ncXj] Xj = Xjdes[jnt];
  matrix[N, ncXp] Xp = Xpdes[pref];
  matrix[N, ncXd1] Xd1 = Xd1des[d1];
  
}
parameters{
  real mu;
  vector[ncXs] bs;
  vector[ncXj] bj;
  vector[ncXp] bp;
  vector[ncXd1] bd1;
}
transformed parameters{
  // vector[N] eta = mu + Xs*bs + Xj*bj + Xp*bp + Xd1*bd1;
  vector[N] eta = mu + Xs*bs + Xj*bj + Xp*bp + Xd1*bd1;
} 
model{
  target += logistic_lpdf(mu | 0, 1);
  target += normal_lpdf(bs | 0, ss);
  target += normal_lpdf(bj | 0, sj);
  target += normal_lpdf(bp | 0, sp);
  target += normal_lpdf(bd1 | 0, sd1);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
}
