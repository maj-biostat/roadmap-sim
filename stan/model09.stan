// group all your data up first
// instructions:
data{ 
  int N;
  array[N] int y;
  array[N] int n;
  
  // int S;
  array[N] int silo;
  
  int N1;
  int N2;
  int N3;
  array[N1] int ix1;
  array[N2] int ix2;
  array[N3] int ix3;
  
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
  matrix[N1, ncXs] X1s = Xsdes[silo[ix1]];
  matrix[N1, ncXj] X1j = Xjdes[jnt[ix1]];
  matrix[N1, ncXd1] X1d1 = Xd1des[d1[ix1]];
  
  matrix[N2, ncXs] X2s = Xsdes[silo[ix2]];
  matrix[N2, ncXj] X2j = Xjdes[jnt[ix2]];
  matrix[N2, ncXp] X2p = Xpdes[pref[ix2]];
  matrix[N2, ncXd1] X2d1 = Xd1des[d1[ix2]];
  
  matrix[N3, ncXs] X3s = Xsdes[silo[ix3]];
  matrix[N3, ncXj] X3j = Xjdes[jnt[ix3]];
  matrix[N3, ncXd1] X3d1 = Xd1des[d1[ix3]];
  
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
  vector[N1] eta1 = mu + X1s*bs + X1j*bj + X1d1*bd1;
  vector[N2] eta2 = mu + X2s*bs + X2j*bj + X2p*bp + X2d1*bd1;
  vector[N3] eta3 = mu + X3s*bs + X3j*bj + X3d1*bd1;
  vector[N] eta;
  eta[1:N1] = eta1;
  eta[(N1+1):(N1+N2)] = eta2;
  eta[(N1+N2+1):N] = eta3;
  
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
