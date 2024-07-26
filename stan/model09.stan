// group all your data up first
// instructions:
data{ 
  int N;
  int N1;
  int N2;
  int N3;
  array[N1] int ixs1;
  array[N2] int ixs2;
  array[N3] int ixs3;
  
  array[N] int y;
  array[N] int n;
  
  array[N1+N2] int silo_1;
  array[N1+N2] int jnt_1;
  array[N1+N2] int pref_1;
  array[N1+N2] int g1_1;
  
  array[N3] int silo_2;
  array[N3] int jnt_2;
  array[N3] int pref_2;
  array[N3] int d1_2;

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
  // domain 1 (non-randomised surgery)
  int nrXg1;
  int ncXg1;
  matrix[nrXg1, ncXg1] Xg1des;
  vector[ncXg1] sg1;
  // domain 1 randomised surgery (rev portions selected)
  int nrXd1;
  int ncXd1;
  matrix[nrXd1, ncXd1] Xd1des;
  vector[ncXd1] sd1;
 
  int prior_only;
}
transformed data{
  
  // build full design matrices
  matrix[N1, ncXs] X1s = Xsdes[silo_1[ixs1]];
  matrix[N1, ncXj] X1j = Xjdes[jnt_1[ixs1]];
  matrix[N1, ncXj] X1p = Xjdes[pref_1[ixs1]];
  matrix[N1, ncXg1] X1g1 = Xg1des[g1_1[ixs1]];
  
  matrix[N2, ncXs] X2s = Xsdes[silo_1[ixs2]];
  matrix[N2, ncXj] X2j = Xjdes[jnt_1[ixs2]];
  matrix[N2, ncXj] X2p = Xjdes[pref_1[ixs2]];
  matrix[N2, ncXg1] X2g1 = Xg1des[g1_1[ixs2]];
  
  matrix[N3, ncXs] X3s = Xsdes[silo_2[ixs3]];
  matrix[N3, ncXj] X3j = Xjdes[jnt_2[ixs3]];
  matrix[N3, ncXj] X3p = Xjdes[pref_2[ixs3]];
  matrix[N3, ncXd1] X3d1 = Xd1des[d1_2[ixs3]];

}
parameters{
  real mu;
  vector[ncXs] bs;
  vector[ncXj] bj;
  vector[ncXp] bp;
  vector[ncXg1] bg1;
  vector[ncXd1] bd1;
}
transformed parameters{
  vector[N1] eta1 = mu + X1s*bs + X1j*bj + X1p*bp + X1g1*bg1;
  vector[N2] eta2 = mu + X2s*bs + X2j*bj + X2p*bp + X2g1*bg1;
  vector[N3] eta3 = mu + X3s*bs + X3j*bj + X3p*bp + X3d1*bd1;
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
  target += normal_lpdf(bg1 | 0, sg1);
  target += normal_lpdf(bd1 | 0, sd1);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
}
