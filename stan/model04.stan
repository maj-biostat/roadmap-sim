// group all your data up first
// instructions:
data{ 
  int N; // total number of records
  int S; // number of silos
  int J; // number of joints
  int O; // number of surgery types
  int P; // number of abdu types
  int Q; // number of abty types
  
  array[N] int y;
  array[N] int n;

  array[N] int silo;
  array[N] int joint;
  array[N] int rand_surg;
  array[N] int surg;
  array[N] int rand_abdu;
  array[N] int abdu;
  array[N] int rand_abty;
  array[N] int abty;

  // prior only
  int prior_only;
}
transformed data{
  matrix [O,O] XS = diag_matrix(rep_vector(1,O));
  matrix [O,O-1] XS_qr;
  matrix [P,P] XD = diag_matrix(rep_vector(1,P));
  matrix [P,P-1] XD_qr;
  matrix [Q,Q] XT = diag_matrix(rep_vector(1,Q));
  matrix [Q,Q-1] XT_qr;

  for(i in 1:(O-1)){ XS[O,i] = -1; }
  for(i in 1:(P-1)){ XD[P,i] = -1; }
  for(i in 1:(Q-1)){ XT[Q,i] = -1; }

  XS[O,O] = 0;
  XD[P,P] = 0;
  XT[Q,Q] = 0;
  // The orthogonal matrix in the fat QR decomposition of A,
  // which implies that the resulting matrix is square with
  // the same number of rows as A
  XS_qr = qr_Q(XS)[,1:(O-1)];
  XD_qr = qr_Q(XD)[,1:(P-1)];
  XT_qr = qr_Q(XT)[,1:(Q-1)];

}
parameters{
  vector[S] b_silo;
  real b_joint_raw;
  real b_not_rand_surg_raw;
  real b_not_rand_abdu_raw;
  real b_not_rand_abty_raw;
  vector[O-1] b_surg_raw;
  vector[P-1] b_abdu_raw;
  vector[Q-1] b_abty_raw;
}
transformed parameters{
  // the 3 just scales the prior sd, purely for demonstration
  vector[N] eta;
  vector[J] b_joint;
  vector[O] b_surg = 3 * XS_qr * b_surg_raw;
  vector[P] b_abdu = 3 * XD_qr * b_abdu_raw;
  vector[Q] b_abty = 3 * XT_qr * b_abty_raw;
  vector[2] b_not_rand_surg;
  vector[2] b_not_rand_abdu;
  vector[2] b_not_rand_abty;

  b_joint[1] = 0.0;
  b_joint[2] = b_joint_raw;
  b_not_rand_surg[1] = 0.0;
  b_not_rand_abdu[1] = 0.0;
  b_not_rand_abty[1] = 0.0;
  b_not_rand_surg[2] = b_not_rand_surg_raw; 
  b_not_rand_abdu[2] = b_not_rand_abdu_raw; 
  b_not_rand_abty[2] = b_not_rand_abty_raw; 

  for(i in 1:N){
    eta[i] =  
        b_silo[silo[i]] + b_joint[joint[i]] + 
        b_not_rand_surg[rand_surg[i]] + (rand_surg[i] == 1 ? b_surg[surg[i]] : 0.0) + 
        b_not_rand_abdu[rand_abdu[i]] + (rand_abdu[i] == 1 ? b_abdu[abdu[i]] : 0.0) + 
        b_not_rand_abty[rand_abty[i]] + (rand_abty[i] == 1 ? b_abty[abty[i]] : 0.0);
  }
} 
model{
  target += normal_lpdf(b_silo | 0, 1.5);
  target += normal_lpdf(b_joint_raw | 0, 1);
  target += normal_lpdf(b_not_rand_surg_raw | 0, 1);
  target += normal_lpdf(b_not_rand_abdu_raw | 0, 1);
  target += normal_lpdf(b_not_rand_abty_raw | 0, 1);
  target += normal_lpdf(b_surg_raw | 0, inv(sqrt(1 - inv(O))) );
  target += normal_lpdf(b_abdu_raw | 0, inv(sqrt(1 - inv(P))) );
  target += normal_lpdf(b_abty_raw | 0, inv(sqrt(1 - inv(Q))) );

  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);
  }
}
generated quantities{
  vector[N] p = inv_logit(eta);
}

