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
  array[N] int surg;
  array[N] int abdu;
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
  vector[O-1] b_surg_raw;
  vector[P-1] b_abdu_raw;
  vector[Q-1] b_abty_raw;
}
transformed parameters{
  // the 3 just scales the prior sd, purely for demonstration
  vector[J] b_joint;
  vector[O] b_surg = 3 * XS_qr * b_surg_raw;
  vector[P] b_abdu = 3 * XD_qr * b_abdu_raw;
  vector[Q] b_abty = 3 * XT_qr * b_abty_raw;

  b_joint[1] = 0.0;
  b_joint[2] = b_joint_raw;
} 
model{
  target += normal_lpdf(b_silo | 0, 1.5);
  target += normal_lpdf(b_joint_raw | 0, 1);
  target += normal_lpdf(b_surg_raw | 0, inv(sqrt(1 - inv(O))) );
  target += normal_lpdf(b_abdu_raw | 0, inv(sqrt(1 - inv(P))) );
  target += normal_lpdf(b_abty_raw | 0, inv(sqrt(1 - inv(Q))) );

  if(!prior_only){
    for(i in 1:N){
      // leave as binomial now so negates the possibility of subject level effects
      target += binomial_logit_lpmf(y[i] | n[i], 
        b_silo[silo[i]] + b_joint[joint[i]] + 
        b_surg[surg[i]] +
        b_abdu[abdu[i]] +
        b_abty[abty[i]]);
    }
} }
generated quantities{
  vector[N] p;
  for(i in 1:N){
    p[i] = inv_logit(b_silo[silo[i]] + 
      b_joint[joint[i]] + 
      b_surg[surg[i]] +
      b_abdu[abdu[i]] +
      b_abty[abty[i]] 
      );
  }
}

