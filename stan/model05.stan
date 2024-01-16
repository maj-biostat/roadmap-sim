// group all your data up first
// instructions:
data{ 
  int N; // total number of records
  int S; // number of silos
  int J; 
  int O; // number of surgery types
  int P; // number of abdu types
  
  array[N] int y;
  array[N] int n;

  array[N] int silo;
  array[N] int joint;
  array[N] int surg;
//  array[N] int abdu;

  // prior only
  int prior_only;
}
transformed data{
  matrix [O,O] XS = diag_matrix(rep_vector(1,O));
  matrix [O,O-1] XS_qr;
  matrix [P,P] XD = diag_matrix(rep_vector(1,P));
  matrix [P,P-1] XD_qr;

  for(i in 1:(O-1)){ XS[O,i] = -1; }
  for(i in 1:(P-1)){ XD[P,i] = -1; }

  XS[O,O] = 0;
  XD[P,P] = 0;
  XS_qr = qr_Q(XS)[,1:(O-1)];
  XD_qr = qr_Q(XD)[,1:(P-1)];
}
parameters{
  vector[S] b_silo;
  real b_joint_raw;
  vector[O-1] b_surg_raw;
//  vector[P-1] b_abdu_raw;
}
transformed parameters{
  // the 3 just scales the prior sd, purely for demonstration
  vector[N] eta;
  vector[J] b_joint;
  vector[O] b_surg = 3 * XS_qr * b_surg_raw;
//  vector[P] b_abdu = 3 * XD_qr * b_abdu_raw;
//
  b_joint[1] = 0.0;
  b_joint[2] = b_joint_raw;

  for(i in 1:N){
    //eta[i] =  
    //    b_silo[silo[i]] + b_joint[joint[i]] + 
    //    b_surg[surg[i]] + b_abdu[abdu[i]];
    eta[i] =  
        b_silo[silo[i]] + b_joint[joint[i]] +
        b_surg[surg[i]] ; 
  }
} 
model{
  target += normal_lpdf(b_silo | 0, 1.5);
  target += normal_lpdf(b_joint_raw | 0, 1);
  target += normal_lpdf(b_surg_raw | 0, inv(sqrt(1 - inv(O))) );
 // target += normal_lpdf(b_abdu_raw | 0, inv(sqrt(1 - inv(P))) );

  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);
  }
}
generated quantities{
  vector[N] p = inv_logit(eta);
}

