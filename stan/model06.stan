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
}
parameters{
  real b0;
  real b_silo_raw;
  real b_joint_raw;
  vector[O] z_surg;
  real<lower=0> s_surg;
}
transformed parameters{
  // the 3 just scales the prior sd, purely for demonstration
  vector[N] eta;
  vector[S] b_silo;
  vector[J] b_joint;
  vector[O] b_surg;

  b_silo[1] = 0.0;
  b_silo[2] = b_silo_raw;
  b_joint[1] = 0.0;
  b_joint[2] = b_joint_raw;

  b_surg = z_surg * s_surg;

  for(i in 1:N){
    //eta[i] =  
    //    b_silo[silo[i]] + b_joint[joint[i]] + 
    //    b_surg[surg[i]] + b_abdu[abdu[i]];
    eta[i] =  
        b0 + b_silo[silo[i]] + b_joint[joint[i]] +
        b_surg[surg[i]] ; 
  }
} 
model{
  target += normal_lpdf(b0 | 0, 1.5);
  target += std_normal_lpdf(b_silo_raw);
  target += std_normal_lpdf(b_joint_raw);
  target += std_normal_lpdf(z_surg);
  target += exponential_lpdf(s_surg | 1);

  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);
  }
}
generated quantities{
  vector[N] p = inv_logit(eta);
}

