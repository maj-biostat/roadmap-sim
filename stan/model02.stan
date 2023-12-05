functions {
   // if(in_it(3,{1,2,3,4})) will evaluate as 1
  int in_it(int pos,int[] pos_var) {
    for (p in 1:(size(pos_var))) {
       if (pos_var[p]==pos) {
       // can return immediately, as soon as find a match
          return 1;
       } 
    }
    return 0;
  }
}
data{
  int N;
  array[N] int y;
  array[N] int n;
  // 1 = early, 2 = late, 3 = chronic
  array[N] int silo;
  // 1 = knee, 2 = hip
  array[N] int joint;
  // 1 = dair, 2 = one-stage, 3 = two-stage
  array[N] int surg;
  // 1 = short, 2 = long
  array[N] int abdu;
  // 1 = soc, 2 = soc + rif
  array[N] int abty;
  array[N] int abty_elig;
}
parameters{
  vector[3] b_silo;
  real b_joint_raw;
  real b_all_abty_raw;
  real b_late_surg_rev;
  real b_late_abdu_surg_one_raw;
  real b_late_abdu_surg_two_raw;
  real b_chro_surg;
  real b_chro_abdu_surg_one_raw;
  real b_chro_abdu_surg_two_raw;
}
transformed parameters{
  vector[2] b_joint;
  vector[2] b_all_abty;
  vector[2] b_late_abdu_surg_one;
  vector[2] b_late_abdu_surg_two;
  vector[2] b_chro_abdu_surg_one;
  vector[2] b_chro_abdu_surg_two;
  b_joint[1] = 0.0;
  b_joint[2] = b_joint_raw;
  b_all_abty[1] = 0.0;
  b_all_abty[2] = b_all_abty_raw;
  b_late_abdu_surg_one[1] = 0.0;
  b_late_abdu_surg_one[2] = b_late_abdu_surg_one_raw;
  b_late_abdu_surg_two[1] = 0.0;
  b_late_abdu_surg_two[2] = b_late_abdu_surg_two_raw;
  b_chro_abdu_surg_one[1] = 0.0;
  b_chro_abdu_surg_one[2] = b_chro_abdu_surg_one_raw;
  b_chro_abdu_surg_two[1] = 0.0;
  b_chro_abdu_surg_two[2] = b_chro_abdu_surg_two_raw;
}
model{
  vector[N] eta;
  for(i in 1:N){
    eta[i] = b_silo[silo[i]] + b_joint[joint[i]];
   
    // for late acute silo - revision surgery and ab duration
    if(silo[i] == 2 & surg[i] == 2)){
      eta[i] = eta[i] + b_late_surg_rev + b_late_abdu_surg_one[abdu[i]];
    }
    if(silo[i] == 2 & surg[i] == 3)){
      eta[i] = eta[i] + b_late_surg_rev + b_late_abdu_surg_two[abdu[i]];
    }
    // for chronic silo - one vs two stage surgery and ab duration
    if(silo[i] == 3 & surg[i] == 2)){
      eta[i] = eta[i] + b_chro_abdu_surg_one[abdu[i]];
    }
    if(silo[i] == 3 & surg[i] == 3)){
      eta[i] = eta[i] + b_chro_surg + b_late_abdu_surg_two[abdu[i]];
    }
    // all domains include an ab choice effect if the pt is gram pos
    if(abty_elig[i] == 1){
      eta[i] = eta[i] + b_abty[abty[i]];
    }
  }

  target += normal_lpdf(b_silo | 0, 1.5);
  target += std_normal_lpdf(b_joint_raw);
  target += std_normal_lpdf(b_all_abty_raw);
  target += std_normal_lpdf(b_late_surg_rev);
  target += std_normal_lpdf(b_late_abdu_surg_one_raw);
  target += std_normal_lpdf(b_late_abdu_surg_two_raw);
  target += std_normal_lpdf(b_chro_surg);
  target += std_normal_lpdf(b_chro_abdu_surg_one_raw);
  target += std_normal_lpdf(b_chro_abdu_surg_two_raw);

  target += binomial_logit_lpmf(y | n, eta);

}
