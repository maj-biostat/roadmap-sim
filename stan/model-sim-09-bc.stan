data{ 
  int N;
  array[N] int y;
  array[N] int n;
  int K_reg;
  // lookup: for reg cell k, which silo column (1..4)
  array[K_reg] int reg_silo_idx;  
   // lookup: which d1 row (1=dair, 2=r1, 3=r2)
  array[K_reg] int reg_d1_idx;   
  array[K_reg] int reg_d2_idx;
  array[K_reg] int reg_d3_idx;
  
  array[N] int silo;
  
  // each dair, r1, r2
  array[N] int  d1;
  
  // nad2, wk12, wk6
  array[N] int d2;
  // nad3, none, wk12
  array[N] int d3;
  // nad4, norif, rif
  array[N] int d4;
  
  // priors
  vector[2] pri_b_0;
  vector[2] pri_b_reg; 
  vector[2] pri_b_d4;
  
  int prior_only;
}
transformed data{
}
parameters{
  real b_0;
  vector[2] b_d1_l_raw;
  vector[3] b_d1_lnr;
  vector[3] b_d1_enr;
  vector[3] b_d1_cnr;
  
  vector[2] b_d2_raw;
  vector[2] b_d3_raw;
  vector[2] b_d4_raw;
}
transformed parameters{

  // dair, r1, r2 x l, lnr, enr, cnr
  matrix[3, 4] b_d1;
  
  b_d1[1, 1] = 0.0;
  b_d1[2:3, 1] = b_d1_l_raw;
  b_d1[, 2] = b_d1_lnr;
  b_d1[, 3] = b_d1_enr;
  b_d1[, 4] = b_d1_cnr;
  
  vector[3] b_d2;
  b_d2[1] = 0.0;
  b_d2[2:3] = b_d2_raw;
  
  vector[3] b_d3;
  b_d3[1] = 0.0;
  b_d3[2:3] = b_d3_raw;
  
  vector[3] b_d4;
  b_d4[1] = 0.0;
  b_d4[2:3] = b_d4_raw;
} 
model{
  target += logistic_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  
  target += normal_lpdf(b_d1_l_raw | 0, 2);
  target += normal_lpdf(b_d1_lnr | 0, 2);
  target += normal_lpdf(b_d1_enr | 0, 2);
  target += normal_lpdf(b_d1_cnr | 0, 2);
  
  target += normal_lpdf(b_d2_raw | 0, 2);
  target += normal_lpdf(b_d3_raw | 0, 2);
  target += normal_lpdf(b_d4_raw | 0, 2);
  
  if(!prior_only){
    for(i in 1:N){
      target += binomial_logit_lpmf(
        y[i] | n[i], b_0 + b_d1[d1[i], silo[i]] + b_d2[d2[i]] + b_d3[d3[i]] + b_d4[d4[i]]);  
    }
  }
}
generated quantities{
  vector[K_reg] b_reg;
  for (k in 1:K_reg) {
    b_reg[k] = b_d1[reg_d1_idx[k], reg_silo_idx[k]] + b_d2[reg_d2_idx[k]] + b_d3[reg_d3_idx[k]];
  }
}
