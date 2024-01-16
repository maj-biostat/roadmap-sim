// group all your data up first
// instructions:
data{ 
  int N; // total number of records
  int S; // number of silos
  int N_silo; // total number of records by silo
  int J; // number of joints
  int D; // number domains
  // number trt for surg, abdu, abty
  array[D] int P_trt;
  
  array[N] int y;
  array[N] int n;

  array[N] int silo;
  array[N] int idx_early;
  
  array[N] int joint;
  
  array[N] int joint;

  // Design matrices:
  
  // silo, joint and 
  matrix[N, (S-1)] X_s;
  matrix[N, (J-1)] X_j;
  // non-rand indicator
  matrix[N, D] X_nr;
  
  // Treatments indexes
  array[N] int surg;
  array[N] int abdu;
  array[N] int abty;
  
  // prior only
  int prior_only;
}
transformed data{
  int EARLY = 1; int LATE = 2; int CHRONIC = 3;
  int SURG = 1; int ABDU = 2; int ABTY = 3;
  
}
parameters{
  real a; 
  // silo, joint, non-rand
  vector[S-1] b_silo; 
  vector[J-1] b_joint; 
  vector[D] b_nr; 

  //vector[PE_S-1] be_s_raw; 
  vector[PL_S-1] bl_s_raw; 
  vector[PC_S-1] bc_s_raw; 
  //vector[PE_D-1] be_d_raw; 
  vector[PL_D-1] bl_d_raw; 
  vector[PC_D-1] bc_d_raw;
  vector[PE_T-1] be_t_raw;
  vector[PL_T-1] bl_t_raw;
  vector[PC_T-1] bc_t_raw;
  
}
transformed parameters{
  vector[N] eta;

  //vector[PE_S] be_s; // early, surg
  vector[PL_S] bl_s; // late, surg (dair, rev = one or two)
  vector[PC_S] bc_s; // chronic, surg (one, two)
  //vector[PE_D] be_d; 
  vector[PL_D] bl_d; // late, duration (1, 6, 12 wks)  
  vector[PC_D] bc_d; // chronic, duration (1, 6, 12 wks) 
  vector[PE_T] be_t; // early, type - nonrif, rif 
  vector[PL_T] bl_t; 
  vector[PC_T] bc_t; 

  bl_s[1] = 0.0; bc_s[1] = 0.0;
  bl_d[1] = 0.0; bc_d[1] = 0.0;
  be_t[1] = 0.0; bl_t[1] = 0.0; bc_t = 0.0;

  eta = a + X_s * b_silo + 
    X_j * b_joint + 
    X_nr[,SURG] * b_nr[SURG] + 
    X_nr[,ABDU] * b_nr[ABDU] + 
    X_nr[,ABTY] * b_nr[ABTY] ; 

  for(i in 1:N){
    // early only has abty domain (and only if pt eligible).
    if(silo[i] == EARLY){
      eta[i] = eta[i] + X_nr[i,ABTY] * be_t[abty[i]];
    }

    if(silo[i] == LATE){
      // assumes surg = 2 is one stage and surg = 3 is two stage
      // bl_s[1] dair, bl_s[2] rev
      eta[i] = eta[i] + X_nr[i,SURG] * bl_s[surg[i] == 2 | surg[i] == 3 ? 2 : 1]; 

      eta[i] = eta[i] + X_nr[i,ABDU] * bl_d[surg[i] == 2 ? 2 : 1]; 
      eta[i] = eta[i] + X_nr[i,ABDU] * bl_d[surg[i] == 2 ? 2 : 1]; 
      eta[i] = eta[i] + X_nr[i,ABDU] * bl_d[surg[i] == 2 ? 2 : 1]; 
    }


  }

} 
model{
  target += normal_lpdf(a | 0, 1.5);
  target += std_normal_lpdf(b_silo);
  target += std_normal_lpdf(b_joint);
  target += std_normal_lpdf(b_nr);




}
generated quantities{
}

