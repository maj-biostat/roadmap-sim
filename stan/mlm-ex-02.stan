// Model for:
// independent estimates of treatment arm means
// shared within group variance for deflections from trt arm means 
// due to subgroup membership
// goal being to shrink the subgroup deflections down to the treatment
// arm means if the subgroup variation is low
data {
  int N;
  array[N] int y;
  // intervention arms (one less that total num arms)
  int J;
  array[N] int j;
  // subgroups of interest creating effect heterogeneity
  int K;
  array[N] int k;
  real s;
  int prior_only;
}
transformed data {
  // QR decom setup - sum to zero constraint
  // matrix [J,J] A = diag_matrix(rep_vector(1,J));
  // matrix [J,J-1] A_qr;
  // for(i in 1:J-1){
  //   A[J,i] = -1;
  // }
  // A[J,J] = 0;
  // // The orthogonal matrix in the fat QR decomposition of A,
  // // which implies that the resulting matrix is square with
  // // the same number of rows as A
  // A_qr = qr_Q(A)[,1:(J-1)];
}
parameters{
  // real mu;
  // vector[J-1] z_j_raw;
  
  // convenience to adjust to group mean rather than an offset
  vector[J] mu_j;
  
  // parameters to construct offsets from trt arm means
  real<lower=0> s_j;
  matrix[J,K] z_j_k;
}
transformed parameters{
  // zero centred normal prior on each of the J pars
  // vector[J] z_j;
  
  matrix[J, K] mu_j_k;
  vector[N] eta;
  
  // z_j = s * A_qr * z_j_raw;
  // mu_j = mu + z_j;
  
  for(i in 1:J){
    mu_j_k[i, ] = to_row_vector(mu_j[i] + z_j_k[i, ] .* s_j);
  }
  for(i in 1:N){
    eta[i] = mu_j_k[j[i], k[i]];
  }
}
model{

  // target += logistic_lpdf(mu | 0, 1);
  // target += normal_lpdf(z_j_raw | 0, inv( sqrt(1 - inv(J))));
  
  // target += logistic_lpdf(mu_j | 0, 1);
  
  target += normal_lpdf(mu_j | 0, s);
  
  target += std_normal_lpdf(to_vector(z_j_k));
  target += student_t_lpdf(s_j | 3, 0, 2) -
     1 * student_t_lccdf(0 | 3, 0, 2);
  
  if(!prior_only){
    target += bernoulli_logit_lpmf(y | eta);  
  }
}
generated quantities{
}
