data{
  int P;
  array[3] int idx;
}
parameters{
  vector[P ? P : 0] mu;
  real sig;
}
model{
  for(i in idx){print(i);}

  target += std_normal_lpdf(mu);
  target += std_normal_lpdf(sig);
}


