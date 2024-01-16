functions{
  matrix qr1(int O){
    matrix [O,O] XS = diag_matrix(rep_vector(1,O));
    matrix [O,O] XS_qr;

    for(i in 1:(O-1)){ XS[O,i] = -1; }

    XS[O,O] = 0;
    // The orthogonal matrix in the fat QR decomposition of A,
    // which implies that the resulting matrix is square with
    // the same number of rows as A
    //XS_qr = qr_Q(XS)[,1:(O-1)];
    XS_qr = qr_Q(XS);

    return XS_qr;
  }

  matrix qr2(int O){
    matrix [O,O] XS = diag_matrix(rep_vector(1,O));
    matrix [O,O] RR;
    for(i in 1:(O-1)){ XS[O,i] = -1; }
    XS[O,O] = 0;
    RR = qr_R(XS);
    //RR_inv = inverse(RR);
    return RR;
  } 
}
data{
  int O;
}
transformed data{
}
parameters{
}
model{
}
