real uniuni(real A, real P, real V1, real V2, real Ka, real Keq){
  real num = V1*V2*(A-P/Keq);
  real denom = Ka*V2 + V2*A + V1*P/Keq;
  return num / denom;
}
real ordered_unibi(real A, real P, real Q,
                   real V1, real V2,
                   real Ka, real Kp, real Kq,
                   real Kia, real Kip, real Kiq,
                   real Keq){
  real num = V1*V2*(A-P*Q/Keq); 
  real denom =
    Ka*V2
    + V2*A
    + Kq*V1*P/Keq
    + Kp*V1*Q/Keq
    + V1*P*Q/Keq
    + V2*A*P/Kip;
  return num / denom;
}
real irr_mass_action(real A, real V1){
  return(A*V1);
}
real fixed_flux(real f){
  return f;
}