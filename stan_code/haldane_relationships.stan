real get_Kip_ordered_unibi(real Keq, real Kia, real Kq, real V1, real V2){
  return Keq*Kia*V2 / (Kq*V1);
}
real get_Kiq_ordered_unibi(real Keq, real Ka, real Kp, real V1, real V2){
  return Keq*V2*Ka / (V1*Kp);
}
