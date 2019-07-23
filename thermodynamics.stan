real get_Keq(real delta_g, real T, real R){
  return exp(delta_g / (-R * T));
}
real get_delta_g(real Keq, real T, real R){
  return log(Keq) * (-R * T);
}
