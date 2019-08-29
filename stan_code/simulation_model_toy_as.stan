functions {
#include rate_equations.stan
#include haldane_relationships.stan
#include allostery.stan

  vector get_fluxes(real[] m, real[] k, real[] xr){  
    real r2_Kip = get_Kip_ordered_unibi(k[5], k[11], k[10], k[6], k[7]);
    real r2_Kiq = get_Kiq_ordered_unibi(k[5], k[8], k[9], k[6], k[7]);
    return [uniuni(m[2], m[1], xr[1]*k[2], xr[1]*k[3], k[4], k[1]),
            ordered_unibi(m[1], m[3], m[4], xr[2]*k[6], xr[2]*k[7], k[8], k[9], k[10], k[11], r2_Kip, r2_Kiq, k[5]),
            uniuni(m[4], m[5], xr[3]*k[13], xr[3]*k[14], k[15], k[12]),
            uniuni(m[5], m[6], xr[4]*k[17], xr[4]*k[18], k[19], k[16]),
            uniuni(m[3], m[7], xr[5]*k[21], xr[5]*k[22], k[23], k[20]),
            uniuni(m[3], m[4], xr[6]*k[25], xr[6]*k[26], k[27], k[24])]';
  }

  real[] ode_func(real t, real[] concentration, real[] kinetic_parameters, real[] xr, int[] xi){
    vector[6] fluxes = get_fluxes(concentration, kinetic_parameters, xr);
    return {1*fluxes[1]-1*fluxes[2],
            0,
            1*fluxes[2]-1*fluxes[5]-1*fluxes[6],
            1*fluxes[2]-1*fluxes[3]+1*fluxes[6],
            1*fluxes[3]-1*fluxes[4],
            0,
            0};
  }

  vector steady_state_system(vector balanced, vector theta, real[] xr, int[] xi){
    int N_unbalanced = 3;
    int N_balanced = 4;
    real initial_time = 0;
    real time_step = 0.05;
    real conc[N_balanced + N_unbalanced];
    real balanced_new[N_balanced];
    conc[{1, 3, 4, 5}] = to_array_1d(balanced);
    conc[{2, 6, 7}] = to_array_1d(theta[1:N_unbalanced]);
    balanced_new = integrate_ode_bdf(ode_func,
                                     conc,
                                     initial_time,
                                     rep_array(time_step, 1),
                                     to_array_1d(theta[N_unbalanced+1:]),
                                     xr,
                                     rep_array(0, 1))[1, {1, 3, 4, 5}]; 
    return to_vector(balanced_new) - balanced;
  }
}
data {
  // dimensions
  int<lower=1> N_balanced;    // 'Balanced' metabolites must have constant concentration at steady state
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_kinetic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_experiment;
  int<lower=1> N_known_real;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_conc_measurement;
  // measurements
  int<lower=1,upper=N_experiment> experiment_yconc[N_conc_measurement];
  int<lower=1,upper=N_balanced+N_unbalanced> metabolite_yconc[N_conc_measurement];
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  vector<lower=0>[N_flux_measurement] sigma_flux;
  // hardcoded numbers
  real xr[N_experiment, N_known_real];
  vector<lower=0>[N_balanced] balanced_guess;
  vector<lower=0>[N_kinetic_parameter] kinetic_parameter;
  vector<lower=0>[N_unbalanced] unbalanced[N_experiment];
  // ode configuration
  real abs_tol;
  real f_tol;
  int max_steps;
}
transformed data {
  int xi[0];
}
parameters {
  real z;
}
model {
  z ~ std_normal();
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_flux_measurement] yflux_sim;
  vector<lower=0>[N_balanced+N_unbalanced] conc[N_experiment];
  vector[N_reaction] flux[N_experiment];
  for (e in 1:N_experiment){
    vector[N_unbalanced+N_kinetic_parameter] theta = append_row(unbalanced[e], kinetic_parameter);
    conc[e, {1, 3, 4, 5}] = algebra_solver(steady_state_system, balanced_guess, theta, xr[e], xi, abs_tol, f_tol, max_steps);
    conc[e, {2, 6, 7}] = unbalanced[e];
    flux[e] = get_fluxes(to_array_1d(conc[e]), to_array_1d(kinetic_parameter), xr[e]);
  }
  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], metabolite_yconc[c]]), sigma_conc[c]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
