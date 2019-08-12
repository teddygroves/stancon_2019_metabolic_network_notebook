functions {
#include rate_equations.stan
#include haldane_relationships.stan
#include allostery.stan
real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){  
    real r2_Kip = get_Kip_ordered_unibi(params[5], params[11], params[10], params[6], params[7]);
    real r2_Kiq = get_Kiq_ordered_unibi(params[5], params[8], params[9], params[6], params[7]);

  

  return{   
      uniuni(metabolites[2], metabolites[1], known_reals[1]*params[2], known_reals[1]*params[3], params[4], params[1]),
      ordered_unibi(metabolites[1], metabolites[3], metabolites[4], known_reals[2]*params[6], known_reals[2]*params[7], params[8], params[9], params[10], params[11], r2_Kip, r2_Kiq, params[5]),
      uniuni(metabolites[4], metabolites[5], known_reals[3]*params[13], known_reals[3]*params[14], params[15], params[12]),
      uniuni(metabolites[5], metabolites[6], known_reals[4]*params[17], known_reals[4]*params[18], params[19], params[16]),
      uniuni(metabolites[3], metabolites[7], known_reals[5]*params[21], known_reals[5]*params[22], params[23], params[20]),
      uniuni(metabolites[3], metabolites[4], known_reals[6]*params[25], known_reals[6]*params[26], params[27], params[24])
    };
  }
real[] get_odes(real[] fluxes){
  return {  
      1*fluxes[1]-1*fluxes[2],
      0,
      1*fluxes[2]-1*fluxes[5]-1*fluxes[6],
      1*fluxes[2]-1*fluxes[3]+1*fluxes[6],
      1*fluxes[3]-1*fluxes[4],
      0,
      0
  };
}
real[] steady_state_equation(
          real t,
          real[] metabolites,
          real[] params,
          real[] known_reals,
          int[] known_ints
        ){
        for (m in 1:size(metabolites)){
          if (metabolites[m] < 0){
            reject("Metabolite ", m, " is ", metabolites[m], " but should be greater than zero.");
          }
        }
        return get_odes(get_fluxes(metabolites, params, known_reals));
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
  int<lower=1> N_concentration_measurement;
  // position of balanced and unbalanced metabolites in overall metabolite array 
  int<lower=1,upper=N_balanced+N_unbalanced> pos_balanced[N_balanced];      
  int<lower=1,upper=N_balanced+N_unbalanced> pos_unbalanced[N_unbalanced];
  // which measurement goes with which experiment
  int<lower=1,upper=N_experiment> ix_experiment_concentration_measurement[N_concentration_measurement];
  int<lower=1,upper=N_experiment> ix_experiment_flux_measurement[N_flux_measurement];
  // which measurement goes with which reaction or metabolite
  int<lower=1,upper=N_balanced+N_unbalanced> ix_metabolite_concentration_measurement[N_concentration_measurement];
  int<lower=1,upper=N_reaction> ix_reaction_flux_measurement[N_flux_measurement];
  // measurement scales, i.e. how much the measured values deviate from the true values
  vector<lower=0>[N_flux_measurement] flux_measurement_scale;
  vector<lower=0>[N_concentration_measurement] concentration_measurement_scale;
  // hardcoded numbers
  real known_reals[N_known_real, N_experiment];
  real<lower=0> kinetic_parameter[N_kinetic_parameter];
  real<lower=0> concentration_unbalanced[N_unbalanced, N_experiment];
  // ode stuff
  real initial_time;
  real steady_time;
  real rel_tol;
  real abs_tol;
  int max_steps;
}
transformed data {
  int known_ints[0];
}
parameters {
  real z;
}
model {
  z ~ normal(0, 1);
}
generated quantities {
  vector[N_concentration_measurement] simulated_concentration_measurement;
  vector[N_flux_measurement] simulated_flux_measurement;
  real concentration[N_balanced+N_unbalanced, N_experiment];
  real flux[N_reaction, N_experiment];
  real balanced_metabolite_rate_of_change[N_balanced, N_experiment];
  for (e in 1:N_experiment){
    real initial_concentration[N_balanced+N_unbalanced];
    initial_concentration[pos_balanced] = rep_array(1.0, N_balanced);
    initial_concentration[pos_unbalanced] = concentration_unbalanced[,e];
    concentration[,e] = integrate_ode_bdf(steady_state_equation,
                                          initial_concentration,
                                          initial_time,
                                          {steady_time},
                                          kinetic_parameter,
                                          known_reals[,e],
                                          known_ints,
                                          rel_tol, abs_tol, max_steps)[1];
    flux[,e] = get_fluxes(concentration[,e], kinetic_parameter, known_reals[,e]);
    balanced_metabolite_rate_of_change[, e] = get_odes(flux[, e])[pos_balanced];
  }
  for (mc in 1:N_concentration_measurement){
    simulated_concentration_measurement[mc] =
      lognormal_rng(log(concentration[ix_metabolite_concentration_measurement[mc],
                                      ix_experiment_concentration_measurement[mc]]),
                    concentration_measurement_scale[mc]);
  }
  for (mf in 1:N_flux_measurement){
    simulated_flux_measurement[mf] =
      normal_rng(flux[ix_reaction_flux_measurement[mf],
                      ix_experiment_flux_measurement[mf]],
                 flux_measurement_scale[mf]);
  }
}
