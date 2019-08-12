functions {
#include rate_equations.stan
#include allostery.stan
  real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){  
    real empty_array[0];
    real free_enzyme_ratio_r2 = get_free_enzyme_ratio_uniuni(metabolites[2],
                                                             metabolites[3],
                                                             known_reals[2]*params[6],
                                                             known_reals[2]*params[7],
                                                             params[8],
                                                             params[5]);
    return {
      uniuni(metabolites[1], metabolites[2], known_reals[1]*params[2], known_reals[1]*params[3], params[4], params[1]),
      uniuni(metabolites[2], metabolites[3], known_reals[2]*params[6], known_reals[2]*params[7], params[8], params[5])
        * get_regulatory_effect(empty_array, {metabolites[2]}, free_enzyme_ratio_r2, empty_array, {params[9]}, params[10]),
      uniuni(metabolites[3], metabolites[4], known_reals[3]*params[12], known_reals[3]*params[13], params[14], params[11])
    };
  }
  real[] get_odes(real[] fluxes){
    return {  
      0,
      1*fluxes[1]-1*fluxes[2],
      1*fluxes[2]-1*fluxes[3],
      0
    };
  }
  real[] steady_state_equation(real t,
                               real[] metabolites,
                               real[] params,
                               real[] known_reals,
                               int[] known_ints){
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
  // measurements
  vector[N_flux_measurement] flux_measurement;
  vector[N_concentration_measurement] concentration_measurement;
  // measurement scales, i.e. how much the measured values deviate from the true values
  vector<lower=0>[N_flux_measurement] flux_measurement_scale;
  vector<lower=0>[N_concentration_measurement] concentration_measurement_scale;
  // hardcoded priors
  vector[N_kinetic_parameter] prior_location_kinetic_parameter;
  vector<lower=0>[N_kinetic_parameter] prior_scale_kinetic_parameter;
  real prior_location_unbalanced[N_unbalanced, N_experiment];
  real<lower=0> prior_scale_unbalanced[N_unbalanced, N_experiment];
  // other hardcoded numbers
  real known_reals[N_known_real, N_experiment];
  // ode configuration
  real initial_time;
  real steady_time;
  real rel_tol;
  real abs_tol;
  int max_steps;
  // likelihood configuration - set to 0 for priors-only mode
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  int known_ints[0];
}
parameters {
  real<lower=0> kinetic_parameter[N_kinetic_parameter];
  real<lower=0> concentration_unbalanced[N_unbalanced, N_experiment];
}
transformed parameters {
  real concentration[N_balanced+N_unbalanced, N_experiment];
  real flux[N_reaction, N_experiment];
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
  }
}
model {
  kinetic_parameter ~ lognormal(prior_location_kinetic_parameter, prior_scale_kinetic_parameter);
  for (e in 1:N_experiment){
    concentration_unbalanced[,e] ~ lognormal(prior_location_unbalanced[,e], prior_scale_unbalanced[,e]);
  }
  if (LIKELIHOOD == 1){
    real concentration_hat[N_concentration_measurement];
    real flux_hat[N_flux_measurement];
    for (mc in 1:N_concentration_measurement){
      concentration_hat[mc] = concentration[ix_metabolite_concentration_measurement[mc],
                                            ix_experiment_concentration_measurement[mc]];
    }
    for (mf in 1:N_flux_measurement){
      flux_hat[mf] = flux[ix_reaction_flux_measurement[mf], ix_experiment_flux_measurement[mf]];
    }
    concentration_measurement ~ lognormal(log(concentration_hat), concentration_measurement_scale);
    flux_measurement ~ normal(flux_hat, flux_measurement_scale);
  }
}
generated quantities {
  vector[N_concentration_measurement] simulated_concentration_measurement;
  vector[N_flux_measurement] simulated_flux_measurement;
  real balanced_metabolite_rate_of_change[N_balanced, N_experiment];
  for (e in 1:N_experiment){
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