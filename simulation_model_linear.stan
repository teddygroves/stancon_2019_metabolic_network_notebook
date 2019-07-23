functions {
#include rate_equations.stan
#include allostery.stan
  real[] get_fluxes(real[] metabolites, real[] params, real[] known_reals){
    real empty_array[0];
    real free_enzyme_ratio_r2 = get_free_enzyme_ratio_uniuni(metabolites[1],
                                                             metabolites[2],
                                                             known_reals[2]*params[2],
                                                             known_reals[2]*params[3],
                                                             params[4], params[1]);
    return {
      fixed_flux(known_reals[1]),
      uniuni(metabolites[1], metabolites[2], known_reals[2]*params[2], known_reals[2]*params[3], params[4], params[1])
      * get_regulatory_effect(empty_array, {metabolites[1]}, free_enzyme_ratio_r2, empty_array, {params[5]}, params[6]),
      irr_mass_action(metabolites[2], params[7])
    };
  }
  real[] get_odes(real[] fluxes){
    return {
      1*fluxes[1]-1*fluxes[2],
      1*fluxes[2]-1*fluxes[3]
    };
  }
  real[] steady_state_equation(real t, real[] metabolites, real[] params, real[] known_reals, int[] known_ints){
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
  int<lower=1> N_metabolite;
  int<lower=1> N_param;
  int<lower=1> N_reaction;
  int<lower=1> N_experiment;
  int<lower=1> N_known_real;
  int<lower=1> N_measurement_flux;
  int<lower=1> N_measurement_conc;
  // measurements
  int<lower=0,upper=N_metabolite> metabolite_ix[N_measurement_conc];
  int<lower=1,upper=N_experiment> experiment_ix_conc[N_measurement_conc];
  vector<lower=0>[N_measurement_conc] measurement_scale_conc;
  int<lower=0,upper=N_reaction> reaction_ix[N_measurement_flux];
  int<lower=1,upper=N_experiment> experiment_ix_flux[N_measurement_flux];
  vector<lower=0>[N_measurement_flux] measurement_scale_flux;
  // hardcoded
  real known_reals[N_known_real, N_experiment];
  real params[N_param];
  // ode stuff
  real initial_concentration[N_metabolite];
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
  vector[N_measurement_flux] flux_pred;
  vector[N_measurement_conc] conc_pred;
  real metabolite_flux[N_metabolite, N_experiment]; 
  {
    real metabolite_concentration[N_metabolite, N_experiment];
    real flux[N_reaction, N_experiment];
    for (e in 1:N_experiment){
      metabolite_concentration[,e] = integrate_ode_bdf(steady_state_equation,
                                                       initial_concentration,
                                                       initial_time,
                                                       {steady_time},
                                                       params,
                                                       known_reals[,e],
                                                       known_ints,
                                                       rel_tol, abs_tol, max_steps)[1];
      flux[,e] = get_fluxes(metabolite_concentration[,e], params, known_reals[,e]);
      metabolite_flux[, e] = get_odes(flux[, e]);
    }
    for (mc in 1:N_measurement_conc){
      conc_pred[mc] = lognormal_rng(log(metabolite_concentration[metabolite_ix[mc], experiment_ix_conc[mc]]),
                                    measurement_scale_conc[mc]);
    }
    for (mf in 1:N_measurement_flux){
      flux_pred[mf] = normal_rng(flux[reaction_ix[mf], experiment_ix_flux[mf]], measurement_scale_flux[mf]);
    }
  }
}
