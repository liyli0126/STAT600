#include <Rcpp.h>
#include <cmath>
#include <random>
#include <vector>
#include "ReactionDiffusionCpp.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// log normal prior of D and rho
double log_prior(double D, double rho,
                 double D_prior_mean_log, double D_prior_sd_log,
                 double rho_prior_mean_log, double rho_prior_sd_log) {
  if (D <= 0 || rho <= 0) return R_NegInf;
  double lp_D = -0.5 * std::pow((std::log(D) - D_prior_mean_log) / D_prior_sd_log, 2)
    - std::log(D_prior_sd_log) - std::log(D);
  double lp_rho = -0.5 * std::pow((std::log(rho) - rho_prior_mean_log) / rho_prior_sd_log, 2)
    - std::log(rho_prior_sd_log) - std::log(rho);
  return lp_D + lp_rho;
}

// Gaussian likelihood function
double log_likelihood_cn(double D, double rho,
                         const std::vector<double>& u_initial,
                         const std::vector<double>& u_observed,
                         double noise_std,
                         double domain_length,
                         int grid_size,
                         double dt,
                         double t_obs) {
  try {
    // Stability test
    double dx = domain_length / (grid_size - 1);
    double cfl = D * dt / (dx * dx);
    if (cfl > 0.5) {
      Rcpp::warning("CFL condition may be violated");
    }
    
    ReactionDiffusionSolver1D solver(D, rho, domain_length, grid_size, dt);
    std::vector<double> u_sim = solver.solve_crank_nicolson(u_initial, t_obs);
    if (u_sim.size() != u_observed.size()) {
      return R_NegInf;
    }
    double loglik = 0.0;
    double var = noise_std * noise_std;
    double inv_var = 1.0 / var;
    for (size_t i = 0; i < u_observed.size(); ++i) {
      double r = u_sim[i] - u_observed[i];
      loglik -= 0.5 * r * r * inv_var;
    }
    loglik -= 0.5 * u_observed.size() * std::log(2 * M_PI * var);
    return loglik;
  } catch (...) {
    return R_NegInf;
  }
}


// MCMC implementation
//[[Rcpp::export]]
List run_adaptive_mcmc_chain(
  const NumericVector& u_initial_r, const NumericVector u_observed_r, double noise_std, double domain_length,
  double t_obs, double D_prior_mean_log, double D_prior_sd_log, double rho_prior_mean_log, double rho_prior_sd_log,
  int grid_size = 100, double dt = 0.01, int num_iterations = 2000, int burnin = 1000, double initial_logD_prop_sd = 0.1,
  double initial_logrho_prop_sd = 0.1, double target_acceptance = 0.4, unsigned int seed = 0){
  
  if (u_initial_r.size() != u_observed_r.size()) {
    stop("u_initial and u_observed must have the same length.");
  }
  if (grid_size < 3) {
    stop("grid_size must be at least 3.");
  }
  if (dt <= 0) {
    stop("dt must be positive.");
  }
  
  std::vector<double> u_initial = as<std::vector<double>>(u_initial_r);
  std::vector<double> u_observed = as<std::vector<double>>(u_observed_r);
  
  // Generate random number
  std::mt19937 gen(seed == 0 ? std::random_device{}() : seed);
  std::uniform_real_distribution<> unif(0.0, 1.0);
  std::normal_distribution<> norm(0.0, 1.0);
  
  // Adaptive parameters
  double logD_prop_sd = initial_logD_prop_sd;
  double logrho_prop_sd = initial_logrho_prop_sd;
  const int adaptation_interval = 100;
  double adaptation_factor = 1.0;
  
  // Initialization: with middle of prior
  double logD_curr = D_prior_mean_log;
  double logrho_curr = rho_prior_mean_log;
  double D_curr = std::exp(logD_curr);
  double rho_curr = std::exp(logrho_curr);
  
  // Initial posterior
  double log_prior_curr = log_prior(D_curr, rho_curr, D_prior_mean_log, D_prior_sd_log, rho_prior_mean_log, rho_prior_sd_log);
  double log_lik_curr = log_likelihood_cn(D_curr, rho_curr, u_initial, u_observed, noise_std, domain_length, grid_size, dt, t_obs);
  double log_post_curr = log_prior_curr + log_lik_curr;
  
  if (!std::isfinite(log_post_curr)) {
    // If not finite, add disturbation
    warning("Initial posterior is not finite. Trying small perturbation.");
    for (int attempt = 0; attempt < 10; ++attempt) {
      logD_curr = D_prior_mean_log + 0.1 * norm(gen);
      logrho_curr = rho_prior_mean_log + 0.1 * norm(gen);
      D_curr = std::exp(logD_curr);
      rho_curr = std::exp(logrho_curr);
      
      log_prior_curr = log_prior(D_curr, rho_curr, D_prior_mean_log, D_prior_sd_log,
                                 rho_prior_mean_log, rho_prior_sd_log);
      log_lik_curr = log_likelihood_cn(D_curr, rho_curr, u_initial, u_observed,
                                       noise_std, domain_length, grid_size, dt, t_obs);
      log_post_curr = log_prior_curr + log_lik_curr;
      
      if (std::isfinite(log_post_curr)) break;
    }
    if (!std::isfinite(log_post_curr)) {
      stop("Failed to initialize with a finite posterior after 10 attempts.");
    }
  }

  std::vector<double> D_samples, rho_samples;
  std::vector<double> acceptance_rates;
  int total_accept = 0;
  int window_accept = 0;
  int total_iters = num_iterations + burnin;
  
  for (int iter = 0; iter < total_iters; ++iter) {
    double logD_prop = logD_curr + norm(gen) * logD_prop_sd;
    double logrho_prop = logrho_curr + norm(gen) * logrho_prop_sd;
    double D_prop = std::exp(logD_prop);
    double rho_prop = std::exp(logrho_prop);
    
    double log_prior_prop = log_prior(D_prop, rho_prop, D_prior_mean_log, D_prior_sd_log,
                                      rho_prior_mean_log, rho_prior_sd_log);
    double log_lik_prop = log_likelihood_cn(D_prop, rho_prop, u_initial, u_observed,
                                            noise_std, domain_length, grid_size, dt, t_obs);
    double log_post_prop = log_prior_prop + log_lik_prop;
    
    double log_alpha = log_post_prop - log_post_curr;
    if (std::isfinite(log_alpha)) {
      if (log_alpha > 0 || unif(gen) < std::exp(log_alpha)) {
        logD_curr = logD_prop;
        logrho_curr = logrho_prop;
        D_curr = D_prop;
        rho_curr = rho_prop;
        log_post_curr = log_post_prop;
        total_accept++;
        window_accept++;
      }
    }
    // Adaptive adjustion, only within burin region
    if (iter < burnin && iter > 0 && iter % adaptation_interval == 0) {
      double current_acceptance = static_cast<double>(window_accept) / adaptation_interval;
      if (current_acceptance > target_acceptance) {
        adaptation_factor = 1.0 + 0.1 * (current_acceptance - target_acceptance) / target_acceptance;
      } else {
        adaptation_factor = 1.0 / (1.0 + 0.1 * (target_acceptance - current_acceptance) / target_acceptance);
      }
      
      // Limit adjust range
      adaptation_factor = std::max(0.5, std::min(2.0, adaptation_factor));
      logD_prop_sd *= adaptation_factor;
      logrho_prop_sd *= adaptation_factor;
      
      logD_prop_sd = std::max(0.01, std::min(1.0, logD_prop_sd));
      logrho_prop_sd = std::max(0.01, std::min(1.0, logrho_prop_sd));
      
      acceptance_rates.push_back(current_acceptance);
      window_accept = 0;
    }
     
     // save samples
     if (iter >= burnin) {
       D_samples.push_back(D_curr);
       rho_samples.push_back(rho_curr);
     }
  }
  double final_accept_rate = static_cast<double>(total_accept) / total_iters;
  
  return List::create(
    _["D"] = wrap(D_samples),
    _["rho"] = wrap(rho_samples),
    _["accept_rate"] = final_accept_rate,
    _["final_logD"] = logD_curr,
    _["final_logrho"] = logrho_curr,
    _["adaptation_history"] = List::create(
      _["acceptance_rates"] = acceptance_rates,
      _["final_logD_prop_sd"] = logD_prop_sd,
      _["final_logrho_prop_sd"] = logrho_prop_sd
    )
  );
}


// PDE
// [[Rcpp::export]]
NumericVector simulate_pde_cpp(NumericVector u0, double D, double rho, double domain_length, double dt, double t_end, int grid_size = -1){
  if (grid_size == -1) {
    grid_size = u0.size();
  }
  ReactionDiffusionSolver1D solver(D, rho, domain_length, grid_size, dt);
  std::vector<double> u0_std = as<std::vector<double>>(u0);
  std::vector<double> result = solver.solve_crank_nicolson(u0_std, t_end);
  return wrap(result);
}
