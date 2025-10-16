// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <Rcpp.h>
using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// Five MC method implementation
// returning estimation
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// Target density function
// [[Rcpp::export]]
double q_target(double x){
  return exp(-pow(fabs(x), 3) / 3.0);
}

// Variance function
//[[Rcpp::export]]
NumericVector colVars(NumericMatrix mat){
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  NumericVector vars(ncol);
  
  for (int j = 0; j < ncol; j++){
    double mean = 0.0;
    double sum_sq = 0.0;
    
    for(int i = 0; i < nrow; i++){
      mean += mat(i,j);
    }
    mean /= nrow;
    
    for(int i = 0; i < nrow; i++){
      double diff = mat(i, j) - mean;
      sum_sq += diff*diff;
    }
    vars[j] = sum_sq / (nrow - 1.0);
  }
  return vars;
}

// Importance sampling - with standard weight
// [[Rcpp::export]]
List importance_sampling_std(int n){
  // Use standard normal distribution as envelope
  NumericVector samples = rnorm(n);
  NumericVector weights(n);
  
  for (int i = 0; i < n; i++){
    double x = samples[i];
    weights[i] = q_target(x) / R::dnorm(x, 0.0, 1.0, false);
  }
  
  // Standardized weight
  double sum_weights = sum(weights);
  NumericVector std_weights = weights / sum_weights;
  
  // Estimate E[X^2]
  double estimate = 0.0;
  for (int i = 0; i < n; i++){
    estimate += std_weights[i] * pow(samples[i], 2);
  }
  
  return List::create(
    _["estimate"] = estimate,
    _["samples"] = samples,
    _["weights"] = std_weights
  );
}

// Rejection sampling
// [[Rcpp::export]]
List rejection_sampling(int n){
  // Use normal distribution as the envelope
  double M = 2; // M >=1
  
  NumericVector accepted_samples;
  int attempts = 0;
  
  while (accepted_samples.size() < n){
    attempts++;
    double x = rnorm(1)[0];
    double u = runif(1)[0]; // uniform distribution
    
    double ratio = q_target(x) / (R::dnorm(x, 0.0, 1.0, false) * M);
    
    if (u <= ratio){
      accepted_samples.push_back(x);
    }
  }  
    // Estimate E[X^2]
    double estimate = 0.0;
    for (int i = 0; i < n; i++){
      estimate += pow(accepted_samples[i], 2);
    }
    estimate /= n;

  return List::create(
    _["estimate"] = estimate,
    _["samples"] = accepted_samples,
    _["acceptance_rate"] = double(n) / attempts
  );
}


// Sampling importance resampling (SIR)
// [[Rcpp::export]]
List sir_sampling(int m, int n_resample){
  // Define g(x) as envelope
  NumericVector candidates = rnorm(m);
  NumericVector weights(m);
  
  for(int i = 0; i < m; i++){
    double x = candidates[i];
    weights[i] = q_target(x) / R::dnorm(x, 0.0, 1.0, false);
  }
  
  // Standardized weights
  NumericVector std_weights = weights / sum(weights);
  
  // Resampling based on weights
  NumericVector resampled = sample(candidates, n_resample, true, std_weights);
  
  // Estimate E[X^2]
  double estimate = 0.0;
  for (int i = 0; i < n_resample; i++){
    estimate += pow(resampled[i], 2);
  }
  estimate /= n_resample;
  
  return List::create(
    _["estimate"] = estimate,
    _["samples"] = resampled,
    _["effective_sampling_size"] = 1.0 / sum(pow(std_weights, 2))
  );
}


// Riemann sum strategy
//[[Rcpp::export]]
double riemann_estimator(NumericVector samples){
  int n = samples.size();
  NumericVector sorted_samples = clone(samples);
  std::sort(sorted_samples.begin(), sorted_samples.end()); // ordered sample
  
  double numerator = 0.0;
  double denumerator = 0.0;
  
  for (int i = 0; i < n-1; i++){
    double delta_x = sorted_samples[i+1] - sorted_samples[i]; // X_[i+1] - X_[i]
    double h_val = pow(sorted_samples[i], 2); // h(X_[i]) =X_[i]^2
    double q_val = q_target(sorted_samples[i]); // q(X_[i])
    
    // Sum
    numerator += delta_x * h_val * q_val;
    denumerator += delta_x * q_val;
  }
  
  return numerator / denumerator;
}

// Comparison
//[[Rcpp::export]]
List comparison(int n_reps, int n_samples){
  NumericMatrix results(n_reps, 4);
  NumericVector times(4);
  for(int rep = 0; rep < n_reps; rep++){
    // Importance sampling
    clock_t start = clock();
    List is_result = importance_sampling_std(n_samples);
    times[0] += double(clock() - start) / CLOCKS_PER_SEC;
    results(rep, 0) = as<double>(is_result["estimate"]);
    
    // Rejection sampling
    start = clock();
    List rs_result = rejection_sampling(n_samples);
    times[1] += double(clock() - start) / CLOCKS_PER_SEC;
    results(rep, 1) = as<double>(rs_result["estimate"]);
    
    // SIR sampling
    start = clock();
    List sir_result = sir_sampling(10*n_samples, n_samples); // m = 10*n
    times[2] += double(clock() - start) / CLOCKS_PER_SEC;
    results(rep, 2) = as<double>(sir_result["estimate"]);
    
    // Riemann sum strategy with rejection sampling as output
    start = clock();
    NumericVector rs_sample = rs_result["samples"];
    double riemann_est = riemann_estimator(rs_sample);
    times[3] += double(clock() - start) / CLOCKS_PER_SEC;
    results(rep, 3) = riemann_est;
  }
  
  times = times/n_reps;
  
  // Compute mean and variance
  NumericVector mean_estimates = colMeans(results);
  NumericVector variance_estimates = colVars(results);
  
  return List::create(
    _["estimate"] = results,
    _["mean_times"] = times,
    _["mean_estimates"] = mean_estimates,
    _["variance_estimates"] = colVars(results)
  );
}




