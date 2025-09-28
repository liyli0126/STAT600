#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// EM optimization for a mixture of exponentials
// [[Rcpp::export]]
List EM_mix(arma::vec y, double p0, double lambda0, double mu0, double tol = 1e-6, int max_iter = 1000){
  // Initialization
  int n = y.n_elem;
  double p = p0;
  double lambda = lambda0;
  double mu = mu0;
  double diff = 1.0;  // For convergence criterion
  double iter = 0;
  
  while(diff > tol && iter < max_iter){
    // E step: compute delta_{i}
    arma::vec delta(n);
    for(int i=0; i < n; i++){
      double num = p * lambda* exp(-lambda* y(i));
      double denom = p * lambda* exp(-lambda* y(i)) + (1-p) * mu * exp(-mu * y(i));
      delta(i) = num / denom;
    }
    
    // M step: update parameters
    double p_new = arma::mean(delta);
    double lambda_new = arma::sum(delta) / arma::sum(delta % y);
    double mu_new =  arma::sum(1 - delta) / arma::sum((1 - delta) % y);
    
    // Compute difference of parameters: consider L1 norm
    diff = fabs(p_new - p) + fabs(lambda_new - lambda) + fabs(mu_new - mu);
    
    // Update parameters
    p = p_new;
    lambda = lambda_new;
    mu = mu_new;
    iter++;
  }
  
  return List::create(
    _["p"] = p,
    _["lambda"] = lambda,
    _["mu"] = mu,
    _["iter"] = iter,
    _["converged"] = (iter < max_iter)
  );
}