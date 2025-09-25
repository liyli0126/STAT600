#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// Function to compute l'(theta)
double l_prime(double theta, const arma::vec& x){
  double sum = 0.0;
  for(unsigned int i = 0; i < x.n_elem; i++){
    double diff = x(i) - theta;
    sum += 2.0 * diff/ (1.0 +  diff * diff);
  }
  return sum;
}


// Function to compute l''(theta)
double l_double_prime(double theta, const arma::vec& x){
  double sum = 0.0;
  for(unsigned int i = 0; i < x.n_elem; i++){
    double diff = x(i) - theta;
    double denom = 1.0 + diff * diff;
    sum += -2.0 * ( 1 - diff * diff) / (denom * denom);
  }
  return sum;
}


// Bisection method
// [[Rcpp::export]]
List bisection(double a, double b, const arma::vec& x, double tol = 1e-6, int max_iter = 1000){
  double ga = l_prime(a, x);
  double gb = l_prime(b, x);
  if (ga * gb > 0){
    stop("Bad inputs: Initial values a and b must have opposite signs");
  }
  double c = 0.0;
  int iter = 0 ;
  for(iter = 0; iter < max_iter; iter++){
    c = (a + b) / 2.0;
    double gc = l_prime(c, x);
    if (std::abs(gc) < tol|| (b - a) / 2.0 < tol){
      break;
    }
    if (ga * gc < 0){
      b = c; 
      gb = gc;
    } else{
      a = c;
      ga = gc;
    }
  }
  return List::create(_["theta"] = c, _["iteration"] = iter);
}

// Newton-Raphson method
// [[Rcpp::export]]
List newton(double theta0, const arma::vec& x, double tol = 1e-6, int max_iter = 1000){
  double theta = theta0;
  int iter = 0;
  double change = 1.0;
  
  while(iter < max_iter && change > tol){
    double gd = l_prime(theta, x);
    double gdd = l_double_prime(theta, x);
    if (gdd == 0){
      stop("Bad: Divided by zero in Newton's method" );
    }
    double theta_new = theta - gd / gdd;
    change = std::abs(theta_new - theta);
    theta = theta_new;
    iter++;
  }
  return List::create(_["theta"] = theta, _["iteration"] = iter);
}

// Fisher scoring
// [[Rcpp::export]]
List fisher_scoring(double theta0, const arma::vec& x, double tol = 1e-6, int max_iter = 1000){
  double theta = theta0;
  int iter = 0;
  double change = 1.0;
  double n = x.n_elem;
  double I = n / 2.0; // Fisher information of Cauchy distribution is n/2
  while(iter < max_iter && change > tol){
    double gd = l_prime(theta, x);
    if(std::abs(gd) < tol){
      break;
    }
    double theta_new = theta + gd / I;
    change = std::abs(theta_new - theta);
    theta = theta_new;
    iter++;
  }
  return List::create(_["theta"] = theta, _["iteration"] = iter);
}

// Secant method
// [[Rcpp::export]]
List secant(double theta0, double theta1, const arma::vec& x, double tol = 1e-6, int max_iter = 1000){
  double theta_prev = theta0;
  double theta_curr = theta1;
  double g_prev = l_prime(theta_prev, x);
  double g_curr = l_prime(theta_curr, x);
  int iter = 0;
  double change = 1.0;
  while(iter < max_iter && change > tol){
    if (g_prev == g_curr){
      break;
    }
    double theta_next = theta_curr - g_curr * (theta_curr - theta_prev) / (g_curr - g_prev);
    change = std::abs(theta_next - theta_curr);
    theta_prev = theta_curr;
    theta_curr = theta_next;
    g_prev = g_curr;
    g_curr = l_prime(theta_curr, x);
    iter++;
  }
  return List::create(_["theta"] = theta_curr, _["iteration"] = iter);
}


// Predicted probability function
// [[Rcpp::export]]
arma::vec logistic_predict(const arma::mat& X, const arma::vec& beta){
  vec eta = X * beta;
  return 1.0 / (1.0 + exp(-eta));
}

// Predicted log MLE function
// [[Rcpp::export]]
double logistic_lik(const arma::mat& X, const arma::vec& response, const arma::vec& beta){
  vec eta = X * beta;
  double ll = 0.0;
  for(unsigned int i = 0; i < response.n_elem; i++){
    double exp_neg_eta = exp(-eta(i));
    ll += response(i) * eta(i) - log(1.0 + exp_neg_eta); // response ={0,1}
  }
  return ll;
}

//  Multivaribale Newton's method
// [[Rcpp::export]]
List multi_newton(const arma::mat& X, const arma::vec& y, double tol = 1e-6, int max_iter = 1000){
  // int n_obs = X.n_rows;
  int n_params = X.n_cols;  // Number of parameters
  
  // Initialize beta to zeros
  vec beta = zeros<vec>(n_params);
  int iter = 0;
  double change = 1;
  bool converged = false;
  
  while(iter < max_iter && change > tol){
    // Compute linear predictor and probabilities
    vec eta = X * beta;
    vec p = 1.0 / (1.0 + exp(-eta));
    
    // Compute gradient
    vec gradient = X.t() * (y - p);
    
    // Compute Hessian
    mat W = diagmat(p % (1.0 - p));  // use % for element-wise multiplication
    mat hessian = -X.t() * W * X;
    
    // Newton-Raphson update
    vec beta_new;
    try{
      beta_new  = beta - solve(hessian, gradient);
    }catch(const std::exception& e){
      beta_new = beta - 0.01*gradient;
      Rcpp::warning("Hessian matrix was singular, used regularization");
    }
    
    // Check convergence
    change = norm(beta_new - beta, 2);
    beta = beta_new;
    iter++;
    
    // Save the last results
    converged = (change <=tol);
  }
  
  // Compute output results
  vec eta_final = X * beta;
  vec p_final = 1.0 / (1.0 + exp(-eta_final));
  vec w_final = p_final % (1.0 - p_final);
  mat hessian_final = -X.t() * diagmat(w_final) * X;
  
  // Compute standard errors and log-likelihood
  mat cov_matrix;
  try {
    cov_matrix = -inv(hessian_final);
  } catch(const std::exception& e) {
    // If failure to compute inverse of hessian matrix
    mat hessian_reg = hessian_final + 1e-6 * eye<mat>(n_params, n_params);
    cov_matrix = -inv(hessian_reg);
    Rcpp::warning("Hessian matrix was singular, used regularization for inverse");
  }
  vec se = sqrt(cov_matrix.diag());
  
  // Compute log-likelihood
  double log_lik_final = logistic_lik(X, y, beta);
  
  
  return List::create(
    _["coefficients"] = beta,
    _["standard_errors"] = se,
    _["log_likelihood"] = log_lik_final,
    _["iterations"] = iter,
    _["cov_matrix"] = cov_matrix,
    _["converged"] = converged,
    _["fitted_values"] = p_final
  );
}

// Help function
// Compute score function (logistic gradient function)
// [[Rcpp::export]]
arma::vec score(const arma::vec& beta, const arma::mat& X, const arma::vec& y){
  vec eta = X * beta;
  vec p = 1.0 / (1.0 + exp(-eta));
  return X.t() * (y - p);
}

// Compute Hessian matrix
// [[Rcpp::export]]
arma::mat hessian(const arma::vec& beta, const arma::mat& X, const arma::vec& y){
  vec eta = X * beta;
  vec p = 1.0 / (1.0 + exp(-eta));
  vec w = p % (1.0 - p);  
  mat W = diagmat(w);
  return -X.t() * W * X;
}