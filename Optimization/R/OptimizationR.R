#' Univariable optimization via Rcpp
#' 
#' @useDynLib Optimization, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name OptimizationR
NULL


#' Cauchy Distribution Log-Likelihood and Derivatives
#' 
#' Helper functions for Cauchy distribution optimization
#' 
#' @param theta Location parameter
#' @param x Numeric vector of data
#' @return Function value
#' @export

cauchy_log_likelihood <- function(theta, x) {
  if(!is.numeric(theta) || !is.numeric(x)){
    stop("Bad output: theta and x must be numeric scalars")
  }
  output <- cauchy_log_likelihood_cpp(theta, x)
  return(output)
}

#' @rdname cauchy_log_likelihood
#' @export
cauchy_score <- function(theta, x) {
  if(!is.numeric(theta) || !is.numeric(x)){
    stop("Bad output: theta and x must be numeric scalars")
  }
  output <- cauchy_score_cpp(theta, x)
  return(output)
}

#' @rdname cauchy_log_likelihood
#' @export
cauchy_observed_info <- function(theta, x) {
  if(!is.numeric(theta) || !is.numeric(x)){
    stop("Bad output: theta and x must be numeric scalars")
  }
  output <- cauchy_observed_info_cpp(theta, x)
  return(output)
}

# Helper functions for input validation
validate_numeric_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1) {
    stop(paste(name, "must be a numeric scalar"))
  }
  if (is.na(value)) {
    stop(paste(name, "cannot be NA"))
  }
}

validate_numeric_vector <- function(value, name) {
  if (!is.numeric(value)) {
    stop(paste(name, "must be a numeric vector"))
  }
  if (length(value) == 0) {
    stop(paste(name, "cannot be empty"))
  }
  if (any(is.na(value))) {
    stop(paste(name, "cannot contain NA values"))
  }
}

#' Bisection Method for Root Finding
#' 
#' @param a Lower bound of interval
#' @param b Upper bound of interval
#' @param x Data vector
#' @param tol Tolerance for convergence
#' @param max_iter Maximum number of iterations
#' @return List containing the root and number of iterations
#' @export
bisection <- function(a, b, x, tol = 1e-6, max_iter =1000){
  # Check inputs
  if(!is.numeric(a) || !is.numeric(b)){
    stop("Bad output: a and b must be numeric scalars")
  }
  if (length(a) != 1 || length(b)!= 1){
    stop("Bad output: a and b must be scalars")
  }
  if (a >= b){
    stop("Bad output: a must be less than b")
  }
  if(!is.numeric(x) ){
    stop("Bad output: x must be vectors")
  }
  if(length(x) == 0){
    stop("Bad input: x must not be empty")
  }
  
  # Call the C++ function
  output <- bisectionCpp(a, b, x, tol, max_iter)
  
  return(output)
}

#' Newton-Raphson Method for Root Finding
#' 
#' @param theta0 Initial guess
#' @param x Data vector
#' @param tol Tolerance for convergence
#' @param max_iter Maximum number of iterations
#' @return List containing the root and number of iterations
#' @export
newton <- function(theta0, x, tol = 1e-6, max_iter =1000){
  # Check inputs
  if(!is.numeric(theta0)){
    stop("Bad output: theta0 must be numeric scalars")
  }
  if (length(theta0) != 1 ){
    stop("Bad output: theta0 must be scalars")
  }
  if(!is.numeric(x) ){
    stop("Bad output: x must be vectors")
  }
  if(length(x) == 0){
    stop("Bad input: x must not be empty")
  }
  
  # Call the C++ function
  output <- newtonCpp(theta0, x, tol, max_iter)
  
  return(output)
}


#' Fisher Scoring Method for Root Finding
#' 
#' @param theta0 Initial guess
#' @param x Data vector
#' @param tol Tolerance for convergence
#' @param max_iter Maximum number of iterations
#' @return List containing the root and number of iterations
#' @export
fisher_scoring <- function(theta0, x, tol = 1e-6, max_iter =1000){
  # Check inputs
  if(!is.numeric(theta0)){
    stop("Bad output: theta0 must be numeric scalars")
  }
  if (length(theta0) != 1 ){
    stop("Bad output: theta0 must be scalars")
  }
  if(!is.numeric(x) ){
    stop("Bad output: x must be vectors")
  }
  if(length(x) == 0){
    stop("Bad input: x must not be empty")
  }
  
  # Call the C++ function
  output <- fisher_scoringCpp(theta0, x, tol, max_iter)
  
  return(output)
}


#' Secant Method for Root Finding
#' 
#' @param theta0 First Initial guess
#' @param theta1 Second Initial guess
#' @param x Data vector
#' @param tol Tolerance for convergence
#' @param max_iter Maximum number of iterations
#' @return List containing the root and number of iterations
#' @export
secant <- function(theta0, x, tol = 1e-6, max_iter =1000){
  # Check inputs
  if(!is.numeric(theta0)||!is.numeric(theta1)){
    stop("Bad output: theta0 and theta1 must be numeric scalars")
  }
  if (length(theta0) != 1 || length(theta1) != 1){
    stop("Bad output: theta0 and theta1 must be scalars")
  }
  if(!is.numeric(x) ){
    stop("Bad output: x must be vectors")
  }
  if(length(x) == 0){
    stop("Bad input: x must not be empty")
  }
  
  # Call the C++ function
  output <- secantCpp(theta0, theta1, x, tol, max_iter)
  
  return(output)
}


#' Multivaribale Newton-Raphson Method for Root Finding
#' 
#' @param beta Initial guess
#' @param X Data matrix
#' @param tol Tolerance for convergence
#' @param max_iter Maximum number of iterations
#' @return List containing the root and number of iterations
#' @export
multi_newton <- function(beta, X, tol = 1e-6, max_iter =1000){
  # Check inputs
  if(!is.numeric(beta)){
    stop("Bad output: beta must be vectors")
  }
  if(length(beta) == 0){
    stop("Bad input: beta must not be empty")
  }
  if(!is.matrix(X)){
    stop("Bad output: X must be vectors")
  }
  if(dim(X) == 0){
    stop("Bad input: X must not be empty")
  }
  if (X.n_rows != y.n_elem) {
    stop("X and y must have the same number of observations");
  }
  # Call the C++ function
  output <- multi_newtonCpp(beta, X, tol, max_iter)
  
  return(output)
}

