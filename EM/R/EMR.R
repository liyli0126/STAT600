#' EM optimization via Rcpp
#' 
#' @useDynLib EM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name EMR
NULL


#' Mixture of exponentials with EM optimization
#' 
#' @param p0 The first parameter
#' @param lambda0 The second parameter
#' @param mu0 The third parameter
#' @param y Observed data vector
#' @param tol Tolerance for convergence
#' @param max_iter Maximum number of iterations
#' @return List containing the updated parameters and number of iterations
#' @export
EM_mix <- function(y, p0, lambda0, mu0, tol = 1e-6, max_iter =1000){
# Check inputs
if(!is.numeric(p0) || !is.numeric(lambda0) || !is.numeric(mu0)){
  stop("Bad output: p0, lambda0 and mu0 must be numeric scalars")
}
if(!is.numeric(y) ){
  stop("Bad output: y must be vectors")
}
if(length(y) == 0){
  stop("Bad input: y must not be empty")
}

# Call the C++ function
output <- EMCpp(y, p0, lambda0, mu0, tol, max_iter)

return(output)
}