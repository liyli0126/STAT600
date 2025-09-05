#' Simple Linear Regression via Rcpp
#'
#' @param x numeric vector
#' @param y numeric vector
#' @return A list of regression results
#' @export
SimpLinR <- function(x, y){
  # Check inputs
  if(!is.numeric(x) || !is.numeric(y)){
    stop("Bad output: x and y must be vectors")
  }
  if (length(x) != length(y)){
    stop("Bad output: x and y must have the same length")
  }

  # Compute t_val
  n <- length(x)
  t_val <- qt(0.975, df = n - 2)  # 95% confidence interval

  # Call  the C++ function
  output <- SimpLinCpp(x, y, t_val)

  # Return the list
  return(output)
}
