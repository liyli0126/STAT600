#' Monte Carlo method via Rcpp
#'
#' This function provides various Monte Carlo methods for estimating E[X^2] where
#' X has density proportional to exp{-|x|^3/3}.
#'
#' @param method Character string specifying the Monte Carlo method to use.
#'               Options: "importance" (importance sampling), "rejection" (rejection sampling),
#'               "sir" (sampling importance resampling), "riemann" (Riemann strategy).
#' @param n_samples Integer, number of samples to generate (for SIR, this is the resample size).
#' @param n_reps Integer, number of replications for comparison (only used when method = "comparison").
#' @param m Integer, initial sample size for SIR method (default = 10 * n_samples).
#'
#' @return A list containing the estimate and method-specific information.
#'
#' @export
MCR <- function(method = c("importance", "rejection", "sir", "riemann", "comparison"),
                n_samples = 1000, n_reps = 100, m = NULL) {
  
  # Input validation
  method <- match.arg(method)
  
  if (!is.numeric(n_samples) || length(n_samples) != 1 || n_samples <= 0) {
    stop("n_samples must be a positive integer")
  }
  
  if (!is.numeric(n_reps) || length(n_reps) != 1 || n_reps <= 0) {
    stop("n_reps must be a positive integer")
  }
  
  n_samples <- as.integer(n_samples)
  n_reps <- as.integer(n_reps)
  
  # Method-specific warnings and adjustments
  if (method == "sir") {
    if (is.null(m)) {
      m <- 10 * n_samples
      warning("For SIR method, using default m = 10 * n_samples. ",
              "Consider adjusting this ratio based on your specific needs.")
    }
    if (m < n_samples) {
      warning("For SIR method, m should be larger than n_samples. ",
              "Current ratio m/n_samples = ", round(m/n_samples, 2))
    }
  }
  
  if (method == "riemann") {
    warning("Riemann strategy requires samples from another method. ",
            "Using rejection sampling internally to generate samples.")
  }
  if (n_samples < 100 && method != "comparison") {
    warning("Small sample size (n_samples = ", n_samples, ") may lead to inaccurate estimates. ",
            "Consider increasing n_samples for better precision.")
  }
  
  # Call the appropriate C++ function
  result <- switch(method,
                   "importance" = {
                     importance_sampling_std(n_samples)
                   },
                   "rejection" = {
                     rejection_sampling(n_samples)
                   },
                   "sir" = {
                     if (is.null(m)) m <- 10 * n_samples
                     sir_sampling(m, n_samples)
                   },
                   "riemann" = {
                     # For Riemann strategy, we need samples first
                     rs_result <- rejection_sampling(n_samples)
                     samples <- rs_result$samples
                     riemann_est <- riemann_estimator(samples)
                     list(
                       estimate = riemann_est,
                       samples = samples,
                       method = "riemann"
                     )
                   },
                   "comparison" = {
                     if (n_reps < 10) {
                       warning("Small number of replications (n_reps = ", n_reps, 
                               ") may not provide reliable performance comparisons.")
                     }
                     comparison_experiment(n_reps, n_samples)
                   }
  )
  
  # Add method information to result
  if (method != "comparison") {
    result$method <- method
    result$n_samples <- n_samples
  }
  
  return(result)
}

#' Extended Monte Carlo Comparison with Different Sample Sizes
#'
#' Compare Monte Carlo methods across different sample sizes.
#'
#' @param sample_sizes Numeric vector of sample sizes to test.
#' @param n_reps Integer, number of replications for each sample size.
#' @param methods Character vector of methods to compare.
#'
#' @return A data frame with comparison results.
#'
#' @examples
#' \dontrun{
#' # Compare methods with different sample sizes
#' results <- mc_comparison_study(
#'   sample_sizes = c(100, 500, 1000, 5000),
#'   n_reps = 50,
#'   methods = c("importance", "rejection", "sir")
#' )
#' }
#'
#' @export
mc_comparison_study <- function(sample_sizes = c(100, 500, 1000, 5000),
                                n_reps = 50,
                                methods = c("importance", "rejection", "sir", "riemann")) {
  
  # Input validation
  if (!all(sample_sizes > 0)) {
    stop("Bad input: All sample sizes must be positive")
  }
  
  if (!all(methods %in% c("importance", "rejection", "sir", "riemann"))) {
    stop("Bad input: Invalid method specified. Choose from: 'importance', 'rejection', 'sir', 'riemann'")
  }
  
  if (n_reps < 5) {
    warning("Bad input: Very few replications (n_reps = ", n_reps, 
            "). Results may not be reliable.")
  }
  
  # Initialize results data frame
  results <- data.frame()
  
  # Progress message
  total_iterations <- length(sample_sizes) * length(methods) * n_reps
  current_iteration <- 0
  
  for (n in sample_sizes) {
    for (method in methods) {
      method_results <- numeric(n_reps)
      
      for (rep in 1:n_reps) {
        current_iteration <- current_iteration + 1
        
        # Progress update
        if (current_iteration %% 100 == 0) {
          message("Progress: ", round(current_iteration/total_iterations * 100, 1), "%")
        }
        
        tryCatch({
          result <- monte_carlo_estimation(method, n_samples = n, n_reps = 1)
          method_results[rep] <- result$estimate
        }, error = function(e) {
          warning("Error in ", method, " with n=", n, ", rep=", rep, ": ", e$message)
          method_results[rep] <- NA
        })
      }
      
      # Store results
      temp_df <- data.frame(
        SampleSize = n,
        Method = method,
        Estimate = method_results,
        Variance = var(method_results, na.rm = TRUE),
        Mean = mean(method_results, na.rm = TRUE)
      )
      
      results <- rbind(results, temp_df)
    }
  }
  
  return(results)
}

