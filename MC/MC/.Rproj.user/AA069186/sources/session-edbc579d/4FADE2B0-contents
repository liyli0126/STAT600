# Load necessary libraries
library(EM)

library(ggplot2)
library(dplyr) # Summary data
library(kableExtra)

# Simulate data function
simulate_data <- function(n, p, lambda, mu){
  # Generate missing variable
  delta <- rbinom(n, 1, p)
  
  # Compute observed data y
  y <- numeric(n)
  for (i in 1:n){
    if (delta[i] == 1){
      y[i] <- rexp(1, lambda)
    }else{
      y[i] <- rexp(1, mu)
    }
  }
  return(y)
}

# Set parameters
n <- 100
p_true <- 0.25
lambda_true <- 1
mu_true <- 2
n_datasets <- 100

# Simulate 100 datasets
set.seed(26)
datasets <- lapply(1:n_datasets, function(i) simulate_data(n, p_true, lambda_true, mu_true))

df <- data.frame(y = datasets[[1]])

# Check if correct
p <- ggplot(df, aes(x = y)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightgray", color = "black") +
  geom_density(color = "blue", linewidth = 1.2) +
  ggtitle("A simulated dataset") +
  xlab("y")

# Ensure if a folder exist to save figures, if not, build one
setwd("/Users/yinglili/Desktop/STAT600/")  
dir.create("Figure", showWarnings = FALSE)

ggsave("Figure/simulated_dataset.png", plot = p, width = 8, height = 6, dpi = 300)



# Initialize storage data
results <- matrix(NA, nrow = n_datasets, ncol = 3)
colnames(results) <- c("p", "lambda", "mu")
iterations <- numeric(n_datasets)
converged <- logical(n_datasets)

# Loop each dataset
for(i in 1:n_datasets){
  y <- datasets[[i]]
  output <- EM_mix(y, p0 = 0.3, lambda0 = 1.0, mu0 = 2.5)
  if(!is.null(output) && output$converged){
    results[i,] <- c(output$p, output$lambda, output$mu)
    iterations[i] <- output$iter
    converged[i] <- output$converged
  } else {
    converged[i] <- FALSE
  }
}

# Output the mean value of estimations
summary_df <- data.frame(
  Convergence_rate    = mean(converged, na.rm = TRUE),
  Average_iterations  = mean(iterations, na.rm = TRUE),
  Estimate_p          = mean(output$p, na.rm = TRUE),
  Estimate_lambda     = mean(output$lambda, na.rm = TRUE),
  Estimate_mu         = mean(output$mu, na.rm = TRUE)
)

saveRDS(summary_df, file = "~/Desktop/STAT600/Figure/hw3_summary_df.rds") 


# Convert results to a data frame (will filter after computing SE)
df_results <- as.data.frame(na.omit(results))

# p estimates 
p1 <- ggplot(df_results, aes(x = p)) +
  geom_histogram(bins = 18, fill = "skyblue", alpha = 0.8, color = "black") +
  geom_vline(xintercept = p_true, color = "red", linewidth = 1) +
  geom_vline(xintercept = mean(df_results$p, na.rm = TRUE), 
             color = "blue", linetype = "dashed", linewidth = 1) +
  labs(
    title = expression(paste("Distribution of ", hat(p))),
    x = expression(hat(p)),
    y = "Count"
  ) +
  theme_minimal(base_size = 10)

# lambda estimates
p2 <- ggplot(df_results, aes(x = lambda)) +
  geom_histogram(bins = 18, fill = "lightgreen", alpha = 0.8, color = "black") +
  geom_vline(xintercept = lambda_true, color = "red", linewidth = 1) +
  geom_vline(xintercept = mean(df_results$lambda, na.rm = TRUE), 
             color = "blue", linetype = "dashed", linewidth = 1) +
  labs(
    title = expression(paste("Distribution of ", hat(lambda))),
    x = expression(hat(lambda)),
    y = "Count"
  ) +
  theme_minimal(base_size = 10)

# mu estimates
p3 <- ggplot(df_results, aes(x = mu)) +
  geom_histogram(bins = 18, fill = "lightyellow", alpha = 0.8, color = "black") +
  geom_vline(xintercept = mu_true, color = "red", linewidth = 1) +
  geom_vline(xintercept = mean(df_results$mu, na.rm = TRUE), 
             color = "blue", linetype = "dashed", linewidth = 1) +
  labs(
    title = expression(paste("Distribution of ", hat(mu))),
    x = expression(hat(mu)),
    y = "Count"
  ) +
  theme_minimal(base_size = 10)

# Combine the three plots horizontally
p_combined <- (p1) / (p2) / p3 +
  plot_annotation(
    title = "Parameter Estimation Distributions",
    subtitle = "Red: True value   |   Blue dashed: Mean estimate")

# Show combined plot
p_combined

# Save figures
ggsave("Figure/hw3_parameter_estimation_distributions.png", 
       plot = p_combined, width = 5, height = 8, dpi = 300)


# ============================================================================
# Three Methods for Standard Error Estimation
# ============================================================================

# Method 1: Bootstrap Method
# ============================================================================
bootstrap_se <- function(y, B = 50){
  n <- length(y)
  boot_results <- matrix(NA, nrow = B, ncol = 3)
  colnames(boot_results) <- c("p", "lambda", "mu")
  
  for (b in 1:B) {
    # Sampling with replacement from the original data
    idx <- sample(1:n, n, replace = TRUE)
    y_boot <- y[idx] 
    
    # Apply EM algorithm
    output <- tryCatch({
      EM_mix(y_boot, p0 = 0.28, lambda0 = 1.2, mu0 = 1.8)
    }, error = function(e){
      NULL
    })
    
    if(!is.null(output) && !any(is.na(unlist(output[c("p", "lambda", "mu")])))){
      boot_results[b, ] <- unlist(output[c("p", "lambda", "mu")])
    }
  }
  
  # Remove failure bootstrap sampling
  boot_results <- boot_results[complete.cases(boot_results), ]
  
  if (nrow(boot_results) < 2){
    return(rep(NA, 3))
  }
  
  # Compute standard error of bootstrap (equals to the estimate se of bootstrap)
  se <- apply(boot_results, 2, sd, na.rm = TRUE)
  return(se)
}

# Method 2: Louis Method (Missing Information Principle)
# ============================================================================
louis_se <- function(y, p, lambda, mu) {
  n <- length(y)
  
  # Calculate conditional expectation of delta_i
  delta_hat <- sapply(y, function(y_i) {
    num <- p * lambda * exp(-lambda * y_i)
    denom <- num + (1 - p) * mu * exp(-mu * y_i)
    if (denom == 0) return(0.5)
    num / denom
  })
  
  delta_hat <- pmin(pmax(delta_hat, 1e-10), 1 - 1e-10)
  
  # Calculate complete information matrix I_C 
  I_pp <- sum(delta_hat) / p^2 + sum(1 - delta_hat) / (1 - p)^2
  I_lambda_lambda <- sum(delta_hat) / lambda^2
  I_mu_mu <- sum(1 - delta_hat) / mu^2
  I_C <- diag(c(I_pp, I_lambda_lambda, I_mu_mu))
  
  # Calculate missing information matrix components
  H_pp <- 0
  H_plambda <- 0
  H_pmu <- 0
  H_lambdalambda <- 0
  H_mumu <- 0
  H_lambdamu <- 0
  
  for(i in 1:n) {
    y_i <- y[i]
    d_i <- delta_hat[i]
    term <- p * lambda * exp(-lambda * y_i) + (1-p) * mu * exp(-mu * y_i)
    
    # Derivatives of delta_hat with respect to parameters
    d_dp <- (lambda * exp(-lambda * y_i) - mu * exp(-mu * y_i)) / term - 
      d_i * (lambda * exp(-lambda * y_i) - mu * exp(-mu * y_i)) / term
    
    d_dlambda <- (p * exp(-lambda * y_i) * (1 - lambda * y_i)) / term -
      d_i * (p * exp(-lambda * y_i) * (1 - lambda * y_i)) / term
    
    d_dmu <- ((1-p) * exp(-mu * y_i) * (1 - mu * y_i)) / term -
      d_i * ((1-p) * exp(-mu * y_i) * (1 - mu * y_i)) / term
    
    # Accumulate second derivatives
    H_pp <- H_pp + (d_dp^2) / (d_i * (1 - d_i))
    H_lambdalambda <- H_lambdalambda + (d_dlambda^2) / (d_i * (1 - d_i))
    H_mumu <- H_mumu + (d_dmu^2) / (d_i * (1 - d_i))
    H_plambda <- H_plambda + (d_dp * d_dlambda) / (d_i * (1 - d_i))
    H_pmu <- H_pmu + (d_dp * d_dmu) / (d_i * (1 - d_i))
    H_lambdamu <- H_lambdamu + (d_dlambda * d_dmu) / (d_i * (1 - d_i))
  }
  
  # Construct missing information matrix
  I_M <- matrix(c(H_pp, H_plambda, H_pmu,
                  H_plambda, H_lambdalambda, H_lambdamu,
                  H_pmu, H_lambdamu, H_mumu), nrow = 3, ncol = 3)
  
  # Observed information matrix
  I_Y <- I_C - I_M
  
  # Calculate standard errors with error handling
  se <- tryCatch({
    I_Y_inv <- solve(I_Y)
    diag_vals <- diag(I_Y_inv)
    if (any(diag_vals <= 0)) {
      return(rep(NA, 3))
    }
    sqrt(diag_vals)
  }, error = function(e) {
    rep(NA, 3)
  })
  
  names(se) <- c("p", "lambda", "mu")
  return(se)
}

# Method 3: SEM (Supplemented EM) Method  
# ============================================================================
SEM_se_optimized <- function(y, theta_hat, max_iter_sem = 50, tol_sem = 1e-6){
  p <- length(theta_hat)
  names(theta_hat) <- c("p", "lambda", "mu")
  
  set.seed(26)
  theta_start <- theta_hat * (1 + runif(p, -0.05, 0.05))
  
  R <- matrix(0, p, p)
  R_old <- matrix(Inf, p, p)
  theta_current <- theta_start
  iter <- 0
  converged <- FALSE
  
  while(iter < max_iter_sem && !converged){
    for (j in 1:p) {
      theta_t <- theta_hat
      theta_t[j] <- theta_current[j]
      
      output <- EM_mix(y, p0 = theta_t[1], lambda0 = theta_t[2], mu0 = theta_t[3])
      Psi_theta <- c(output$p, output$lambda, output$mu)
      
      if(abs(theta_current[j] - theta_hat[j]) > 1e-10){
        for (i in 1:p) {
          R[i, j] <- (Psi_theta[i] - theta_hat[i]) / (theta_current[j] - theta_hat[j])
        }
      }
    }
    
    output_current <- EM_mix(y, p0 = theta_current[1], lambda0 = theta_current[2], mu0 = theta_current[3])
    theta_current <- c(output_current$p, output_current$lambda, output_current$mu)
    
    if(max(abs(R - R_old)) < tol_sem){
      converged <- TRUE
    }
    R_old <- R
    iter <- iter + 1
  }
  
  Psi_prime <- R
  
  delta_hat <- E_step_mix(y, theta_hat[1], theta_hat[2], theta_hat[3])
  n <- length(y)
  
  I_complete <- matrix(0, p, p)
  I_complete[1,1] <- n / (theta_hat["p"] * (1 - theta_hat["p"]))
  I_complete[2,2] <- sum(delta_hat) / (theta_hat["lambda"]^2)
  I_complete[3,3] <- sum(1 - delta_hat) / (theta_hat["mu"]^2)
  
  I_Y <- t(diag(p) - Psi_prime) %*% I_complete %*% (diag(p) - Psi_prime)
  var_theta <- solve(I_Y)
  se <- sqrt(diag(var_theta))
  names(se) <- c("p", "lambda", "mu")
  return(se)
}

# ============================================================================
# Compute Standard Errors using all Three Methods
# ============================================================================

# Initialize result matrices
bootstrap_se_results <- matrix(NA, nrow = n_datasets, ncol = 3)
louis_se_results <- matrix(NA, nrow = n_datasets, ncol = 3)
SEM_se_results <- matrix(NA, nrow = n_datasets, ncol = 3)
colnames(bootstrap_se_results) <- c("p", "lambda", "mu")
colnames(louis_se_results) <- c("p", "lambda", "mu")
colnames(SEM_se_results) <- c("p", "lambda", "mu")

# Compute standard errors for each dataset (use already computed EM results)
for (i in 1:n_datasets) {
  y <- datasets[[i]]
  
  # Bootstrap method (works for all datasets)
  bootstrap_se_results[i, ] <- tryCatch({
    bootstrap_se(y, B = 100)
  }, error = function(e) {
    rep(NA, 3)
  })
  
  # Louis and SEM methods (only for converged datasets)
  if (converged[i]) {
    # Get parameter estimates from already computed results
    theta_hat <- results[i, ]
    names(theta_hat) <- c("p", "lambda", "mu")
    
    # Louis method
    louis_se_results[i, ] <- tryCatch({
      louis_se(y, p = theta_hat["p"], lambda = theta_hat["lambda"], mu = theta_hat["mu"])
    }, error = function(e) {
      rep(NA, 3)
    })
    
    # SEM method
    SEM_se_results[i, ] <- tryCatch({
      SEM_se_optimized(y, theta_hat)
    }, error = function(e) {
      rep(NA, 3)
    })
  }
}

# ============================================================================
# Filter All Results to Keep Only Converged Datasets
# ============================================================================

valid_idx <- which(converged)
results <- results[valid_idx, ]
bootstrap_se_results <- bootstrap_se_results[valid_idx, ]
louis_se_results <- louis_se_results[valid_idx, ]
SEM_se_results <- SEM_se_results[valid_idx, ]

cat("Kept", nrow(results), "converged datasets out of", n_datasets, "total\n\n")

# ============================================================================
# Individual Method Results and Summaries
# ============================================================================

# Method 1: Bootstrap Results
# ----------------------------------------------------------------------------
bootstrap_summary <- data.frame(
  Parameter = c("p", "lambda", "mu"), 
  Method = "Bootstrap",
  Mean_SE = c(
    mean(bootstrap_se_results[, "p"], na.rm = TRUE),
    mean(bootstrap_se_results[, "lambda"], na.rm = TRUE),
    mean(bootstrap_se_results[, "mu"], na.rm = TRUE)
  ),
  SD_SE = c(
    sd(bootstrap_se_results[, "p"], na.rm = TRUE),
    sd(bootstrap_se_results[, "lambda"], na.rm = TRUE),
    sd(bootstrap_se_results[, "mu"], na.rm = TRUE)
  ),
  Median_SE = c(
    median(bootstrap_se_results[, "p"], na.rm = TRUE),
    median(bootstrap_se_results[, "lambda"], na.rm = TRUE),
    median(bootstrap_se_results[, "mu"], na.rm = TRUE)
  )
)
cat("\nSummary statistics:\n")
print(bootstrap_summary, row.names = FALSE)
saveRDS(bootstrap_summary, file = "~/Desktop/STAT600/Figure/hw3_bootstrap_summary.rds")

# Compare with empirical SE
empirical_se_temp <- apply(results, 2, sd, na.rm = TRUE)
bootstrap_comparison <- data.frame(
  Parameter = c("p", "lambda", "mu"),
  Empirical_SE = empirical_se_temp,
  Bootstrap_SE = colMeans(bootstrap_se_results, na.rm = TRUE),
  Ratio = colMeans(bootstrap_se_results, na.rm = TRUE) / empirical_se_temp
)
cat("\nComparison with Empirical SE:\n")
print(bootstrap_comparison, row.names = FALSE)
saveRDS(bootstrap_comparison, file = "~/Desktop/STAT600/Figure/hw3_bootstrap_comparison.rds")
cat("\n\n")

# Method 2: Louis Results
# ----------------------------------------------------------------------------
louis_summary <- data.frame(
  Parameter = c("p", "lambda", "mu"), 
  Method = "Louis",
  Mean_SE = c(
    mean(louis_se_results[, "p"], na.rm = TRUE),
    mean(louis_se_results[, "lambda"], na.rm = TRUE),
    mean(louis_se_results[, "mu"], na.rm = TRUE)
  ),
  SD_SE = c(
    sd(louis_se_results[, "p"], na.rm = TRUE),
    sd(louis_se_results[, "lambda"], na.rm = TRUE),
    sd(louis_se_results[, "mu"], na.rm = TRUE)
  ),
  Median_SE = c(
    median(louis_se_results[, "p"], na.rm = TRUE),
    median(louis_se_results[, "lambda"], na.rm = TRUE),
    median(louis_se_results[, "mu"], na.rm = TRUE)
  )
)
cat("\nSummary statistics:\n")
print(louis_summary, row.names = FALSE)
saveRDS(louis_summary, file = "~/Desktop/STAT600/Figure/hw3_louis_summary.rds")

# Compare with empirical SE
louis_comparison <- data.frame(
  Parameter = c("p", "lambda", "mu"),
  Empirical_SE = empirical_se_temp,
  Louis_SE = colMeans(louis_se_results, na.rm = TRUE),
  Ratio = colMeans(louis_se_results, na.rm = TRUE) / empirical_se_temp
)
cat("\nComparison with Empirical SE:\n")
print(louis_comparison, row.names = FALSE)
saveRDS(louis_comparison, file = "~/Desktop/STAT600/Figure/hw3_louis_comparison.rds")
cat("\n\n")

# Method 3: SEM Results
# ----------------------------------------------------------------------------
SEM_summary <- data.frame(
  Parameter = c("p", "lambda", "mu"), 
  Method = "SEM",
  Mean_SE = c(
    mean(SEM_se_results[, "p"], na.rm = TRUE),
    mean(SEM_se_results[, "lambda"], na.rm = TRUE),
    mean(SEM_se_results[, "mu"], na.rm = TRUE)
  ),
  SD_SE = c(
    sd(SEM_se_results[, "p"], na.rm = TRUE),
    sd(SEM_se_results[, "lambda"], na.rm = TRUE),
    sd(SEM_se_results[, "mu"], na.rm = TRUE)
  ),
  Median_SE = c(
    median(SEM_se_results[, "p"], na.rm = TRUE),
    median(SEM_se_results[, "lambda"], na.rm = TRUE),
    median(SEM_se_results[, "mu"], na.rm = TRUE)
  )
)
cat("\nSummary statistics:\n")
print(SEM_summary, row.names = FALSE)
saveRDS(SEM_summary, file = "~/Desktop/STAT600/Figure/hw3_SEM_summary.rds")

# Compare with empirical SE
SEM_comparison <- data.frame(
  Parameter = c("p", "lambda", "mu"),
  Empirical_SE = empirical_se_temp,
  SEM_SE = colMeans(SEM_se_results, na.rm = TRUE),
  Ratio = colMeans(SEM_se_results, na.rm = TRUE) / empirical_se_temp
)
cat("\nComparison with Empirical SE:\n")
print(SEM_comparison, row.names = FALSE)
saveRDS(SEM_comparison, file = "~/Desktop/STAT600/Figure/hw3_SEM_comparison.rds")
cat("\n\n")

# ============================================================================
# Standard Error Distribution Plots for Each Method
# ============================================================================

# Plot 1: Bootstrap SE Distribution
par(mfrow = c(1, 3))
for (j in 1:3) {
  param_name <- colnames(bootstrap_se_results)[j]
  hist(bootstrap_se_results[, j], 
       main = paste("Bootstrap SE Distribution\n", param_name),
       xlab = "Standard Error", 
       col = "#3498db", 
       breaks = 20,
       border = "white")
  abline(v = mean(bootstrap_se_results[, j], na.rm = TRUE), col = "red", lwd = 2)
  abline(v = empirical_se_temp[j], col = "darkgreen", lty = 2, lwd = 2)
  legend("topright", 
         legend = c("Mean Bootstrap SE", "Empirical SE"), 
         col = c("red", "darkgreen"), 
         lty = c(1, 2), 
         lwd = 2,
         cex = 0.8)
}

# Plot 2: Louis SE Distribution
par(mfrow = c(1, 3))
for (j in 1:3) {
  param_name <- colnames(louis_se_results)[j]
  hist(louis_se_results[, j], 
       main = paste("Louis SE Distribution\n", param_name),
       xlab = "Standard Error", 
       col = "#e74c3c", 
       breaks = 20,
       border = "white")
  abline(v = mean(louis_se_results[, j], na.rm = TRUE), col = "red", lwd = 2)
  abline(v = empirical_se_temp[j], col = "darkgreen", lty = 2, lwd = 2)
  legend("topright", 
         legend = c("Mean Louis SE", "Empirical SE"), 
         col = c("red", "darkgreen"), 
         lty = c(1, 2), 
         lwd = 2,
         cex = 0.8)
}

# Plot 3: SEM SE Distribution
par(mfrow = c(1, 3))
for (j in 1:3) {
  param_name <- colnames(SEM_se_results)[j]
  hist(SEM_se_results[, j], 
       main = paste("SEM SE Distribution\n", param_name),
       xlab = "Standard Error", 
       col = "#2ecc71", 
       breaks = 20,
       border = "white")
  abline(v = mean(SEM_se_results[, j], na.rm = TRUE), col = "red", lwd = 2)
  abline(v = empirical_se_temp[j], col = "darkgreen", lty = 2, lwd = 2)
  legend("topright", 
         legend = c("Mean SEM SE", "Empirical SE"), 
         col = c("red", "darkgreen"), 
         lty = c(1, 2), 
         lwd = 2,
         cex = 0.8)
}

# ============================================================================
# Combined Three-Method Summary Table
# ============================================================================
cat("=" , rep("=", 70), "\n", sep = "")
cat("Combined Summary: All Three Methods\n")
cat("=" , rep("=", 70), "\n", sep = "")

# Combine all three summary tables
all_methods_summary <- rbind(bootstrap_summary, louis_summary, SEM_summary)
cat("\nCombined Standard Error Summary:\n")
print(all_methods_summary, row.names = FALSE)
saveRDS(all_methods_summary, file = "~/Desktop/STAT600/Figure/hw3_all_methods_summary.rds")

# Success rates summary
cat("\n\nSuccess Rates (Effective Calculation Ratios):\n")
success_table <- data.frame(
  Method = c("Bootstrap", "Louis", "SEM"),
  Success_Rate = c(
    mean(!is.na(bootstrap_se_results[,1])),
    mean(!is.na(louis_se_results[,1])),
    mean(!is.na(SEM_se_results[,1]))
  )
)
print(success_table, row.names = FALSE)
saveRDS(success_table, file = "~/Desktop/STAT600/Figure/hw3_success_rates.rds")

# Combined comparison with empirical SE
cat("\n\nAll Methods Comparison with Empirical SE:\n")
empirical_se_temp <- apply(results, 2, sd, na.rm = TRUE)
all_methods_comparison <- data.frame(
  Parameter = rep(c("p", "lambda", "mu"), 3),
  Method = rep(c("Bootstrap", "Louis", "SEM"), each = 3),
  Method_SE = c(
    colMeans(bootstrap_se_results, na.rm = TRUE),
    colMeans(louis_se_results, na.rm = TRUE),
    colMeans(SEM_se_results, na.rm = TRUE)
  ),
  Empirical_SE = rep(empirical_se_temp, 3),
  Ratio = c(
    colMeans(bootstrap_se_results, na.rm = TRUE) / empirical_se_temp,
    colMeans(louis_se_results, na.rm = TRUE) / empirical_se_temp,
    colMeans(SEM_se_results, na.rm = TRUE) / empirical_se_temp
  )
)
print(all_methods_comparison, row.names = FALSE)
saveRDS(all_methods_comparison, file = "~/Desktop/STAT600/Figure/hw3_all_methods_comparison.rds")
cat("\n\n")

# ============================================================================
# Summary Statistics and Comparisons
# ============================================================================

# Calculate empirical statistics from EM estimates
true_params <- c(p = p_true, lambda = lambda_true, mu = mu_true)
avg_estimates <- colMeans(results, na.rm = TRUE)
bias <- avg_estimates - true_params
empirical_se <- apply(results, 2, sd, na.rm = TRUE)

# Average standard errors from each method
avg_bootstrap_se <- colMeans(bootstrap_se_results, na.rm = TRUE)
avg_louis_se <- colMeans(louis_se_results, na.rm = TRUE)
avg_SEM_se <- colMeans(SEM_se_results, na.rm = TRUE)

# Success rates for each method
success_bootstrap <- mean(!is.na(bootstrap_se_results[,1]))
success_louis <- mean(!is.na(louis_se_results[,1]))
success_SEM <- mean(!is.na(SEM_se_results[,1]))

# ============================================================================
# Calculate Coverage Probabilities (95% CI)
# ============================================================================

alpha <- 0.05
z_value <- qnorm(1 - alpha/2)

# Calculate coverage for each method (results already filtered)
# Bootstrap coverage
ci_lower_boot <- results - z_value * bootstrap_se_results
ci_upper_boot <- results + z_value * bootstrap_se_results
coverage_boot <- colMeans((ci_lower_boot <= matrix(rep(true_params, nrow(results)), 
                                                   nrow = nrow(results), byrow = TRUE)) & 
                            (matrix(rep(true_params, nrow(results)), 
                                    nrow = nrow(results), byrow = TRUE) <= ci_upper_boot), na.rm = TRUE)

# Louis coverage
ci_lower_louis <- results - z_value * louis_se_results
ci_upper_louis <- results + z_value * louis_se_results
coverage_louis <- colMeans((ci_lower_louis <= matrix(rep(true_params, nrow(results)), 
                                                     nrow = nrow(results), byrow = TRUE)) & 
                             (matrix(rep(true_params, nrow(results)), 
                                     nrow = nrow(results), byrow = TRUE) <= ci_upper_louis), na.rm = TRUE)

# SEM coverage
ci_lower_SEM <- results - z_value * SEM_se_results
ci_upper_SEM <- results + z_value * SEM_se_results
coverage_SEM <- colMeans((ci_lower_SEM <= matrix(rep(true_params, nrow(results)), 
                                                 nrow = nrow(results), byrow = TRUE)) & 
                           (matrix(rep(true_params, nrow(results)), 
                                   nrow = nrow(results), byrow = TRUE) <= ci_upper_SEM), na.rm = TRUE)

# ============================================================================
# Comprehensive Comparison Table
# ============================================================================

comparison_table <- data.frame(
  Parameter = rep(names(true_params), 1),
  True_Value = true_params,
  Avg_Estimate = avg_estimates,
  Bias = bias,
  Empirical_SE = empirical_se,
  Bootstrap_SE = avg_bootstrap_se,
  Louis_SE = avg_louis_se,
  SEM_SE = avg_SEM_se,
  Coverage_Bootstrap = coverage_boot,
  Coverage_Louis = coverage_louis,
  Coverage_SEM = coverage_SEM
)

cat("Comprehensive Comparison Table:\n")
comparison_table_print <- comparison_table
comparison_table_print[, -1] <- round(comparison_table_print[, -1], 4)
print(comparison_table_print, row.names = FALSE)
cat("\n")

# Save comparison table
saveRDS(comparison_table, file = "~/Desktop/STAT600/Figure/hw3_comparison_table.rds")

# ============================================================================
# Visualization 1: Standard Error Estimates by Method
# ============================================================================

# Prepare data for plotting
se_data <- data.frame(
  Parameter = rep(c("p", "lambda", "mu"), 3),
  Method = rep(c("Bootstrap", "Louis", "SEM"), each = 3),
  SE = c(avg_bootstrap_se, avg_louis_se, avg_SEM_se)
)

p_se_comparison <- ggplot(se_data, aes(x = Parameter, y = SE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "Standard Error Estimates by Method",
    subtitle = "Comparison of Bootstrap, Louis, and SEM Methods",
    x = "Parameter",
    y = "Average Standard Error",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("Bootstrap" = "#3498db", "Louis" = "#e74c3c", "SEM" = "#2ecc71")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(p_se_comparison)

# Save figure
ggsave("Figure/hw3_se_comparison.png", plot = p_se_comparison, width = 10, height = 6, dpi = 300)

# ============================================================================
# Visualization 2: Comprehensive Comparison
# ============================================================================

# Create long format data for comprehensive comparison
comp_data_estimates <- data.frame(
  Parameter = names(true_params),
  Metric = "Avg_Estimate",
  Value = avg_estimates
)

comp_data_bias <- data.frame(
  Parameter = names(true_params),
  Metric = "Bias",
  Value = bias
)

comp_data_se <- rbind(
  data.frame(Parameter = names(true_params), Metric = "Bootstrap_SE", Value = avg_bootstrap_se),
  data.frame(Parameter = names(true_params), Metric = "Louis_SE", Value = avg_louis_se),
  data.frame(Parameter = names(true_params), Metric = "SEM_SE", Value = avg_SEM_se),
  data.frame(Parameter = names(true_params), Metric = "Empirical_SE", Value = empirical_se)
)

comp_data_coverage <- rbind(
  data.frame(Parameter = names(true_params), Metric = "Coverage_Bootstrap", Value = coverage_boot),
  data.frame(Parameter = names(true_params), Metric = "Coverage_Louis", Value = coverage_louis),
  data.frame(Parameter = names(true_params), Metric = "Coverage_SEM", Value = coverage_SEM)
)

# Plot 1: Bias
p_bias <- ggplot(comp_data_bias, aes(x = Parameter, y = Value, fill = Parameter)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Bias", y = "Bias") +
  scale_fill_manual(values = c("p" = "#3498db", "lambda" = "#e74c3c", "mu" = "#2ecc71")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot 2: Standard Errors
p_se_all <- ggplot(comp_data_se, aes(x = Parameter, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Standard Errors", y = "Standard Error", fill = "Method") +
  scale_fill_manual(values = c("Bootstrap_SE" = "#3498db", "Louis_SE" = "#e74c3c", 
                               "SEM_SE" = "#2ecc71", "Empirical_SE" = "#95a5a6")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

# Plot 3: Coverage
p_coverage <- ggplot(comp_data_coverage, aes(x = Parameter, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Coverage Probability (95% CI)", y = "Coverage", fill = "Method") +
  scale_fill_manual(values = c("Coverage_Bootstrap" = "#3498db", "Coverage_Louis" = "#e74c3c", 
                               "Coverage_SEM" = "#2ecc71")) +
  ylim(0, 1) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

# Combine all plots
library(patchwork)
p_comprehensive <- (p_bias | p_se_all) / p_coverage +
  plot_annotation(
    title = "Comprehensive Comparison: Estimates, Bias, Standard Errors, and Coverage",
    subtitle = paste("True values: p =", p_true, ", λ =", lambda_true, ", μ =", mu_true),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5, size = 12))
  )

print(p_comprehensive)

# Save comprehensive figure
ggsave("Figure/hw3_comprehensive_comparison.png", plot = p_comprehensive, 
       width = 14, height = 10, dpi = 300)

