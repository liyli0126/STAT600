# Implement MCMC sampling with half-normal prior for different sigma values
run_hn_model <- function(sigma_prior, tau = 0.3, iter = 10000, burnin = 1000) {
  # Initialize parameters
  lambda1 <- 3
  lambda2 <- 1
  theta <- 56  # Fixed theta to focus on lambda1/lambda2 (as per problem scope)
  
  # Storage for samples (only lambda1, lambda2)
  samples <- matrix(0, nrow = iter, ncol = 2)
  colnames(samples) <- c("lambda1", "lambda2")
  accept_lambda1 <- 0
  accept_lambda2 <- 0
  
  for (i in 1:iter) {
    # --- Update lambda1 using Metropolis-Hastings ---
    lambda1_star <- rlnorm(1, meanlog = log(lambda1), sdlog = tau)
    S1 <- sum(disasters[1:theta])
    
    # Log acceptance ratio for lambda1
    log_A <- (S1 + 1) * (log(lambda1_star) - log(lambda1)) - 
      theta * (lambda1_star - lambda1) - 
      (lambda1_star^2 - lambda1^2) / (2 * sigma_prior^2)
    
    if (log(runif(1)) < log_A) {
      lambda1 <- lambda1_star
      accept_lambda1 <- accept_lambda1 + 1
    }
    
    # --- Update lambda2 using Metropolis-Hastings ---
    lambda2_star <- rlnorm(1, meanlog = log(lambda2), sdlog = tau)
    S2 <- sum(disasters[(theta + 1):n_years])
    
    # Log acceptance ratio for lambda2
    log_A <- (S2 + 1) * (log(lambda2_star) - log(lambda2)) - 
      (n_years - theta) * (lambda2_star - lambda2) - 
      (lambda2_star^2 - lambda2^2) / (2 * sigma_prior^2)
    
    if (log(runif(1)) < log_A) {
      lambda2 <- lambda2_star
      accept_lambda2 <- accept_lambda2 + 1
    }
    
    # Store samples
    samples[i, ] <- c(lambda1, lambda2)
  }
  
  # Remove burn-in
  samples_post <- samples[-(1:burnin), ]
  
  return(list(
    samples = samples_post,
    accept_rate = c(accept_lambda1 / iter, accept_lambda2 / iter),
    sigma_prior = sigma_prior
  ))
}

# Run models for different sigma values
sigma_values <- c(1, 2, 3, 4, 5)
results <- list()

for (sigma in sigma_values) {
  results[[as.character(sigma)]] <- run_hn_model(sigma_prior = sigma)
}

# Generate comparison table
generate_comparison_table <- function(results) {
  comparison <- data.frame()
  
  for (sigma in names(results)) {
    samples <- results[[sigma]]$samples
    accept_rates <- results[[sigma]]$accept_rate
    
    post_summary <- data.frame(
      Sigma = as.numeric(sigma),
      Lambda1_mean = mean(samples[, "lambda1"]),
      Lambda1_sd = sd(samples[, "lambda1"]),
     # Lambda1_q2.5 = quantile(samples[, "lambda1"], 0.025),
      #Lambda1_q97.5 = quantile(samples[, "lambda1"], 0.975),
      Lambda2_mean = mean(samples[, "lambda2"]),
      Lambda2_sd = sd(samples[, "lambda2"]),
     # Lambda2_q2.5 = quantile(samples[, "lambda2"], 0.025),
      #Lambda2_q97.5 = quantile(samples[, "lambda2"], 0.975),
      Accept_rate_lambda1 = accept_rates[1],
      Accept_rate_lambda2 = accept_rates[2]
    )
    
    comparison <- rbind(comparison, post_summary)
  }
  
  return(comparison)
}

# Generate and print comparison table
comparison_table <- generate_comparison_table(results)
print(comparison_table)


latex_table <- xtable(
  comparison_table,
  caption = "Comparison table for \\(\\sigma = 1, \\dots, 5\\)",
  label = "tab:comparison",
  digits = 3  
)

print(latex_table, 
      type = "latex",
      file = "~/Desktop/STAT600/Table/hw5_comparison_table.tex",
      sanitize.text.function = identity)


