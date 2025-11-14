library(coda)
library(xtable)
library(ggplot2)
library(patchwork)
library(HDInterval)
# Load data
disasters_data <- read.table("~/Desktop/STAT600/HW/HW5/coal.dat", head = TRUE)
disasters <- disasters_data$disasters


n_years <- length(disasters)
iter <- 10000
burnin <- 1000

# Define function to run a single chain
run_single_chain <- function(seed, disasters, iter = 10000, burnin = 1000) {
  set.seed(seed)
  n_years <- length(disasters)
  
  # Initialize parameters INSIDE the function
  lambda1 <- 1
  lambda2 <- 1
  alpha <- 1
  theta <- 56
  
  # Storage samples INSIDE the function
  samples <- matrix(0, nrow = iter, ncol = 4)
  colnames(samples) <- c("lambda1", "lambda2", "alpha", "theta")
  
  # Gibbs sampling
  for (i in 1:iter) {
    # Sample lambda1: ~ Gamma(S1 + 3, theta + alpha)
    sum_x1 <- sum(disasters[1:theta])
    shape1 <- sum_x1 + 3
    rate1 <- theta + alpha
    lambda1 <- rgamma(1, shape = shape1, rate = rate1)
    
    # Sample lambda2
    sum_x2 <- sum(disasters[(theta + 1):n_years])
    shape2 <- sum_x2 + 3
    rate2 <- (n_years - theta) + alpha
    lambda2 <- rgamma(1, shape = shape2, rate = rate2)
    
    # Sample alpha
    shape_alpha <- 16
    rate_alpha <- lambda1 + lambda2 + 10
    alpha <- rgamma(1, shape = shape_alpha, rate = rate_alpha)  
    
    # Sample theta
    log_prob <- numeric(111)
    for (j in 1:111) {
      sum1 <- sum(disasters[1:j])
      sum2 <- sum(disasters[(j + 1):n_years])
      log_prob[j] <- sum1 * log(lambda1) - j * lambda1 + 
        sum2 * log(lambda2) - (n_years - j) * lambda2
    }
    prob <- exp(log_prob - max(log_prob))  # Enhanced stability
    prob <- prob / sum(prob)
    theta <- sample(1:111, size = 1, prob = prob)
    
    samples[i, ] <- c(lambda1, lambda2, alpha, theta)
  }
  
  # Remove burn-in
  samples_post <- samples[-(1:burnin), ]
  
  return(as.mcmc(samples_post))
}

# Run several chains independently to compute Gelman-Rubin
chain1 <- run_single_chain(seed = 26, disasters = disasters, iter = iter, burnin = burnin)
chain2 <- run_single_chain(seed = 123, disasters = disasters, iter = iter, burnin = burnin)
chain3 <- run_single_chain(seed = 234, disasters = disasters, iter = iter, burnin = burnin)

# Create mcmc list
mcmc_list <- mcmc.list(chain1, chain2, chain3)

# Execute convergence diagnostics
#gelman.diag(mcmc_list)
gr_diag <- gelman.diag(mcmc_list, multivariate = FALSE)
gr_table <- as.data.frame(gr_diag$psrf)
rownames(gr_table) <- c("$\\lambda_1$", "$\\lambda_2$", "$\\alpha$", "$\\theta$")

latex_table <- xtable(
  gr_table,
  caption = "Gelman-Rubin convergence diagnostics ($\\hat{R}$ statistics)",
  label = "tab:gelman_rubin",
  digits = 3  
)

print(latex_table, 
      type = "latex",
      file = "~/Desktop/STAT600/Table/hw5_gelman_rubin.tex",
      sanitize.text.function = identity)

# Combine all chains for posterior analysis
combined_chains <- do.call(rbind, mcmc_list)

#plot(mcmc_list)
png("~/Desktop/STAT600/Figure/hw5_trace_plots.png", 
    width = 2000, height = 1600, res = 300, type = "cairo")
par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(mcmc_list, trace = TRUE)
dev.off()


#autocorr.plot(mcmc_list)
png("~/Desktop/STAT600/Figure/hw5_autocorr_plots.png", 
    width = 2000, height = 1600, res = 300, type = "cairo")
par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
autocorr.plot(mcmc_list, lag.max = 20)  
dev.off()

effectiveSize(mcmc_list)

# Optional: Check summary of combined chains
summary(mcmc_list)

p1 <- ggplot(data.frame(x = combined_chains[,"lambda1"]), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins=30, fill="grey", color="black") +
  ggtitle("Lambda1 Posterior") + theme_minimal()

p2 <- ggplot(data.frame(x = combined_chains[,"lambda2"]), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins=30, fill="grey", color="black") +
  ggtitle("Lambda2 Posterior") + theme_minimal()

p3 <- ggplot(data.frame(x = combined_chains[,"alpha"]), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins=30, fill="grey", color="black") +
  ggtitle("Alpha Posterior") + theme_minimal()

p4 <- ggplot(data.frame(x = combined_chains[,"theta"]), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins=30, fill="grey", color="black") +
  ggtitle("Theta Posterior") + theme_minimal()

combined_plot <- (p1 | p2) / (p3 | p4)

ggsave("~/Desktop/STAT600/Figure/hw5_posterior_hist.png", plot = combined_plot, width = 12, height = 8)


summary_stats <- data.frame(
  Parameter = c("Lambda1", "Lambda2", "Alpha", "Theta"),
  Mean = apply(combined_chains, 2, mean),
  Median = apply(combined_chains, 2, median),
  SD = apply(combined_chains, 2, sd)
)

# Compute HPD intervals using coda
hpd_matrix <- HPDinterval(as.mcmc(combined_chains), prob = 0.95)
summary_stats$HPD_lower <- hpd_matrix[,1]
summary_stats$HPD_upper <- hpd_matrix[,2]

latex_tab <- xtable(
  summary_stats,
  caption = "Posterior summary statistics with 95\\% HPD intervals",
  label = "tab:posterior_summary",
  digits = 3
)

print(
  latex_tab,
  type = "latex",
  file = "~/Desktop/STAT600/Table/hw5_posterior_summary.tex",
  sanitize.text.function = identity,
  include.rownames = FALSE
)
hpd_matrix <- HPDinterval(as.mcmc(combined_chains), prob = 0.95)
summary_stats$HPD_lower <- hpd_matrix[,1]
summary_stats$HPD_upper <- hpd_matrix[,2]

print(summary_stats)

