# Load necessary libraries
library(Optimization)
library(ggplot2)
library(dplyr) # Summary data


# Dataset
data1 <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 
           1.09, 1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 6.83)

l_prime <- function(theta, x) {
  sum(2 * (x - theta) / (1 + (x - theta)^2))
}

# Use functions defined in the package
bisection_result <- bisection(0, 2, data1)
newton_result <- newton(median(data1), data1)
fisher_result <- fisher_scoring(median(data1), data1)
secant_result <- secant(0, 2, data1)

c(bisection_result$theta,
  newton_result$theta,
  fisher_result$theta,
  secant_result$theta)

# Create a table
results_table <- data.frame(
  Method = c("Bisection", "Newton-Raphson", "Fisher Scoring", "Secant"),
  Theta = c(
    bisection_result$theta,
    newton_result$theta,
    fisher_result$theta,
    secant_result$theta
  ),
  Iterations =c(
    bisection_result$iteration,
    newton_result$iteration,
    fisher_result$iteration,
    secant_result$iteration
    
  )
)

# Print table
print(results_table)

# Ensure if a folder exist to save figures, if not, build one
if(!dir.exists("Figure")){
  dir.create("Figure")
}
saveRDS(results_table, file = "Figure/result_table.rds")
saveRDS(bisection_result, file = "Figure/bisection_result.rds")
saveRDS(newton_result, file = "Figure/newton_result.rds")
saveRDS(fisher_result, file = "Figure/fisher_result.rds")
saveRDS(secant_result, file = "Figure/secant_result.rds")

# Draw the figure of l'(theta)
theta_seq <- seq(-2, 4, by = 0.01)
l_prime_val <- sapply(theta_seq, l_prime, x = data1)

plot(theta_seq, l_prime_val, type = "l", 
     xlab = expression(theta), ylab = expression(l * "'"(theta)),
     main = "Derivative of Log Likelihood")
abline(h = 0, col = "red", lty = 2)
abline(v = c(0, 2), col = "blue", lty = 2)
points(c(0, 2), c(l_prime(0, data1), l_prime(2, data1)), 
       col = "blue", pch = 19)


# Compute the standard error of I(\theta)
l_double_prime <- function(theta, x) {
  sum(2 *((x - theta)^2 - 1) / (1 + (x - theta)^2)^2)
}
theta_hat <- fisher_result$theta
se <- 1/sqrt(-l_double_prime(theta_hat, data))
print(se)



data2 <- c(-8.34,  -1.73,  -0.40,  -0.24,   0.60,   0.94,   1.05,   1.06,   1.45,  1.50, 
           1.54,   1.72,   1.74,   1.88,   2.04,   2.16,   2.39,   3.01,   3.01,  3.08,
           4.66,   4.99,   6.01,   7.06,  25.45)

data_full <- c(data1, data2)
data_full_vec <- as.numeric(data_full)

# Run optimization on full data
 bisection_full <- bisection(0, 2, data_full_vec)
 newton_full <- newton(median(data_full), data_full_vec)
 fisher_full <- fisher_scoring(median(data_full), data_full_vec)
 secant_full <- secant(0, 2, data_full_vec)
 
 # Compute standard error
 theta_hat_full <- fisher_full$theta
 se_full <- 1/sqrt(-l_double_prime(theta_hat, data_full))
 print(se_full)
 
 # Create a table
 full_data_table <- data.frame(
   Method = c("Bisection", "Newton-Raphson", "Fisher Scoring", "Secant"),
   Theta = c(
     bisection_result$theta,
     newton_result$theta,
     fisher_result$theta,
     secant_result$theta
   ),
   Iterations =c(
     bisection_full$iteration,
     newton_full$iteration,
     fisher_full$iteration,
     secant_full$iteration
     
   )
 )
 
 # Print table
 print(full_data_table)
 
 saveRDS(full_data_table, file = "Figure/full_data_table.rds")
 saveRDS(bisection_result, file = "Figure/bisection_full.rds")
 saveRDS(newton_result, file = "Figure/newton_full.rds")
 saveRDS(fisher_result, file = "Figure/fisher_full.rds")
 saveRDS(secant_result, file = "Figure/secant_full.rds")
 saveRDS(se_full, file = "Figure/se_full.rds")
 print(se_full)
 
 