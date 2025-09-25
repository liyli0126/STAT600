# Load necessary libraries
library(Optimization)
library(ggplot2)
library(dplyr) # Summary data
library(knitr)
library(kableExtra)
library(webshot2)

data <- data.frame(
  cups = c(0, 2, 4, 5, 0, 2, 4, 5),
  gender = c(1, 1, 1, 1, 0, 0, 0, 0), # male=1, female=0
  r = c(9, 94, 53, 60, 11, 59, 53, 28),
  n = c(41, 213, 127, 142, 67, 211, 133, 76)
)

# Calculate observed probabilities
data$p_obs <- data$r / data$n

# Design matrix with intercept
X <- model.matrix(~ cups + gender, data = data)
y <- data$p_obs

# Run Newton's method
output <- multi_newton(X, y)

# Output results
coef <- c(output$coefficients)
se <- c(output$standard_errors)
param_names <- colnames(X)

print(se)
# Output as a table
results_df <- data.frame(
  Parameter = param_names,
  Estimate = round(coef, 3),
  Std_error = round(se, 3),
  Z_value = round(coef/se, 3),
  P_value = round(2 * pnorm(-abs(coef/se)), 4),
  stringsAsFactors = FALSE
)

print(results_df, row.names = FALSE)

# Create directory for figures if it doesn't exist
data_dir <- "~/Desktop/STAT600/STAT600_Code_test/Optimization/Figure" 
if(!dir.exists("Figure")){
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
}

# Save results
saveRDS(results_df, file = paste0(data_dir, "/multi_newton.rds"))

# Define predicted probability function
predict_prob <- function(X, beta){
  eta <- X %*% beta
  1 / (1 + exp(-eta))
}

# Create data for plotting
coffee_seq <- seq(0, 5, length.out = 100)

# Predictions for males (gender 1)
X_male <- cbind(1, coffee_seq, 1)
prob_male <- predict_prob(X_male, output$coefficients)
male_data <- data.frame(cups = coffee_seq, probability = prob_male, gender = "Male")

# Predictions for females (gender 0)
X_female <- cbind(1, coffee_seq, 0)
prob_female <- predict_prob(X_female, output$coefficients)
female_data <- data.frame(cups = coffee_seq, probability = prob_female, gender = "Female")

# Combine prediction data
pred_data <- rbind(male_data, female_data)

# Prepare observed data for plotting
obs_data <- data %>%
  mutate(
    probability = r / n,
    gender = ifelse(gender == 1, "Male", "Female"),
    se = sqrt(probability * (1 - probability) / n) # standard error
  )

# Determine predicted range
pred_range <- range(c(pred_data$probability, 
                      obs_data$probability - 1.96 * obs_data$se,
                      obs_data$probability + 1.96 * obs_data$se),
                    na.rm = TRUE)

# Set sutiable y range
y_limit <- c(max(0, pred_range[1] - 0.05), min(1, pred_range[2] + 0.05))

cat(sprintf("Using y-axis limits: %.3f to %.3f\n", y_limit[1], y_limit[2]))


# Create the plot
p <- ggplot() +
  # Prediction lines
  geom_line(data = pred_data, aes(x = cups, y = probability, color = gender, linetype = gender),
            linewidth = 1.2) + 
  # Observed data points
  geom_point(data = obs_data, aes(x = cups, y = probability, color = gender, shape = gender), 
             size = 3, alpha = 0.8) + 
  # Error bars for observed proportions (binormal standard errors)
  geom_errorbar(data = obs_data,
                aes(x = cups, y = probability,
                    ymin = pmax(probability - 1.96 * se, 0),   # Ensure greater than 0
                    ymax = pmin(probability + 1.96 * se, 1),   # Ensure less than 1
                    color = gender),
                width = 0.1) +
  # Labels and theme
  labs(title = "Coffee Consumption vs Probabiltiy of Pancreatic Cancer",
       x = "Coffee Consumption (cups per day)",
       y = "Probability of Cancer",
       color = "Gender",
       linetype = "Gender",
       shape = "Gender") +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold") 
  ) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed")) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 15)) +
  scale_y_continuous(limit = y_limit, breaks = seq(0, 1, by = 0.1))

# Save the plot
ggsave("~/Desktop/STAT600/STAT600_Code_test/Optimization/Figure/coffee_cancer_plot.png", p, width = 8, height = 6, dpi = 300)

# Print the plot
print(p)


# Wald tests for significance
z_scores <- output$coefficients / output$standard_errors
p_values <- 2 * pnorm(-abs(z_scores))

significance <- ifelse(p_values < 0.05, "Significant", "Not Significant")

significance_table <- data.frame(
  Parameter = c("Intercept", "Coffee", "Gender"),
  Estimate = output$coefficients,
  Std_Error = output$standard_errors,
  Z_Score = z_scores,
  P_Value = p_values,
  Significance = significance
)

cat("Wald test results \n")
print(significance_table, row.names = FALSE)

# Save results
saveRDS(significance_table, file = paste0(data_dir, "/significance_table.rds"))

cat("Likelihood Ratio Test:\n")
cat("Test Statistic:", lr_statistic, "\n")
cat("P-value:", lr_p_value, "\n")


# Likelihood ratio test (comparing with null model)
X_null <- model.matrix(~ 1, data = data) # Null model (intercept only)
X_no_coffee <- model.matrix(~ gender, data = data) # No coffee
X_no_gender <- model.matrix(~ cups, data = data) # No gender

# Simulate model
output_null <- multi_newton(X_null, y)
output_no_coffee <- multi_newton(X_no_coffee, y)
output_no_gender <- multi_newton(X_no_gender, y)

# Likelihood ratio test
perform_lrt_test <- function(full_model, reduced_model, test_name, df){
  lrt_stat <- -2 * (reduced_model$log_likelihood - full_model$log_likelihood)
  p_value <- pchisq(lrt_stat, df = df, lower.tail = FALSE)
  significance <- ifelse(p_value < 0.05, "Significant", "Not Significant")
  
  return(data.frame(
    Test = test_name,
    LRT_Statistic = round(lrt_stat, 4),
    DF = df,
    P_value = round(p_value, 4),
    Significance = significance
  ))
}

# Output the result
lrt_results <- rbind(
  perform_lrt_test(output, output_null, "Full vs Null Model", 2),
  perform_lrt_test(output, output_no_coffee, "Coffee Effect", 1),
  perform_lrt_test(output, output_no_gender, "Gender Effect", 1)
)

cat("Likelihood Ratio Test Results \n")
print(lrt_results, row.names = FALSE)

# Save results
saveRDS(lrt_results, file = paste0(data_dir, "/lrt_results.rds"))


# Comparison models
model_comparison <- data.frame(
  Model = c("Full Model", "Null model", "No Coffee", "No Gender"),
  Parameters = c(3, 1, 2, 2),
  Log_likelihood = c(
    output$log_likelihood,
    output_null$log_likelihood,
    output_no_coffee$log_likelihood,
    output_no_gender$log_likelihood
  ),
  AIC = c(
    -2 * output$log_likelihood + 2 * 3,
    -2 * output_null$log_likelihood + 2 * 1,
    -2 * output_no_coffee$log_likelihood + 2 * 2,
    -2 * output_no_gender$log_likelihood + 2 * 2
  )
)

cat("Model Comparison \ n")
print(model_comparison, row.names = FALSE)

# Save results
saveRDS(model_comparison, file = paste0(data_dir, "/model_comparison.rds"))


