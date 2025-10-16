# Load necessary libraries
library(MC)

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

set.seed(26)
n_reps <- 100
n_samples <- 1000

results <- comparison(n_reps, n_samples)

# Create a data framework
results_df <- data.frame(
  Method = rep(c("Importance Sampling", "Rejection Sampling", 
                 "SIR", "Riemann Strategy"), each = n_reps),
  Estimate = c(results$estimate[,1], results$estimate[,2],
               results$estimate[,3], results$estimate[,4])
)

# Summary performance
performance_stats <- data.frame(
  Method = c("Importance Sampling", "Rejection Sampling", "SIR", "Riemann Strategy"),
  Mean_Estimate = results$mean_estimates,
  Variance = results$variance_estimates,
  MSE = (results$mean_estimates - mean(results_df$Estimate))^2 + results$variance_estimates,
  Time_per_Rep = results$mean_times
)

print(performance_stats)
saveRDS(performance_stats, file = "~/Desktop/STAT600/Figure/hw4_performance_stats.rds") 


p1 <- ggplot(results_df, aes(x = Method, y = Estimate, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "A: Distribution of Estimates by Method",
       subtitle = paste("n =", n_samples, "samples,", n_reps, "replications"),
       y = "Estimate of σ²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p2 <- ggplot(performance_stats, aes(x = Method, y = Variance, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  labs(title = "B: Variance of Estimates by Method",
       y = "Variance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p3 <- ggplot(performance_stats, aes(x = Method, y = Time_per_Rep, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  labs(title = "C: Computation Time by Method",
       y = "Time per Replication (seconds)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Combine figures
combined_plot <- (p1 | p2 | p3) + 
  plot_annotation(title = "Comparison of Monte Carlo Methods for Estimating \u03C3\u00B2",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

# Save the combined plot
ggsave("~/Desktop/STAT600/Figure/hw4_combined_results.png", 
       plot = combined_plot, width = 16, height = 8, dpi = 300)



# Compare with different number of sampling 
sample_sizes <- c(100, 500, 1000, 5000)

# Initialize a storage list 
size_results <- list(results = list())
for( i in 1:length(sample_sizes)){
  size_results$results[[i]] <- comparison(50, sample_sizes[i])
}

# Tackle the results of different sampling 
variance_data <- data.frame()
for (i in 1:length(sample_sizes)) {
  n <- sample_sizes[i]
  result <- size_results$results[[i]]
  variances <- result$variance_estimates
  
  temp_df <- data.frame(
    SampleSize = n,
    Method = c("Importance Sampling", "Rejection Sampling", "SIR", "Riemann Strategy"),
    Variance = variances
  )
  variance_data <- rbind(variance_data, temp_df)
}


# Create a comparison figure
p4 <- ggplot(variance_data, aes(x = SampleSize, y = Variance, color = Method, shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Variance Reduction with Increasing Sample Size",
       x = "Sample Size (log scale)", 
       y = "Variance (log scale)") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Create method to compare tables
performance_table <- performance_stats %>%
  mutate(across(where(is.numeric), ~round(., 6)))

p_table <- ggplot() +
  annotation_custom(gridExtra::tableGrob(performance_table, rows = NULL, 
                                         theme = gridExtra::ttheme_minimal())) +
  labs(title = "Performance Comparison Table") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Combine the figure together
final_combined <- (combined_plot / (p4 | p_table)) + 
  plot_annotation(title = "Comprehensive Monte Carlo Methods Analysis",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 18)))
print(combined_plot)
print(p4)

# Save figures
ggsave("~/Desktop/STAT600/Figure/hw4_combined_plot_for_performance.png", 
       plot = final_combined, width = 16, height = 8, dpi = 300)



