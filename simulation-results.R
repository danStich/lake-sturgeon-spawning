# Libraries ----
library(tidyverse)


# . Data files used for simulation ----
# Load abundance estimates from empirical study
load(file = "results/grid_means.rda")

# Load detection estimates from empirical study
load(file = "results/p_summary.rda")

# Load capture histories from empirical study
# to get days sampled each year
load(file = "results/caps.rda")


# . Simulation output ----
# Find files with `simulation_output` in the results directory
files <- dir("results")[grep('simulation_output', dir("results"))]

# Make an empty list to hold the contents of those files
result <- vector(mode='list', length=length(files))

# For each file with `simulation_output` in name
for(i in 1:length(files)){
  # Load the file
  load(paste0("results/", files[i]))
  
  # Add the res dataframe from that simulation to the list
  result[[i]] <- res
}

# Rbind all of the simulation output dataframes into one big one
res <- do.call(rbind, result)


# Results ----
# . Bias in estimated abundance ----
# Calculations
median(res$bias)
quantile(res$bias, c(0.025, 0.975))

# Plot the bias in estimated abundance
n_bias_plot <- ggplot(res, aes(x = bias)) +
  geom_histogram(bins = 30) +
  xlab(expression(paste(hat(italic("N")), " - ", italic("N")))) +
  ylab("Frequency")

# Show the plot
n_bias_plot
# 
# # Save it to a file
# jpeg(file = "results/n_bias_plot.jpg",
#      height = 1800,
#      width = 2000,
#      res = 300)
# n_bias_plot
# dev.off()

# . Bias in estimated detection probability ----
# Calculations
# Calculate bias in estimated detection probability
res$p_bias = res$p_est - mean(p_summary$fit)

# Get mean and 95% CI for bias
mean(res$p_bias)
quantile(res$p_bias, c(0.025, 0.975))

# Plot the bias in detection probability
p_bias_plot <- ggplot(res, aes(x = p_bias)) +
  geom_histogram(bins = 30) +
  xlab(expression(paste(hat(italic("p")), " - ", italic("p")))) +
  ylab("Frequency")

# Show the plot
p_bias_plot

# # Save it to a file
# jpeg(file = "results/p_bias_plot.jpg",
#      height = 1800,
#      width = 2000,
#      res = 300)
# p_bias_plot
# dev.off()

# . Correlations between bias in p and N ----
# Calculations
# Calculate the Pearson correlation coefficient 
# between bias in N and bias in p
cor(res$bias, res$p_bias)

# Plot the correlation between bias in N and bias in p
cor_plot <- ggplot(res, aes(y = bias, x = p_bias)) +
  geom_point() +
  geom_smooth(method = "lm", color = 'black') +
  xlab(expression(paste(hat(italic("p")), " - ", italic("p")))) +
  ylab(expression(paste(hat(italic("N")), " - ", italic("N"))))

# Show the plot  
cor_plot

# Save it to a file
# jpeg(file = "results/cor_plot.jpg",
#      height = 1800,
#      width = 2000,
#      res = 300)
# cor_plot
# dev.off()



# Composite figure for manuscript -----
# Save it to a file
jpeg(file = "results/Figure5.jpg",
     height = 2400,
     width = 2000,
     res = 300)
gridExtra::grid.arrange(n_bias_plot, p_bias_plot, cor_plot)
dev.off()


