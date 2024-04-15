# Libraries----
library(tidyverse)
library(lubridate)
library(reshape)
library(R2jags)
library(mgcv)
library(snowfall)

# Data read ----
# Load abundance estimates from empirical study
load(file = "results/grid_means.rda")

# Load detection estimates from empirical study
load(file = "results/p_summary.rda")

# Load capture histories from empirical study
# to get days sampled each year
load(file = "results/caps.rda")

# Set up a parallel back-end ----
# Initialize snowfall
sfInit(parallel = TRUE, cpus = 10, type="SOCK")


# Wrapper function ----
wrapper <- function(x){
  
# . Design considerations ----
p_true <- p_summary$fit #rep(0.5, length(p_summary$fit)) 
# p_true <- p_summary$fit
n_true <- ceiling(grid_means$fit)
  
nsites <- nrow(grid_means)
nreps <- 5
nyears <- 13
ndays <- length(p_true)

# . Simulated abundance ----

# Container to hold observed counts by sites, days, years, reps
y <- array(NA, dim = c(nsites,           # Sites
                       ndays,            # Days
                       nyears,           # Years
                       nreps))           # Replicates

# Simulate observed capture histories from N and p
for(i in 1:nsites){
  for(t in 1:ndays){
    for(j in 1:nyears){
      for(k in 1:nreps){
        y[i, t, j, k] <- rbinom(1, size = n_true[i], prob = p_true[t])
      }
    }
  }
}

y[is.na(caps)] <- NA

# Model Likelihood
# for(i in 1:nsite){
#   for(t in 1:ndays){
#     for(j in 1:nyears){
#       for(k in 1:nreps){
#         # Response
#         y[i, t, j, k] ~ dbin(p[t], N[i, j])
#       }
#     }
#   }
# }

# . Spatial GAM ----
# . Get unique sites with coords for GAM ----
my_sites <- unique(grid_means$id)
site_info <- unique(cbind(grid_means$Easting, grid_means$Northing, grid_means$Bed))
site_coords <- unique(cbind(grid_means$Easting, grid_means$Northing))
sites <- data.frame(my_sites, site_info)
names(sites) <- c('id', 'Easting', 'Northing', 'Bed')

# . Write spatial GAM parameter structure to a file ----
knots = 15
tmp_jags_spatial <- mgcv::jagam(
  response ~ s(x, y, k = knots, bs = "ds", m = c(1, 0.5)),
  data = data.frame(
    response = rep(1, length(my_sites)),
    x = site_coords[, 1],
    y = site_coords[, 2]
  ),
  family = "poisson",
  file = "models/tmp_sim.jags")

# . JAGS data ----
data_list <- list(
  y = y,
  X = tmp_jags_spatial$jags.data$X,           # Spatial GAM Coordinates
  S1 = tmp_jags_spatial$jags.data$S1,         # Spatial GAM parameters             
  zero = tmp_jags_spatial$jags.data$zero,     # Spatial GAM parameters
  nsite = nsites,                             # Site ID
  ndays = ndays,                              # Number of days in data
  nyears = nyears,                            # Number of years in data
  nreps = nreps,                              # Number of replicates
  knots = knots)                              # Knots in the GAM

# . Initial values ----
inits <- function(){list(
  N = matrix(max(y, na.rm = TRUE), 
             nsites,
             nyears))
}


# . Parameters to save ----
params <- c("N", "p")


# . Compile the model ----
multi_year_fit <- jags(
  model.file = "models/multi-year-sim-lightweight.jags",
  parameters.to.save = params,
  data = data_list,
  inits = inits,
  n.chains = 3,
  n.iter = 100000,
  n.burnin = 25000,
  n.thin = 5,
  refresh = 1)

# Print summary
print(multi_year_fit, digits = 3)

# . Results ----
# Extract posterior estimates
posts <- multi_year_fit$BUGSoutput$sims.list

# Get estimates of N
n_posts <- melt(posts$N)
glimpse(n_posts)
names(n_posts) <- c("iteration", "id", "year", "n_est")

n_est <- n_posts %>% 
  group_by(id) %>% 
  summarize(n_est = mean(n_est)) %>% 
  select(n_est) %>% 
  unlist()

bias_sum <- sum(n_est) - sum(n_true)

# . Write results to list (UPDATE) -----
sim <- list(
  nreps = nreps,
  nyears = nyears,
  ndays = ndays,
  # p = p_test,
  p_est = mean(posts$p),
  n_true = median(n_true),
  n_est = median(n_est),
  n_true_mu = mean(n_true),
  n_est_mu = mean(n_est),  
  bias = median(n_est) - median(n_true))

return(sim)

}

# Simulation settings ----
# Necessary libraries
sfLibrary(tidyverse)
sfLibrary(lubridate)
sfLibrary(reshape)
sfLibrary(R2jags)
sfLibrary(mgcv)
sfLibrary(snowfall)

sfExport("grid_means")
sfExport("p_summary")
sfExport("caps")

# Number of simulations
n_iterations <- 100

# Start timer
start_time <- Sys.time()

# Run simulation ----
# Run the simulation n_iterations times
result <- sfLapply(1:n_iterations, wrapper) 

# Stop snowfall
sfStop()

# Stop timer
total_time <- Sys.time()- start_time
total_time

# Process output ----
out <- lapply(result, data.frame)
res <- do.call(rbind, out)

save(res, file = "results/simulation_output_02.rda")

# Results ----
# . Bias in estimated abundance ----
# Calculations
median(res$bias)
quantile(res$bias, c(0.025, 0.975))

# Plot the bias in estimated abundance
n_bias_plot <- ggplot(res, aes(x = bias)) +
  geom_histogram(bins = 10) +
  xlab(expression(paste(hat(italic("N")), " - ", italic("N")))) +
  ylab("Frequency")

# Show the plot
n_bias_plot

# . Bias in estimated detection probability ----
# Calculations
# Calculate bias in estimated detection probability
res$p_bias = res$p_est - mean(p_summary$fit)

# Get mean and 95% CI for bias
mean(res$p_bias)
quantile(res$p_bias, c(0.025, 0.975))

# Plot the bias in detection probability
p_bias_plot <- ggplot(res, aes(x = p_bias)) +
  geom_histogram(bins = 10) +
  xlab(expression(paste(hat(italic("p")), " - ", italic("p")))) +
  ylab("Frequency")

# Show the plot
p_bias_plot

# . Correlations between bias in p and N ----
# Calculations
# Calculate the Pearson correlation coefficient 
# between bias in N and bias in p
cor(res$bias, res$p_bias)

# Plot the correlation between bias in N and bias in p
cor_plot <- ggplot(res, aes(x = bias, y = p_bias)) +
  geom_point() +
  geom_smooth(method = "lm", color = 'black') +
  ylab(expression(paste("Bias in ", hat(italic("p"))))) +
  xlab(expression(paste("Bias in ", hat(italic("N")))))

# Show the plot  
cor_plot
