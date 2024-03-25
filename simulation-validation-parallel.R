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

# Set up a parallel back-end ----
# Initialize snowfall
sfInit(parallel=TRUE, cpus=10, type="SOCK")

# Wrapper function ----
wrapper <- function(x){
  
# . Design considerations ----
nreps <- 10
nyears <- 1
ndays <- 1

# . Simulated abundance ----
p_test <- 0.038

n_sim <- matrix(NA, nrow = nrow(grid_means), ncol = nreps)

for(i in 1:ncol(n_sim)){
  n_sim[, i] <- rbinom(nrow(n_sim), size = round(grid_means$fit), prob = p_test)
}

# . Capture histories ----
sturgeon_sim_data <- data.frame(grid_means, n_sim)
sturgeon_sim_data$day <- 1
sturgeon_sim_data$year <- 1

caps_prep <- sturgeon_sim_data %>% 
  pivot_longer(cols = 8:(ncol(sturgeon_sim_data)-2),
               names_to = "rep",
               values_to = "count"
               )

caps <- cast(caps_prep, formula = id ~ day ~ year ~ rep,
             value = "count",
             fun.aggregate = "mean",
             add.missing = TRUE,
             fill = NA)

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
  y = caps,
  X = tmp_jags_spatial$jags.data$X,           # Spatial GAM Coordinates
  S1 = tmp_jags_spatial$jags.data$S1,         # Spatial GAM parameters             
  zero = tmp_jags_spatial$jags.data$zero,     # Spatial GAM parameters
  nsite = length(unique(grid_means$id)),      # Site ID
  ndays = 1,                                  # Number of days in data
  nyears = 1,                                 # Number of years in data
  knots = knots,
  nreps = nreps)     

# . Initial values ----
inits <- function(){list(
  N = matrix(max(caps_prep$count, na.rm = TRUE), 
             length(unique(caps_prep$id)),
             length(unique(caps_prep$year))))
}


# . Parameters to save ----
params <- c("b", "rho", "log_nlambda", 
            "nlambda", "N", "mu_p", "p")


# . Compile the model ----
multi_year_fit <- jags(
  model.file = "models/multi-year.jags",
  parameters.to.save = params,
  data = data_list,
  inits = inits,
  n.chains = 3,
  n.iter = 500000,
  n.burnin = 250000,
  n.thin = 5)

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

# Smash it all back together
sim_preds <- data.frame(grid_means, n_est)

# Calculate bias
sim_preds$bias <- round(sim_preds$n_est) - round(sim_preds$fit)


# . Write results to list (ADD) -----
sim <- list(
  nreps = nreps,
  nyears = nyears,
  ndays = ndays,
  p = p_test,
  bias = mean(sim_preds$bias)
)

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

# Process results ----
out <- lapply(result, data.frame)
res <- do.call(rbind, out)

ggplot(res, aes(x = bias)) +
  geom_histogram()

median(res$bias)

save(out, file = "results/simulation_output_01.rda")
