# Libraries----
library(tidyverse)
library(lubridate)
library(reshape)
library(R2jags)

# Dataset----
# Read in the data file
fish <- read.csv("data/sturgeon_grid_data.csv")
str(fish)

# Date manipulation----
# . QA/QC ----

# Remove rows in which there was an invalid count in Excel
fish <- fish %>% 
  filter(count != "#REF!")

# Remove spaces in the date column
fish$date <- gsub(pattern = " ", replacement = "",
                  x = fish$date)
# Remove spaces from bed name 
fish$bed <- gsub(pattern = " ", replacement = "",
                 x = fish$bed)
fish$bed <- gsub(pattern = "Upper", replacement = "Upstream",
                 x = fish$bed)

# Get average Easting and Northing per grid cell to 
# calculate the center point of the cell
fish$Easting <- rowMeans(fish[,c('left', 'right')])
fish$Northing <- rowMeans(fish[,c('top', 'bottom')])

# . Summarize by grid ----
# Remove an erroneous data point, then add counts together in 
# each grid cell
fish <- fish %>% 
  filter(!(bed == "Downstream" & Northing < 4983600)) %>% 
  group_by(date, bed, id, Easting, Northing) %>% 
  summarize(count = sum(count)) %>% 
  mutate(rep = row_number(), site = paste0(bed, "-", id))


# . Pivot to long form to add zeroes and NAs (for unsampled days) ----
# Convert date column to a Date object in a way that R can use it
fish$date <- as.Date(fish$date, format = "%m/%d/%Y")

# Extract day of year from the newly formatted date column
fish$day <- yday(fish$date)

# Extract year from the newly formatted date column
fish$year <- year(fish$date)

# Make a 3-dimensional array (matrix of matrices) to get the 
# grid cells that had zero fish counted that day
zeroes <- reshape::cast(fish, 
                        formula = bed ~ site ~ day ~ year,
                        value = "count",
                        fun.aggregate = "sum",
                        fill = 0,
                        drop = FALSE,
                        add.missing = TRUE)

# For each spawning bed in each year on each day, if the sum of all
# the counts for that day were zero, then it was unsampled, so 
# assign that grid cell a value of NA
for(i in 1:2){ # Spawning bed
  for(t in 1:12){ # Year
    for(d in 1:length(unique(fish$day))){ # Day
      if(sum(zeroes[i, , d, t]) == 0){
        zeroes[i, , d, t] <- NA
      }      
    }
  }
}


# Collapse the data back into long-form to make
# capture histories with or without site, date, bed, year
sturgeon <- reshape2::melt(zeroes)
names(sturgeon)[5] <- "count"

# Get rid of erroneous combos for site and bed that resulted from how
# we casted the zeroes array
sturgeon$check <- 0
sturgeon$check[sturgeon$bed == "Downstream" & grepl("Downstream", sturgeon$site)] <- 1
sturgeon$check[sturgeon$bed == "Upstream" & grepl("Upstream", sturgeon$site)] <- 1

sturgeon <- sturgeon %>% filter(check == 1)

# Get Eastings and Northings from original data for each grid cell
sites <- data.frame(unique(cbind(fish$site, fish$Easting, fish$Northing)))
names(sites) <- c("site", "Easting", "Northing")
sites[,2:3] <- apply(sites[, 2:3], 2, as.numeric)

# Merge the Eastings and Northings with the sturgeon data that contains zeroes
sturgeon <- merge(sturgeon, sites)

# Set aside data (for subsetting if needed)
test_data <- sturgeon

# .. Get sum of counts ----
caps <- cast(test_data, formula = site ~ day ~ year,
             value = "count",
             fun.aggregate = "sum")

# Spatial GAM for abundance ----
# Code adapted from:
# https://masonfidino.com/generalized_additive_occupancy_model/
# Following methods used by:


# . Spatial GAM for counts ----
# .. Get unique sites with coords for GAM ----
my_sites <- unique(test_data$site)
site_coords <- unique(cbind(test_data$Easting, test_data$Northing))

# .. Write spatial GAM parameter structure to a file ----
tmp_jags_spatial <- mgcv::jagam(
  response ~ s(x, y, k = 15, bs = "ds", m = c(1, 0.5)),
  data = data.frame(
    response = rep(1, length(my_sites)),
    x = site_coords[, 1],
    y = site_coords[, 2]
  ),
  family = "poisson",
  file = "models/tmp.jags")

# . JAGS data ----
data_list <- list(
  y = caps,
  X = tmp_jags_spatial$jags.data$X,                 # Spatial GAM Coordinates
  S1 = tmp_jags_spatial$jags.data$S1,               # Spatial GAM parameters             
  zero = tmp_jags_spatial$jags.data$zero,           # Spatial GAM parameters
  nsite = length(unique(test_data$site)),           # Site ID
  ndays = length(unique(test_data$day)),            # Number of days in data
  nyears = length(unique(test_data$year)))          # Number of years in data


# . Initial values ----
inits <- function(){list(
  N = matrix(max(test_data$count, na.rm = TRUE), 
          length(unique(test_data$site)),
          length(unique(test_data$year))
          ))
}


# . Parameters to save ----
params <- c("b", "rho", "log_nlambda", 
            "nlambda", "N", "mu_p", "p",
            "fit", "fit.rep", "lambda_pop")


# Compile the model ----
multi_year_fit <- jags(
  model.file = "models/multi-year.jags",
  parameters.to.save = params,
  data = data_list,
  inits = inits,
  n.chains = 3,
  n.iter = 500,
  n.burnin = 250,
  n.thin = 5)

# Print summary
# print(multi_year_fit, digits = 3)

# Save results to .rda file
# save(multi_year_fit, file = "results/multi_year_fit.rda")
load("results/multi_year_fit.rda")

# Results ----
# Posteriors -----
posts <- multi_year_fit$BUGSoutput$sims.list

# Figure S1 ----
# Make a dataframe to plot fitted vs predicted residuals
resids <- data.frame(
  Fitted = posts$fit,
  Predicted = posts$fit.rep
)

# Plot predicted vs fitted with a 1:1 line
ggplot(resids, aes(x = Fitted, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1)

# Calculate Bayesian p-value
mean(posts$fit.rep > posts$fit) 


# Figure 2 ----
n_posts <- melt(posts$N)
names(n_posts) <- c("iteration", "site", "year", "N")
n_posts$site <- unique(as.character(test_data$site))[n_posts$site]
n_posts$year <- unique(test_data$year)[n_posts$year]
n_posts <- left_join(n_posts, sites, by = "site")
n_posts$bed <- gsub("-.*","", n_posts$site)

# . Plot count by Northing and Easting ----
# Both spawning beds
sums <- n_posts %>% 
  group_by(bed, site, year, Easting, Northing) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) 

sums %>%
  arrange(fit) %>%
  ggplot(aes(x = Easting, y = Northing, color = fit, fill = fit)) +
  geom_tile() +
  scale_color_gradient(low = "gray87", high = "red") +
  scale_fill_gradient(low = "gray87", high = "red") +
  # coord_sf(default_crs = sf::st_crs(26914)) +
  scale_x_continuous(breaks = c(2000, 2100)) +
  labs(color = "Count", fill = "Count") +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        panel.spacing.x=unit(0, "lines"),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3)) +
  facet_wrap(~year, nrow = 1)


# Figure 3 ----
# . Max density per site per year ----
maxes <- sums %>% 
  group_by(bed, year) %>% 
  summarize(
    fit = max(fit),
    lwr = max(lwr),
    upr = max(upr)
    )

max_plot <- ggplot(maxes, aes(x = year, y = fit, color = bed, fill = bed)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr, color = NULL), 
              alpha = 0.25) +
  ylab(expression(paste("max(", italic("N"["i,j"]), ")"))) +
  xlab("") +
  scale_color_manual(labels = c("Downstream", "Upstream", "Both"),
                     values = c("gray60", "black", "gray40")) +
  scale_fill_manual(labels = c("Downstream", "Upstream", "Both"),
                     values = c("gray60", "black", "gray40")) +
  scale_x_continuous(breaks = seq(2012, 2022, 2)) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.title.y = element_text(vjust = 3)
  )

# . Max abundance per bed per year across sites ----
overalls <- sums %>% 
  group_by(bed, year) %>% 
  summarize(
    fit = sum(fit),
    lwr = sum(lwr),
    upr = sum(upr)) 

overalls_plot <- ggplot(overalls, aes(x = year, y = fit, color = bed, fill = bed)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr, color = NULL), 
              alpha = 0.25) +
  scale_color_manual(labels = c("Downstream", "Upstream"),
                     values = c("gray60", "black")) +
  scale_fill_manual(labels = c("Downstream", "Upstream"),
                    values = c("gray60", "black")) +
  xlab("") +
  ylab(expression(paste(
    {Sigma^italic(I)}[italic("i"), "= 1"], "(",italic("N"["i, j"]),")"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(vjust = 3)
        )
  

# . Max abundance at whole study area ----
total <- sums %>% 
  group_by(year) %>% 
  summarize(
    fit = sum(fit),
    lwr = sum(lwr),
    upr = sum(upr)) 

total_plot <- ggplot(total, aes(x = year, y = fit)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr), 
              alpha = 0.25) +
  xlab("Year") +
  ylab(expression(paste(
    Sigma, {Sigma^italic(I)}[italic("i"), "= 1"], "(",italic("N"["i, j"]),")"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(vjust = 3)
  )

jpeg("results/Figure3.jpg",
     height = 2400,
     width = 1800,
     res = 300
     )
  gridExtra::grid.arrange(max_plot, overalls_plot, total_plot)
dev.off()

# Figure 4 ----
# . Population growth rate ----
head(total)

# Calculate population growth rate as:
# r = ln(lambda); lambda = (Nt/Nt-1)/t
r_mat <- total %>%  mutate(
  fit = log(fit/lag(fit)),
  lwr = log(lwr/lag(lwr)),
  upr = log(upr/lag(upr))
  )

mean(r_mat$fit, na.rm = TRUE)
mean(r_mat$lwr, na.rm = TRUE)
mean(r_mat$upr, na.rm = TRUE)

r_plot <- r_mat %>% filter(year >= 2012) %>% 
ggplot(aes(x = year, y = fit)) +
  geom_hline(yintercept = a_fit, linetype = 2) +
  geom_ribbon(aes(xmax = year, ymin = a_lwr, ymax = a_upr), alpha = 0.25) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmax = fit, ymin = lwr, ymax = upr), width = 0) +
  xlab("Year") +
  ylab(expression(paste("Population growth rate (", italic(r), ")"))) +
  scale_x_continuous(breaks = seq(2012, 2022, 2)) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3))

jpeg("results/Figure4.jpg",
     width = 2400,
     height = 1800,
     res = 300
     )
r_plot
dev.off()


# Figure S2 Detection probability ----
p_posts <- melt(posts$p)
names(p_posts) <- c("iteration", "day", "p")
p_posts$day <- sort(unique(test_data$day))[p_posts$day]

p_summary <- p_posts %>% 
  group_by(day) %>% 
  summarize(fit = mean(p),
            lwr = quantile(p, 0.025),
            upr = quantile(p, 0.975))

p_plot <- ggplot(p_summary, aes(x = day, y = fit)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmax = day, ymin = lwr, ymax = upr),
              alpha = 0.5, width = 0) +
  xlab("Day of year") +
  ylab(expression(paste("Individual detection probability", italic("p"["t"])))) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3))  

jpeg("results/FigureS2.jpg",
     width = 2400,
     height = 1800,
     res = 300
)
p_plot
dev.off()



# Summary statistics ----
# .. Average number per site by bed ----
avgs <- n_posts %>% 
  group_by(bed) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) 
avgs

# .. Differences in max abundance between sites within year ----
max_diffs <- maxes %>% 
  select(c(bed, year, fit)) %>% 
  spread(bed, fit) %>% 
  mutate(diff = Downstream - Upstream)
mean(max_diffs$diff)

max_diffs <- maxes %>% 
  select(c(bed, year, lwr)) %>% 
  spread(bed, lwr) %>% 
  mutate(diff = Downstream - Upstream)
mean(max_diffs$diff)

max_diffs <- maxes %>% 
  select(c(bed, year, upr)) %>% 
  spread(bed, upr) %>% 
  mutate(diff = Downstream - Upstream)
mean(max_diffs$diff)

# .. Differences in sum of abundance between beds ----
overall_diffs <- overalls %>% 
  select(c(bed, year, fit)) %>% 
  spread(bed, fit) %>% 
  mutate(diff = Downstream - Upstream)
mean(overall_diffs$diff)

overall_diffs <- overalls %>% 
  select(c(bed, year, lwr)) %>% 
  spread(bed, lwr) %>% 
  mutate(diff = Downstream - Upstream)
mean(overall_diffs$diff)

overall_diffs <- overalls %>% 
  select(c(bed, year, upr)) %>% 
  spread(bed, upr) %>% 
  mutate(diff = Downstream - Upstream)
mean(overall_diffs$diff)
