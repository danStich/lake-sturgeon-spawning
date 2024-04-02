# Libraries----
library(tidyverse)
library(lubridate)
library(reshape)
library(R2jags)
library(mgcv)
library(cowplot)

# Datasets----
# . Fish counts ----
fish <- read.csv("data/sturgeon_grid_data.csv")
str(fish)

fish <- fish %>% 
  filter(count != "#REF!" & !is.na(id)) %>% 
  mutate(count = as.numeric(count))

fish <- fish[, c('date', 'id', 'count', 'file_name')]

# Sum within grid cells for a single transect (file_name)
fish <- fish %>%
  group_by(date, id, file_name) %>%
  summarize(count = sum(count))

# .Transect data ----
transects <- read.csv("data/all_transects_2011_2023.csv")
names(transects) <- gsub("\\.", "_", names(transects))

transects <- transects[, c('id', 'date', 'bed', 
                           'left', 'top', 'right', 'bottom')]


# Get average Easting and Northing per grid cell
transects$Easting <- rowMeans(transects[,c('left', 'right')])
transects$Northing <- rowMeans(transects[,c('top', 'bottom')])

# Drop any cells missing coords
transects <- transects[!is.na(transects$Easting), ]
transects <- transects[!duplicated(transects[1:ncol(transects)]), ]


# . Combined data ----
# Merge fish data and transect data
fish_merge <- merge(x = fish, y = transects, 
                    by = c('id', 'date'),
                    all.y = TRUE)

nrow(fish_merge)

fish_merge$count[is.na(fish_merge$count)] <- 0


# . Date management ----
fish_merge$date <- as.Date(fish_merge$date, format = "%m/%d/%Y")
fish_merge$day <- yday(fish_merge$date)
fish_merge$year <- year(fish_merge$date)


# Date manipulation----
# . Summarize by grid ----
sturgeon_check <- fish_merge %>% 
  filter(!(bed == "downstream" & Northing < 4983950)) %>% 
  group_by(date, year, day, bed, id, Easting, Northing) %>% 
  mutate(rep = row_number(), site = as.character(id)) %>% 
  filter(!is.na(Easting) & !is.na(Northing))


# . Plot by longitude and latitude ----
# Both spawning beds
sturgeon_check %>%
  # filter(bed == 'upstream' & year == 2018) %>%
  arrange(count) %>%
  ggplot(aes(x = Easting, y = Northing, color = count)) +
  geom_point() +
  scale_color_gradient(low = "black", high = "red") +
  coord_sf(default_crs = sf::st_crs(26918)) +
  labs(color = "Count", fill = "Count") +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.spacing.x=unit(2.5, "lines"),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3))

# . Subsetting ----
sturgeon <- sturgeon_check
sturgeon <- sturgeon %>%
  filter(day >= 120)

# .. Capture history ----
caps <- cast(sturgeon, formula = id ~ day ~ year ~ rep,
             value = "count",
             fun.aggregate = "mean",
             add.missing = TRUE,
             fill = NA)

# Spatial GAM for abundance ----
# Code adapted from:
# https://masonfidino.com/generalized_additive_occupancy_model/
# Following methods used by:


# . Spatial GAM for counts ----
# .. Get unique sites with coords for GAM ----
my_sites <- unique(sturgeon$id)
site_info <- unique(cbind(sturgeon$Easting, sturgeon$Northing, sturgeon$bed))
site_coords <- unique(cbind(sturgeon$Easting, sturgeon$Northing))
sites <- data.frame(my_sites, site_info)
names(sites) <- c('id', 'Easting', 'Northing', 'Bed')

# .. Write spatial GAM parameter structure to a file ----
knots = 15
tmp_jags_spatial <- mgcv::jagam(
  response ~ s(x, y, k = knots, bs = "ds", m = c(1, 0.5)),
  data = data.frame(
    response = rep(1, length(my_sites)),
    x = site_coords[, 1],
    y = site_coords[, 2]
  ),
  family = "poisson",
  file = "models/tmp_ds.jags")

# . JAGS data ----
data_list <- list(
  y = caps,
  X = tmp_jags_spatial$jags.data$X,           # Spatial GAM Coordinates
  S1 = tmp_jags_spatial$jags.data$S1,         # Spatial GAM parameters             
  zero = tmp_jags_spatial$jags.data$zero,     # Spatial GAM parameters
  nsite = length(unique(sturgeon$site)),      # Site ID
  ndays = length(unique(sturgeon$day)),       # Number of days in data
  nyears = length(unique(sturgeon$year)),     # Number of years in data
  knots = knots,
  nreps = length(unique(sturgeon$rep)))     

# . Initial values ----
inits <- function(){list(
  N = matrix(max(sturgeon$count, na.rm = TRUE), 
          length(unique(sturgeon$id)),
          length(unique(sturgeon$year))))
}


# . Parameters to save ----
params <- c("b", "rho", "log_nlambda", 
            "nlambda", "N", "mu_p", "p",
            "fit", "fit.rep")


# . Compile the model ----
multi_year_fit <- jags(
  model.file = "models/multi-year.jags",
  parameters.to.save = params,
  data = data_list,
  inits = inits,
  n.chains = 3,
  n.iter = 100000,
  n.burnin = 50000,
  n.thin = 5)

# Print summary
print(multi_year_fit, digits = 3)

# Results ----
# . Save results to .rda file ----
save(multi_year_fit, file = "results/multi_year_fit.rda")
# load("results/multi_year_fit.rda")

# Posteriors -----
posts <- multi_year_fit$BUGSoutput$sims.list

# Figure S1 (Posterior predictive check) ----
# Make a dataframe to plot fitted vs predicted residuals
resids <- data.frame(
  Fitted = posts$fit,
  Predicted = posts$fit.rep
)

# Calculate Bayesian p-value
mean(posts$fit.rep > posts$fit) 

# Figure
jpeg("results/FigureS1.jpg",
     height = 1800,
     width = 1800,
     res = 300)
  # Plot predicted vs fitted with a 1:1 line
  ggplot(resids, aes(x = Fitted, y = Predicted)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1)
dev.off()

# Figure 2 (Spatial plot) ----
n_posts <- melt(posts$N)
names(n_posts) <- c("iteration", "id_num", "year", "N")
n_posts$id <- as.numeric(unique(as.character(sturgeon$id))[n_posts$id])
n_posts$year <- unique(sturgeon$year)[n_posts$year]
n_posts <- left_join(n_posts, sites, by = "id")

# . Plot count by Northing and Easting ----
# Both spawning beds
sums <- n_posts %>% 
  group_by(id, year, Easting, Northing) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) %>% 
  mutate(Easting = as.numeric(Easting),
         Northing = as.numeric(Northing))

Figure2 <- sums %>%
  arrange(fit) %>%
  ggplot(aes(x = Easting, y = Northing, color = fit, fill = fit)) +
  geom_tile() +
  scale_color_gradient(low = "gray87", high = "red") +
  scale_fill_gradient(low = "gray87", high = "red") +
  coord_sf(default_crs = sf::st_crs(26914)) +
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

jpeg("results/Figure2.jpg",
     height = 1800,
     width = 2400,
     res = 300
)
Figure2
dev.off()



# Figure 4 (Abundance plots) ----
# . Mean density per site per year ----
means <- n_posts %>% 
  group_by(id, year, Easting, Northing, Bed) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) %>% 
  mutate(Easting = as.numeric(Easting),
         Northing = as.numeric(Northing)) %>% 
  group_by(year, Bed) %>% 
  summarize(fit = mean(fit),
            lwr = mean(lwr),
            upr = mean(upr))

means_plot <- ggplot(means, aes(x = year, y = fit, color = Bed, fill = Bed)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr, color = NULL), 
              alpha = 0.25) +
  ylab(expression(paste("Mean density per 900 m"^2))) +
  xlab("") +
  labs(color = "Spawning Bed", fill = "Spawning Bed") +
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


# . Sum of abundance per bed per year across sites ----
overalls <- n_posts %>% 
  group_by(id, year, Easting, Northing, Bed) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) %>% 
  mutate(Easting = as.numeric(Easting),
         Northing = as.numeric(Northing)) %>% 
  group_by(year, Bed) %>% 
  summarize(fit = sum(fit),
            lwr = sum(lwr),
            upr = sum(upr))

overalls_plot <- ggplot(overalls, aes(x = year, y = fit, color = Bed, fill = Bed)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr, color = NULL), 
              alpha = 0.25) +
  scale_color_manual(labels = c("Downstream", "Upstream"),
                     values = c("gray60", "black")) +
  scale_fill_manual(labels = c("Downstream", "Upstream"),
                    values = c("gray60", "black")) +
  scale_x_continuous(breaks = seq(2012, 2022, 2)) +
  xlab("") +
  # labs(color = "Spawning Bed", fill = "Spawning Bed") +
  ylab("Total abundance per site") +
  theme_bw() +
  theme(legend.position = "NULL",
        legend.direction = "horizontal",
        axis.title.y = element_text(vjust = 3)
        )
  
# jpeg("results/ppt_beds.jpg",
#      height = 1800,
#      width = 2400,
#      res = 300
# )
# overalls_plot
# dev.off()

# . Sum of abundance at whole study area ----
total <- n_posts %>% 
  group_by(id, year, Easting, Northing, Bed) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) %>% 
  mutate(Easting = as.numeric(Easting),
         Northing = as.numeric(Northing)) %>% 
  group_by(year) %>% 
  summarize(fit = sum(fit),
            lwr = sum(lwr),
            upr = sum(upr))

total_plot <- ggplot(total, aes(x = year, y = fit)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr), 
              alpha = 0.25) +
  scale_x_continuous(breaks = seq(2012, 2022, 2)) +
  xlab("Year") +
  labs(color = "Spawning Bed", fill = "Spawning Bed") +
  ylab("Total abundance across sites") +
  theme_bw() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        axis.title.y = element_text(vjust = 3)
  )

# jpeg("results/ppt_total.jpg",
#      height = 1800,
#      width = 2400,
#      res = 300
# )
# total_plot
# dev.off()


# .. Figure ----
jpeg("results/Figure4.jpg",
     height = 2400,
     width = 1800,
     res = 300
     )
  ggdraw() +
    draw_plot(means_plot, x = 0, y = .66, width = 1, height = .35) +
    draw_plot(overalls_plot, x = 0, y = .33, width = 1, height = .3) +
    draw_plot(total_plot, x = 0, y = 0, width = 1, height = 0.3)
  
dev.off()

# Figure 5 (Population growth) ----
# . Population growth rate ----
head(total)

# Calculate population growth rate as:
# r = ln(lambda); lambda = (Nt/Nt-1)/t
r_mat <- overalls %>%  
  group_by(Bed) %>% 
  mutate(
    fit = log(fit/lag(fit)),
    lwr = log(lwr/lag(lwr)),
    upr = log(upr/lag(upr))
  )

r_mat %>% 
  group_by(Bed) %>% 
  summarize(fit = mean(fit, na.rm = TRUE),
            lwr = mean(lwr, na.rm = TRUE),
            upr = mean(upr, na.rm = TRUE))

mean(r_mat$fit, na.rm = TRUE)
mean(r_mat$lwr, na.rm = TRUE)
mean(r_mat$upr, na.rm = TRUE)


r_plot <- r_mat %>% filter(year >= 2012) %>% 
ggplot(aes(x = year, y = fit)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmax = fit, ymin = lwr, ymax = upr), width = 0) +
  xlab("Year") +
  ylab(expression(paste("Population growth rate (", italic(r), ")"))) +
  scale_x_continuous(breaks = seq(2012, 2022, 2)) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3)) +
  facet_wrap(~Bed)

jpeg("results/Figure5.jpg",
     width = 2400,
     height = 1800,
     res = 300
     )
r_plot
dev.off()


# Figure 3 (Detection probability) ----
p_posts <- melt(posts$p)
names(p_posts) <- c("iteration", "day", "p")
p_posts$day <- sort(unique(sturgeon$day))[p_posts$day]

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
  ylab(expression(paste("Individual detection probability (", italic("p"["t"]), ")"))) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3))  

jpeg("results/Figure3.jpg",
     width = 2400,
     height = 1800,
     res = 300
)
p_plot
dev.off()



# Summary statistics and post-processing files ----
# . Abundance ----
# .. Average number per site by bed ----
avgs <- n_posts %>% 
  group_by(Bed) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) 
avgs

# .. Differences in mean abundance between sites within year ----
mean_diffs <- n_posts %>% 
  select(iteration, Bed, year, N) %>% 
  group_by(iteration, Bed, year) %>%   
  summarize(N = mean(N)) %>%   
  group_by(iteration, Bed) %>% 
  summarize(N = mean(N)) %>% 
  pivot_wider(names_from = Bed, values_from = N) %>% 
  mutate(diff = downstream - upstream)

mean(mean_diffs$diff)
quantile(mean_diffs$diff, c(0.025, 0.975))

# .. Differences in sum of abundance between beds ----
overall_diffs <- n_posts %>% 
  select(iteration, Bed, year, N) %>% 
  group_by(iteration, Bed, year) %>%   
  summarize(N = sum(N)) %>%   
  group_by(iteration, Bed) %>% 
  summarize(N = mean(N)) %>% 
  pivot_wider(names_from = Bed, values_from = N) %>% 
  mutate(diff = downstream - upstream)

mean(overall_diffs$diff)
quantile(overall_diffs$diff, c(0.025, 0.975))

# .. Data files for simulating capture histories ----
grid_means <- n_posts %>% 
  group_by(id, Easting, Northing, Bed) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) %>% 
  mutate(Easting = as.numeric(Easting),
         Northing = as.numeric(Northing))

grid_means_year <- n_posts %>% 
  group_by(year, id, Easting, Northing, Bed) %>% 
  summarize(
    fit = mean(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) %>% 
  mutate(Easting = as.numeric(Easting),
         Northing = as.numeric(Northing))

save(grid_means, file = "results/grid_means.rda")
save(grid_means_year, file = "results/grid_means_year.rda")

# . Detection ----
# .. Mean detection ----
mean(posts$p)
quantile(posts$p, c(0.025, 0.975))

# .. Highest detection ----
p_summary[p_summary$fit == max(p_summary$fit), ]
