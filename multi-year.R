# Libraries----
library(tidyverse)
library(lubridate)
library(reshape)
library(R2jags)

# Dataset----
fish <- read.csv("data/sturgeon_grid_data.csv")
str(fish)

fish <- fish %>% 
  filter(count != "#REF!")

# Date manipulation----
# . QA/QC ----
# Remove spaces
fish$date <- gsub(pattern = " ", replacement = "",
                  x = fish$date)
# Remove spaces from bed name 
fish$bed <- gsub(pattern = " ", replacement = "",
                 x = fish$bed)
fish$bed <- gsub(pattern = "Upper", replacement = "Upstream",
                 x = fish$bed)

# Get average Easting and Northing per grid cell
fish$Easting <- rowMeans(fish[,c('left', 'right')])
fish$Northing <- rowMeans(fish[,c('top', 'bottom')])

# . Summarize by grid ----
fish <- fish %>% 
  filter(!(bed == "Downstream" & Northing < 4983600)) %>% 
  group_by(date, bed, id, Easting, Northing) %>% 
  summarize(count = sum(count)) %>% 
  mutate(rep = row_number(), site = paste0(bed, "-", id))


# . Pivot to long form to add zeroes and NAs (for unsampled days) ----
fish$date <- as.Date(fish$date, format = "%m/%d/%Y")
fish$day <- yday(fish$date)
fish$year <- year(fish$date)

zeroes <- reshape::cast(fish, 
                        formula = bed ~ site ~ day ~ year,
                        value = "count",
                        fun.aggregate = "sum",
                        fill = 0,
                        drop = FALSE,
                        add.missing = TRUE)

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
# capture histories with or without site, date, bed
sturgeon <- reshape2::melt(zeroes)
names(sturgeon)[5] <- "count"

# Get rid of erroneous combos for site/bed
sturgeon$check <- 0
sturgeon$check[sturgeon$bed == "Downstream" & grepl("Downstream", sturgeon$site)] <- 1
sturgeon$check[sturgeon$bed == "Upstream" & grepl("Upstream", sturgeon$site)] <- 1

sturgeon <- sturgeon %>% filter(check == 1)

# Get Eastings and Northings from original data
sites <- data.frame(unique(cbind(fish$site, fish$Easting, fish$Northing)))
names(sites) <- c("site", "Easting", "Northing")
sites[,2:3] <- apply(sites[, 2:3], 2, as.numeric)

# Merge with sturgeon data
sturgeon <- merge(sturgeon, sites)

# Set aside data (for subsetting if needed)
test_data <- sturgeon


# Spatial occupancy test ----
# . The temporary GAM we will take apart ----
# https://masonfidino.com/generalized_additive_occupancy_model/
# .. Get year, site, day combos ----
caps <- cast(test_data, formula = site ~ day ~ year,
             value = "count",
             fun.aggregate = "sum")


# .. Get unique sites with coords for GAM ----
my_sites <- unique(test_data$site)
site_coords <- unique(cbind(test_data$Easting, test_data$Northing))

# .. Write gam parameter structure to a file ----
tmp_jags <- mgcv::jagam(
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
  X = tmp_jags$jags.data$X,                 # GAM Northings/Eastings
  S1 = tmp_jags$jags.data$S1,               # GAM parameters             
  zero = tmp_jags$jags.data$zero,           # GAM parameters
  nsite = length(unique(test_data$site)),   # Site ID
  ndays = length(unique(test_data$day)),    # Number of days in data set
  nyears = length(unique(test_data$year)))

# . Initial values ----
inits <- function(){list(
  N = matrix(max(test_data$count, na.rm = TRUE), 
          length(unique(test_data$site)),
          length(unique(test_data$year))
          ))#,
  # z = rep(1, length(unique(test_data$day))))
}

# Compile the model ----
my_mod <- jags(
  model.file = "models/multi-year.jags",
  parameters.to.save = c("b", "rho", "log_nlambda", 
                         "nlambda", "N", "mu_p", "p"),
  data = data_list,
  inits = inits,
  n.chains = 3,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 5)

# Print summary
print(my_mod, digits = 3)

save(my_mod, file = "results/multi-year.rda")

# Posteriors -----
posts <- my_mod$BUGSoutput$sims.list


# . Local abundance ----
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
    fit = median(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) 

sums %>%
  arrange(fit) %>%
  ggplot(aes(x = Easting, y = Northing, color = fit)) +
  geom_point(shape = 15, size = 2) +
  scale_color_gradient(low = "yellow", high = "red") +
  coord_sf(default_crs = sf::st_crs(26914)) +
  labs(color = "Count", fill = "Count") +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.spacing.x=unit(2.5, "lines"),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3)) +
  facet_wrap(~year)


# .. Summary statistics ----
# ... Max density per site per year ----
maxes <- sums %>% 
  group_by(bed, year) %>% 
  summarize(
    fit = max(fit),
    lwr = max(lwr),
    upr = max(upr)
    ) 

ggplot(maxes, aes(x = year, y = fit, color = bed, fill = bed)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr, color = NULL), 
              alpha = 0.25)

# ... Max abundance per bed per year across sites ----
overalls <- sums %>% 
  group_by(bed, year) %>% 
  summarize(
    fit = sum(fit),
    lwr = sum(lwr),
    upr = sum(upr)) 

ggplot(overalls, aes(x = year, y = fit, color = bed, fill = bed)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr, color = NULL), 
              alpha = 0.25)

# ... Max abundance at whole study area ----
overalls <- sums %>% 
  group_by(year) %>% 
  summarize(
    fit = sum(fit),
    lwr = sum(lwr),
    upr = sum(upr)) 

ggplot(overalls, aes(x = year, y = fit)) + 
  geom_line() +
  geom_ribbon(aes(xmax = year, ymin = lwr, ymax = upr, color = NULL), 
              alpha = 0.25)
