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


# . Plot by longitude and latitude ----
# Both spawning beds
# fish %>%
#   # filter(bed == "Upstream") %>%
#   arrange(count) %>%
#   ggplot(aes(x = Easting, y = Northing, color = count)) +
#   geom_point() +
#   scale_color_gradient(low = "yellow", high = "red") +
#   coord_sf(default_crs = sf::st_crs(26914)) +
#   labs(color = "Count", fill = "Count") +
#   theme_bw() +
#   theme(strip.text = element_text(size = 8),
#         axis.text = element_text(size = 8),
#         panel.spacing.x=unit(2.5, "lines"),
#         axis.title.x = element_text(vjust = -1),
#         axis.title.y = element_text(vjust = 3))

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

for(i in 1:2){ # Site
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


# Spatial occupancy test ----
# . Get a single season of data ----
test_data <- sturgeon %>% 
  filter(year == 2022) %>%
  filter(grepl("Downstream", site)) %>%
  ungroup()


# . The temporary GAM we will take apart ----
# https://masonfidino.com/generalized_additive_occupancy_model/
# .. Get year, site, day combos ----
caps <- cast(test_data, formula = site ~ day,
             value = "count",
             fun.aggregate = "sum")


# .. Get unique sites with coords for GAM ----
my_sites <- unique(test_data$site)
site_coords <- unique(cbind(test_data$Easting, test_data$Northing))

tmp_jags <- mgcv::jagam(
  response ~ s(x, y, k = 15, bs = "ds", m = c(1, 0.5)),
  data = data.frame(
    response = rep(1, length(my_sites)),
    x = site_coords[, 1],
    y = site_coords[, 2]
  ),
  family = "poisson",
  file = "tmp.jags")

# . Single-season n-mixture model with spatial GAM ----
model_string <- 
  "model{
  for(i in 1:nsite){
    # Latent state, b is the smoothing term and X is site coords
    log_nlambda[i] <- inprod(b, X[i,])
    log(nlambda[i]) <- log_nlambda[i]
    N[i] ~ dpois(nlambda[i])
  }
  
  for(i in 1:nsite){
   for(t in 1:ndays){
     # the priors
     lp[i, t] ~ dnorm(mu_lp, tau_lp)
     logit(p[i, t]) <- lp[i, t]
     y[i, t] ~ dbin(p[i, t], N[i])
   }
  }
  
  # Priors for detection
  logit(mu_p) <- mu_lp
  mu_lp ~ dnorm(0, 1)
  tau_lp <- pow(sigma_lp, -2)
  sigma_lp ~ dunif(0, 10)
  
  ## Parametric effect priors for abundance
  b[1] ~ dnorm(0, 1/0.95^2)
  ## prior for s(x,y)
  K1 <- S1[1:14,1:14] * lambda[1] 
  b[2:15] ~ dmnorm(zero[2:15], K1) 
  ## smoothing parameter priors CHECK...
  lambda ~ dgamma(.05, .005)
  rho <- log(lambda)
}"


# . JAGS data ----
ch <- as.matrix(caps[, 2:ncol(caps)])
# ch[is.na(ch)] <- 0

data_list <- list(
  y = ch,
  X = tmp_jags$jags.data$X,                 # GAM Northings/Eastings
  S1 = tmp_jags$jags.data$S1,               # GAM parameters             
  zero = tmp_jags$jags.data$zero,           # GAM parameters
  nsite = length(unique(test_data$site)),   # Site ID
  ndays = length(unique(test_data$day))     # Number of days in data set

)

inits <- function(){list(
  N = rep(max(test_data$count, na.rm = TRUE), length(unique(test_data$site))))
}

# Compile the model ----
my_mod <- jags(
  model = textConnection(model_string),
  parameters.to.save = c(
    "b", "rho", "log_nlambda", "nlambda", "N", "mu_p"),
  data = data_list,
  inits = inits,
  n.chains = 3,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 5)

# Print summary
print(my_mod, digits = 3)

save(my_mod, file = "results/my_mod_test.rda")

# Posteriors -----
posts <- my_mod$BUGSoutput$sims.list


# . Local abundance ----
n_posts <- melt(posts$N)
glimpse(n_posts)
names(n_posts) <- c("iteration", "site", "N")
n_posts$site <- unique(as.character(test_data$site))[n_posts$site]
n_posts <- left_join(n_posts, sites, by = "site")


# . Plot count by Northing and Easting ----
# Both spawning beds
sums <- n_posts %>% 
  group_by(site, Easting, Northing) %>% 
  summarize(
    fit = median(N),
    lwr = quantile(N, 0.025),
    upr = quantile(N, 0.975)) 

sums %>% 
  arrange(fit) %>% 
  ggplot(aes(x = Easting, y = Northing, color = fit)) +
  geom_point(shape = 15, size = 10) +
  scale_color_gradient(low = "yellow", high = "red") +
  coord_sf(default_crs = sf::st_crs(26914)) +
  labs(color = "Count", fill = "Count") +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.spacing.x=unit(2.5, "lines"),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3))


# .. Summary statistics ----
sum(sums$fit)
max(sums$fit)
max(test_data$count, na.rm = TRUE)


# . Detection (totally haven't done this yet) ----
# p_posts <- data.frame(melt(posts$mu_p))
# names(p_posts) <- c("iteration", "group", "p")
# # ggplot(p_posts, aes(p)) +
# #   geom_histogram()
# 
# hist(p_posts$p)
