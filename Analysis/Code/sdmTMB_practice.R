library(sdmTMB)
library(dplyr)
library(visreg)
library(ggplot2)
library(sf)
library(blockCV)

# PLOT MAPS ----------------------------------------------------
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    coord_fixed()
}

# make a mesh ------------------------------------------------
m <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10) # cutuoff = minimum allowed distance between points in units of x and y

#plot(mesh)

# RUNNING BASIC MODELS ----------------------------------------
# fit spatial model with smoother for depth
fit <- sdmTMB(density ~ s(depth, k = 3), data = pcod, mesh = mesh, family = tweedie(link = "log"), spatial = "on")

fit

# Testing NaN error messages and singularity issues
m <- sdmTMB(
  data = pcod,
  formula = density ~ s(depth_scaled, k = 3) + s(depth_scaled2, k = 3) + as.factor(year),
  mesh = make_mesh(pcod, c("X", "Y"), cutoff = 10),
  family = tweedie(link = "log"),
  spatial = "on"
)

# Cross validation
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)

pcod.sf <- st_as_sf(pcod, coords = c("lon", "lat")) %>% 
  st_set_crs(4326) %>% st_transform(5070)

# Calculate appropriate block size
sac.ln1 <- cv_spatial_autocor(x = pcod.sf,  column = "density", plot = T)
blocksize <- sac.ln1$range * 2

# Several ways of determining folds - spatial, buffer, nearest neighbor
folds.spatial <- cv_spatial(x = pcod.sf, size = blocksize/1000, 
                                      column = "density",
                                      k = 5,
                                      selection = "random",
                                      iteration = 100)

# Assign fold #s
pcod$fold = folds.spatial$folds_ids

m_cv <- sdmTMB_cv(
  density ~ s(depth, k = 5) + as.factor(year),
  data = pcod, mesh = mesh, fold_ids = "fold",
  spatial = "on",
  family = tweedie(link = "log"), 
  k_folds = 5
)

# extract parameters as a dataframe
# range: A derived parameter that defines the distance at which 2 points are effectively independent (actually about 13% correlated). If the share_range argument is changed to FALSE then the spatial and spatiotemporal ranges will be unique, otherwise the default is for both to share the same range.
# phi: Observation error scale parameter (e.g., SD in Gaussian).
# sigma_O: SD of the spatial process ("Omega").
# sigma_E: SD of the spatiotemporal process ("Epsilon").
# tweedie_p: Tweedie p (power) parameter; between 1 and 2.

tidy(fit, effects = "ran_pars", conf.int = TRUE)

# run basic model checks
sanity(fit)

# plot smooother on depth in link space with randomized quantile partial residuals (??)
visreg(fit, xvar = "depth", xlim = c(50,500))

# predict on new data
# qcs_grid is an example 2x2 km prediction grid for Queen Charlotte sound
p <- predict(fit, newdata = qcs_grid)
head(p)

# plot predictions
# "est" is the column in p, exponentiated because of the link? Same with transform in viridis?
ggplot(p, aes(X, Y, fill = exp(est))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "sqrt")

# presence absence model
# use binomial distribution obviously
# "present" is a column in the dataset
fit <-
  sdmTMB(
    present ~ s(depth),
    data = pcod,
    mesh = mesh,
    family = binomial(link = "logit")
  )

# hurdle/delta model
fit <-
  sdmTMB(
    density ~ s(depth),
    data = pcod,
    mesh = mesh,
    family = delta_gamma(link1 = "logit", link2 = "log")
  )

# spatiotemporal model
# takes slightly longer (~30s total)
# choose whether spatial random fields are "IID", 1st order autoregressive "AR1"
# or as a random walk "RW"
# Turning "spatial" to on/off only slightly changes results
fit_spatiotemporal <-
  sdmTMB(
    density ~ s(depth, k = 5),
    data = pcod,
    mesh = mesh,
    time = "year",
    family = tweedie(link = "log"),
    spatial = "on",
    spatiotemporal = "ar1"
  )

# 'filling in' time slices with "extra_time" argument
# takes longer (~2 min total)
fit_spatiotemporal <-
  sdmTMB(
    density ~ s(depth, k = 5),
    data = pcod,
    mesh = mesh,
    time = "year",
    family = tweedie(link = "log"),
    spatial = "off",
    spatiotemporal = "ar1",
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016)
  )

# Time-varying effects
fit_timevarying <- sdmTMB(density ~ 0 + as.factor(year),
                          data = pcod,
                          time_varying = ~ 0 + depth_scaled + depth_scaled2,
                          mesh = mesh,
                          family = tweedie(link = "log"),
                          spatial = "on",
                          time = "year", 
                          spatiotemporal = "IID")

fit_timevarying
sanity(fit_timevarying)

# MODEL DIAGNOSTICS ----------------------------------------
# Randomized quantile residuals
pcod$resids <- residuals(fit_spatiotemporal)
qqnorm(pcod$resids)
qqline(pcod$resids)

# Would want more warmup and iterations in practice
# This causes R to crash every time!
# r <- residuals(fit, "mle-mcmc", mcmc_warmup = 100, mcmc_iter = 101)

# PREDICTIONS ------------------------------------------

# area weighted standardized population index 
# predict on a grid covering the entire survey (qcs_grid) with grid cell area 4 and pass predictions to get_index()
# qcs_grid also has data on all covariates used to predict the model
# get_index is an sdmTMB function that extracts relative biomass/abundance or center of gravity
# area is the grid cell area, a vector of length newdata or a value of length 1 which will be repeated
predictions_st <- predict(fit_spatiotemporal, newdata = qcs_grid)

index = get_index(obj = predictions_st, area = rep(4, nrow(qcs_grid)))

# Year versus biomass
ggplot(index, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  labs(x = "Year", y = "Biomass (kg)") +
  theme_classic()

# Show various effects
plot_map(predictions_st, "exp(est)") +
  scale_fill_viridis_c(
    trans = "sqrt",
    na.value = "yellow", limits = c(0, quantile(exp(predictions_st$est), 0.995))
  ) +
  facet_wrap(~year) + 
  ggtitle("Prediction (fixed effects + random effects)",
          subtitle = paste("maximum estimated biomass density =", round (max(exp(predictions_st$est)))))

# time-varying effects
timevarying.df <- expand.grid(depth_scaled = seq(min(pcod$depth_scaled) + 0.2,
                                                 max(pcod$depth_scaled) - 0.2,
                                                 length.out = 50),
                              year = unique(pcod$year))

timevarying.df$depth_scaled2 = timevarying.df$depth_scaled^2

predict <- predict(fit_timevarying, newdata = timevarying.df, se_fit = TRUE, re_form = NA)

ggplot(predict, aes(depth_scaled, exp(est),
                    ymin = exp(est - 1.96 + est_se),
                    ymax = exp(est + 1.96 * est_se),
                    group = as.factor(year))) +
  geom_line(aes(color = year), lwd = 1) + 
  geom_ribbon(aes(fill = year), alpha = 0.1) +
  theme_classic() +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  scale_x_continuous(labels = function(x) round(exp(x * pcod$depth_sd[1] + pcod$depth_mean[1])))+
  coord_cartesian(expand = FALSE) +
  labs(x = "Depth (m)", y = "Biomass density (kg/km2)")

# FURTHER ANALYSIS -------------------------------------
# Center of gravity
cog = get_cog(p_st, format = "wide")
ggplot(cog, aes(est_x, est_y, colour = year)) +
  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y)) +
  scale_colour_viridis_c()

# time-varying coefficients
pcod$year_factor <- as.factor(pcod$year)
fit <- sdmTMB(density ~ 0 + s(depth, k = 5) + (1 | year_factor),
              data = pcod, mesh = mesh,
              time = "year", family = tweedie(link = "log"), 
              silent = FALSE)
fit
tidy(fit, effects = "ran_pars", conf.int = TRUE)

# MY DATA------------------------------------------------
species = "Trachurus symmetricus"

source("Analysis/Code/getSpeciesData.R")

# Get species data
data = getspeciesData(species, speciesRangeSubset = "speciesRange", allgear = F)
data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>% st_transform(5070)
data$year_scaled = scale(data$year)

# Make mesh - very simple version
mesh <- make_mesh(data, xy_cols = c("X",  "Y"), n_knots = 200, type= "cutoff_search")

formula = as.formula(paste0("abundance_logN1_scaled ~ s(sst_scaled, k = 3) + s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3) + s(bottom_depth_scaled, k = 3) + s(month, bs = 'cc', k = 12)"))

fit <- sdmTMB(data = data,
         # Different versions of this formula (different covariates, no month or gear term) all sometimes result in singularity issue
         formula = formula,
         mesh = mesh,
         family = tweedie(link = "log"),
         spatial = "on", # singularity issue occurs with both on and off
         spatiotemporal = "off",
         silent = F
  )


write.csv(data, file = "Data/Tarletonbeania_crenularis.csv")
