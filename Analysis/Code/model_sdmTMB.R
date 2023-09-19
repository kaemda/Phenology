# LOAD PACKAGES ---------------------------------------------------------------
# Data management
library(xlsx)
library(readxl)
library(tidyverse)
library(dplyr)

# Analysis
library(sdmTMB)
library(ape)

# Computational logistics
library(parallel)
library(future)

# Mapping
library(sf)

# Visualization
library(visreg)
library(ggplot2)
library(RColorBrewer)

setwd("C://KDale/Projects/Phenology/")

all_tows_roms <- read.csv("Data/AllTows_200nm.csv") %>%
  subset(., year >= 1995 & year <= 2019) %>% 
  subset(gearGeneral == "Bongo/Ring" | gearGeneral == "Cobb MWT" | gearGeneral == "Manta")

northAmerica <- read_sf("C://KDale/GIS/North_South_America.shp")
northAmerica <- sf::st_transform(northAmerica, crs = 5070)

#TRADEOFF MODELS-----------------------------------------------------------------
tradeoffModels <- function(formula = NA, mesh, mainFormula, data, shortFormula, species, speciesRange, varying, programSubset = "allPrograms", dependentVar, gearTerm) {
  
  if (varying == "pheno") {
    # VARYING PHENOLOGY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock)) +", gearTerm))
    
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1),
                  spatial = "on",
                  spatiotemporal = "off",
                  silent = FALSE)
    
  } else if (varying == "geo") {
    # VARYING GEOGRAPHY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12) +", gearTerm))
    
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1),
                  spatial = "on",
                  spatiotemporal = "iid",
                  time = "timeblock",
                  silent = FALSE)
    
  } else if (varying == "both") {
    # BOTH VARYING GEOGRAPHY AND VARYING PHENOLOGY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock)) +", gearTerm))
    
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1),
                  spatial = "on",
                  spatiotemporal = "iid",
                  time = "timeblock",
                  silent = FALSE)
  } else {
    # BASE
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12) +", gearTerm))
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1),
                  spatial = "on",
                  spatiotemporal = "off",
                  silent = FALSE)
  }
  
  fit.estimates <- bind_rows(tidy(fit, effects = "ran_par", conf.int = TRUE) %>%
                               filter(term %in% c("sigma_O", "sigma_E")),
                             tidy(fit, effects = "fixed", conf.int = TRUE)  %>% 
                               filter(!grepl('year', term))) %>%
    mutate(term = factor(term))
  
  # Save model fit and other data
  programs <- paste(programSubset, collapse = "-")
  save(fit, data, programs, speciesRange, mesh, file = paste0("Results/", species, "/Models/", dependentVar,"_", speciesRange, "_", programs, "_", varying, "_", shortFormula, "_", gearTerm, ".rdata"))
  write.csv(x = fit.estimates, file = paste0("Results/",  species, "/", species, "_covariate_coefficients_", dependentVar,"_", speciesRange, "_", programs, "_", varying, "_", shortFormula, "_", gearTerm, ".csv"))
  
  # Return model fit
  return(fit)
  
}

# MODEL CHECKS -----------------------------------------------------------------
modelChecks <- function(fit) {
  
  # Basic sanity checks
  sdmTMB::sanity(fit)
  
  # Randomized quantile residuals - quick check
  data$resids <- residuals(fit)
  qqnorm(data$resids)
  qqline(data$resids)
  
  # View confidence intervals and extract parameters as a dataframe
  # range: A derived parameter that defines the distance at which 2 points are effectively independent
  # (actually about 13% correlated). If the share_range argument is changed to FALSE then the spatial
  # and spatiotemporal ranges will be unique, otherwise the default is for both to share the same range.
  # phi: Observation error scale parameter (e.g., SD in Gaussian).
  # sigma_O: SD of the spatial process ("Omega").
  # sigma_E: SD of the spatiotemporal process ("Epsilon").
  # tweedie_p: Tweedie p (power) parameter; between 1 and 2.
  tidy(fit, effects = "ran_pars", conf.int = TRUE)

}

# RUN FUNCTIONS -----------------------------------------------------------------

# Set species and models
species = "Glyptocephalus zachirus"
tradeoffs = c("both", "geo", "pheno", "base")
tradeoffs = c("geo")

# Get species data
source("Analysis/Code/getSpeciesData.R")
speciesRange = "speciesRange"
data <- getspeciesData(species = species, speciesRangeSubset = speciesRange) %>% 
  mutate(year_scaled = scale(year), gear_factor = as.numeric(as.factor(gearGeneral)), timeblock_factor = as.numeric(timeblock)) 

dependentVar = "logN1"
shortFormula ="sst+ssh+salinity+dfs+bd+month"
mainFormula = "s(sst_scaled, k = 3) + s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3) + s(distance_from_shore_scaled, k = 3) + s(bottom_depth_scaled, k = 3)"
# shortFormula ="ssh+dfs+month"
# mainFormula = "s(ssh_scaled, basis = 'cv')+s(distance_from_shore_scaled, )"
gearTerm = "as.factor(gearGeneral)"

for (i in 1:length(tradeoffs)) {
  
  # Make mesh
  mesh <- make_mesh(data, xy_cols = c("X",  "Y"), n_knots = 200, type= "cutoff_search")
  
  fit <- tradeoffModels(
    #formula = "logN1 ~ s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled) + as.factor(gearGeneral) + s(month, bs = 'cc', k = 12)",
    data = data,
    mesh = mesh,
    mainFormula = mainFormula , 
    shortFormula = shortFormula,
    species = species,
    speciesRange = speciesRange,
    varying = tradeoffs[i],
    programSubset = c("allPrograms"), # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm) # "as.factor(gearGeneral)" or ""

modelChecks(fit)

fit %>% summary()

}
#-----------------------------------------------------------------------------#
# Reload model
modelType = c("logN1_speciesRange_allPrograms_base_sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)")
load(file = paste0("Results/", species, "/Models/", modelType, ".rdata"))

# Plot smooother on variables (include "scale = response" to plot on response scale)
# Doesn't work on geo model (error about time column missing?)
visreg(fit, xvar = "sst_scaled")
visreg(fit, xvar = "ssh_scaled")
visreg(fit, xvar = "salinity_scaled")
visreg(fit, xvar = "month", scale = "response")

# SANDBOX --------------------------------------------------------------------

# Trying to match models that are in the original NSF proposal
formula = as.formula(paste0("logN1 ~", mainFormula, "+ s(month, bs = 'cc', k = 12) +", gearTerm))

fit <- sdmTMB(formula = "logN1 ~ s(ssh_scaled) + s(salinity_scaled) + s(bottom_depth_scaled) + s(month, bs = 'cc', k = 12, by = sst_scaled) + as.factor(gearGeneral)",
              data = data,
              mesh = mesh,
              spatial_varying = ~ 0 + timeblock, time = "timeblock",
              family = tweedie(link = "log"),
              control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1),
              spatial = "on",
              spatiotemporal = "off",
              silent = FALSE)

# ASSUMPTIONS -----------------------------------------------------------------
# Test for normality

## Test if a spatial model should be used (Moran's I) ----------
inv_dists <- as.matrix(dist(speciesData[,c("X","Y")]))
diag(inv_dists) <- 0
Moran.I(speciesData$logN1, inv_dists)

## Test for covariance
glm <- glm(data = data, formula = logN1 ~ ssh_scaled + sst_scaled + salinity_scaled + bottom_depth_scaled + distance_from_shore_scaled)

summary(glm)
car::vif(glm)

#create vector of VIF values
vif_values <- car::vif(glm)

#create horizontal bar chart to display each VIF value
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")
abline(v = 5, lwd = 3, lty = 2)

vars <- data[ , c("ssh_scaled", "sst_scaled", "salinity_scaled", "bottom_depth_scaled", "distance_from_shore_scaled")]

#create correlation matrix
cor(vars)

## Latitude model --------
gam.lat <- mgcv::gam(data = data.latModel, latitude ~ s(day_of_year_weighted) + as.factor(gearGeneral))
summary(gam.lat)
plot(gam.lat)
p.gam = predict(gam.lat, newdata = data.latModel)
plot(x = data.latModel$day_of_year, y = p.gam)

