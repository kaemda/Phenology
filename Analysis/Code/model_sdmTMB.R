# LOAD PACKAGES ---------------------------------------------------------------
# Data management
library(readxl)
library(xlsx)
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
library(ggplot2)

all_tows_roms <- read.csv("Data/AllTows_200nm_ROMS.csv") %>%
  subset(., year >= 1995 & year <= 2019) %>% 
  subset(gear == "Bongo/Ring" | gear == "Manta")

northAmerica <- st_read("Data/Shapefiles/NorthAmerica/boundary_p_v2.shp") %>%
  subset(COUNTRY != "water/agua/d'eau" & COUNTRY != "FN") %>% 
  group_by(COUNTRY) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  st_transform(5070)

# ADD ERROR MESSAGES------------------------------------------------------
addErrorMessages <- function(fit) {
  
  if("message" %in% colnames(fit)) {
    fit$error = fit$message
  } else {
    fit$error = "None"
  }
  return(fit)
}

#TRADEOFF MODELS-----------------------------------------------------------------
tradeoffModels <- function(formula = NA, mesh, mainFormula, data, shortFormula, species, spatialTerm, speciesRange, modelType, dependentVar) {
  
  if (modelType == "pheno") {
    # VARYING PHENOLOGY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock))"))
    
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                  spatial = spatialTerm,
                  spatiotemporal = "off",
                  silent = T)
    
  } else if (modelType == "geo") {
    # VARYING GEOGRAPHY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12)"))
    
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                  spatial = "on",
                  spatiotemporal = "iid",
                  time = "timeblock",
                  silent = T)
    
  } else if (modelType == "both") {
    # BOTH VARYING GEOGRAPHY AND VARYING PHENOLOGY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock))"))
    
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                  spatial = "on",
                  spatiotemporal = "iid",
                  time = "timeblock",
                  silent = T)
  } else {
    # BASE
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12)"))
    fit <- sdmTMB(formula = formula,
                  data = data,
                  mesh = mesh,
                  family = tweedie(link = "log"),
                  control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                  spatial = spatialTerm,
                  spatiotemporal = "off",
                  silent = T)
  }
  
  fit.estimates <- bind_rows(tidy(fit, effects = "ran_par", conf.int = TRUE) %>%
                               filter(term %in% c("sigma_O", "sigma_E")),
                             tidy(fit, effects = "fixed", conf.int = TRUE)  %>% 
                               filter(!grepl('year', term))) %>%
    mutate(term = factor(term))
  
  write.csv(x = fit.estimates, file = paste0("Results/",  species, "/", species, "_covariate_coefficients_", dependentVar,"_", speciesRange, "_", modelType, "_", shortFormula, ".csv"))
  
  
  # Return model fit
  return(fit)
  
}

# MODEL CHECKS -----------------------------------------------------------------
modelChecks <- function(fit, data) {
  
  # Basic sanity checks
  sanity <- sdmTMB::sanity(fit)
  
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
  
  return(sanity)
  
}

# RUN FUNCTIONS -----------------------------------------------------------------
source("Analysis/Code/getSpeciesData.R")
source("Analysis/Code/getMesh.R")

# Set species and models
speciesLevels = read_xlsx("Data/Species_Info.xlsx", sheet = 1)$species

covariates <- read_xlsx("Analysis/Covariates.xlsx", sheet = 1)
tradeoffs <- read_xlsx("Results/02_Tradeoff Comparison/Tradeoffs_summary.xlsx", sheet = 1) %>% group_by_at(c("species", "main_formula")) %>%
  filter(., delta_sum_loglik == min(delta_sum_loglik)) %>% 
  merge(., covariates[c("main_formula", "short_formula")]) %>% 
  mutate(species = factor(species, levels = speciesLevels)) %>% 
  arrange(species)

# SINGLE SPECIES
# tradeoffs <- subset(tradeoffs, species == "Stenobrachius leucopsarus")

speciesList = tradeoffs$species
modelTypes = c("base", "pheno", "geo", "both")

dependentVar = "abundance_logN1_scaled"
shortFormulas = tradeoffs$short_formula
mainFormulas = tradeoffs$main_formula

# Set up progress bar
pb <- txtProgressBar(min = 0, max = nrow(tradeoffs), char = "=", style = 3)

# Run sdmTMB models
for (i in 1:nrow(tradeoffs)) {
  
  species = as.character(speciesList[i])
  print(species)
  
  # Get species data
  speciesRange = "speciesRange"
  data <- getspeciesData(species = species, speciesRangeSubset = speciesRange) %>% 
    mutate(year_scaled = scale(year), gear_factor = as.numeric(as.factor(gear)), timeblock_factor = as.numeric(timeblock)) 
  
  # Create spatial version
  data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
    st_set_crs(4326) %>% st_transform(5070)
  
  # Make mesh
  mesh <- getMesh(data, data.sf)
  
  for (j in 1:length(modelTypes)) {

    model = modelTypes[j]
    
    fit <- tryCatch({tradeoffModels(
      #formula = "logN1 ~ s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled) + as.factor(gear) + s(month, bs = 'cc', k = 12)",
      data = data,
      mesh = mesh,
      mainFormula = mainFormulas[i], 
      shortFormula = shortFormulas[i],
      species = species,
      speciesRange = speciesRange,
      modelType = model,
      spatialTerm = "on",
      dependentVar = dependentVar)
      
    },  error = function(err) {
      return(data.frame(aic = NA, converged = FALSE, error = err[1]))
    })
    
    fit <- addErrorMessages(fit)
    
    # Save model fit and other data
    save(fit, data, mesh, file = paste0("Results/", species, "/Models/", dependentVar,"_", speciesRange, "_", model, ".rdata"))
    gc()
    
    setTxtProgressBar(pb, i)
  }
}

close(pb) # Close progress bar


