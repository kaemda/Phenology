# Data management
library(xlsx)
library(readxl)
library(tidyverse)
library(dplyr)
library(sf)

# Analysis
library(sdmTMB)
library(ape)
library(blockCV)

# Computational logistics
library(parallel)
library(future)

setwd("C://KDale/Projects/Phenology/")
all_tows_roms <- read.csv("Data/AllTows_200nm.csv") %>%
  subset(., year >= 1995 & year <= 2019) 

# RUN MODEL COMPARISON --------------------------------------------------------------------------------
runModelComparison <- function(mainFormula, formula = NA, data, programSubset = "allPrograms", dependentVar, spatiotemporalType, gearTerm, timeTerm) {
  
  # Construct formula
  if (!is.na(formula)) {
    formula = as.formula(formula)
  } else if (timeTerm == "" & gearTerm == "") {
    formula = as.formula(paste0(dependentVar,"~ ", mainFormula, "+ s(month, bs = 'cc', k = 12)"))
  } else if (gearTerm == "") {
    formula = as.formula(paste0(dependentVar, "~ ", timeTerm, " +", mainFormula, "+ s(month, bs = 'cc', k = 12)"))
  } else if (timeTerm == "") {
    formula = as.formula(paste0(dependentVar, "~ ", gearTerm, " +", mainFormula, "+ s(month, bs = 'cc', k = 12)"))
  } else {
    formula = as.formula(paste0(dependentVar, "~", mainFormula))
  }
  
  if (programSubset != "allPrograms") {
    data <- data[data$program %in% programSubset,]
  }
  
  # Make mesh
  mesh <- make_mesh(data, xy_cols = c("X",  "Y"), n_knots = 200, type= "cutoff_search")
  
  if(dependentVar == "presence") {
    
    fit <- sdmTMB_cv(formula = formula,
                     data = data,
                     mesh = mesh,
                     k_folds = max(data$fold),
                     fold_ids = data$fold,
                     family = binomial(link = "logit"),
                     spatial = "on",
                     spatiotemporal = spatiotemporalType,
                     time = timeTerm,
                     silent = FALSE)
    
  } else if (dependentVar == "logN1" | dependentVar == "catch_anomaly_positive") {
    
    fit <- sdmTMB_cv(formula = formula,
                     data = data,
                     mesh = mesh,
                     k_folds = max(data$fold),
                     fold_ids = data$fold,
                     family = tweedie(link = "log"),
                     spatial = "on",
                     spatiotemporal = spatiotemporalType,
                     time = timeTerm,
                     silent = FALSE)
  }
  return(fit)
}

# TRADEOFF MODEL COMPARISON ----------------------------------------------------
tradeoffModelComparison <- function(formula = NA, mainFormula, data, shortFormula, species, speciesRange, varying, programSubset = "allPrograms", dependentVar, gearTerm) {
  
  if (programSubset != "allPrograms") {
    data <- data[data$program %in% programSubset,]
  }  
  
  # Make mesh
  mesh <- make_mesh(data, xy_cols = c("X",  "Y"), n_knots = 200, type= "cutoff_search")
  
  if (varying == "pheno") {
    # VARYING PHENOLOGY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock)) +", gearTerm))
    
    fit <- sdmTMB_cv(formula = formula,
                     data = data,
                     mesh = mesh,
                     k_folds = max(data$fold),
                     fold_ids = data$fold,
                     family = tweedie(link = "log"),
                     control = sdmTMBcontrol(nlminb_loops = 2),
                     spatial = "on",
                     spatiotemporal = "off",
                     silent = FALSE)
    
  } else if (varying == "geo") {
    # VARYING GEOGRAPHY
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12) +", gearTerm))
    
    fit <- sdmTMB_cv(formula = formula,
                     data = data,
                     mesh = mesh,
                     k_folds = max(data$fold),
                     fold_ids = data$fold,
                     family = tweedie(link = "log"),
                     control = sdmTMBcontrol(nlminb_loops = 2),
                     spatial = "on",
                     spatiotemporal = "iid",
                     time = "timeblock",
                     silent = FALSE)
    
  } else if (varying == "both") {
    
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock)) +", gearTerm))
    
    fit <- sdmTMB_cv(formula = formula,
                     data = data,
                     mesh = mesh,
                     k_folds = max(data$fold),
                     fold_ids = data$fold,
                     family = tweedie(link = "log"),
                     control = sdmTMBcontrol(nlminb_loops = 2),
                     spatial = "on",
                     spatiotemporal = "iid",
                     time = "timeblock",
                     silent = FALSE)
  } else {
    # BASE
    formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12) +", gearTerm))
    fit <- sdmTMB_cv(formula = formula,
                     data = data,
                     mesh = mesh,
                     k_folds = max(data$fold),
                     fold_ids = data$fold,
                     family = tweedie(link = "log"),
                     control = sdmTMBcontrol(nlminb_loops = 2),
                     spatial = "on",
                     spatiotemporal = "off",
                     silent = FALSE)
  } 
  
  # Return model fit
  return(fit)
}


# INDIVIDUAL MODEL FOR COMPARISON-----------------------------------------------------------------------
# species = "Sardinops sagax"
# 
# source("Analysis/Code/getSpeciesData.R")
# data = getspeciesData(species, speciesRangeSubset = "speciesRange")
# data$year_scaled = scale(data$year)
# 
# for (i in 1:length(species)) {
#   response = "logN1"
#   
#  fit.cv <- runModelComparison(
#     #formula = "presence ~ s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled) + (1 | gearGeneral_factor)",
#     mainFormula = "s(sst_scaled) + s(ssh_scaled) +s(salinity_scaled) +s(distance_from_shore_scaled)", 
#     programSubset = "allPrograms", # Default = "AllPrograms"
#     data = data,
#     dependentVar = response, # "presence" or "catch_anomaly_positive"
#     spatiotemporalType = "ar1", # "ar1" or "rw" or "iid" or "off"
#     timeTerm = "", # "as.factor(timeblock)" or ""
#     varying = "none", # "time", "spatial", "none"
#     gearTerm = "as.factor(gearGeneral)" # "as.factor(gearGeneral)" or ""
#   )
#   
# }
# fit.cv$sum_loglik
# fit.cv$elpd

# TRADEOFF MODEL COMPARISON ----------------------------

species = "Anoplopoma fimbria" # "Sebastes jordani", "Glyptocephalus zachirus")

source("Analysis/Code/getSpeciesData.R")

data = getspeciesData(species, speciesRangeSubset = "speciesRange")
data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>% st_transform(5070)
data$year_scaled = scale(data$year)

# Figure out appropriate block size
sac.ln1 <- cv_spatial_autocor(x = data.sf, column = "abundance", plot = T)
# plot(sac.ln1$variograms[[1]])
blocksize <- sac.ln1$range

# Several ways of determining folds - spatial, buffer, nearest neighbor
folds.spatial <- cv_spatial(x = data.sf, size = blocksize, 
                    column = "abundance",
                    k = 5,
                    selection = "random",
                    iteration = 100)

data$fold = folds.spatial$folds_ids

for (i in 1:length(species)) {
  
  dependentVar = "logN1"
  shortFormula ="0 + sst+ssh+salinity+dfs+month"
  mainFormula = "0 + (sst_scaled) + (ssh_scaled) + (salinity_scaled) + (distance_from_shore_scaled)"
  #shortFormula ="sst+salinity+dfs+month"
  #mainFormula = "(sst_scaled) + (salinity_scaled) + (distance_from_shore_scaled)"
  gearTerm = "as.factor(gearGeneral)"
  
  base.cv <- tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula , 
    shortFormula = shortFormula,
    species = species[i],
    speciesRange = speciesRange,
    varying = "base", #"geo", "pheno", or "base"
    programSubset = c("allPrograms"), # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm # "as.factor(gearGeneral)" or ""
  )
  geo.cv <- tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula, 
    shortFormula = shortFormula,
    species = species[i],
    speciesRange = speciesRange,
    varying = "geo", #"geo", "pheno", or "base"
    programSubset = c("allPrograms"), # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm # "as.factor(gearGeneral)" or ""
  )
  pheno.cv <- tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula, 
    shortFormula = shortFormula,
    species = species[i],
    speciesRange = speciesRange,
    varying = "pheno", #"geo", "pheno", or "base"
    programSubset = c("allPrograms"), # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm # "as.factor(gearGeneral)" or ""
  )
  
  both.cv <- tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula,
    shortFormula = shortFormula,
    species = species[i],
    speciesRange = speciesRange,
    varying = "both",
    programSubset = c("allPrograms"), # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm
  )
}
base.cv$elpd
base.cv$sum_loglik
geo.cv$elpd
geo.cv$sum_loglik
pheno.cv$elpd
pheno.cv$sum_loglik
both.cv$elpd
both.cv$sum_loglik

# ENVIRO COVARIATE COMPARISON -------------------------------------------------------------------------------

species = "Anoplopoma fimbria"

# Set column names
colnames <- c("species", "responses", "gearTerm","programs","tows", "spatiotemporal", "covariates", "sum_loglik", "elpd")

# Create data frame with 0 rows to hold comparison results
modelComparisonResults <- data.frame(matrix(ncol = length(colnames), nrow = 40))
colnames(modelComparisonResults) = colnames

# Specify model terms to test
responses = c("logN1") # "logN1", "presence" 
speciesRangeSubset = c("speciesRange") #"allTows"
programSubset = c("allPrograms") # Default = "AllPrograms"
spatiotemporalTypes = c("iid")
gearTerms = c("as.factor(gearGeneral)") # ""
timeTerms = "timeblock"
covariates = c("s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled) + s(distance_from_shore_scaled)",
               "s(ssh_scaled) + s(salinity_scaled) + s(distance_from_shore_scaled)",
               "s(ssh_scaled) + s(sst_scaled) + s(distance_from_shore_scaled)",
               "s(sst_scaled) + s(salinity_scaled) + s(distance_from_shore_scaled)",
               "s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled)",
               "s(ssh_scaled) + s(sst_scaled)",
               "s(salinity_scaled) + s(distance_from_shore_scaled)",
               "s(ssh_scaled) + s(salinity_scaled)",
               "s(sst_scaled) + s(distance_from_shore_scaled)",
               "s(ssh_scaled) + s(distance_from_shore_scaled)",
               "s(salinity_scaled)",
               "s(distance_from_shore_scaled)",
               "s(ssh_scaled)",
               "s(sst_scaled)") # "s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled)"

# Run loop
no_cores <- detectCores() - 1
plan(multicore, workers = no_cores)
# plan(multisession)

source("Analysis/Code/getSpeciesData.R")
row = 1 # necessary for adding to overall results dataframe

for(response in 1:length(responses)) {
  for (range in 1:length(speciesRangeSubset)) {
    
    speciesData <- getspeciesData(species, speciesRangeSubset = speciesRangeSubset[range]) %>% 
      mutate(year_scaled = scale(year)[,1])
    
    for(programs in 1:length(programSubset)) {
      for(time in 1:length(timeTerms)) {
        for(gear in 1:length(gearTerms)) {
          for (spatiotemporal in 1:length(spatiotemporalTypes)) {
            for (covariate in 1:length(covariates)) {
              
              fit <- runModelComparison(
                data = speciesData,
                mainFormula = covariates[covariate], 
                programSubset = programSubset[programs],
                dependentVar = responses[response], 
                spatiotemporalType = spatiotemporalTypes[spatiotemporal], 
                gearTerm = gearTerms[gear],
                timeTerm = timeTerms[time]
              )
              
              modelComparisonResults$species[row] = species
              modelComparisonResults$responses[row] = responses[response]
              modelComparisonResults$tows[row] = speciesRangeSubset[range]
              modelComparisonResults$programs[row] =   paste(programSubset, collapse = "-")
              modelComparisonResults$gearTerm[row] = gearTerms[gear]
              modelComparisonResults$spatiotemporal[row] = spatiotemporalTypes[spatiotemporal]
              modelComparisonResults$covariates[row] = covariates[covariate]
              modelComparisonResults$sum_loglik[row] = fit$sum_loglik
              modelComparisonResults$elpd[row] = fit$elpd
              
              row = row + 1
              
              print(paste("Model complete:", species, speciesRangeSubset[range], programSubset[programs], gearTerms[gear], timeTerms[time], sep = " "))
              print(row)
            } 
          }
        }
      }
    }
  }
}


# Save comparison results dataframe
write.xlsx(x = modelComparisonResults, file = paste0("Results/ModelComparisonResults_", species, ".xlsx"))

corrplot(cor(data[c("sst_scaled", "ssh_scaled", "salinity_scaled", "distance_from_shore_scaled")]), method = "number")
