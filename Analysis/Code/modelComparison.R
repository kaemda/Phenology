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
library(corrplot)

# Computational logistics
library(parallel)
library(future)
library(retry)

setwd("C://KDale/Projects/Phenology/")
all_tows_roms <- read.csv("Data/AllTows_200nm.csv") %>%
  subset(., year >= 1995 & year <= 2019) 

# ENVIRO COVARIATE COMPARISON --------------------------------------------------------------------------------
enviroCovariateComparison <- function(mainFormula, formula = NA, data, programSubset = "allPrograms", dependentVar, spatiotemporalType, gearTerm, timeTerm, mesh) {
  
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
    
  } else if (dependentVar == "logN1" | dependentVar == "abundance_logN1_scaled") {
    
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
  return(fit)
}

# TRADEOFF MODEL COMPARISON ----------------------------------------------------
tradeoffModelComparison <- function(formula = NA, mainFormula, data, shortFormula, species, speciesRange, varying, programSubset = "allPrograms", dependentVar, gearTerm, mesh) {
  
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

# RUN TRADEOFF MODEL COMPARISON  ----------------------------

speciesModels <- read_xlsx(path = "Results/BestModels.xlsx", sheet = 1)

speciesList = speciesModels$species

source("Analysis/Code/getSpeciesData.R")

for (i in 1:length(speciesList)) {
  
  species = speciesList[i]
  
  # Get species data
  data = getspeciesData(species, speciesRangeSubset = "speciesRange")
  data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
    st_set_crs(4326) %>% st_transform(5070)
  data$year_scaled = scale(data$year)
  
  programSubset = c("allPrograms")
  
  if (programSubset != "allPrograms") {
    data <- data[data$program %in% programSubset,]
  }  
  
  # Make mesh
  mesh <- make_mesh(data, xy_cols = c("X",  "Y"), n_knots = 100, type= "cutoff_search")
  
  # Calculate appropriate block size
  sac.ln1 <- cv_spatial_autocor(x = data.sf,  column = "abundance", plot = T)
  # plot(sac.ln1$variograms[[1]])
  blocksize <- sac.ln1$range * 2
  
  # Several ways of determining folds - spatial, buffer, nearest neighbor
  folds.spatial <- tryCatch({cv_spatial(x = data.sf, size = blocksize, 
                                        column = "abundance_scaled",
                                        k = 5,
                                        selection = "random",
                                        iteration = 100)},
                            error = function(err) {
                              folds.spatial.try <- cv_spatial(x = data.sf, size = blocksize/1000, 
                                                              column = "abundance_scaled",
                                                              k = 5,
                                                              selection = "random",
                                                              iteration = 100)
                              return(folds.spatial.try)
                            })
  
  # Assign fold #s
  data$fold = folds.spatial$folds_ids
  
  # Set function arguments
  dependentVar = "abundance_logN1_scaled"
  shortFormula = speciesModels$short_formula[i]
  mainFormula = speciesModels$main_formula[i]
  #shortFormula ="sst+ssh+salinity+bd+month"
  #mainFormula = "s(sst_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile) + s(salinity_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)"
  gearTerm = "as.factor(gearGeneral)"
  
  # Run models
  base.cv <- tryCatch({tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula, 
    shortFormula = shortFormula,
    species = species,
    speciesRange = speciesRange,
    varying = "base", #"geo", "pheno", or "base"
    programSubset = programSubset, # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm, # "as.factor(gearGeneral)" or ""
    mesh = mesh
  )},  error = function(err) {
    return(data.frame(elpd = NA, sum_loglik = NA, converged = FALSE))
  })
  
  geo.cv <- tryCatch({tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula, 
    shortFormula = shortFormula,
    species = species,
    speciesRange = speciesRange,
    varying = "geo", #"geo", "pheno", or "base"
    programSubset = programSubset, # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm, # "as.factor(gearGeneral)" or ""
    mesh = mesh
  )},  error = function(err) {
    return(data.frame(elpd = NA, sum_loglik = NA, converged = FALSE))
  })
  
  pheno.cv <- tryCatch({ tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula, 
    shortFormula = shortFormula,
    species = species,
    speciesRange = speciesRange,
    varying = "pheno", #"geo", "pheno", or "base"
    programSubset = programSubset, # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm, # "as.factor(gearGeneral)" or ""
    mesh = mesh
  )},  error = function(err) {
    return(data.frame(elpd = NA, sum_loglik = NA, converged = FALSE))
  })
  
  both.cv <- tryCatch({ tradeoffModelComparison(
    data = data,
    mainFormula = mainFormula,
    shortFormula = shortFormula,
    species = species,
    speciesRange = speciesRange,
    varying = "both",
    programSubset = programSubset, # some combination of programs or "allPrograms" (default)
    dependentVar = dependentVar, # "presence", "logN1", or "catch_anomaly_positive"
    gearTerm = gearTerm, # "as.factor(gearGeneral)" or ""
    mesh = mesh
  )},  error = function(err) {
    return(data.frame(elpd = NA, sum_loglik = NA, converged = FALSE))
  })
  
  modelComparisonOutput <- data.frame(model = c("base", "geo", "pheno", "both"), elpd = NA, elpd_se = NA, sum_loglik = NA)
  
  modelComparisonOutput$elpd = c(base.cv$elpd, geo.cv$elpd, pheno.cv$elpd, both.cv$elpd)
  modelComparisonOutput$sum_loglik = c(base.cv$sum_loglik, geo.cv$sum_loglik, pheno.cv$sum_loglik, both.cv$sum_loglik)
  modelComparisonOutput$converged = c(base.cv$converged, geo.cv$converged, pheno.cv$converged, both.cv$converged)
  
  write.csv(modelComparisonOutput, file = paste0("Results/", species, "/", species, "_tradeoffComparisonOutput.csv"))
  
  # clear memory
  gc()
}

# ENVIRO COVARIATE COMPARISON -------------------------------------------------------------------------------
speciesModels <- read_xlsx(path = "Results/ModelComparisonResults_all_species.xlsx", sheet = 1) %>% subset(., best_model == "x")

speciesList = speciesModels$species

# Specify model terms to test
responses = c("abundance_logN1_scaled")
speciesRangeSubset = c("speciesRange") #"allTows"
programSubset = c("allPrograms") # Default = "AllPrograms"
spatiotemporalTypes = c("iid")
gearTerms = c("as.factor(gearGeneral)") # ""
timeTerms = "timeblock"

covariates = c("s(sst_scaled, k = 3) + s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3, ) + s(bottom_depth_scaled, k = 3)",
               "s(sst_scaled, k = 3)  + s(salinity_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(sst_scaled, k = 3) + s(ssh_scaled, k = 3) +  s(bottom_depth_scaled, k = 3)",
               "s(sst_scaled, k = 3) + s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3)",
               "s(spice_scaled, k = 3) + s(ssh_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(spice_scaled, k = 3) + s(ssh_scaled, k = 3)",
               "s(spice_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(sst_scaled, k = 3)",
               "s(spice_scaled, k = 3)",
               "s(ssh_scaled, k = 3)",
               "s(salinity_scaled, k = 3)",
               "s(bottom_depth_scaled, k = 3)",
               "s(sst_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(ssh_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(salinity_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(sst_scaled, k = 3) + s(ssh_scaled, k = 3)",
               "s(salinity_scaled, k = 3) + s(ssh_scaled, k = 3)",
               "s(sst_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile) + s(salinity_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(ssh_scaled, k = 3, by = rangePercentile) + s(salinity_scaled, k = 3, by = rangePercentile, ) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(sst_scaled, k = 3, by = rangePercentile)  + s(salinity_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(sst_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile) +  s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(sst_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile) + s(salinity_scaled, k = 3, by = rangePercentile)",
               "s(spice_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(spice_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile)",
               "s(spice_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(sst_scaled, k = 3, by = rangePercentile)",
               "s(spice_scaled, k = 3, by = rangePercentile)",
               "s(ssh_scaled, k = 3, by = rangePercentile)",
               "s(salinity_scaled, k = 3, by = rangePercentile)",
               "s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(sst_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(ssh_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(salinity_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
               "s(sst_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile)",
               "s(salinity_scaled, k = 3, by = rangePercentile) + s(ssh_scaled, k = 3, by = rangePercentile)",
               "s(salinity_scaled, k = 3, by = rangePercentile) + s(sst_scaled, k = 3, by = rangePercentile)")



# speciesList = "Glyptocephalus zachirus"
# covariates = "s(sst_scaled, k = 3, by = rangePercentile)"

# Create matrix to hold results
colnames <- c("species", "responses", "gearTerm","programs","tows", "spatiotemporal", "main_formula", "converged", "sum_loglik", "elpd")
modelComparisonResults <- data.frame(matrix(ncol = length(colnames), nrow = length(covariates)*length(speciesList)))
colnames(modelComparisonResults) = colnames

# Run loop
no_cores <- detectCores() - 1
plan(multicore, workers = no_cores)
# plan(multisession)

source("Analysis/Code/getSpeciesData.R")
row = 1 # necessary for adding to overall results dataframe
ntries = 5

for(s in 8:length(speciesList)) {
  for(response in 1:length(responses)) {
    for (range in 1:length(speciesRangeSubset)) {
      for(programs in 1:length(programSubset)) {
        
        species = speciesList[s]
        
        # Get species data
        data = getspeciesData(species, speciesRangeSubset = "speciesRange", allgear = F)
        data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
          st_set_crs(4326) %>% st_transform(5070)
        data$year_scaled = scale(data$year)
        
        programSubset = programSubset[programs]
        
        if (programSubset != "allPrograms") {
          data <- data[data$program %in% programSubset,]
        }  
        
        # Make mesh
        mesh <- make_mesh(data, xy_cols = c("X",  "Y"), n_knots = 100, type= "cutoff_search")
        
        # Calculate appropriate block size
        sac.ln1 <- cv_spatial_autocor(x = data.sf,  column = "abundance", plot = T)
        # plot(sac.ln1$variograms[[1]])
        blocksize <- sac.ln1$range * 2
        
        # Several ways of determining folds - spatial, buffer, nearest neighbor
        folds.spatial <- tryCatch({cv_spatial(x = data.sf, size = blocksize, 
                                              column = "abundance_scaled",
                                              k = 5,
                                              selection = "random",
                                              iteration = 100)},
                                  error = function(err) {
                                    folds.spatial.revised <- cv_spatial(x = data.sf, size = blocksize/1000, 
                                                                        column = "abundance_scaled",
                                                                        k = 5,
                                                                        selection = "random",
                                                                        iteration = 100)
                                    return(folds.spatial.revised)
                                  })
        
        # Assign fold #s
        data$fold = folds.spatial$folds_ids
        
        for(time in 1:length(timeTerms)) {
          for(gear in 1:length(gearTerms)) {
            for (spatiotemporal in 1:length(spatiotemporalTypes)) {
              for (covariate in 1:length(covariates)) {
                
                # fit <- retry::retry(enviroCovariateComparison(
                #   data = data,
                #   mainFormula = covariates[covariate],
                #   programSubset = programSubset[programs],
                #   dependentVar = responses[response],
                #   spatiotemporalType = spatiotemporalTypes[spatiotemporal],
                #   gearTerm = gearTerms[gear],
                #   timeTerm = timeTerms[time]
                # ), max_tries = ntries)
                
                fit <- tryCatch({
                  expr = enviroCovariateComparison(
                    data = data,
                    mainFormula = covariates[covariate],
                    programSubset = programSubset[programs],
                    dependentVar = responses[response],
                    spatiotemporalType = spatiotemporalTypes[spatiotemporal],
                    gearTerm = gearTerms[gear],
                    timeTerm = timeTerms[time],
                    mesh = mesh)
                }, error = function(e) {
                  print(e)
                  return(data.frame(sum_loglik = NA, elpd = NA, converged = FALSE))
                }
                )
                
                # testing
                # fit <- enviroCovariateComparison(
                #   data = data,
                #   mainFormula = covariates[covariate],
                #   programSubset = programSubset[programs],
                #   dependentVar = responses[response],
                #   spatiotemporalType = spatiotemporalTypes[spatiotemporal],
                #   gearTerm = gearTerms[gear],
                #   timeTerm = timeTerms[time],
                #   mesh = mesh
                # )
                
                modelComparisonResults$species[row] = species
                modelComparisonResults$responses[row] = responses[response]
                modelComparisonResults$tows[row] = speciesRangeSubset[range]
                modelComparisonResults$programs[row] =   paste(programSubset, collapse = "-")
                modelComparisonResults$gearTerm[row] = gearTerms[gear]
                modelComparisonResults$spatiotemporal[row] = spatiotemporalTypes[spatiotemporal]
                modelComparisonResults$main_formula[row] = covariates[covariate]
                
                modelComparisonResults$converged[row] = fit$converged
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
  
}
# Save comparison results dataframe
write.xlsx(x = modelComparisonResults, file = paste0("Results/ModelComparisonResults_multispecies.xlsx"))
modelComparisonResults <- read.xlsx(file = "Results/ModelComparisonResults_multispecies_notes.xlsx", sheetIndex = 1)

modelComparison_filtered = subset(modelComparisonResults, !is.na(sum_loglik) & converged == T & rangePercentile == "N")
bestModels_sumloglik = modelComparison_filtered %>% group_by(species) %>% mutate(delta_elpd = max(elpd) - elpd, delta_sumloglik = max(sum_loglik) - sum_loglik) %>%  top_n(n = 1, wt = sum_loglik) %>% ungroup()
bestModels_elpd = modelComparison_filtered %>% group_by(species) %>% mutate(delta_elpd = max(elpd) - elpd, delta_sumloglik = max(sum_loglik) - sum_loglik) %>% top_n(n = 1, wt = elpd) %>% ungroup()

write.xlsx(x = bestModels_sumloglik, file = "Results/BestModels.xlsx")


# RUN FUNCTIONS ----------------------------------------------------------------

# SANDBOX --------------------------------------------------------------------

# Check for collinearity
corrplot(cor(data[c("sst_scaled", "ssh_scaled", "salinity_scaled", "distance_from_shore_scaled", "bottom_depth_scaled", "spice_scaled")]), method = "number")

species = "Trachurus symmetricus"
data = getspeciesData(species = species, speciesRangeSubset = "speciesRange", allgear = T)

programComparisonData = subset(data, program == "PRS_juveniles" | program == "PRS_larvae") %>% 
  group_by_at(c("year", "date",  "program")) %>% summarize(abundance = mean(abundance_logN1_scaled, na.rm = T)) %>% 
  pivot_wider(names_from = program, values_from = abundance)

ggplot(programComparisonData, mapping = aes(x = PRS_juveniles, y = PRS_larvae)) +
  geom_point(color = "goldenrod2") +
  theme_classic(base_size = 14) +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = 2) +
  labs(x = "Average scaled\nabundance (trawl data)", y = "Average scaled\nabundance (bongo data)") +
  ggtitle(species) +
  geom_smooth(method = "lm", color = "gold", fill = "gold")

lm(data = programComparisonData, PRS_juveniles ~ PRS_larvae) %>% summary()
