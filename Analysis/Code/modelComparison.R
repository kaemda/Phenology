# LOAD LIBRARIES ------------------------------
# Data management
library(readxl)
library(dplyr)
library(sf)

# Analysis & models
library(sdmTMB)
library(blockCV)
library(corrplot)

# Computational logistics
library(parallel)
library(future)

# Plotting
library(Cairo)

# GET FOLDS --------------------------------------------------------------------
getFolds <- function(data.sf) {
  
  column = "presence"
  selection = "random"
  smallBlocks = 10
  smallerBlocks = 100
  smallestBlocks = 1000
  
  # Get spatial autocorrelation distance and examine variogram plot
  sac.dist <- cv_spatial_autocor(x = data.sf,  column = "abundance_logN1_scaled", plot = F)
  plot(sac.dist$variograms[[1]])
  
  # Calculate block size (recommended to be substantially larger than SAC)
  blocksize <- sac.dist$range * 2
  
  # Several ways of determining folds - spatial, buffer, nearest neighbor
  folds.spatial <- tryCatch({
    
    cv_spatial(x = data.sf, size = blocksize,
               column = column,
               k = 5,
               selection = selection,
               iteration = 100, report = F)
    
    # print(paste(data.sf$scientific_name[1],"Normal folds"))
  },
  error = function(err) {
    tryCatch({folds.spatial.revised <- cv_spatial(x = data.sf, size = blocksize/smallBlocks, 
                                                  column = column,
                                                  k = 5,
                                                  selection = selection,
                                                  iteration = 100, report = F)
    print(paste(data.sf$scientific_name[1], paste0("Folds/", smallBlocks)))
    return(folds.spatial.revised )
    },
    error = function(err) {
      folds.spatial.revised  <- cv_spatial(x = data.sf, size = blocksize/smallerBlocks, 
                                           column = column,
                                           k = 5,
                                           selection = selection,
                                           iteration = 100, report = F)
      print(paste(data.sf$scientific_name[1], paste0("Folds/", smallerBlocks)))
      return(folds.spatial.revised )
    },
    error = function(err) {
      folds.spatial.revised  <- cv_spatial(x = data.sf, size = blocksize/smallestBlocks, 
                                           column = column,
                                           k = 5,
                                           selection = selection,
                                           iteration = 100, report = F)
      print(paste(data.sf$scientific_name[1], paste0("Folds/", smallestBlocks)))
      return(folds.spatial.revised )
    }
    
    )
    
  })
  
  cairo_pdf(paste0("Figures/", data.sf$scientific_name[i], "/Folds.pdf")) 
  print(cv_plot(cv = folds.spatial, x = data.sf))
  dev.off()
  
  return(folds.spatial)
}

# GET COVARIATES ----------------------------------------------------------------
getCovariates <- function(covariates) {
  # Create a list of combinations for 1 to length(covariates)
  combinations_list <- lapply(1:length(covariates), function(r) {
    gtools::combinations(n = length(covariates), r = r, v = covariates, repeats.allowed = FALSE) %>% 
      as.data.frame()
  })
  
  # Combine all combinations into one data frame
  vars <- bind_rows(combinations_list)
  
  # Unite all columns into a single formula string, removing NA values
  vars <- vars %>%
    mutate(main_formula = apply(., 1, function(row) paste(na.omit(row), collapse = " + ")))
  
  return(vars)
}

# ADD ERROR MESSAGES------------------------------------------------------
addErrorMessages <- function(fit) {
  
  if("message" %in% colnames(fit)) {
    fit$error = fit$message
  } else {
    fit$error = "None"
  }
  return(fit)
}

# RUN ENVIRO COVARIATE COMPARISON -------------------------------------------------------------------------------
runEnviroCovariateComparison <- function(seeds) {
  
  # # Select species
  # species_changed = c("Engraulis mordax", "Sebastes jordani")
  # speciesList = speciesInfo[speciesInfo$species %in% species_changed,]
  
  # Specify model terms to test
  responses = "abundance_logN1_scaled"
  speciesRangeSubset = "speciesRange"
  programSubset = "allPrograms"
  spatialTerms = "on"
  timeTerm = "timeblock"
  
  # Run model for every non-collinear combination of covariates N times
  for (seed in 1:length(seeds)) {
    
    # Create matrix to hold results
    colnames <- c("species", "responses", "programs","tows", "spatial", "main_formula", "converged", "sum_loglik", "error", "seed")
    modelComparisonResults <- data.frame(matrix(ncol = length(colnames), nrow = nrow(covariates)*nrow(speciesList)*length(spatialTerms)))
    colnames(modelComparisonResults) = colnames
    
    row = 1 # necessary for adding to overall results dataframe
    
    # Set seed
    set.seed(seeds[seed])
    
    for(i in 1:nrow(speciesList)) {
      
      # Get species data
      species = speciesInfo$species[i]
      data = getspeciesData(species, speciesRangeSubset = speciesRangeSubset, allgear = F)
      
      # Create spatial version
      data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
        st_set_crs(4326) %>% st_transform(5070)
      
      # Make mesh
      mesh <- getMesh(data, data.sf)
      
      # Assign folds
      data$fold = getFolds(data.sf = data.sf)$folds_ids
      
      for(response in 1:length(responses)) {
        for (j in 1:length(spatialTerms)) {
          for (covariate in 1:nrow(covariates)) {
            
            # Construct formula
            formula = as.formula(paste0(responses[response], "~", covariates$main_formula[covariate], "+ s(month, bs = 'cc', k = 12)"))
            
            no_cores <- detectCores() - 1
            
            # Run model
            fit <- tryCatch({
              
              plan(multisession, workers = no_cores)
              
              expr = sdmTMB_cv(data = data,
                               formula = formula,
                               mesh = mesh,
                               k_folds = max(data$fold),
                               fold_ids = data$fold,
                               family = tweedie(link = "log"),
                               control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                               spatial = spatialTerms[j],
                               spatiotemporal = "off",
                               silent = TRUE,
                               parallel = T)
            }, error = function(err) {
              print(err)
              return(data.frame(sum_loglik = NA, converged = FALSE, error = err[1]))
            }
            )
            
            # Populate results dataframe
            modelComparisonResults$species[row] = species
            modelComparisonResults$responses[row] = responses[response]
            modelComparisonResults$tows[row] = speciesRangeSubset
            modelComparisonResults$programs[row] =   paste(programSubset, collapse = "-")
            modelComparisonResults$spatial[row] = spatialTerms[j]
            modelComparisonResults$main_formula[row] = covariates$main_formula[covariate]
            modelComparisonResults$short_formula[row] = covariates$short_formula[covariate]
            modelComparisonResults$converged[row] = fit$converged
            modelComparisonResults$sum_loglik[row] = fit$sum_loglik
            modelComparisonResults$seed[row] = seeds[seed]
            
            # Extract error messages, if any
            if("message" %in% colnames(fit)) {
              modelComparisonResults$error[row] = fit$message
            } else {
              modelComparisonResults$error[row] = "None"
            }
            
            print(paste("Model complete:", species, covariates$main_formula[covariate], sep = " "))
            print(row)
            
            row = row + 1
            
            gc() # clear unused memory
            
          }
        } 
      }
    }
    # Save comparison results dataframe
    write.xlsx2(x = modelComparisonResults, file = paste0("Results/01_EnviroComparison/EnviroComparisonResults_seed", seeds[seed], ".xlsx"), rowNames = F)
    
    gc()
  }
}

# EVALUATE ENV COVARIATES ----------------------------------------------------
evaluateEnviroCovariates <- function(seeds) {
  
  # Storage dataframe for results
  env_var_all = data.frame()
  
  # For each seed, extract results (includes all species)
  for (i in 1:length(seeds)) {
    
    file = paste0("Results/01_EnviroComparison/EnviroComparisonResults_seed", seeds[i], ".xlsx")
    
    if(file.exists(file)) {
      env_var_seed <- read_xlsx(file) 
      # Add to overall dataframe
      env_var_all <- bind_rows(env_var_all, env_var_seed) %>% subset(., !is.na(species))
    } else {
      # if that seed hadn't been run for that species
      next
    }
  }
  
  # Retain only models that ran successfully
  modelComparison_filtered = subset(env_var_all, !is.na(sum_loglik))
  
  # Select models with the lowest sum log likelihood that converged
  bestModels_sumloglik_converged = modelComparison_filtered %>%
    subset(main_formula %in% covariates$main_formula) %>% 
    merge(covariates) %>% 
    subset(converged == T) %>%
    group_by(species) %>%
    mutate(delta_sumloglik = max(sum_loglik) - sum_loglik) %>% 
    top_n(n = 1, wt = sum_loglik) %>% 
    ungroup() %>% 
    distinct(across(-seed), .keep_all = T) %>% 
    subset()
  
  write.xlsx2(x = modelComparison_filtered, "Results/01_EnviroComparison/allEnviroModels.xlsx", row.names = F)
  write.xlsx2(x = ungroup(bestModels_sumloglik_converged), file = paste0("Results/01_EnviroComparison/BestModels.xlsx"), rowNames = F)
}

# RUN TRADEOFF MODEL COMPARISON  ------------------------------------------------------------
runTradeoffModelComparison <- function() {  
  
  # Load most-predictive models
  speciesModels <- read_xlsx(path = "Results/01_EnviroComparison/BestModels.xlsx", sheet = 1) %>% 
    merge(covariates[c("main_formula", "short_formula")]) %>% 
    arrange(., species)
  
  # species_changed = c("Engraulis mordax")
  # speciesModels = speciesModels[speciesModels$species %in% species_changed, ]
  
  # Set basics
  dependentVar = "abundance_logN1_scaled"
  speciesRangeSubset = "speciesRange"
  models <- c("base", "geo", "pheno", "both")

  spatialTerm = "on"
  
  no_cores <- detectCores() - 1
  
  for (i in 1:nrow(speciesModels)) {
    
    # Get species data
    species = speciesModels$species[i]
    data = getspeciesData(species, speciesRangeSubset = speciesRangeSubset, allgear = F)
    
    # Create spatial version
    data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
      st_set_crs(4326) %>% st_transform(5070)
    
    # Get mesh
    mesh <- getMesh(data, data.sf)
    
    # Set function arguments
    shortFormula = speciesModels$short_formula[i]
    mainFormula = speciesModels$main_formula[i]
    
    # Construct storage dataframe for results
    modelComparisonOutput <- data.frame(species = rep(species, length(models)), model = NA, sum_loglik = NA, aic = NA, converged.cv = NA, converged = NA, error.cv = NA, error = NA)
    
    # Set parallel session
    plan(list(tweak(multisession, workers = no_cores), multisession))
    
    for (seed in 1:length(seeds)) {
      
      set.seed(seeds[seed])
      
      # Assign folds
      data$fold =  getFolds(data.sf)$folds_ids
      
      # Create list to hold model objects and a unique iterating value
      modelList <- list()
      modelListTypes <- data.frame(modelType = c(NA, NA, NA, NA))
      listIndex = 1
      
      for (j in 1:length(models)) {
        
        tradeoff <- models[j]
        
        # Run regular model
        fit <- tryCatch({tradeoffModels(
          func = "standard",
          data = data,
          mesh = mesh, 
          mainFormula = mainFormula, 
          shortFormula = shortFormula,
          species = species,
          speciesRange = speciesRange,
          tradeoff = models[j],
          spatialTerm = spatialTerm,
          dependentVar = dependentVar,
        )},  error = function(err) {
          return(data.frame(aic = NA, converged = FALSE, error = err[1]))
        })
        
        fit <- addErrorMessages(fit)
        
        # Save regular model object
        save(fit, data, mesh, file = paste0("Results/", species, "/Models/", dependentVar,"_", speciesRangeSubset, "_", models[j], ".rdata"))
        
        # Run cross-validation models
        fit.cv <- tryCatch({tradeoffModels(
          func = "cv",
          data = data,
          mainFormula = mainFormula, 
          shortFormula = shortFormula,
          species = species,
          speciesRange = speciesRange,
          tradeoff = models[j],
          spatialTerm = spatialTerm,
          dependentVar = dependentVar,
          mesh = mesh
        )},  error = function(err) {
          return(data.frame(sum_loglik = NA, converged = FALSE, error = err[1]))
        })
        
        fit.cv <- addErrorMessages(fit.cv)
        
        # Add any working models to the model list
        if (fit.cv$error == "None") {
          modelList[[listIndex]] <- fit.cv
          modelListTypes$modelType[listIndex] <- models[j]
          listIndex = listIndex + 1  
        }
        
        gc()
        
        # Fill in output dataframe
        modelComparisonOutput$model = models
        modelComparisonOutput$sum_loglik[j] = fit.cv$sum_loglik
        modelComparisonOutput$converged.cv[j] = fit.cv$converged
        modelComparisonOutput$main_formula[j] = mainFormula
        modelComparisonOutput$short_formula[j] = shortFormula
        modelComparisonOutput$error.cv[j] = fit.cv$error
        modelComparisonOutput$error[j] = fit$error
        modelComparisonOutput$spatial[j] = spatialTerm
        
        if(fit$error == "None") {
          modelComparisonOutput$aic[j] = AIC(fit)
        }
        print(paste("Models for", species, "under", models[j], "complete"))
      }
      
      # Retrieve and save model weights
      if(!is_empty(modelList) & length(modelList) > 1) {
        fit.weights <- data.frame(modelType = subset(modelListTypes, !is.na(modelType)), weights = sdmTMB_stacking(model_list = modelList)) 
        write.csv(fit.weights, file = paste0("Results/", species, "/", species, "_modelWeights_", seeds[seed], ".csv"), row.names = F)
      }
      
      # Write files
      write.csv(modelComparisonOutput, file = paste0("Results/", species, "/", species, "_tradeoffComparisonOutput_", seeds[seed], ".csv"), row.names = F)
      
    }
  }
}
# TRADEOFF MODELS ----------------------------------------------------
tradeoffModels <- function(func, data, mesh, mainFormula, shortFormula, species, speciesRange, tradeoff, spatialTerm, dependentVar) {
  
  if (func == "cv") {
    if (tradeoff == "pheno") {
      # VARYING PHENOLOGY
      formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock))"))
      
      fit <- sdmTMB_cv(formula = formula,
                       data = data,
                       mesh = mesh,
                       k_folds = max(data$fold),
                       fold_ids = data$fold,
                       family = tweedie(link = "log"),
                       control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                       spatial = spatialTerm,
                       spatiotemporal = "off",
                       silent = TRUE)
      
    } else if (tradeoff == "geo") {
      # VARYING GEOGRAPHY
      formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12)"))
      
      fit <- sdmTMB_cv(formula = formula,
                       data = data,
                       mesh = mesh,
                       k_folds = max(data$fold),
                       fold_ids = data$fold,
                       family = tweedie(link = "log"),
                       control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                       spatial = "on",
                       spatiotemporal = "iid",
                       time = "timeblock",
                       silent = TRUE)
      
    } else if (tradeoff == "both") {
      
      formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock))"))
      
      fit <- sdmTMB_cv(formula = formula,
                       data = data,
                       mesh = mesh,
                       k_folds = max(data$fold),
                       fold_ids = data$fold,
                       family = tweedie(link = "log"),
                       control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                       spatial = "on",
                       spatiotemporal = "iid",
                       time = "timeblock",
                       silent = T)
    } else {
      # BASE
      formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12)"))
      fit <- sdmTMB_cv(formula = formula,
                       data = data,
                       mesh = mesh,
                       k_folds = max(data$fold),
                       fold_ids = data$fold,
                       family = tweedie(link = "log"),
                       control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                       spatial = spatialTerm,
                       spatiotemporal = "off",
                       silent = T)
      
    } 
  } else if (func == "standard") {
    if (tradeoff == "pheno") {
      # VARYING PHENOLOGY
      formula = as.formula(paste0(dependentVar, "~", mainFormula, "+ s(month, bs = 'cc', k = 12, by = as.factor(timeblock))"))
      
      fit <- sdmTMB(formula = formula,
                    data = data,
                    mesh = mesh,
                    family = tweedie(link = "log"),
                    control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                    spatial = spatialTerm,
                    spatiotemporal = "off",
                    silent = TRUE)
      
    } else if (tradeoff == "geo") {
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
                    silent = TRUE)
      
    } else if (tradeoff == "both") {
      
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
    
  }
  # Return model fit
  return(fit)
}
# EVALUATE TRADEOFFS -----------------------------------------------------------
evaluateTradeoffs <- function(seeds) {
  
  # Extract tradeoff datasheets and combine
  for (j in 1:length(seeds)) {
    tradeoffs_seed = data.frame(species = NA, model = NA,  sum_loglik = NA, converged = NA)
    
    for (i in 1:length(speciesList)) {
      species = speciesList[i]
      
      # File name for seed result
      file = paste0("Results/", species, "/", species, "_tradeoffComparisonOutput_", seeds[j],  ".csv")
      
      if(file.exists(file)) {
        # Extract csv with tradeoffs results for each species
        modelComparisonOutput <- read.csv(file = file) %>%
          subset(., converged.cv == TRUE) %>% # subset to models that ran successfully
          mutate(., species = speciesList[i]) %>%
          mutate(., across("spatial", as.character)) 
        
        # If there were no models that functioned for the species, skip
        if(nrow(modelComparisonOutput) == 0) next
        
        # Load weights
        weightFile = paste0("Results/", species, "/", species, "_modelWeights_", seeds[j],  ".csv")
        if(file.exists(weightFile)) {
          weights <- read.csv(weightFile) %>% 
            dplyr::select(-X)
          modelComparisonOutput <- merge(modelComparisonOutput, weights, by.x = "model", by.y = "modelType")
        }
        
        # Add to overall dataframe
        tradeoffs_seed <- bind_rows(tradeoffs_seed, modelComparisonOutput) %>% subset(., !is.na(species))
        
        # Factor and remove unnecessary column (artifact of previous functions to save xlsx files)
        tradeoffs_seed$species = factor(tradeoffs_seed$species, levels = speciesLevels)
        tradeoffs_seed$model = factor(tradeoffs_seed$model, levels = c("base", "geo", "pheno", "both"))
        tradeoffs_seed <- tradeoffs_seed %>%
          dplyr::select(-X)
        
      } else { # if that seed hadn't been run for that species, skip
        next
      }
    }
    
    # Save file
    write.xlsx2(x = tradeoffs_seed, file = paste0("Results/02_Tradeoff Comparison/Tradeoffs" , seeds[j], ".xlsx"), row.names = F)
  }
  
  # Put all runs from all seeds from all species together
  for (i in 1:length(seeds)) {
    if (i == 1) {
      tradeoffs_all_seeds = read.xlsx(paste0("Results/02_Tradeoff Comparison/Tradeoffs", seeds[1], ".xlsx"), sheetIndex = 1)
    } else {
      tradeoffs_all_seeds <- bind_rows(tradeoffs_all_seeds, read.xlsx(paste0("Results/02_Tradeoff Comparison/Tradeoffs", seeds[i], ".xlsx"), sheetIndex = 1)) %>% 
        subset(!is.na(species))
    }
  }
  
  # Add species information and factor
  tradeoffs_all_seeds <- tradeoffs_all_seeds %>% 
    merge(speciesInfo[c("species", "common_name")]) %>%
    mutate(common_name = factor(
      common_name,
      levels = speciesLevelsCommonName))

  # Summarize tradeoff results
  tradeoffs_summary <- tradeoffs_all_seeds %>%
    subset(error == "None") %>%  # retain only models that ran successfully
    ungroup() %>% group_by_at(c("species", "model", "main_formula", "converged.cv")) %>%
    summarize(mean_sum_loglik = mean(sum_loglik, na.rm = T),
              sd_sum_loglik = sd(sum_loglik, na.rm = T),
              sd_weight = sd(weights, na.rm = T),
              weight = mean(weights, na.rm = T),
              mean_aic = mean(aic)) %>% 
    group_by(species) %>% 
    mutate(.after = "mean_sum_loglik",
           delta_sum_loglik = mean_sum_loglik - max(mean_sum_loglik, na.rm = T),
           delta_aic = mean_aic - min(mean_aic, na.rm = T)) %>%
    merge(speciesInfo[c("species", "common_name")]) %>%
    ungroup() %>% 
    mutate(common_name = factor(
      common_name,
      levels = speciesLevelsCommonName)) %>% 
    mutate(model = factor(model, levels = c("base", "geo", "pheno", "both")))
  
  # Create draft table for manuscript
  table2 <- tradeoffs_summary %>% 
    mutate(weight = round(weight, digits = 3),
           delta_sum_loglik = round(delta_sum_loglik, digits = 2),
           mean_sum_loglik = round(mean_sum_loglik, digits = 2)) %>% 
    arrange(., common_name, model) %>% 
    rename(., Species = species, "Common name" = common_name,
           Model = model, "Converged" = converged.cv,
           "Mean model weight" = weight,
           "Main formula" = main_formula,
           "Mean sum log-likelihood" = mean_sum_loglik,
           "\u394 Sum log-likelihood" = delta_sum_loglik)
  
  # Write files
  write.xlsx2(x = table2, file = "Manuscript/Table2_ModelingResults.xlsx", row.names = F)
  write.xlsx2(x = tradeoffs_summary, file = "Results/02_Tradeoff Comparison/Tradeoffs_summary.xlsx", row.names = F)
  write.xlsx2(x = tradeoffs_all_seeds, file = "Results/02_Tradeoff Comparison/Tradeoffs_allSeeds.xlsx", row.names = F)
}

# EXAMINE FOLDS --------------------------------------------------------------------------
examineFolds <- function() {
  
  for (i in 1:length(speciesList)) {
    
    data = getspeciesData(species = speciesList[i], speciesRangeSubset = "speciesRange", allgear = F)
    
    # Create spatial version
    data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
      st_set_crs(4326) %>% st_transform(5070)
    
    getFolds(data.sf = data.sf)
    
    gc()
    
  }
}

# LOAD DATA & SET GLOBAL VARIABLES ---------------------------------------------------------------
source("Analysis/Code/getSpeciesData.R")
source("Analysis/Code/getMesh.R")

all_tows_roms <- read.csv("Data/AllTows_200nm_ROMS.csv") %>%
  subset(., year >= 1995 & year <= 2019) 

speciesInfo <- read_xlsx(path = "Data/Species_Info.xlsx", sheet = 1)
speciesList = speciesInfo$species
speciesLevels = arrange(speciesInfo, speciesInfo$phylogenetic_order)$species
speciesLevelsCommonName = arrange(speciesInfo, speciesInfo$phylogenetic_order)$common_name

northAmerica <- st_read("C://KDale/GIS/NorthAmerica_New2/boundary_p_v2.shp") %>%
  subset(COUNTRY != "water/agua/d'eau" & COUNTRY != "FN") %>% 
  group_by(COUNTRY) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  st_transform(5070) # Albers equal area conic (aka conic albers) - projected

# Evaluate collinearity
corrMatrix = cor(all_tows_roms[c("sst_roms", "spice_roms", "ssh_roms", "salinity_roms", "bottom_depth")])
corrplot(corrMatrix)

# Retrieve covariates -- this must be manually edited to remove collinear formulas!
# covariates = getCovariates(c("s(sst_scaled, k = 3)", "s(ssh_scaled, k = 3)", "s(salinity_scaled, k = 3)", "s(bottom_depth_scaled, k = 3)", "s(spice_scaled, k = 3)"))
# write.xlsx2(vars, file = "Analysis/Covariates_full.xlsx", rowNames = F)
covariates = read.xlsx("Analysis/Covariates.xlsx", sheetIndex = 1)

# RUN FUNCTIONS -----
seeds = seq(1, 10, by = 1)

# Run comparison loops
runEnviroCovariateComparison(seeds) # Takes a long time (~2.5 days for 16 spp and 10 seeds)
evaluateEnviroCovariates(seeds)
runTradeoffModelComparison() # Takes a moderate amount of time (~18 hr for 16 spp and 10 seeds)
evaluateTradeoffs(seeds)

# Optional to just run fold creation process
examineFolds()
