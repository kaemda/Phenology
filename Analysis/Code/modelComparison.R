# Data management
library(xlsx)
library(readxl)
library(tidyverse)
library(dplyr)
library(sf)

# Analysis & models
library(fmesher)
library(sdmTMB)
library(ape)
library(blockCV)
library(corrplot)

# Computational logistics
library(parallel)
library(future)
library(retry)

# Plotting
library(Cairo)

setwd("C://KDale/Projects/Phenology/")
all_tows_roms <- read.csv("Data/AllTows_200nm.csv") %>%
  subset(., year >= 1995 & year <= 2019) 

# # ENVIRO COVARIATE COMPARISON --------------------------------------------------------------------------------
# enviroCovariateComparison <- function(mainFormula, formula = NA, data, programSubset = "allPrograms", dependentVar, spatiotemporalType, gearTerm, timeTerm, mesh) {
#   
#   # Construct formula
# 
#   return(fit)
# }

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
                     spatial = "off",
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

# ENVIRO COVARIATE COMPARISON -------------------------------------------------------------------------------
speciesList <- read_xlsx(path = "Data/species_of_interest.xlsx", sheet = 1) %>% subset(., selected == "x")
speciesList <- speciesList$scientific_name

# Specify model terms to test
responses = "abundance_logN1_scaled"
speciesRangeSubset = "speciesRange"
programSubset = "allPrograms"
spatialTerms = c("off")
gearTerm = "as.factor(gearGeneral)"
timeTerm = "timeblock"

covariates = c("s(sst_scaled, k = 3) + s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
               "s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3) + s(bottom_depth_scaled, k = 3)",
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
               "s(ssh_scaled, k = 3, by = rangePercentile) + s(salinity_scaled, k = 3, by = rangePercentile) + s(bottom_depth_scaled, k = 3, by = rangePercentile)",
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



 speciesList = "Trachurus symmetricus"
# covariates = "s(sst_scaled, k = 3) + s(ssh_scaled, k = 3) + s(salinity_scaled, k = 3) + s(bottom_depth_scaled, k = 3)"

# Create matrix to hold results
ntries = 15
colnames <- c("species", "responses", "gearTerm","programs","tows", "spatial", "main_formula", "converged", "sum_loglik", "error")
modelComparisonResults <- data.frame(matrix(ncol = length(colnames), nrow = length(covariates)*length(speciesList)*length(spatialTerms)))
colnames(modelComparisonResults) = colnames

# Run loop
no_cores <- detectCores() - 1
plan(multicore, workers = no_cores)
# plan(multisession)

northAmerica <- read_sf("C://KDale/GIS/NorthAmerica.shp") %>% sf::st_union()
northAmerica <- sf::st_transform(northAmerica, crs = 5070)

source("Analysis/Code/getSpeciesData.R")
row = 1 # necessary for adding to overall results dataframe

#for (iter in 1:ntries) {
  for(s in 13:length(speciesList)) {
    
    #set.seed(4)
    
    species = speciesList[s]
    
    # Get species data
    data = getspeciesData(species, speciesRangeSubset = speciesRangeSubset, allgear = F)
    
    data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>% 
      st_set_crs(4326) %>% st_transform(5070)
    
    # Make mesh
    #mesh <- make_mesh(data, xy_cols = c("X","Y"), n_knots = 200, type = "kmeans")
    
    # Coastline mesh
    domain <- fmesher::fm_nonconvex_hull_inla(as.matrix(st_coordinates(data.sf)),
                                  concave = -.07,convex = -0.05, resolution=c(100,100))
    
    # Coastline boundary
    #northAmericaBoundary <- as(northAmerica, "Spatial") %>% INLA::inla.sp2segment()
    
    # Calculate max edge. The spatial range is about 1/3 of the study area;
    # the max.edge is about 1/5 of that. 
    spatial.range = diff(range(data$Y))/3
    max.edge = spatial.range/5
    
    mesh <- make_mesh(data, xy_cols = c("X", "Y"), fmesher_func = fmesher::fm_mesh_2d_inla,
                      max.edge = c(1,2)*max.edge, # inner and outer max triangle lengths
                      offset = c(1,2)*max.edge,  # inner and outer border widths
                      cutoff = max.edge/5 # shortest allowed distance between points
    )
    
    # mesh <- make_mesh(data, xy_cols = c("X", "Y"), fmesher_func = fmesher::fm_mesh_2d_inla, 
    #                   max.edge = c(100,300), # inner and outer max triangle lengths
    #                   offset = c(100,300),  # inner and outer border widths
    #                   cutoff = 100 # shortest allowed distance between points
    # )
    #plot(mesh)
    
    # Calculate appropriate block size
    sac.ln1 <- cv_spatial_autocor(x = data.sf,  column = "abundance_scaled", plot = T)
    # plot(sac.ln1$variograms[[1]])
    blocksize <- sac.ln1$range * 2
    
    # Several ways of determining folds - spatial, buffer, nearest neighbor
    folds.spatial <- tryCatch({cv_spatial(x = data.sf, size = blocksize, 
                                          column = "abundance_scaled",
                                          k = 5,
                                          selection = "random",
                                          iteration = 100)},
                              error = function(err) {
                                folds.spatial.revised <- cv_spatial(x = data.sf, size = blocksize/100, 
                                                                    column = "abundance_scaled",
                                                                    k = 5,
                                                                    selection = "random",
                                                                    iteration = 100)
                                return(folds.spatial.revised)
                              })
    
    # Assign fold #s
    data$fold = folds.spatial$folds_ids
    
    for(response in 1:length(responses)) {
      for (i in 1:length(spatialTerms)) {
        for (covariate in 1:length(covariates)) {
          
          formula = as.formula(paste0(responses[response], "~", covariates[covariate], "+ s(month, bs = 'cc', k = 12)"))
          
          #TESTING
 
          fit <- tryCatch({
            expr = sdmTMB_cv(data = data,
                             formula = formula,
                             #formula = abundance_scaled ~ s(ssh_scaled, k = 3) +s(bottom_depth_scaled, k = 3) + s(salinity_scaled, k = 3) + s(sst_scaled, k = 3) + s(month, bs = 'cc', k = 12) + as.factor(gearGeneral),
                             mesh = mesh,
                             k_folds = max(data$fold),
                             fold_ids = data$fold,
                             family = tweedie(link = "log"),
                             control = sdmTMBcontrol(newton_loops = 2, nlminb_loops = 2),
                             spatial = spatialTerms[i],
                             spatiotemporal = "off",
                             silent = FALSE)
          }, error = function(err) {
            print(err)
            return(data.frame(sum_loglik = NA, converged = FALSE, error = err[1]))
          }
          )
          
          modelComparisonResults$species[row] = species
          modelComparisonResults$responses[row] = responses[response]
          modelComparisonResults$tows[row] = speciesRangeSubset
          modelComparisonResults$programs[row] =   paste(programSubset, collapse = "-")
          modelComparisonResults$gearTerm[row] = gearTerm
          modelComparisonResults$spatial[row] = spatialTerms[i]
          modelComparisonResults$main_formula[row] = covariates[covariate]
          
          modelComparisonResults$converged[row] = fit$converged
          modelComparisonResults$sum_loglik[row] = fit$sum_loglik
          
          if("message" %in% colnames(fit)) {
            modelComparisonResults$error[row] = fit$message
          } else {
            modelComparisonResults$error[row] = "None"
          }
          
          row = row + 1
          
          print(paste("Model complete:", species, speciesRangeSubset, programSubset, gearTerm, timeTerm, sep = " "))
          print(row)
          
          gc() # clear unused memory
          
        }
      } 
    }
  }
#}
# Save comparison results dataframe
write.xlsx(x = modelComparisonResults, file = paste0("Results/ModelComparisonResults_Trachurus.xlsx"))
modelComparisonResults <- read.xlsx(file = "Results/ModelComparisonResults_multispecies_notes.xlsx", sheetIndex = 1)

# TESTING
run1 <- read.xlsx(file = "Results/ModelComparisonResults_run1.xlsx", sheetIndex = 1)
run2 <- read.xlsx(file = "Results/ModelComparisonResults_run2.xlsx", sheetIndex = 1)
runs <- merge(run1, run2, by = c("species", "main_formula"))
write.xlsx(x = runs, file = paste0("Results/run_comparison.xlsx"))

modelComparison_filtered = subset(modelComparisonResults, !is.na(sum_loglik))
bestModels_sumloglik_converged = modelComparison_filtered %>% subset(converged == T) %>% group_by(species) %>% mutate(delta_sumloglik = max(sum_loglik) - sum_loglik) %>%  top_n(n = 1, wt = sum_loglik) %>% ungroup()
bestModels_sumloglik = modelComparison_filtered %>% group_by(species) %>% mutate(delta_sumloglik = max(sum_loglik) - sum_loglik) %>%  top_n(n = 1, wt = sum_loglik) %>% ungroup()

write.xlsx(x = bestModels_sumloglik, file = "Results/BestModels.xlsx")

# RUN TRADEOFF MODEL COMPARISON  ----------------------------

speciesModels <- read_xlsx(path = "Results/BestModels.xlsx", sheet = 1)
#speciesModels = subset(speciesModels, speciesModels$species == "Stenobrachius leucopsarus")

speciesList = speciesModels$species 

source("Analysis/Code/getSpeciesData.R")

for (i in 1:length(speciesList)) {
  
  set.seed(4)
  
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
                              folds.spatial.revised <- cv_spatial(x = data.sf, size = blocksize/1000, 
                                                                  column = "abundance_scaled",
                                                                  k = 5,
                                                                  selection = "random",
                                                                  iteration = 100)
                              return(folds.spatial.revised)
                            })
  
  # Assign fold #s
  data$fold = folds.spatial$folds_ids
  
  # Set function arguments
  dependentVar = "abundance_logN1_scaled"
  shortFormula = speciesModels$short_formula[i]
  mainFormula = speciesModels$main_formula[i]
  gearTerm = "as.factor(gearGeneral)"
  
  base.cv <- tradeoffModelComparison(
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
  )
  
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
    return(data.frame(sum_loglik = NA, converged = FALSE))
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
    return(data.frame(sum_loglik = NA, converged = FALSE))
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
    return(data.frame(sum_loglik = NA, converged = FALSE))
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
    return(data.frame(sum_loglik = NA, converged = FALSE))
  })
  
  modelComparisonOutput <- data.frame(species = species, model = c("base", "geo", "pheno", "both"), sum_loglik = NA, converged = NA)
  
  # Fill in output sheet
  modelComparisonOutput$sum_loglik = c(base.cv$sum_loglik, geo.cv$sum_loglik, pheno.cv$sum_loglik, both.cv$sum_loglik)
  modelComparisonOutput$converged = c(base.cv$converged, geo.cv$converged, pheno.cv$converged, both.cv$converged)
  modelComparisonOutput$main_formula = speciesModels$main_formula[i]
  
  write.csv(modelComparisonOutput, file = paste0("Results/", species, "/", species, "_tradeoffComparisonOutput.csv"))
  
  # clear memory
  gc()
}


# EVALUATE TRADEOFFS -----------------------------------------------------------
speciesModels <- read_xlsx(path = "Results/BestModels.xlsx", sheet = 1)
speciesList = speciesModels$species

tradeoffSummary = data.frame(species = NA, model = NA,  sum_loglik = NA, delta_sum_loglik = NA, converged = NA)

for (i in 1:length(speciesList)) {
  
  species = speciesList[i]
  
  # Extract csv with tradeoffs results for each species
  modelComparisonOutput <- read.csv(file = paste0("Results/", species, "/", species, "_tradeoffComparisonOutput.csv")) %>%
    subset(., !is.na(sum_loglik)) %>% # subset to models that ran successfully
    mutate(., species = speciesList[i], delta_sum_loglik = max(sum_loglik) - sum_loglik) # calculate deltas
  
  # If there were no models that functioned for the species
  if(nrow(modelComparisonOutput) == 0) next
  
  # Add to overall dataframe
  tradeoffSummary <- bind_rows(tradeoffSummary, modelComparisonOutput) %>% subset(., !is.na(species))
  
}

tradeoffSummary$species = factor(tradeoffSummary$species, levels = c("Trachurus symmetricus", "Engraulis mordax", "Sardinops sagax", "Citharichthys stigmaeus", "Citharichthys sordidus", "Sebastes jordani", "Sebastes paucispinis", "Glyptocephalus zachirus", "Parophrys vetulus", "Merluccius productus", "Stenobrachius leucopsarus", "Tarletonbeania crenularis", "Triphoturus mexicanus", "Vinciguerria lucetia"))

# Save file
write.xlsx(x = tradeoffSummary, "Results/Tradeoffs.xlsx")

# Plot tradeoffs 
CairoPDF("Figures/ModelComparison.pdf", width = 7, height = 6)
ggplot(tradeoffSummary, mapping = aes(x = model, y = species, fill = delta_sum_loglik)) +
  geom_tile() +
  theme_classic(base_size = 16) +
  scale_fill_gradient("\u2206 sum of\nlog-likelihood", low = "goldenrod3", high = "lightyellow2") +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(y = "")
dev.off()

# SANDBOX ----------------------------------------------------------------------

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
