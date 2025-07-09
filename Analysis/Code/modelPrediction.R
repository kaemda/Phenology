# LOAD PACKAGES ---------
# Data management
library(tidyr)
library(dplyr)
library(readxl)
library(writexl)
library(stringr)

# Analysis
library(boot)
library(glmnet)
library(mgcv)
library(mgcViz)
library(sdmTMB)
library(ggeffects)

# Mapping
library(sf)

# Visualization
library(ggplot2)
library(cowplot)
library(ggnewscale)
library(ggpubr)
library(scales)
library(viridis)
library(Cairo)
library(colorspace)

# Computational management
library(furrr)

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# CREATE PREDICTION OBJECTS---------------------------------------------------------------
createPredictionObjects <- function(species, fitName, predictionObjectName, gridFilename, makeNewGrid, makeNewPrediction, stacked) {
  
  # If the prediction grid hasn't already been created, create it (or override)
  if (!file.exists(gridFilename) | makeNewGrid == T) { # Load prediction grid
    source("Analysis/Code/CreatePredictionGrid.R")
    createPredictionGrid(species = species, path = paste0("Analysis/PredictionGrids/", species, "_grid.rdata"))
  }
  
  # Then load prediction grid (and grid.df) from source
  load(gridFilename) 
  
  # Create factored versions of year and timeblock
  prediction_grid_roms <- merge(prediction_grid_roms, timeblocks) %>% 
    mutate(timeblock = factor(timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")))
  prediction_grid_roms$year_scaled = as.vector(scale(prediction_grid_roms$year, center = T, scale = T)[,1])
  
  weights = list(modelType = NA, weights = NA)
  
  # MODEL STACKING
  if (stacked == TRUE) {
    
    for (i in 1:10) {
      # Compile individual weight results
      weightFile = paste0("Results/", species, "/", species, "_modelWeights_", i, ".csv")
      if(! file.exists(weightFile)) next
      weight <- read.csv(weightFile)
      weight <- merge(modelNames, weight, all.x = T) %>% 
        tidyr::replace_na(., list(weights = 0))
      
      if(is.null(nrow(weights))) { # If the storage dataset is empty
        weights = weight
      } else {
        weights = bind_rows(weights, weight)
      }
    }
    
    # Find mean weight across all seeds
    weights = weights %>% group_by(modelType) %>% 
      summarize(weights = mean(weights))
    
    # Write output
    write.csv(weights, paste0("Results/", species, "/Mean_weights.csv"), row.names = F)
    
    models <- list()
    
    # Load model results one at a time, calculate estimates
    # and original stations (p.original)
    for (i in 1:length(weights$modelType)) {
      
      load(paste0("Results/", species, "/Models/abundance_logN1_scaled_speciesRange_", weights$modelType[i],".rdata"))   
      if(fit$error == "None") {
        models[[i]] = fit
      } else {
        models[[i]] = NA
      }
    }
    
    models <- subset(models, !is.na(models))
    
    # Create storage dataframes
    estWeighted.p.original = data.frame(matrix(ncol = length(weights$modelType), nrow = nrow(data)))
    estWeighted.p = data.frame(matrix(ncol = length(weights$modelType), nrow = nrow(prediction_grid_roms)))
    
    # Create prediction objects
    for (i in 1:length(models)) {
      if(typeof(models[[i]]) != "list") next
      p.original <- predict(object = models[[i]], newdata = data)
      p <- predict(object = models[[i]], newdata = prediction_grid_roms) %>% 
        mutate(towID = paste(latitude, longitude, year, month, sep = "_")) 
      
      # Calculate estimated weighted abundance for the model and add to dataframe
      estWeighted.p.original[,i] = models[[i]]$family$linkinv(p.original$est) * weights$weights[i]
      estWeighted.p[,i] = models[[i]]$family$linkinv(p$est) * weights$weights[i]
    }
    
    # Sum across weighted abundances
    p.original$est_retransform_weighted = rowSums(estWeighted.p.original, na.rm = T)
    p$est_retransform_weighted = rowSums(estWeighted.p, na.rm = T)
    
    # Calculate p.obj
    p.obj <- predict(object = models[[i]], newdata = prediction_grid_roms, return_tmb_object = T)
    p.obj$data$est_retransform_weighted = rowSums(estWeighted.p, na.rm = T)
    
  } else if (stacked == FALSE) {
    # Load model results
    if(file.exists(fitName)) {
      load(file = fitName)
    } else {
      print("Model object does not exist for", species, "and", fitName, "-- run modelsdmTMB.R")
      break
    }
    
    # Predict on both original data and on the prediction grid
    p.original <- predict(object = fit, newdata = data) %>%  
      mutate(est_retransform_weighted = fit$family$linkinv(est)) 
    
    # Prediction object
    p.obj <- predict(object = fit, newdata = prediction_grid_roms, return_tmb_object = T)
    
    # Prediction
    p <- predict(object = fit, newdata = prediction_grid_roms) %>% 
      mutate(towID = paste(latitude, longitude, year, month, sep = "_")) %>% 
      mutate(est_retransform_weighted = fit$family$linkinv(est)) 
  }
  
  # Add Albers equal area to p.original
  coords= cbind(p.original$longitude, p.original$latitude)
  scale = 1000
  albert_equal_area <- sf::sf_project(from = "EPSG:4326", to = 'EPSG:5070', pts = coords)/scale
  p.original$X = albert_equal_area[,1]
  p.original$Y = albert_equal_area[,2]
  
  # Add Albers equal area to p object
  coords= cbind(p$longitude, p$latitude)
  scale = 1000
  albert_equal_area <- sf::sf_project(from = "EPSG:4326", to = 'EPSG:5070', pts = coords)/scale
  p$X = albert_equal_area[,1]
  p$Y = albert_equal_area[,2]
  
  # Save prediction objects
  save(p, p.obj, p.original, file = predictionObjectName)
  
}

# MAKE HEX ----------------------------------------------------------------
makeHex <- function(data) {
  
  # Cast species dataset to a spatial object
  # Not casting/combining causes sf functions to stall
  data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>%
    st_cast(., to = "MULTIPOINT") %>%
    st_combine(.) %>%
    st_set_crs(4326) %>% # WGS 84 - geographic reference system
    st_transform(5070) # Albers equal area conic (aka conic albers) - projected
  
  # Create study area, with an offshore buffer of 10 km (equal to grid)
  studyArea <-
    st_convex_hull(data.sf) %>%
    st_make_valid() %>%
    st_difference(., northAmerica_forHex) %>%
    st_buffer(., dist = -10000) # buffer distance in m
  
  # Create grid that covers just the study area. Cell size in EPSG 5070 is in meters. 50000 cell size = 50 km (same size as prediction grid)
  hex <- st_make_grid(studyArea, crs = st_crs(northAmerica_forHex), cellsize = 50000, square = FALSE, what = "polygons") %>%
    st_make_valid() %>%
    st_intersection(., studyArea) %>% 
    st_difference(., gulfOfCalifornia) %>% 
    st_as_sf() %>% 
    mutate(., grid_id = seq(1:nrow(.)))
  
  return(hex)
}

# ROTATE MONTHS-----------------------------------------------
rotateMonths <- function(data, abundanceTerm, months, direction) {
  
  # Find month with lowest empirical abundance
  summary.month <- data %>% ungroup() %>%  group_by_at("month") %>% 
    summarize(mean_abund = mean(!!sym(abundanceTerm)))
  lowestMonth <- which(summary.month$mean_abund == min(summary.month$mean_abund))
  
  # Ensure rotation_amount is within the range 1 to 12
  rotation_amount <- lowestMonth %% 12
  
  # Perform the rotation
  if (direction == "rotate") {
    rotated_months <- (((months + 1 - rotation_amount) %% 12) - 12) %% 12
  } else if (direction == "reverse") {
    rotated_months <- (((months - 1 + rotation_amount) %% 12) + 12) %% 12
  }
  
  # Replace 0s with 12
  rotated_months[rotated_months == 0] <- 12
  return(rotated_months)
  
}


# CENTER OF GRAVITY -------------------------------------------------------

calculateCOG <- function(p, species) {
  
  # Create dataframe to store results
  years = unique(p$year)
  cog.year = data.frame(matrix(ncol = 0, nrow = length(years)), year = NA, cog.lat = NA, cog.month = NA, species = species)
  
  # Loop through each year and calculate COG for that year
  for (i in 1:length(years)) {
    
    # Subset to each year
    yearSubset = subset(p, year == years[i])
    
    # Calculate center of gravities
    cog.year$cog.lat[i] = sum(yearSubset$est_retransform_weighted*yearSubset$latitude)/sum(yearSubset$est_retransform_weighted)
    cog.year$cog.long[i] = sum(yearSubset$est_retransform_weighted*yearSubset$longitude)/sum(yearSubset$est_retransform_weighted)
    cog.year$cog.month_rotated[i] = sum(yearSubset$est_retransform_weighted*yearSubset$month_rotated)/sum(yearSubset$est_retransform_weighted)
    
    # Include sum of larval abundance
    cog.year$sum_est[i] = sum(yearSubset$est_retransform_weighted)
    
    # Include year
    cog.year$year[i] = years[i]
  }
  
  # Reverse rotation
  cog.year$cog.month = rotateMonths(data = p, abundanceTerm = "est_retransform_weighted", months = cog.year$cog.month_rotated, direction = "reverse")
  
  return(cog.year)
  
}
#---------------------------------------------------------------------------#
gridCOG <- function(species, commonName, modelType) {
  
  # Create dataframe to store results
  gridids = unique(p$gridid)
  cog.grid = data.frame(matrix(ncol = 0, nrow = length(gridids)), gridid = NA, cog.month = NA, species = species)
  
  # For each grid cell, calculate COG
  for (i in 1:length(gridids)) {
    gridSubset = subset(p, gridid == gridids[i])
    
    # Calculate COG on rotated month
    cog.grid$cog.month_rotated[i] = sum(gridSubset$est_retransform_weighted*gridSubset$month_rotated)/sum(gridSubset$est_retransform_weighted)
    
    # Include sums
    cog.grid$sum_est[i] = sum(gridSubset$est_retransform_weighted)
    cog.grid$mean_est[i] = mean(gridSubset$est_retransform_weighted)
    
    # Include gridid and species info
    cog.grid$gridid[i] = gridids[i]
  }
  
  # Reverse months
  cog.grid$cog.month = rotateMonths(data = p, abundanceTerm = "est_retransform_weighted", months = cog.grid$cog.month_rotated, direction = "reverse")
  
  # Add lat and long
  cog.grid <- merge(cog.grid, grid.df[c("gridid", "latitude", "longitude")])
  
  # Make sf object
  cog.grid.sf <- cog.grid %>% 
    st_as_sf(., coords = c("longitude", "latitude")) %>% 
    st_set_crs(4326) %>% # WGS 84 - geographic reference system
    st_transform(5070)
  
  species <- as.character(species)
  
  ## Create map of grid points colored by seasonal COG
  cog_grid_plot <- ggplot() +
    ggtitle(label = species) +
    geom_sf(data = northAmerica) +
    geom_sf(data = cog.grid.sf, aes(col = cog.month, size = mean_est)) +
    xlim(min(data$longitude)-1, max(data$longitude)+1) +
    ylim(min(data$latitude)-1, max(data$latitude)+1) +
    scale_color_gradientn(limits = c(0,12), "Seasonal COG \n(month)", colors = monthColors) +
    scale_size_continuous(name = "Estimated \nabundance [Log(N+1)]", limits = c(0,1.5), range = c(0.2, 5)) +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          axis.line = element_line(color = NA),
          axis.text.x = element_text(angle = 60,hjust=1)) +
    labs(x = "", y = "") +
    ggtitle(label = commonName, subtitle = bquote(~italic(.(species))))
  
  # Save individual map
  pdf(paste0("Figures/",species, "/CentralTendency.map.overall_", modelType, ".pdf"), width = 6, height = 6)
  print(cog_grid_plot)
  dev.off()
  
  # Center of gravity per year in relation to environmental variables
  cog.year <- calculateCOG(p, species)
  
  return(cog_grid_plot)
} 
#------------------------------------------------------------------------------#
cogMultispecies <- function(tradeoffs) {
  
  speciesList = tradeoffs$species
  
  # Set up storage dataframes
  cog.all = data.frame(year = integer(), cog.lat = double(), cog.long = double(), cog.month = double())
  cog.month.means = data.frame(cog.month.mean = NA, species = speciesList)
  lms <- data.frame(species = speciesList, y = "cog.lat", x = "cog.month")
  
  for (s in 1:length(speciesList)) {
    
    species = as.character(speciesList[s])
    
    # Load prediction object for each species
    predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", tradeoffs$short_formula[s], ".rdata")
    load(predictionObjectName)
    
    # Rotate months
    p$month_rotated = rotateMonths(p,abundanceTerm = "est_retransform_weighted",  p$month, direction = "rotate")
    
    # Calculate cog for each species based on rotated month
    cog <- calculateCOG(p, species)
    
    # Row-bind
    cog.all <- bind_rows(cog.all, cog)
    
    # Calculate mean COG for each month and rotate months back
    cog.month.means$cog.month.mean[s] = mean(cog$cog.month_rotated) %>% 
      rotateMonths(data = p, abundanceTerm = "est_retransform_weighted", months = ., direction = "reverse")
    
  }
  
  # Summarize results
  cog.all.summary <- cog.all %>% group_by_at("species") %>%
    summarize(cog.month.var = var(cog.month_rotated),
              cog.lat.mean = mean(cog.lat),
              cog.lat.var = var(cog.lat),
              cog.long.mean = mean(cog.long),
              cog.long.var = var(cog.long)) %>% 
    merge(cog.month.means)
  
  tradeoffs <- merge(tradeoffs, cog.all.summary)
  cog.all.tradeoffs <- merge(cog.all, tradeoffs, by = "species")
  
  return(list(tradeoffs, cog.all.tradeoffs))
  
}

# NICHE HYPERVOLUME ----------------------------------------------------------------------------------------------
hypervolumeMultiSpecies <- function(tradeoffs, recalculate) {
  
  speciesList = as.character(tradeoffs$species)
  modelTypes = tradeoffs$model
  lifeHistories = tradeoffs$lifeHistory
  
  colnames = c("species", "smith.lower", "smith", "smith.upper", "lifeHistory", "model")
  nicheBreadthCombined = data.frame(matrix(ncol = length(colnames), nrow = length(speciesList)))
  colnames(nicheBreadthCombined) = colnames
  
  # Set up progress bar
  pb <- txtProgressBar(min = 0, max = nrow(tradeoffs), char = "=", style = 3)
  
  for (i in 1:length(speciesLevels)) {
    
    # Calculate hypervolume
    shortFormula = tradeoffs$short_formula[i]
    modelType = modelTypes[i]
    species = speciesList[i]
    hypervolume_results = calculateHypervolume(species = species, modelType = modelType, shortFormula = shortFormula)
    
    # Bind niche breadth datasets together
    nicheBreadthCombined$species[i] = species
    nicheBreadthCombined$common_name[i] = tradeoffs$common_name[i]
    nicheBreadthCombined$smith.lower[i] = hypervolume_results[[1]]
    nicheBreadthCombined$smith[i] = hypervolume_results[[2]]
    nicheBreadthCombined$smith.upper[i] = hypervolume_results[[3]]
    nicheBreadthCombined$lifeHistory[i] = lifeHistories[i]
    nicheBreadthCombined$model[i] = modelType
    
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  
  nicheBreadthCombined$species = factor(nicheBreadthCombined$species, levels = speciesLevels)
  nicheBreadthCombined$common_name = factor(nicheBreadthCombined$common_name, levels = speciesInfo$common_name)
  nicheBreadthCombined$model = factor(nicheBreadthCombined$model, levels = c("base", "geo", "pheno", "both"))
  
  pdf(paste0("Figures/Fig#_NicheBreadth/AllLifeHistories_NicheBreadth.pdf"), width = 7, height = 6)
  print(ggplot(data = nicheBreadthCombined, mapping = aes(x = common_name, y = smith, ymin = smith.lower, ymax = smith.upper, color = model)) +
          geom_point(size = 3) +
          geom_linerange(linewidth = 1.2) +
          geom_vline(xintercept = 5.5, linewidth = 0.75, color = "gray80", linetype = 2) +
          geom_vline(xintercept = 11.5, linewidth = 0.75, color = "gray80", linetype = 2) +
          labs(x = "", y = "Smith's measure") +
          theme_classic(base_size = 16) +
          coord_cartesian(ylim = c(0,1)) +
          theme(panel.border = element_rect(fill = NA, color = "black"),
                axis.line = element_line(color = NA),
                legend.position.inside = c(0.13,0.2), legend.background = element_rect(fill = NA, color = NA), 
                axis.text.x = element_text(angle = 60,hjust=1))  +
          scale_color_manual("Model", values = modelColors, labels = c("base" ="Base", "geo" = "Geography", "pheno" = "Phenology", "both" = "Both")))
  dev.off()
  
  tradeoffs$smith = nicheBreadthCombined$smith
  
  return(tradeoffs)
}
#--------------#
calculateHypervolume <- function(species, modelType, shortFormula) {
  
  # Filenames
  fitName = paste0("Results/", species, "/Models/abundance_logN1_scaled_speciesRange_allPrograms_", modelType, ".rdata")
  predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", shortFormula, ".rdata")
  gridFilename = paste0("Analysis/PredictionGrids/", species, "_grid.rdata")
  
  # If the prediction object doesn't exist, create it
  if(!file.exists(predictionObjectName)) {
    createPredictionObjects(species = species, makeNewGrid = F, makeNewPrediction = T, fitName = fitName, predictionObjectName = predictionObjectName, gridFilename = gridFilename, stacked = T)
  } 
  
  load(predictionObjectName)
  
  # Summarize data by bins
  # Breaks: sst = 1 deg, ssh = 0.2 m, salinity = 0.5 ppt, bottom depth = 50, 200
  # Spiciness not included due to collinearity
  p_niche <- p %>% ungroup() %>%
    mutate(., sst_binned = cut(sst_roms, breaks = seq(min(p$sst_roms), max(p$sst_roms), by = 1), dig.lab = 0, include.lowest = T)) %>% 
    mutate(ssh_binned = cut(ssh_roms, breaks = seq(min(p$ssh_roms), max(p$ssh_roms), by = 0.2), include.lowest = T, dig.lab = 1)) %>% 
    mutate(salinity_binned = cut(salinity_roms, breaks = seq(min(p$salinity_roms), max(p$salinity_roms), by = 0.5), dig.lab = 0, include.lowest = T)) %>% 
    mutate(bottom_depth_binned = cut(bottom_depth,  breaks = c(max(p$bottom_depth), -50, -100, -200, min(p$bottom_depth)), dig.lab = 0, include.lowest = T)) %>%
    group_by_at(c("month", "gridid", "sst_binned", "salinity_binned", "ssh_binned", "bottom_depth_binned")) %>%
    summarize(sum_est = sum(est_retransform_weighted), n = n()) %>% 
    ungroup() %>% mutate(est_prop = sum_est/sum(sum_est), n_prop = n/sum(n))
  
  if (sum(p_niche$sum_est) == 0) {
    stop("Sum of estimates is zero; cannot compute hypervolume bounds.")
  }
  
  # Calculate Smith's hypervolume
  hypervolume = sum(sqrt(p_niche$est_prop*p_niche$n_prop))
  
  # Calculate upper and lower bounds
  hypervolume_upper = sin(asin(hypervolume)+1.96/(2*sqrt(sum(p_niche$sum_est))))
  hypervolume_lower = sin(asin(hypervolume)-1.96/(2*sqrt(sum(p_niche$sum_est))))
  
  return(list(hypervolume_lower, hypervolume, hypervolume_upper))
}


# SPAWNING LATITUDINAL RANGES -------------------
spawningRange <- function(species = species, commonName = commonName, shortFormula = shortFormula, modelType = modelType) {
  
  # Filenames
  fitName = paste0("Results/", species, "/Models/abundance_logN1_scaled_speciesRange_", modelType, ".rdata")
  predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", shortFormula, ".rdata")
  gridFilename = paste0("Analysis/PredictionGrids/", species, "_grid.rdata")
  
  # Load objects
  load(fitName) # data and model
  load(predictionObjectName) # prediction objects
  load(gridFilename) # grid.df
  
  years <- unique(p$year)
  
  # Create storage dataframe
  spawningRangeSummary <- data.frame(
    species = species,
    year = years, 
    trailingLat = NA, 
    leadingLat = NA, 
    spawningRange = NA, 
    trailingLon = NA, 
    leadingLon = NA, 
    spawningRangeLon = NA
  )
  
  # Function to calculate leading/trailing thresholds
  getPercentiles <- function(data, coord_col, low, high) {
    data <- data %>%
      group_by(.data[[coord_col]]) %>%
      summarize(meanAbundance = mean(est_retransform_weighted, na.rm = TRUE)) %>%
      ungroup() %>%
      arrange(coord_col) %>%
      mutate(cum_sum = cumsum(meanAbundance))
    
    trailingThresh <- max(data$cum_sum, na.rm = TRUE) * low
    leadingThresh <- max(data$cum_sum, na.rm = TRUE) * high
    
    # Find longitude/longitude that corresponds to the cumulative abundance threshold
    trailingLat <- max(data[[coord_col]][data$cum_sum < trailingThresh], na.rm = TRUE)
    leadingLat <- max(data[[coord_col]][data$cum_sum < leadingThresh], na.rm = TRUE)
    
    return(c(trailingLat, leadingLat))
  }
  
  # Loop over years
  for (i in 1:length(years)) {
    
    yearSubset <- subset(p, year == years[i])
    
    leading = 0.85
    trailing = 0.15
    
    # Latitude
    lat_percentiles <- getPercentiles(data = yearSubset, coord_col = "latitude", low = trailing, high = leading)
    spawningRangeSummary$trailingLat[i] <- lat_percentiles[1]
    spawningRangeSummary$leadingLat[i] <- lat_percentiles[2]
    spawningRangeSummary$spawningRange[i] <- lat_percentiles[2] - lat_percentiles[1]
    
    # Longitude
    lon_percentiles <- getPercentiles(yearSubset, coord_col = "longitude", low = trailing, high = leading)
    spawningRangeSummary$trailingLon[i] <- lon_percentiles[1]
    spawningRangeSummary$leadingLon[i] <- lon_percentiles[2]
    spawningRangeSummary$spawningRangeLon[i] <- lon_percentiles[2] - lon_percentiles[1]
  }
  
  # Create plot for each species
  plot <- ggplot(spawningRangeSummary) +
    geom_errorbar(mapping = aes(x = years, ymin = trailingLat, ymax = leadingLat)) +
    #geom_point(mapping = aes(x = year, y = spawningRange)) +
    labs(y = "Latitude", x = "Year") +
    ggtitle(label = commonName, subtitle = bquote(~italic(.(species))))
  
  return(list(plot, spawningRangeSummary))
  
}


# CUMULATIVE CATCH ------------------------------------------
cumulativeCatch <- function(species) {
  
  # Summarize modeled larval abundance by month, year, and timeblock
  predictedCatch.time <- p %>% group_by_at(., c("month", "year", "timeblock")) %>%
    dplyr::summarize(., mean_est_retransform = mean(est_retransform_weighted)) %>% # summarize across grid cells
    mutate(., est_x_month = mean_est_retransform * month) %>%
    mutate(., timeblock = factor(timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")))
  
  # Rotate months
  predictedCatch.time$rotated_month <- rotateMonths(data = predictedCatch.time, abundanceTerm = "mean_est_retransform", months = predictedCatch.time$month, direction = "rotate")
  
  # Calculate cumulative sum of mean abundance
  predictedCatch.time.cumulative <- predictedCatch.time %>%
    ungroup() %>% group_by_at(c("year", "rotated_month")) %>%
    summarize(month = max(month), mean_est_retransform = mean(mean_est_retransform)) %>%
    mutate(cum_catch = cumsum(mean_est_retransform))
  
  # Find 15th and 85th percentiles of cumulative abundance
  cumulativeAbundance <- predictedCatch.time.cumulative %>% 
    summarize(fifteen = max(cum_catch)*0.15, fifty = max(cum_catch)*0.5, eightyfive = max(cum_catch)*0.85) %>% 
    mutate(fifteen.month.rot = NA, fifty.month.rot = NA, eightyfive.month.rot = NA, fifteen.month = NA, fifty.month = NA, eightyfive.month = NA)
  
  timesteps = unique(cumulativeAbundance$year)
  
  for (i in 1:length(timesteps)) {
    
    # Subset timesteps
    timeSubset = subset(predictedCatch.time.cumulative, year == timesteps[i])
    
    # Find month at which 15%, 50%, and 85% are passed (ROTATED YEAR)
    fifteen.month.rot = max(timeSubset$rotated_month[timeSubset$cum_catch < cumulativeAbundance$fifteen[i]], na.rm = T)
    if(fifteen.month.rot == -Inf) {
      fifteen.month.rot = 1
    } 
    
    slope = (timeSubset$cum_catch[fifteen.month.rot+1]-timeSubset$cum_catch[fifteen.month.rot])/(1) # differences between months are always 1
    b = timeSubset$cum_catch[fifteen.month.rot] - slope * fifteen.month.rot #y-mx
    cumulativeAbundance$fifteen.month.rot[i] = (cumulativeAbundance$fifteen[i] - b) / slope #y-b/m
    
    fifty.month.rot = max(timeSubset$rotated_month[timeSubset$cum_catch < cumulativeAbundance$fifty[i]])
    slope = (timeSubset$cum_catch[fifty.month.rot+1]-timeSubset$cum_catch[fifty.month.rot])/(1) # differences between months are always 1
    b = timeSubset$cum_catch[fifty.month.rot] - slope * fifty.month.rot #y-mx
    cumulativeAbundance$fifty.month.rot[i] = (cumulativeAbundance$fifty[i] - b) / slope #y-b/m
    
    eightyfive.month.rot = max(timeSubset$rotated_month[timeSubset$cum_catch < cumulativeAbundance$eightyfive[i]])
    slope = (timeSubset$cum_catch[eightyfive.month.rot+1]-timeSubset$cum_catch[eightyfive.month.rot])/(1) # differences between months are always 1
    b = timeSubset$cum_catch[eightyfive.month.rot] - slope * eightyfive.month.rot #y-mx
    cumulativeAbundance$eightyfive.month.rot[i] = (cumulativeAbundance$eightyfive[i] - b) / slope #y-b/m
  }
  
  # Reverse rotation on months
  cumulativeAbundance <- cumulativeAbundance %>% 
    mutate(fifteen.month = rotateMonths(data = p, abundanceTerm = "est_retransform_weighted", months = fifteen.month.rot, direction = "reverse")) %>% 
    mutate(fifty.month = rotateMonths(data = p, abundanceTerm = "est_retransform_weighted", months = fifty.month.rot, direction = "reverse")) %>% 
    mutate(eightyfive.month = rotateMonths(data = p,abundanceTerm = "est_retransform_weighted", months = eightyfive.month.rot, direction = "reverse"))
  
  # Correct for limits that cross the month boundary
  for (i in 1:nrow(cumulativeAbundance)) {
    if(cumulativeAbundance$fifteen.month[i] > cumulativeAbundance$fifty.month[i]) {
      cumulativeAbundance$fifteen.month[i] = cumulativeAbundance$fifteen.month[i] - 12
    }
    if(cumulativeAbundance$eightyfive.month[i] < cumulativeAbundance$fifty.month[i]) {
      cumulativeAbundance$eightyfive.month[i] = cumulativeAbundance$eightyfive.month[i] + 12
    }
  }
  
  # Plot cumulative abundance
  plot <- ggplot(cumulativeAbundance) +
    geom_errorbar(mapping = aes(x = year, ymin = eightyfive.month, ymax = fifteen.month))+
    geom_point(mapping = aes(x = year, y = fifty.month)) +
    labs(y = "Month", x = "Year") +
    ggtitle(bquote(~italic(.(species)))) +
    scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = ""))+
    scale_x_continuous(n.breaks = length(timesteps)/2, labels = scales::number_format(accuracy = 1, big.mark = ""))
  
  pdf(file = paste0("Figures/", species, "/", "Cumulative_Thresholds_.pdf"), width = 5, height = 5)
  print(plot)
  dev.off()
  
  write_xlsx(x = cumulativeAbundance,  path = paste0("Results/", species, "/", species, "_cumulativeAbundance.xlsx"), format_headers = F)
  
  return(p)
  
}

# CONDITIONAL EFFECTS ----------------------------------------------------
conditionalEffects <- function(species, data, grid.df, fit, shortFormula, prediction_grid_roms, newPrediction) {
  
  # Get formula
  formula = str_split(string = shortFormula, pattern = "\\+")[[1]]
  
  # Create storage list for plots
  plotlist <- list()
  
  scale_names <- c(
    sst_scaled = "sst_roms",
    ssh_scaled = "ssh_roms",
    bottom_depth_scaled = "bottom_depth",
    spice_scaled = "spice_roms",
    salinity_scaled = "salinity_roms"
  )
  
  labels <- c(
    sst_scaled = "Sea surface temperature (\u00b0)",
    ssh_scaled = "Sea surface height (m)",
    bottom_depth_scaled = "Bottom depth (m)",
    spice_scaled = "Spiciness",
    salinity_scaled = "Salinity (ppt)"
  )
  
  # Calculate center/scale values from empirical data
  scale_params <- lapply(scale_names, function(var) {
    vals <- fit$data[[var]]
    list(center = mean(vals, na.rm = TRUE), scale = sd(vals, na.rm = TRUE))
  })
  
  plan(multisession)
  
  plotlist <- future_map(1:length(formula), function(i) {
    
    # Get all covariates in formula
    covariates <- c(
      formula[i], # Set covariate of interest to first slot
      formula[-i]  # Move other components to the end
    )
    
    year <- seq(from = min(data$year), to = max(data$year), by = 1)
    month <- seq(from = 1, to = 12, by = 1)
    timeblock = factor(x = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019"))
    
    # Get column name of covariate
    covariate = paste0(covariates[1], "_scaled") # The first entry is always the one of interest
    
    covariate_unscaled <- scale_names[[covariate]]
    covariate_label <- labels[[covariate]]
    
    if(newPrediction == T) {
      # If not already created, run conditional effects, and save output (takes a while)
      # Margin = "mean_reference" is What is the predicted/expected larval abundance at meaningful values or levels of for a 'typical' observation?
      pred_response <- predict_response(model = fit, terms = paste0(covariate, " [all]"), type = "fixed", margin = "mean_reference")
      save(pred_response, file = paste0("Results/", species, "/", covariate, "_posterior_response.Rdata"))
    } else {
      load(pred_response, file = paste0("Results/", species, "/", covariate, "_posterior_response.Rdata"))
    }
    
    # Reverse scaling for plotting purposes
    center <- scale_params[[covariate]]$center
    scale <- scale_params[[covariate]]$scale
    pred_response$x_unscaled <- (pred_response$x * scale) + center
    
    # Plot
    plot <- ggplot(pred_response, aes(x_unscaled, predicted)) +
      geom_line(linewidth = 0.75) +
      theme_classic() +
      labs(x = covariate_label, y = "") +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)
    
    # If the first plot, add y axis, otherwise leave axis blank
    if (i == 1) {
      plot = plot + labs(y = "Abundance [log(N+1)]")
    }
    
    # Return to top level
    return(plot)
    
  })
  
  # Plot all covariates together
  conditionalEffectsPlot <- ggarrange(plotlist = plotlist, nrow = 1)
  
  pdf(paste0("Figures/", species, "/EnviroPosteriors.pdf"), width = 7, height = 4)
  print(conditionalEffectsPlot)
  dev.off()
  
  # Return plot
  return(conditionalEffectsPlot)
  
}

# TRADEOFFS (GENERAL) -----------------------------------------------------------------
tradeoffsBarPlots <- function(tradeoffs) {
  
  tradeoffs.summary <- group_by_at(tradeoffs, c("model", "lifeHistory")) %>% summarize(n = n())
  tradeoffs.summary$model <- factor(tradeoffs.summary$model, levels = c("base", "geo", "pheno", "both"))
  
  envCovariates <- read_xlsx("Results/Environmental_covariates_by_life_history.xlsx")
  
  # Fig. 2 ------ 
  # Environmental covariate bar plot
  p1 <- print(ggplot(envCovariates, mapping = aes(x = Covariate, y = N, fill = `Life history`)) +
                geom_col(position = position_dodge2(preserve = c("single")), width = .8) +
                theme_classic(base_size = 16) +
                scale_y_continuous(minor_breaks = waiver()) +
                scale_fill_manual("Adult habitat", values = lifeHistoryColors) +
                theme(axis.text.x= element_text(angle = 45,hjust=1), legend.position = "none",
                      axis.line = element_line(color = NA),
                      panel.border = element_rect(fill = NA, color = "black")) +
                coord_cartesian(ylim = c(0.5, 7), xlim = c(1,5)) +
                labs(y = "N species", x= ""))
  
  # Tradeoff bar plot
  p2 <- print(ggplot(tradeoffs.summary) +
                geom_col(aes(x = model, y = n, fill = lifeHistory), position = position_dodge2(preserve = c("single"))) +
                labs(x = "", y = "") +
                theme_classic(base_size = 16) +
                scale_x_discrete(labels = c("base" = "Base","geo" =  "Shifting geography","pheno" = "Shifting phenology", "both" = "Both")) +
                coord_cartesian(ylim = c(0.5, 7), xlim = c(1.15,3)) +
                theme(axis.text.x= element_text(angle = 45,hjust=1), legend.position.inside = c(0.05,0.76),
                      legend.background = element_rect(fill = NA, color = NA),
                      axis.line = element_line(color = NA),
                      legend.justification = "left",
                      panel.border = element_rect(fill = NA, color = "black")) +
                scale_fill_manual("Adult habitat", values = lifeHistoryColors))
  
  
  pdf(file = "Figures/Fig2_TradeoffSummary/Fig2_Tradeoffs and EnviroCovariates.pdf", width = 9, height = 5)
  print(ggarrange(p1, p2, align = "hv", widths = c(0.58, 0.42), labels = c("A", "B"), font.label = list(size = 20)))
  dev.off()
  
  jpeg("Figures/Fig2_TradeoffSummary/Fig2_Tradeoffs and EnviroCovariates.jpeg", width = 9, height = 5, res = 500, units = "in")
  print(ggarrange(p1, p2, align = "hv", widths = c(0.58, 0.42), labels = c("A", "B"), font.label = list(size = 20)))
  dev.off()
}

# TRADEOFFS (DETAILED) ------------------
tradeoffsDetailed <- function() {
  
  # Load full tradeoff model comparison sheet
  tradeoffs_summary <- read_xlsx("Results/02_Tradeoff Comparison/Tradeoffs_summary.xlsx")
  
  # Factor common name in reverse (necessary for plotting)
  tradeoffs_summary$common_name = factor(
    tradeoffs_summary$common_name,
    levels = rev(speciesInfo$common_name))
  
  # Factor model types
  tradeoffs_summary$model = factor(tradeoffs_summary$model, levels = c("base", "geo", "pheno", "both"))
  
  # Fig. 3 -------
  cairo_pdf(paste0("Figures/Fig3_TradeoffDetailed/ModelComparison_weight.pdf"), width = 6, height = 6)
  print(ggplot(tradeoffs_summary, mapping = aes(x = model, y = common_name, fill = weight)) +
          geom_tile() +
          coord_cartesian(expand = F) +
          theme_classic(base_size = 16) +
          geom_hline(linewidth = 0.5, color = "gray40", yintercept = c(seq(from = 1.5, to = 16.5, by = 1))) +
          labs(x = "") +
          scale_x_discrete(labels = c("base" = "Base", "geo" = "Geography", "pheno" = "Phenology", "both" = "Both")) +
          scale_fill_gradientn("Model\nweight", colors = viridis::mako(n = 5, direction = -1, end = 1, begin = 0)) +
          theme(panel.border = element_rect(fill = NA, linewidth = 0.5),
                axis.line = element_line(color = NA),
                axis.text.x = element_text(angle = 45, hjust=1), plot.margin = margin(l = 2, t = 0.5, unit = "cm")) +
          labs(y = ""))
  dev.off()
  
}
# EMPIRCAL VS PREDICTED --------------------------------------------------------
empiricalVsPredicted <- function(species, p.original, commonName) {
  
  # Create sf object out of prediction on original stations
  predictedCatch.grid.original <- p.original %>%
    st_as_sf(., coords = c("longitude", "latitude")) %>% 
    st_set_crs(4326) %>% # WGS 84 - geographic reference system
    st_transform(5070)
  
  # Join with hex (for plotting purposes)
  hex.predict <- st_join(hex, predictedCatch.grid.original) %>% 
    group_by_at(., c("grid_id")) %>%
    dplyr::summarize(., mean_est = mean(est_retransform_weighted)) %>% 
    st_transform(5070)
  
  # Create sf object out of empirical data
  data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>%
    st_set_crs(4326) %>%  # WGS 84 - geographic reference system
    st_transform(5070) 
  
  # Join with hex (for plotting purposes)
  hex.empirical <- st_join(hex, data.sf) %>% 
    group_by_at(., c("grid_id")) %>%
    dplyr::summarize(., mean_abundance_logN1_scaled = mean(abundance_logN1_scaled)) %>% 
    st_transform(5070)
  
  # Run Pearson correlation test
  pearson <- cor.test(p.original$est_retransform_weighted, p.original$abundance_logN1_scaled, method = "pearson")
  pvalue = pearson$p.value
  if(pvalue == 0) pvalue = "<0.001"
  
  nullStationColor = "#d0e8ed"
  
  # Empirical map
  p1 <- ggplot() +
    geom_sf(data = northAmerica, fill = "gray70") +
    geom_sf(data = subset(hex.empirical, !is.na(mean_abundance_logN1_scaled)), color = NA, mapping = aes(fill = mean_abundance_logN1_scaled)) + 
    geom_sf(data = subset(hex.empirical, mean_abundance_logN1_scaled == 0), color = NA, fill = nullStationColor) + 
    xlim(min(data$longitude), max(data$longitude)) +
    ylim(min(data$latitude)-1, max(data$latitude)) +
    scale_fill_gradient("Average abundance\n[log(N+1)]", limits = c(0, max(hex.empirical$mean_abundance_logN1_scaled, na.rm = T)), na.value = "gray90", low = "lightgoldenrod2", high =  "firebrick") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          legend.text = element_text(size = 10), panel.border = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(angle = 60,hjust=1),
          axis.line = element_line(color = NA),
          # legend.position = "right", legend.key.width = unit(0.8, units = "cm"),
          legend.key.spacing = unit(0.5, units = "cm")) +
    labs(x = "", y = "") 
  
  # Prediction map
  p2 <- ggplot() +
    geom_sf(data = northAmerica, fill = "gray70") +
    geom_sf(data = subset(hex.predict, !is.na(mean_est)), mapping = aes(fill = mean_est), color = NA,) + 
    geom_sf(data = subset(hex.empirical, mean_abundance_logN1_scaled < 0.0001), color = NA, fill = nullStationColor) + 
    xlim(min(data$longitude), max(data$longitude)) +
    ylim(min(data$latitude)-1, max(data$latitude)) +
    scale_fill_gradient("Average\nabundance\n[log(N+1)]", limits = c(0, max(hex.empirical$mean_abundance_logN1_scaled, na.rm = T)), na.value = "gray90", low = "lightgoldenrod2", high =  "firebrick") +
    theme_classic(base_size = 14) +
    theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9),
          legend.background = element_rect(fill = "white", color = NA),
          axis.text.y = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(angle = 60,hjust=1),
          axis.line = element_line(color = NA),
          plot.margin = unit(c(0,0,0,0), units = "cm"), 
          legend.position = "inside",
          legend.position.inside = c(0.77,0.8),
          legend.direction = "vertical") +
    labs(x = "", y = "") 
  
  # Run a GAM to examine nonlinear trends between empirical and modeled data
  gamModel <- mgcv::gam(p.original$est_retransform_weighted ~ s(p.original$abundance_logN1_scaled))
  gamPlotResult <- mgcViz::getViz(gamModel)
  gamPlot <- plot(mgcViz::sm(gamPlotResult, 1) )
  
  
  # Fig. 4 --------
  gamPlot <- gamPlot +
    l_fitLine(colour = "gray30", linewidth = 1) +
    l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_points(colour = "gray30", shape = 16, size = 2, alpha = 0.4) +
    # l_ciLine(colour = "cadetblue", level = 0.95, linetype = 2, linewidth = 1) + 
    l_ciPoly(fill = "lightblue2", level = 0.95, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, lty = 1, linewidth = 1, color = "black") +
    theme_classic(base_size = 14) +
    annotate(hjust = 1, size  = 4.5, geom = "text", label = paste("Pearson coeff. =", round(pearson$estimate, digits = 3), "\np =", pvalue),
             x = max(p.original$abundance_logN1_scaled),
             y = max(p.original$est_retransform_weighted)-max(p.original$est_retransform_weighted)*0.2, color="black") +
    theme(panel.border = element_rect(fill = NA, color = "black"), axis.line = element_line(color = NA)) +
    labs(y = "s(Empirical abundance [log(N+1)])", x = "Empirical abundance [log(N+1)]\n")
  
    # Plot for the supplemental figure (all species together)
  gamPlot_group <- gamPlot  +
    ggtitle(label = commonName, subtitle = bquote(~italic(.(species))))
  
  # Individual species
  pdf(paste0("Figures/", species,"/", species, "_GAM_Empirical_vs_Predicted.pdf"), width = 4, height = 4)
  print(gamPlot_group$ggObj)
  dev.off()
  
  maps <- ggarrange(p1+ggtitle("A. Empirical"),
                    p2+ggtitle("B. Predicted"),
                    nrow = 1, common.legend = F, align = "hv")
  
  
  fig <- ggarrange(maps, gamPlot$ggObj+ggtitle("C. GAM"), widths = c(1,2), nrow = 1)
  
  cairo_pdf(paste0("Figures/", species,"/", species, "_Fig4.pdf"), width = 10, height = 6)
  print(fig)
  dev.off()
  
  # Plots for the supplemental figure (all species together)
  p1_group = p1+annotate(geom = "text", size = 5, label = "Empirical", x = min(data$longitude), y = min(data$latitude)+0.5, hjust = 0)   +
    ggtitle(label = commonName, subtitle = bquote(~italic(.(species))))
  p2_group = p2+annotate(geom = "text", size = 5, label = "Predicted", x = min(data$longitude), y = min(data$latitude)+0.5, hjust = 0)
  maps_group <- ggarrange(p1_group, p2_group, nrow = 1, common.legend = T, align = "hv", legend = "bottom", legend.grob = get_legend(p1))
  
  return(list(gamPlot_group$ggObj, maps_group))
}
# TRENDS WITH YEAR ------------------------------------------------
trendsWithYear <- function(metrics.allYears) {
  
  vars = list(
    list("cog.month_rotated", "Seasonal COG (month)"),
    list("cog.lat", "Latitudinal COG (\u00b0)"),
    list("cog.long","Longitudinal COG (\u00b0)"),
    list("trailingLat", "Trailing latitude (\u00b0)"),
    list("leadingLat", "Leading latitude (\u00b0)")
  )
  
  varLabels = c("cog.month_rotated" = "Seasonal COG (month)",
                "cog.lat" = "Latitudinal COG (\u00b0)",
                "cog.long" = "Longitudinal COG (\u00b0)",
                "trailingLat" = "Trailing latitude (\u00b0)",
                "leadingLat" = "Leading latitude (\u00b0)")
  
  # Create storage dataframe for results
  colnames = c("species", "var", "r2adj", "f", "p", "df", "slope")
  lmSummary <- as.data.frame(matrix(nrow = length(speciesLevels)*length(vars), ncol = length(colnames)))
  colnames(lmSummary) = colnames
  
  index = 1
  plotlist = list()
  
  # Run linear regressions
  for (i in 1:length(speciesLevels)) {
    for (j in 1:length(vars)) {
      
      # Subset to individual species
      subset <- subset(metrics.allYears, species == speciesLevels[i])
      
      # Construct formula
      formula = as.formula(paste0(vars[[j]][1], "~", "year"))
      
      lm <- lm(data = subset, formula) %>% summary()
      lmSummary$species[index] = speciesLevels[i]
      lmSummary$var[index] = vars[[j]][1]
      lmSummary$varname[index] = vars[[j]][2]
      lmSummary$r2adj[index] = lm$adj.r.squared
      lmSummary$f[index] = lm$fstatistic[1]
      lmSummary$p[index] = lm$coefficients[8]
      lmSummary$df[index] = list(lm$df)
      lmSummary$slope[index] = lm$coefficients[2]
      
      index = index + 1
    }
  }
  
  # Correct p values and add useful columns to results dataframe
  lmSummary <- lmSummary %>% 
    mutate(p.adjust = p.adjust(p, method = "BH"), .before = p) %>% 
    mutate(Significance = case_when(
      p.adjust < alpha ~ "*",
      TRUE ~ "")) %>% 
    merge(speciesInfo[c("species", "phylogenetic_order")]) %>% 
    merge(tradeoffs[c("species", "cog.month.mean")])
  
  # Create supplemental table
  supplTable <- lmSummary
  supplTable <- supplTable %>% 
    arrange(phylogenetic_order, p.adjust) %>% 
    select(species, varname, r2adj, f, slope, p.adjust, Significance) %>% 
    mutate(r2adj = round(r2adj, digits = 2), slope = round(slope, digits = 3)) %>% 
    rename("Species" = species,
           "Tradeoff metric" = varname,
           "R2 adj." = r2adj,
           "F" = f,
           "Slope" = slope,
           "Adj. p" = p.adjust,
    ) 
  
  write_xlsx(supplTable, "Manuscript/Suppl_Table4_Changes over time results.xlsx", format_headers = F)
  
  # Extract significant regressions
  lmSummarySignificant <- subset(lmSummary, p.adjust < alpha & p.adjust != "NaN") %>% 
    merge(., speciesInfo[c("species", "lifeHistory")]) %>% 
    mutate(var = factor(var, levels = names(varLabels)))
  
  # Test for trends between slope and seasonal COG
  # "Are species that spawn in different seasons more likely to advance or delay their phenology?"
  lmSummarySeasonal <- subset(lmSummarySignificant, var == "cog.month_rotated")
  lm(lmSummarySeasonal$slope ~ lmSummarySeasonal$cog.month.mean) %>% summary() %>% print()
  
  # Test for trends between slope and latitudinal COG
  # "Are species located at different latitudes more likely to advance or delay their phenology?"
  lmSummarySeasonal <- subset(lmSummarySignificant, var == "cog.lat")
  lm(lmSummarySeasonal$slope ~ lmSummarySeasonal$cog.month.mean) %>% summary() %>% print()
  
  # Summarize number of significant trends per variable
  lmSummarySignificant_vars <- lmSummarySignificant %>% 
    group_by(varname) %>% 
    summarize(n = n(), meanSlope = mean(slope),
              minSlope = min(slope), maxSlope = max(slope),
              sdSlope = sd(slope)) %>% as.data.frame() %>% 
    mutate(varname = factor(varname, levels = c("Seasonal COG (month)",
                                                "Latitudinal COG (\u00b0)",
                                                "Longitudinal COG (\u00b0)",
                                                "Trailing latitude (\u00b0)",
                                                "Leading latitude (\u00b0)")))
  
  # Fig. 5 ------
  # Plot boxplots of slopes for each tradeoff metric
  pdf("Figures/Fig6_TradeoffMetricTrends/Fig6_TradeoffMetricTrends_box_point.pdf", width = 6, height = 5)
  print(ggplot(lmSummarySignificant) +
          geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), color = "gray85") +
          geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
          geom_boxplot(mapping = aes(x = var, y = slope, fill = lifeHistory, color = lifeHistory), outliers = F,  alpha = 0.5) +
          scale_fill_manual("Adult life history", values = lifeHistoryColors) +
          scale_color_manual("Adult life history", values = darken(lifeHistoryColors, amount = 0.3)) +
          new_scale_color() +
          geom_point(mapping = aes(x = var, y = slope, color = lifeHistory), shape = 16, alpha = 1, position = position_dodge(width = 0.75)) +
          scale_color_manual("Adult life history", values = lifeHistoryColors) +
          scale_y_continuous(limits = c(-0.25, 0.33), n.breaks = 10) +
          scale_x_discrete(labels = varLabels) +
          theme(panel.background = element_rect(fill = NA, color = "black"),
                axis.line = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(y = "Slope", x = ""))
  dev.off()
  
  # Test if slopes are significantly different between tradeoff metrics (ANOVA)
  aov <- aov(lmSummarySignificant$slope ~ lmSummarySignificant$var)
  TukeyHSD(x = aov)
  
  write_xlsx(lmSummary, path = "Results/03b_COG Linear Trends with Year/Regressions_Metrics vs year.xlsx", format_headers = F)
  write_xlsx(lmSummarySignificant_vars, path = "Results/03b_COG Linear Trends with Year/Significant regressions by variable.xlsx", format_headers = F)
  
  return(lmSummarySignificant)
}

# TRADEOFF METRIC TRENDS ----------------------------------------------
internalCOGTrends <- function(tradeoffs) {

  # Set dependent and independent variables
  vars = c("cog.month.var", "cog.month.mean",
           "cog.lat.var", "cog.long.var",
           "cog.lat.mean", "cog.long.mean",
           "trailing.var", "leading.var",
           "trailing.mean", "leading.mean")
  
  varNames = c("Seasonal COG\nvariance (month)",
               "Seasonal COG\nmean (month)",
               "Latitudinal COG\nvariance (\u00b0)",
               "Longitudinal COG\nvariance (\u00b0)",
               "Latitudinal COG\nmean (\u00b0)",
               "Longitudinal COG\nmean (\u00b0)",
               "Trailing edge\nvariance (\u00b0)",
               "Leading edge\nvariance (\u00b0)",
               "Trailing edge\nmean (\u00b0)",
               "Leading edge\nmean (\u00b0)")
  
  # Generate all unique pairs of covariates
  formulas <- t(combn(vars, 2, simplify = T)) %>% as.data.frame() %>% 
    mutate(formula = paste0(V1, " ~ ", V2))
  names <- t(combn(varNames, 2, simplify = T)) %>% as.data.frame() %>% 
    rename("varName1" = V1, "varName2" = V2)
  formulas = bind_cols(formulas,names)
  
  colnames = c("var1", "var2", "r2adj", "f", "p", "df", "slope")
  lmSummary <- as.data.frame(matrix(nrow = nrow(formulas), ncol = length(colnames)))
  colnames(lmSummary) = colnames
  
  index = 1
  plotlist = list()
  
  # Run formulas
  for (i in 1:nrow(formulas)) {
    lm <- lm(data = tradeoffs, as.formula(formulas$formula[i])) %>% summary()
    print(lm)
    lmSummary$var1[i] = formulas$V1[i]
    lmSummary$var2[i] = formulas$V2[i]
    lmSummary$varname1[i] = formulas$varName1[i]
    lmSummary$varname2[i] = formulas$varName2[i]
    lmSummary$r2adj[i] = lm$adj.r.squared
    lmSummary$f[i] = lm$fstatistic[1]
    lmSummary$p[i] = lm$coefficients[8]
    lmSummary$df[i] = list(lm$df)
    lmSummary$slope[i] = lm$coefficients[2]
    
    if (lm$coefficients[8] < alpha) {
      plotlist[[index]] <- ggplot(tradeoffs, mapping = aes(y = .data[[formulas$V1[i]]], x = .data[[formulas$V2[i]]])) +
        geom_smooth(method = "lm", color = "black") +
        geom_point(mapping = aes(color = lifeHistory)) +
        scale_color_manual("Life history", values = lifeHistoryColors) +
        theme_classic(base_size = 12) +
        labs(y = formulas$varName1[i], x = formulas$varName2[i]) +
        theme()
      index = index + 1
    }
    
  }
  
  cairo_pdf("Figures/Fig#_COG_Regressions/COG_Regressions.pdf", width = 10, height = 7)
  print(ggarrange(plotlist = plotlist, common.legend = T))
  dev.off()
  
  write_xlsx(x = lmSummary, path = "Results/03a_COG Internal Trends/COG_Season_Lat_Long_Linear Regressions.xlsx", format_headers = F)
  
  lmSummary <- lmSummary %>% 
    mutate(p.adjust = p.adjust(p, method = "bonferroni")) %>% 
    mutate(Significance = case_when(
      p.adjust < alpha ~ "*",
      TRUE ~ ""))
  
  # Create supplementary table
  supplTable <- lmSummary
  
  # supplTable$p <- ifelse(supplTable$p < 0.0001, "p < 0.0001", as.character(supplTable$p))
  
  supplTable <- supplTable %>% 
    select(varname1, varname2, r2adj, f, slope, p, p.adjust, Significance) %>% 
    arrange(p) %>% 
    rename("Metric 1" = varname1,
           "Metric 2" = varname2,
           "R2 adj." = r2adj,
           "F" = f,
           "Adj. p" = p.adjust,
           "Slope" = slope) %>% 
    select(-p)
  
  write_xlsx(x = supplTable, path = "Manuscript/Suppl_Table3_TradeoffRegressionResults.xlsx", format_headers = F)

  # Fig 6 --------------- 
  cairo_pdf("Figures/Fig5_COG/Lat_vs_Long_COG.pdf", width = 6, height = 6)
  print(ggplot(tradeoffs) +
          geom_errorbar(mapping = aes(x = cog.long.mean, ymin = cog.lat.mean-cog.lat.var, ymax = cog.lat.mean+cog.lat.var,color = cog.month.mean)) +
          geom_errorbarh(mapping = aes(y = cog.lat.mean, xmin = cog.long.mean-cog.long.var, xmax = cog.long.mean+cog.long.var, color = cog.month.mean)) +
          geom_point(mapping = aes(x = cog.long.mean, y = cog.lat.mean, color = cog.month.mean, size = cog.month.var)) +
          scale_color_gradientn("Seasonal COG\nmean (month)", limits = c(1,12),  colors = monthColors ) +
          scale_size_continuous("Seasonal COG\nvariance (month)", breaks = waiver()) +
          labs(x = "Longitudinal COG (\u00b0)", y = "Latitudinal COG (\u00b0)") +
          theme_classic(base_size = 14) +
          geom_text(mapping = aes(x = cog.long.mean, y = cog.lat.mean, label = species), hjust=0, vjust=0) +
          theme(panel.background = element_rect(fill = NA, color = "black"),
                plot.background = element_rect(fill = NA, color = "black"),
                legend.background = element_rect(fill = NA, color = NA), legend.axis.line = element_line(color = NA),
                axis.line = element_line(color = NA),
                legend.position.inside = c(0.22, 0.3), legend.spacing.y = unit(0, units = "cm"),
                legend.key.size = unit(0.5, "cm"),
                legend.text = element_text(size = 12), 
                legend.title = element_text(size = 12)))
  dev.off()
}

# SPATIAL/SPATIOTEMPORAL EFFECTS ---------------------------------------------------------------
effectsFigures <- function(species, modelType) {
  
  # Spatial effect per timeblock
  if("epsilon_st" %in% colnames(p)) { # geo/both models
    pdf(file = paste0("Figures/", species, "/Spatiotemporaleffects_map_", modelType,".pdf"), width = 7, height = 4)
    print(ggplot(northAmerica) +
            geom_sf() +
            facet_wrap(~timeblock, nrow = 1) +
            geom_point(data = predictedCatch.grid.timeblock, mapping = aes(x = (X*1000), y = (Y*1000), col = mean_epsilon_st),  size = 2) +
            scale_color_gradient2("Spatiotemporal effects", low = "firebrick", mid = "white", high =  "firebrick", midpoint = 0) +
            xlim(min(data$X)*1000-1000, max(data$X)*1000+1000) +
            ylim(min(data$Y)*1000-1000, max(data$Y)*1000+1000) +
            theme(axis.text.x = element_text(angle = 60,hjust=1))+
            theme_classic(base_size = 8) +
            ggtitle(bquote(~italic(.(species)))) +
            labs(x = "Longitude", y = "Latitude"))
    dev.off()
    
    # Random effects per latitude per timeblock
    pdf(file = paste0("Figures/", species, "/Randomeffects_", modelType,".pdf"), width = 5, height = 4)
    print(ggplot(predictedCatch.lat.timeblock, mapping = aes(x = latitude, y = mean_rf, fill = timeblock, col = timeblock)) +
            geom_point() +
            geom_smooth() +
            labs(x = "Latitude", y = "Mean random effect") +
            scale_color_brewer(palette = "Spectral") +
            theme_classic(base_size = 12) +
            ggtitle(bquote(~italic(.(species)))) +
            scale_fill_brewer(palette = "Spectral"))
    dev.off()
  }
  
}

# MAPPING EMPIRICAL DATA --------------------------------------------------------
mappingEmpiricalData <- function(species) {
  
  # Empirical catch by timeblock map
  pdf(paste0("C://KDale/Projects/Phenology/Figures/", species, "/logN1_Empirical.pdf"), width = 10, height = 6)
  print(ggplot(northAmerica) +
          facet_wrap(~timeblock, nrow = 1) +
          geom_sf(data = subset(data.sf, abundance_logN1_scaled == 0), col = "gray80", size = 0.2, pch = 20, alpha = 0.5, inherit.aes = FALSE) +
          geom_sf(data = subset(data.sf, abundance_logN1_scaled > 0), aes(col = abundance_logN1_scaled), size = 1.5, alpha = 1, pch = 20, inherit.aes = FALSE) +
          geom_sf(fill = "gray20") +
          xlim(min(data$longitude), max(data$longitude)) +
          ylim(min(data$latitude), max(data$latitude)) +
          scale_color_gradient(low = "lightgoldenrod2", high = "firebrick4") +
          theme_classic(base_size = 8) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude") +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
  # Empirical catch overall map
  pdf(paste0("C://KDale/Projects/Phenology/Figures/", species, "/logN1_Empirical_overall.pdf"), width = 7, height = 4)
  print(ggplot(northAmerica) +
          geom_sf(data = subset(data.sf, abundance_logN1_scaled == 0), col = "gray80", size = 0.2, pch = 20, alpha = 0.5, inherit.aes = FALSE) +
          geom_sf(data = subset(data.sf, abundance_logN1_scaled > 0), aes(col = abundance_logN1_scaled), size = 1.5, alpha = 1, pch = 20, inherit.aes = FALSE) +
          geom_sf(fill = "gray20") +
          xlim(min(data$longitude), max(data$longitude)) +
          ylim(min(data$latitude), max(data$latitude)) +
          scale_color_gradient(low = "lightgoldenrod2", high = "firebrick4") +
          theme_classic(base_size = 8) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude") +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
}

# PREDICTED OUPUT-------------------------------------------------------------
mappingPredictedOutput <- function(response, species, modelType) {
  
  #Summarize catch data by grid point, month, year, region, timeblock
  predictedCatch.grid <- p %>% group_by_at(., c("gridid", "month", "year", "timeblock", "sst_roms", "latitude", "longitude","X", "Y", "salinity_roms", "ssh_roms")) %>%
    dplyr::summarize(., mean_est = mean(est_retransform_weighted)) %>% # summarize across grid cells
    mutate(., est_x_month = mean_est * month) %>%
    mutate(., timeblock = factor(timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")))
  
  # Predicted catch abundance for each month, across all years
  predictedCatch.grid.month <- predictedCatch.grid %>%
    st_as_sf(., coords = c("longitude", "latitude")) %>% 
    st_set_crs(4326) %>% # WGS 84 - geographic reference system
    st_transform(5070) %>% 
    st_join(hex, .) %>%
    group_by_at(., c("gridid", "month")) %>% 
    summarize(., mean_est = mean(mean_est)) %>% 
    subset(., !is.na(month))
  
  cairo_pdf(paste0("Figures/", species, "/Predicted_", modelType,"_months.pdf"), width = 10, height = 7)
  print(ggplot(northAmerica) +
          facet_wrap(~month, nrow = 3) +
          geom_sf(data = predictedCatch.grid.month, mapping = aes(fill = mean_est), color = NA) +
          geom_sf(fill = "gray20") +
          scale_fill_gradient(response, low = "lightgoldenrod2", high =  "firebrick", na.value = "gray90") +
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          theme_classic(base_size = 12) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          xlim(min(data$longitude), max(data$longitude)) +
          ylim(min(data$latitude), max(data$latitude)) +
          labs(x = "Longitude", y = "Latitude"))
  dev.off()
  
  # Predicted abundance for each timeblock
  predictedCatch.grid.timeblock.sf <- predictedCatch.grid.timeblock %>%
    st_as_sf(., coords = c("longitude", "latitude")) %>% 
    st_set_crs(4326) %>% # WGS 84 - geographic reference system
    st_transform(5070) %>% 
    st_join(hex, .) %>%
    group_by_at(., c("gridid", "X", "Y", "timeblock")) %>% 
    summarize(., mean_est = mean(mean_est))%>% 
    subset(., !is.na(timeblock))
  
  cairo_pdf(paste0("Figures/", species, "/Predicted_", modelType,"_timeblocks.pdf"), width = 10, height = 6)
  print(ggplot() +
          geom_sf(data = northAmerica, fill = "gray20") +
          geom_sf(data = predictedCatch.grid.timeblock.sf, mapping = aes(fill = mean_est), color = NA) +
          facet_wrap(~timeblock, nrow = 1) +
          scale_fill_gradient(response, low = "lightgoldenrod2", high =  "firebrick", na.value = "gray90") +
          xlim(min(data$longitude), max(data$longitude)) +
          ylim(min(data$latitude), max(data$latitude)) +
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          theme_classic(base_size = 12) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude"))
  dev.off()
  
}

# REGIONAL SEASONALITY --------------------------------------------------
regionalSeasonality <- function() {
  
  speciesList = tradeoffs$species
  
  # Set up storage dataframe
  cog.region = data.frame(species = speciesList, region = NA, cog.month = NA)
  
  # Loop through species, calculating regional seasonal COG for each
  for (i in 1:length(speciesList)) {
    
    species = as.character(speciesList[i])
    
    # Load prediction object for each species
    predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", tradeoffs$short_formula[i], ".rdata")
    load(predictionObjectName)
    
    # Rotate months
    p$month_rotated = rotateMonths(p,abundanceTerm = "est_retransform_weighted",  p$month, direction = "rotate")
    
    # Add latitudinal region
    p$region = 0
    p <- p %>%
      mutate(region = case_when(
        latitude < 34.5 ~ "Southern CCE",
        latitude < 42 ~ "Central CCE",
        latitude < 48.3 ~ "OR/WA",
        latitude < 54.4 & longitude > -141 ~ "British Columbia",
        TRUE ~ "Gulf of Alaska"
      ))
    
    # Calculate cog for each species; add species
    for (j in 1:length(regions)) {
      region <- subset(p, p$region == regions[j]) # subset to each region
      if (nrow(region) == 0) next
      cog <- calculateCOG(region, species) %>% # Calculate COG
        mutate(region = regions[j])
      if (j == 1) {
        cog.region.species = cog
      } else {
        cog.region.species = bind_rows(cog.region.species, cog)
      }
    }
    # Calculate average seasonal COG for each region
    cog.region.species <- cog.region.species %>%
      group_by_at(c("species", "region")) %>% 
      summarize(cog.month = mean(cog.month))
    
    if (i == 1) {
      cog.region = cog.region.species
    } else {
      cog.region = bind_rows(cog.region, cog.region.species)
    }
  }
  
  cog.region$region = factor(cog.region$region, levels = regions)
  cog.region$species = factor(cog.region$species, levels = rev(speciesLevels))
  
  pdf(file = "Figures/Fig#_Regional_COG/Fig#_Regional_COG.pdf", width = 6, height = 6)
  print(ggplot(data = cog.region) +
          geom_tile(mapping = aes(x = region, y = species, fill = cog.month)) +
          theme_classic(base_size = 16) +
          coord_cartesian(expand = F) +
          labs(x = "", y = "") +
          theme(axis.text.y = element_text(face = "italic"),
                axis.text.x = element_text(angle = 45, hjust=1),
                panel.border = element_rect(fill = NA, linewidth = 0.5),
                legend.key.size = unit(0.75, "cm"),
                legend.text = element_text(size = 10), 
                legend.title = element_text(size = 12)) +
          scale_fill_stepsn("Seasonal COG \n(month)", breaks = c(1,2,3,4,5,6,7,8,9,10,11),
                            limits = c(1,12),
                            colors = c("#213ab7", "slateblue",  "salmon", "goldenrod1",  "cadetblue3", "dodgerblue3")))
  
  dev.off()
}

# FREQUENCY OF SAMPLING -------------------------------------------------------
# Plot frequency of stations across years of the study period
frequencyOfSampling <- function(speciesData) {
  
  # Frequency of sampling across months/years by program
  freq_months <- speciesData %>% group_by_at(., c("timeblock", "month", "program")) %>%
    summarize(n = n())
  
  print(ggplot(freq_months) +
          geom_tile(aes(month, timeblock, fill = n)) +
          facet_wrap(vars(program)) +
          scale_fill_gradient2("Number \n of tows", low = "white", mid = "cornflowerblue", high = "darkblue") +
          scale_x_continuous(n.breaks = 12) +
          labs(x = "Month", y = "Program") +
          coord_cartesian(expand = FALSE) +
          theme_classic())
  
  # Frequency of sampling across months/years by region
  freq_months <- speciesData %>% group_by_at(., c("timeblock", "month", "region")) %>%
    summarize(n = n())
  
  pdf(paste0("Figures/",species, "/Frequency_of_sampling_", modelType,".pdf"), width = 4, height = 6)
  print(ggplot(freq_months) +
          geom_tile(aes(month, timeblock, fill = n)) +
          facet_wrap(vars(region)) +
          scale_fill_gradient2("Number \n of tows", low = "white", mid = "cornflowerblue", high = "darkblue") +
          scale_x_continuous(n.breaks = 12) +
          labs(x = "Month", y = "Region") +
          coord_cartesian(expand = FALSE) +
          theme_classic())
  dev.off()
}

# LASSO/RIDGE -------------

lasso_ridge <- function() {
  
  depVars = c(
    "cog.month.var","cog.lat.var", "cog.long.var",  
    "trailing.var","leading.var")
  
  depVarNames = c(
    "A. Seasonal COG\nvariance (month)", "B. Latitudinal\nCOG variance (\u00b0)", "C. Longitudinal\nCOG variance (\u00b0)",  
    "D. Trailing edge\nvariance (\u00b0)","E. Leading edge\nvariance (\u00b0)"
  )
  
  predictorLabels = c("factor(lifeHistory)Coastal pelagic" = "Life history (coastal pelagic)",
                      "factor(lifeHistory)Groundfish" = "Life history (groundfish)",
                      "factor(lifeHistory)Mesopelagic" = "Life history (mesopelagic)",
                      "factor(horz_habitat)Continental shelf" = "Inshore-offshore habitat (continental shelf)",
                      "factor(horz_habitat)Coastal" = "Inshore-offshore habitat (coastal)",
                      "factor(horz_habitat)Oceanic" = "Inshore-offshore habitat (oceanic)",
                      "maxBodySize_cm" = "Maximum body size (cm)",
                      "spawningBreadth_regionalAverage" = "Spawning season length (d)",
                      "maxAge" = "Maximum age (y)",
                      "factor(fisheries_status)Highly commercial" = "Fishery status (highly commercial)",
                      "factor(fisheries_status)Commercial" = "Fishery status (commercial)",
                      "factor(fisheries_status)None" = "Fishery status (none)",
                      "adultLatRange" = "Adult latitudinal range (\u00b0)",
                      "larvalLatRange" = "Larval latitudinal range (\u00b0)",
                      "adultLarvalDiff" = "Adult-larval latitudinal\ndifference (\u00b0)",
                      "maxDepthLove" = "Maximum depth (m)",
                      "trophic_level" = "Trophic level",
                      "smith" = "Niche breadth")
  
  predictorLevels = c("factor(lifeHistory)Coastal pelagic",
                      "factor(lifeHistory)Groundfish",
                      "factor(lifeHistory)Mesopelagic",
                      "factor(horz_habitat)Oceanic",
                      "factor(horz_habitat)Continental shelf",
                      "factor(horz_habitat)Coastal",
                      "spawningBreadth_regionalAverage",
                      "maxAge",
                      "maxBodySize_cm",
                      "adultLatRange",
                      "larvalLatRange",
                      "adultLarvalDiff",
                      "maxDepthLove",
                      "trophic_level",
                      "smith",
                      "factor(fisheries_status)Highly commercial",
                      "factor(fisheries_status)Commercial",
                      "factor(fisheries_status)None")
  
  formula = as.formula(~ factor(lifeHistory) + factor(fisheries_status) +  factor(horz_habitat) +
                         spawningBreadth_regionalAverage + maxAge + maxBodySize_cm + adultLatRange + 
                         + maxDepthLove + trophic_level + smith)
  
  # Convert categorical to factor variables if necessary
  X <- model.matrix(formula, data = tradeoffs)
  
  # Extract independent variables
  X_full <- model.matrix(formula, data = tradeoffs)[, -1]
  
  # Obtain coefficient names
  all_coef_names <- c("(Intercept)", colnames(X_full))  # Ensure intercept is included
  
  # Boot function
  lassoRidge_boot <- function(data, indices, depVar, alpha) {
    
    # Resample the data
    boot_data <- data[indices, ]  
    
    X_boot <- model.matrix(formula, data = tradeoffs)[, -1]
    y_boot <- as.matrix(boot_data[, depVar])  # Ensure correct response variable selection
    
    # cross-validation to estimate best lambda
    cv <- cv.glmnet(X_boot, y_boot, alpha = alpha, grouped = F) 
    
    # Extract lambda with minimum value
    best_lambda <- cv$lambda.min  
    
    # Run model
    model <- glmnet(X_boot, y_boot, alpha = alpha, lambda = best_lambda, standardize = T)
    
    # Extract coefficients
    coefs <- as.numeric(coef(model))
    coef_names <- rownames(coef(model))
    
    # Create a full coefficient vector initialized with zeros
    fixed_coefs <- setNames(rep(0, length(all_coef_names)), all_coef_names)
    
    # Fill in existing coefficients
    common_names <- intersect(all_coef_names, coef_names)  # Ensure matching names
    fixed_coefs[common_names] <- coefs[match(common_names, coef_names)]
    
    return(as.numeric(fixed_coefs))  # Ensure consistent vector length
  }
  
  # Select best value of alpha
  # y <- as.matrix(tradeoffs_scaled[c(depVars[1])])
  # 
  # foldid <- sample(1:3, size = nrow(tradeoffs), replace = TRUE)
  # cv1  <- cv.glmnet(X, y, foldid = foldid, alpha = 1)
  # cv.5 <- cv.glmnet(X, y, foldid = foldid, alpha = 0.5)
  # cv.25 <- cv.glmnet(X, y, foldid = foldid, alpha = 0.25)
  # cv.75 <- cv.glmnet(X, y, foldid = foldid, alpha = 0.75)
  # cv0  <- cv.glmnet(X, y, foldid = foldid, alpha = 0)
  # 
  # par(mfrow = c(2,2))
  # plot(cv1); plot(cv.5); plot(cv0)
  # plot(log(cv1$lambda)   , cv1$cvm , pch = 19, col = "yellow",
  #      xlab = "log(Lambda)", ylab = cv1$name)
  # points(log(cv0$lambda) , cv0$cvm , pch = 19, col = "black")
  # points(log(cv.25$lambda) , cv.25$cvm , pch = 19, col = "orange4")
  #   points(log(cv.5$lambda), cv.5$cvm, pch = 19, col = "orange2")
  # points(log(cv.75$lambda) , cv.75$cvm , pch = 19, col = "goldenrod1")
  # 
  # legend("topleft", legend = c("alpha= 1", "alpha= .5", "alpha 0"),
  #        pch = 19, col = c("red","grey","blue"))
  
  alphas = c(0)
  
  # Loop through dependent variables - LASSO = 1, ridge = 0
  for (j in 1:length(alphas)) {
    for (i in 1:length(depVars)) {
      
      alpha = alphas[j]
      depVar = depVars[i]
      
      print(alpha)
      print(depVar)
      
      # Set up model results matrix
      y <- as.matrix(tradeoffs[c(depVar)])
      
      # Run cross-validation to obtain best lambda
      cv <- cv.glmnet(X, y, alpha = alpha, grouped = F) 
      
      # Extract best lambda with lowest value
      best_lambda <- cv$lambda.min  
      
      # Fit model using best lambda
      model <- glmnet(X, y, alpha = alpha, lambda = best_lambda, standardize = T)
      
      # Extract coefficients from lasso model
      coefs <- coef(model)
      
      # Run bootstrap (adjust R for more accuracy)
      set.seed(456)
      boot_results <- boot(tradeoffs, lassoRidge_boot, depVar = depVar, alpha = alpha, R = 10000)
      
      # Compute standard errors
      coef_se <- apply(boot_results$t, 2, sd, na.rm = TRUE)[-1]
      
      # Remove intercent
      coefs_noInt <- 
        subset(coefs, rownames(coefs) != "(Intercept)")
      
      # Construct dataframe
      coef_df <- data.frame(
        variable = colnames(X_full),
        estimate = as.numeric(coefs_noInt),
        se = coef_se
      )
      
      # Calculate standard errors and add additional information
      coef_df <- coef_df %>% 
        subset(variable != "(Intercept)") %>% 
        mutate(se_upper = estimate + se, se_lower = estimate - se) %>% 
        mutate(se_upper_round = round(se_upper, digits = 2), se_lower_round = round(se_lower, digits = 2)) %>% 
        mutate(significant = case_when(
          estimate == 0 ~ "Not included", 
          se_upper > 0 & se_lower < 0 | se_upper_round == 0 | se_lower_round == 0  ~ "Non-significant",
          TRUE ~ "Significant"
        )) %>% 
        mutate(variable = factor(variable, levels = predictorLevels))
      
      if (alpha == 1) {
        write.csv(coef_df, paste0("Results/04_LASSO_Ridge/", depVar, "_LASSO_coefficients.csv"))
      } else if (alpha == 0) {
        write.csv(coef_df, paste0("Results/04_LASSO_Ridge/", depVar, "_Ridge_coefficients.csv"))
      } else {
        write.csv(coef_df, paste0("Results/04_LASSO_Ridge/", depVar, "_ElasticNet_coefficients.csv"))
      }
    }
  }
  
  plotlist = list()
  
  # Fig. 7 ------
  for (j in 1:length(alphas)) {
    for (i in 1:length(depVars)) {
      
      alpha = alphas[j]
      depVar = depVars[i]
      
      # Access saved dataframes
      if (alpha == 0) {
        coef_df <- read.csv(paste0("Results/04_LASSO_Ridge/", depVar, "_Ridge_coefficients.csv")) %>% 
          mutate(variable = factor(variable, levels = predictorLevels))
        filename = paste0("Figures/Fig7_LassoRidgeCoefficients/Ridge_", depVar, ".pdf")
      } else if (alpha == 1) {
        coef_df <- read.csv(paste0("Results/04_LASSO_Ridge/", depVar, "_LASSO_coefficients.csv")) %>% 
          mutate(variable = factor(variable, levels = predictorLevels))
        filename = paste0("Figures/Fig7_LassoRidgeCoefficients/LASSO_", depVar, ".pdf")
      } else {
        coef_df <- read.csv(paste0("Results/04_LASSO_Ridge/", depVar, "_ElasticNet_coefficients.csv")) %>% 
          mutate(variable = factor(variable, levels = predictorLevels))
        filename = paste0("Figures/Fig7_LassoRidgeCoefficients/ElasticNet_", depVar, ".pdf")
      }
      
      # Create plot
      p <- ggplot() +
        geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "gray70") +
        geom_point(data = coef_df, aes(y = variable, x = estimate, color = significant), size = 2) +
        scale_y_discrete(labels = predictorLabels) +
        geom_errorbar(data = coef_df,aes(y = variable, xmin = estimate - se, xmax = estimate + se, color = significant), width = 0.2, linewidth = 0.5) +
        theme_classic(base_size = 12) +
        scale_color_manual("", values = c("Not included" = "gray80", "Non-significant" = "gray40", "Significant" = "seagreen3")) +
        theme(legend.position = "bottom", axis.line = element_blank(),
              title = element_text(size = 8),
              axis.text = element_text(size = 9),
              axis.title = element_text(size = 10),
              panel.border = element_rect(fill = NA, color = 1)) +
        labs(title = depVarNames[i],
             y = "", x = "Coefficient estimate")
      
      if(i != 1) {
        p = p + theme(axis.text.y = element_blank()) + labs(y = "", x = "")
      }
      
      plotlist[[i]] = p
      
      pdf(filename, width = 8, height = 7)
      print(p)
      dev.off()
    }
    
    if (alpha == 0) {
      filename = paste0("Figures/Fig7_LassoRidgeCoefficients/Ridge_AllPredictors.pdf")
    } else if (alpha == 1) {
      filename = paste0("Figures/Fig7_LassoRidgeCoefficients/LASSO_AllPredictors.pdf")
    } else {
      filename = paste0("Figures/Fig7_LassoRidgeCoefficients/ElasticNet_AllPredictors.pdf")
    }
    
    pdf(file = filename, width = 9, height = 4)
    print(ggarrange(plotlist = plotlist, nrow = 1, common.legend = T, legend = "bottom",
                    legend.grob = get_legend(plotlist[[1]]),
                    widths = c(2.4,1,1,1,1,1)))
    dev.off()
  }
}


# RUN FUNCTIONS -------------------------------------------------------------------------
## LOAD SHAPEFILES --------
all_tows_roms <- read.csv(file = "Data/AllTows_200nm_ROMS.csv")

northAmerica <- st_read("Data/Shapefiles/NorthAmerica/boundary_p_v2.shp") %>%
  subset(COUNTRY != "water/agua/d'eau" & COUNTRY != "FN") %>% 
  group_by(COUNTRY) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  st_transform(4326)

northAmerica_forHex <- read_sf("Data/Shapefiles/North_South_America/North_South_America.shp") %>% st_union() %>%
  st_transform(., crs = "EPSG:5070")

gulfOfCalifornia <- st_read("Data/Shapefiles/World_Seas_IHO_v3/World_Seas_IHO_v3.shp") %>% 
  subset(MRGID == 4314) %>% dplyr::select(., geometry) %>% st_transform(., crs = "EPSG:5070")

## SET GLOBAL VARS ----------
# Colors
regions = c("Southern CCE", "Central CCE", "OR/WA", "British Columbia", "Gulf of Alaska")
lifeHistoryColors = c("Coastal pelagic" = "#92C6EA", "Groundfish" = "#CD7E09", "Mesopelagic" = "gray20")
modelColors = c("base" = "gray80", "geo" = "dodgerblue3", "pheno" = "indianred2", "both" = "orchid4")
monthColors = c("#112040FF", "#263E91FF", "#1E6EA1FF", "#3D9BACFF", "#8AC0BAFF", "#CBDCCD", "#EFE29CFF", "#C2B63AFF", "#759906FF" ,"#237824FF", "#124F2BFF" ,"#172313FF")

# Set ggplot theme
theme_set(theme_classic(base_size = 12))

modelNames = data.frame(modelType =  c("base", "geo", "pheno", "both"))

# Set p value alpha
alpha = 0.05

# Load list of potential covariates and timeblocks
covariates = read_xlsx("Analysis/Covariates.xlsx", sheet = 1)
timeblocks <- read_xlsx("Data/timeblocks.xlsx", sheet = 1)

# Load species info information
speciesInfo <- read_xlsx("Data/Species_Info.xlsx")
speciesLevels = speciesInfo$species

# Load list of best models
tradeoffs <- read_xlsx("Results/02_Tradeoff Comparison/Tradeoffs_summary.xlsx", sheet = 1)

## LOAD TRADEOFF FILE -----------------
# Add species-specific information
tradeoffs <- tradeoffs %>%
  group_by_at(c("species", "main_formula")) %>%
  filter(., delta_sum_loglik == 0) %>% 
  merge(., speciesInfo) %>% 
  merge(., covariates[c("main_formula", "short_formula")]) %>% 
  mutate(species = factor(species, levels = speciesLevels)) %>% 
  arrange(phylogenetic_order)

# Single species
# tradeoffs <- subset(tradeoffs, species == "Vinciguerria lucetia")

## CREATE PREDICTION OBJECTS------------
makeNewPrediction = T
makeNewGrid = F

pb <- txtProgressBar(min = 0, max = nrow(tradeoffs), char = "=", style = 3)

for (i in 1:nrow(tradeoffs)) {
  
  # Set species and model
  species = tradeoffs$species[i]
  modelType = tradeoffs$model[i]
  shortFormula = tradeoffs$short_formula[i]
  
  # Filenames
  fitName = paste0("Results/", species, "/Models/abundance_logN1_scaled_speciesRange_", modelType, ".rdata")
  predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", shortFormula, ".rdata")
  gridFilename = paste0("Analysis/PredictionGrids/", species, "_grid.rdata")
  
  # Create prediction objects
  if(!file.exists(predictionObjectName) | makeNewPrediction == T) {
    createPredictionObjects(
      species = species,
      fitName = fitName,
      predictionObjectName = predictionObjectName,
      gridFilename = gridFilename, 
      makeNewGrid = makeNewGrid, # Change these to force a new grid or new prediction
      makeNewPrediction = makeNewPrediction,
      stacked = TRUE
    )
  }
  setTxtProgressBar(pb, i)
}
close(pb)

## INDIVIDUAL SPECIES FIGURES & ANALYSES ----------------------
# Set up progress bar
pb <- txtProgressBar(min = 0, max = nrow(tradeoffs), char = "=", style = 3)

# Set up storage data structures
cog_grid_plots <- list()
supplFig1Plots <- list()
gam_plots <- list()
cumulativeCatchPlots <- list()
conditionalEffectsPlots <- list()

modelPerformance <- data.frame(species = tradeoffs$species, commonName = tradeoffs$common_name, nPosTows = NA, nTows = NA, pearson = NA, pearson.p = NA)

for (i in 1:nrow(tradeoffs)) {
  
  # Set species and model
  species = as.character(tradeoffs$species[i])
  commonName = tradeoffs$common_name[i]
  modelType = tradeoffs$model[i]
  shortFormula = tradeoffs$short_formula[i]
  
  # Filenames
  fitName = paste0("Results/", species, "/Models/abundance_logN1_scaled_speciesRange_", modelType, ".rdata")
  predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", shortFormula, ".rdata")
  gridFilename = paste0("Analysis/PredictionGrids/", species, "_grid.rdata")
  
  # Load objects
  load(fitName) # data and model
  load(predictionObjectName) # prediction objects
  load(gridFilename) # grid.df
  
  prediction_grid_roms <- merge(prediction_grid_roms, timeblocks) # add timeblocks
  prediction_grid_roms$timeblock = factor(prediction_grid_roms$timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019"))
  
  # Get species dataset, cast to a spatial object
  data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>%
    # st_combine(.) %>%
    st_set_crs(4326) %>% 
    st_transform(5070) %>% # Albers equal area conic (aka conic albers) - projected
    st_as_sf()
  
  # Rotate months
  p$month_rotated = rotateMonths(data = p, abundanceTerm = "est_retransform_weighted", months = p$month, direction = "rotate")
  
  ### Make hex grid ----
  hex <- makeHex(data = data) %>%
    st_transform(., 5070)
  
  # Get response variable
  response = strsplit(as.character(fit$formula), "~")[[1]][1]
  
  ## Predictive performance -------------
  pearson <- cor.test(p.original$est_retransform_weighted, p.original$abundance_logN1_scaled, method = "pearson")
  
  modelPerformance$species[i] = species
  modelPerformance$nPosTows[i]= sum(data$abundance_logN1_scaled > 0)
  modelPerformance$nTows[i] = nrow(data)
  modelPerformance$pearson[i] = pearson$estimate
  modelPerformance$pearson.p[i] = pearson$p.value
  
  if("epsilon_st" %in% colnames(p)) { # For spatiotemporal models
    # Summarize catch data by grid point, month, timeblock
    predictedCatch.grid.timeblock <- p %>% group_by_at(., c("gridid", "timeblock",  "latitude", "longitude","X", "Y")) %>%
      summarize(.,
                mean_est = mean(est_retransform_weighted),
                mean_rf = mean(fit$family$linkinv(est_rf)),
                mean_epsilon_st = mean(fit$family$linkinv(epsilon_st)),
                mean_fixed = mean(fit$family$linkinv(est_non_rf)),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      ) %>%
      mutate(., timeblock = factor(timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")))
    
    # Summarize catch data by latitude, month, timeblock
    predictedCatch.lat.timeblock <- p %>% group_by_at(., c("latitude", "timeblock" )) %>%
      summarize(.,
                mean_est = mean(est_retransform_weighted),
                # mean_rf = mean(fit$family$linkinv(est_rf)),
                mean_fixed = mean(fit$family$linkinv(est_non_rf)),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      )
  } else if ("omega_s" %in% colnames(p)) { # Spatial models
    predictedCatch.grid.timeblock <- p %>% group_by_at(., c("gridid", "timeblock", "latitude", "longitude","X", "Y")) %>%
      summarize(.,
                mean_est = mean(est_retransform_weighted),
                mean_rf = mean(fit$family$linkinv(est_rf)),
                mean_fixed = mean(fit$family$linkinv(est_non_rf)),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      )
    predictedCatch.lat.timeblock <- p %>% group_by_at(., c("latitude", "timeblock")) %>%
      summarize(.,
                mean_est = mean(est_retransform_weighted),
                mean_rf = mean(fit$family$linkinv(est_rf)),
                mean_fixed = mean(fit$family$linkinv(est_non_rf)),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      )
  } else { # Base or geo models (no epsilon or omega)
    predictedCatch.grid.timeblock <- p %>% group_by_at(., c("gridid", "timeblock", "latitude", "longitude","X", "Y")) %>%
      summarize(.,
                mean_est = mean(est_retransform_weighted),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      )
    
    predictedCatch.lat.timeblock <- p %>% group_by_at(., c("latitude", "timeblock")) %>%
      summarize(.,
                mean_est = mean(est_retransform_weighted),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms))
  }
  
  ## Abundances by month -------
  year.month.sum <- group_by_at(p, c("year", "month")) %>%
    summarize(sum_est = sum(est_retransform_weighted))
  
  pdf(paste0("Figures/", species, "/Abundances_by_month.pdf"), width = 6, height = 4)
  print(ggplot(year.month.sum, aes(x = month, y = sum_est, color = year)) +
          geom_point(size = 4))
  dev.off()
  
  ## COG maps --------------
  cog_grid_plots[[i]] <- gridCOG(species = species,commonName = commonName, modelType = modelType)
  
  ## Cumulative curves --------
  cumulativeCatchPlots[[i]] <- cumulativeCatch(species)
  
  ## Empirical and predicted maps ---------
  empiricalPredictedPlots <- empiricalVsPredicted( species = species, commonName = commonName, p.original = p.original)
  mappingEmpiricalData(species = species)
  gam_plots[[i]] <- empiricalPredictedPlots[[1]] # GAM of empirical vs modeled
  supplFig1Plots[[i]] <- empiricalPredictedPlots[[2]] # Maps of empirical vs modeled
  
  ## Env conditional effects figures - takes a long time to run! ----
  if(modelType == "geo") {
    
    load(paste0("Results/", species, "/Models/abundance_logN1_scaled_speciesRange_both.rdata"))
    
    # Change to using the "both" model since predict_response does not work for the geo model
    conditionalEffectsPlots[[i]] <- conditionalEffects(newPrediction = T, species = species,
                                                       fit = fit,
                                                       shortFormula = shortFormula,
                                                       data = data,
                                                       grid.df = grid.df,
                                                       prediction_grid_roms = prediction_grid_roms)
  } else {
    conditionalEffectsPlots[[i]] <- conditionalEffects(newPrediction = T, species = species,
                                                       fit = fit,
                                                       shortFormula = shortFormula,
                                                       data = data,
                                                       grid.df = grid.df,
                                                       prediction_grid_roms = prediction_grid_roms)
  }
  
  ## Mapping predicted output --------
  mappingPredictedOutput(response = response, species = species, modelType = modelType)
  
  ## Frequency of sampling for species----
  # frequencyOfSampling(speciesData = data)
  
  setTxtProgressBar(pb, i)
  
  gc()
}
close(pb)

## MULTISPECIES ANALYSES & FIGURES-----------------
### Model performance ------
suppltable1 = modelPerformance %>% 
  mutate(freqPosTows = nPosTows / nTows) %>% 
  rename("Common name" = commonName,
         "Scientific name" = species, 
         "N total tows" = nTows,
         "Frequency of positive tows" = freqPosTows,
         "Pearson coefficient" = pearson)
suppltable1$`Scientific name` = factor(suppltable1$`Scientific name`, levels = speciesLevels)
suppltable1 = arrange(suppltable1, suppltable1$`Scientific name`)
write_xlsx(x = suppltable1, path = "Manuscript/Suppl_Table2_Model performance.xlsx", format_headers = F)
suppltable1 <- read_xlsx(path = "Manuscript/Suppl_Table2_Model performance.xlsx", sheet = 1)

lm(suppltable1$Pearson.coefficient ~ suppltable1$Frequency.of.positive.tows) %>% 
  summary()

### Plotlists ---------------
pdf("Figures/Suppl_Fig4_COG_Maps/Fig#_COG_Maps_singleLegend.pdf", height = 20, width = 20)
ggarrange(plotlist = cog_grid_plots, nrow = 4, ncol = 4)
dev.off()

pdf("Figures/Suppl_Fig2_GAM plots/SupplFig2_GAM_plots.pdf", height = 20, width = 20)
plot_grid(plotlist = gam_plots, nrow = 4)
dev.off()

pdf("Figures/Suppl_Fig1_Empirical vs modeled abundance/FigSuppl1_Empirical vs modeled abundance.pdf", height = 20, width = 20)
plot_grid(plotlist = supplFig1Plots, nrow = 4, ncol = 4)
dev.off()

pdf("Figures/Fig#_CumulativeCurves/Fig#_CumulativeCurves_AllSpecies.pdf", height = 20, width = 20)
plot_grid(plotlist = cumulativeCatchPlots, nrow = 4)
dev.off()

pdf("Figures/Fig#_EnvCovariates/Fig#_EnvCovariates_AllSpecies.pdf", height = 20, width = 30)
plot_grid(plotlist = conditionalEffectsPlots, nrow = 4, labels = tradeoffs$common_name)
dev.off()

### Tradeoff plots ---------------------
tradeoffsBarPlots(tradeoffs = tradeoffs)
tradeoffsDetailed()

### Leading and trailing variation ------
spawningRangePlots <- list()
spawningRangeSummary <- data.frame()

for (i in 1:nrow(tradeoffs)) {
  species= tradeoffs$species[i]
  commonName = tradeoffs$common_name[i]
  spawningRangeList = spawningRange(commonName = commonName, species = species, modelType = tradeoffs$model[i], shortFormula = tradeoffs$short_formula[i])
  spawningRangePlots[[i]] <- spawningRangeList[[1]]
  spawningRangeSpecies <- spawningRangeList[[2]]
  tradeoffs$leading.var[i] = var(spawningRangeSpecies$leadingLat) # cumulative abundance south to north, so trailing = 85th percentile
  tradeoffs$trailing.var[i] = var(spawningRangeSpecies$trailingLat) # cumulative abundance south to north, so trailing = 15th percentile
  tradeoffs$leading.mean[i] = mean(spawningRangeSpecies$leadingLat) # cumulative abundance south to north, so trailing = 85th percentile
  tradeoffs$trailing.mean[i] = mean(spawningRangeSpecies$trailingLat) # cumulative abundance south to north, so trailing = 15th percentile
  spawningRangeSummary <- bind_rows(spawningRangeSummary, spawningRangeSpecies)
}

pdf("Figures/Fig#_spawningRangePlots/Fig#_spawningRangePlots_AllSpecies.pdf", height = 24, width = 24)
plot_grid(plotlist = spawningRangePlots, nrow = 4)
dev.off()

### Center of gravity --------------
cog <- cogMultispecies(tradeoffs)
tradeoffs <- cog[[1]]
cog.all <- cog[[2]]
internalCOGTrends(tradeoffs = tradeoffs)
regionalSeasonality(speciesList = tradeoffs$species)

### Tradeoff metric regressions with year -----
metrics.allYears <- merge(spawningRangeSummary, cog.all, by = c("species", "year"))
write.csv(metrics.allYears, "Results/Metrics_allYears.csv")
metrics.allYears <- read.csv("Results/Metrics_allYears.csv")
lmSummarySignificant = trendsWithYear(metrics.allYears)

### Hypervolume --------------
tradeoffs_niche <- hypervolumeMultiSpecies(tradeoffs = tradeoffs, recalculate = T)
write_xlsx(tradeoffs_niche, "Results/Tradeoffs_niche_cog.xlsx", format_headers = F)
tradeoffs_niche <- read_xlsx(path = "Results/Tradeoffs_niche_cog.xlsx")
tradeoffs = tradeoffs_niche %>% 
  arrange(., phylogenetic_order)

### LASSO/Ridge --------------
lasso_ridge() 

### Table 1 ----
table1 <- tradeoffs %>% 
  mutate(smith = round(smith, digits = 2)) %>% 
  select("Adult category" = lifeHistory,
         "Common name" = common_name,
         "Scientific name" = species,
         "Niche breadth\n(Smith's hypervolume)" = smith,
         "Trophic level" = trophic_level,
         "Max depth (m)" = maxDepthLove, 
         "Adult latitudinal range (\u00b0)" = adultLatRange,
         "Max body size (cm)" = maxBodySize_cm,
         "Max age (yr)" = maxAge,
         "Spawning season length (d)" = spawningBreadth_regionalAverage, 
         "Fishery status" = fisheries_status,
         "Inshore-offshore\nhabitat" = horz_habitat)

write_xlsx(table1, "Manuscript/Table1_SpeciesInfo.xlsx", format_headers = F)

### Table 3 
table3 <- tradeoffs %>% 
  mutate(across(leading.var:cog.month.mean, ~ round(., digits = 2))) %>% 
  select("Common name" = common_name,
         "Scientific name" = species,
         "Seasonal COG mean (months)" = cog.month.mean,
         "Seasonal COG var. (months)" = cog.month.var,
         "Lat. COG mean (°)" = cog.lat.mean,
         "Lat. COG var. (°)" = cog.lat.var,
         "Long. COG mean (°)" = cog.long.mean ,
         "Long. COG var. (°)" = cog.long.var,
         "Trailing edge mean (°)" = trailing.mean,
         "Trailing edge variance (°)" = trailing.var,
         "Leading edge mean (°)" = leading.mean,
         "Leading edge var. (°)" = leading.var)

write_xlsx(table3, "Manuscript/Table3_TradeoffMetrics.xlsx", format_headers = F)
