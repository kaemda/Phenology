# LOAD PACKAGES ---------
# Data management
library(readxl)
library(xlsx)
library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)

# Analysis
library(sdmTMB)
library(ape)
library(hutils)

# Mapping
library(sf)

# Visualization
library(visreg)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

# Data access
library(ncdf4)

setwd("C://KDale/Projects/Phenology/")

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

# CALCULATE COG  ------------------------------------------------------------
calculateCenterOfGravity <- function(p, species) {
  
  # Center of gravity per year in relation to environmental variables
  years = unique(p$year)
  cog.year = data.frame(matrix(ncol = 0, nrow = length(years)), year = NA, cog.lat = NA, cog.month = NA, avgsst = NA, avgssh = NA, avgsal= NA, avgbd = NA, avgdfs = NA, avgspice = NA)
  
  for (i in 1:length(years)) {
    yearSubset = subset(p, year == years[i])
    cog.year$cog.lat[i] = sum(yearSubset$est_retransform*yearSubset$latitude)/sum(yearSubset$est_retransform)
    cog.year$cog.month[i] = sum(yearSubset$est_retransform*yearSubset$month)/sum(yearSubset$est_retransform)
    cog.year$avgsst[i] = mean(yearSubset$sst_roms)
    cog.year$avgssh[i] = mean(yearSubset$ssh_roms)
    cog.year$avgsal[i] = mean(yearSubset$salinity_roms)
    cog.year$avgdfs[i] = mean(yearSubset$distance_from_shore_m)
    cog.year$avgbd[i] = mean(yearSubset$bottom_depth)
    cog.year$avgspice[i] = mean(yearSubset$spice)
    cog.year$year[i] = years[i]
    cog.year$species[i] = species
  }
  
  return(cog.year)
  
}

# CENTRAL TENDENCY -------------------------------------------------------
centralTendency_figures <- function(response, species, modelType) {
  
  ## Central tendency: Complete grid map (overall)
  pdf(paste0("Figures/",species, "/CentralTendency.map.overall_", modelType, ".pdf"), width = 4, height = 6)
  print(ggplot(northAmerica) +
          ggtitle(label = species, subtitle = str_wrap(paste(fit$formula, speciesRange, programs), 50)) +
          geom_sf() +
          geom_point(data = centralTendency.grid, mapping = aes(x = (X*1000), y = (Y*1000), col = centralTendency, size = mean_est_grid)) +
          scale_color_gradient("Central tendency\n(month)", low = "lightgoldenrodyellow", high =  "darkred") +
          scale_size_continuous(response, range = c(0.2, 3)) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          xlim(min(prediction_grid_roms$X)*1000-1000, max(prediction_grid_roms$X)*1000+1000) +
          ylim(min(prediction_grid_roms$Y)*1000-1000, max(prediction_grid_roms$Y)*1000+1000) +
          labs(x = "Longitude", y = "Latitude") +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
  # Central tendency per region per timeblock: scatterplot
  pdf(paste0("Figures/",species, "/CentralTendency_region.timeblock_", modelType, ".pdf"), width = 10, height = 6)
  print(ggplot(data = centralTendency.region.timeblock, mapping = aes(x = as.numeric(timeblock), y = centralTendency, col = region, fill = region)) +
          ggtitle(label = species, subtitle = str_wrap(paste(fit$formula, speciesRange, programs), 60)) +
          geom_point(size = 2) +
          geom_smooth(method = "lm") +
          scale_color_manual("Region", values = regionColors) +
          scale_fill_manual("Region", values = regionColors) +
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          labs(x = "Time period", y = response) +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
  # Central tendency per region per year: scatterplot
  pdf(paste0("Figures/",species, "/CentralTendency.region.year_", modelType, ".pdf"), width = 10, height = 6)
  print(ggplot(data = centralTendency.region.year, mapping = aes(x = year, y = centralTendency, col = region, fill = region)) +
          geom_point(size = 2) +
          geom_smooth(method = "lm") +
          scale_color_manual("Region", values = regionColors) +
          scale_fill_manual("Region", values = regionColors) +
          theme_classic(base_size = 12) +
          ggtitle(label = species, subtitle = str_wrap(paste(fit$formula, speciesRange, programs), 50)) +
          labs(x = "Year", y = "Central tendency") +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
  # Significant changes in central tendency: map
  pdf(paste0("Figures/",species, "/CentralTendency.significant.changes_", modelType,".pdf"), width = 4, height = 6)
  print(ggplot(northAmerica) +
          ggtitle(label = species, subtitle = str_wrap(paste(fit$formula, speciesRange, programs), 50)) +
          geom_sf() +
          geom_point(data = subset(grid.df, cog_pvalue < 0.05), mapping = aes(x = (X*1000), y = (Y*1000), col = days)) +
          scale_color_gradient2("Shift in central \ntendency (days)",  low = "red", mid = "white", high =  "blue", midpoint = 0) +
          #scale_size_continuous("P", breaks = 0.05, range = c(0, 2), trans = "reverse") +
          xlim(min(prediction_grid_roms$X)*1000-1000, max(prediction_grid_roms$X)*1000+1000) +
          ylim(min(prediction_grid_roms$Y)*1000-1000, max(prediction_grid_roms$Y)*1000+1000) +
          labs(x = "Longitude", y = "Latitude") +
          theme_classic(base_size = 12) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
  # Center of gravity per year in relation to environmental variables
  cog.year <- calculateCenterOfGravity(p, species)
  
  # Linear regressions between central tendency (lat and month) and environmental variables
  xaxes = c("avgsst", "avgssh", "avgsal")
  yaxes = c("cog.lat", "cog.month")
  xaxislabels = c("Average sea surface \n temperature [\u00b0C]", "Average sea\nsurface height [m]", "Average salinity [ppt]")
  yaxislabels = c("Center of gravity:\nlatitude [\u00b0]", "Center of gravity:\nmonth")
  colnames = c("y", "x", "p", "slope", "R2")
  lms <- as.data.frame(matrix(ncol = 0, nrow = length(xaxes)*length(yaxes))) # create dataframe to hold lm results
  index = 1
  for (j in 1:length(yaxes)) {
    for (i in 1:length(xaxes)) {
      lm <- lm(data = cog.year, formula = as.formula(paste0(yaxes[j], "~", xaxes[i]))) %>% summary()
      lms$y[index] = yaxes[j]
      lms$x[index] = xaxes[i]
      lms$p[index] = lm$coefficients[8]
      lms$slope[index] = lm$coefficients[2]
      lms$r2[index] = lm$r.squared
      index = index+1
    }
  }
  
  plotlist <- list()
  index = 1
  for (j in 1:length(yaxes)) {
    for (i in 1:length(xaxes)) {
      
      annotations1 <- data.frame(
        xpos = c(-Inf), ypos =  c(-Inf), #left-bottom
        annotateText = paste0("p = ", round(lms$p[index], digits = 2), " , R2 = ", round(lms$r2[index], digits = 2)),
        hjustvar = c(-.3),   #shifts bottom left 'Text' to the right; make more negative to move it further right),  #shifts top right 'texT' to the left; make more positive to move it further left
        vjustvar = c(-1))    #shifts bottom left 'Text' upward; make more negative to move it further up)
      
      plotlist[[index]] <- print(ggplot(cog.year, mapping = aes_string(x = xaxes[i], y = yaxes[j])) +
                                   geom_point(size = 3) +
                                   geom_smooth(method = "lm", color = "black") +
                                   theme_classic(base_size = 14) +
                                   coord_cartesian(expand = TRUE) +
                                   geom_text(data = annotations1, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText))+
                                   labs(x = xaxislabels[i], y = yaxislabels[j]))
      index = index + 1
    }
  }
  pdf(paste0("Figures/",species, "/COG_month_vs_latitude_", modelType,".pdf"), width = 12, height = 12)
  print(cowplot::plot_grid(plotlist = plotlist, nrow = 4, align = "hv"))
  dev.off()
  
} 

#------------------------------------------------------------------------------#
# CT MULTISPECIES -----------
centralTendencyMultiSpecies <- function(speciesList, modelTypes) {
  
  cog.all = data.frame(year = integer(), cog.lat = double(), cog.month = double(), avgsst = double(), avgssh = double(), avgsal = double(), avgspice = double())
  lms.all <- data.frame(y = character(), x = character(), p = double(), slope = double(), R2 = double(), species = character()) # create dataframe to hold lm results
  
  for (s in 1:length(speciesList)) {
    
    # Load prediction object for each species
    predictionObjectName = paste0("Results/", speciesList[s], "/Models/prediction_objects_", modelTypes[s], ".rdata")
    load(predictionObjectName)
    
    # Calculate cog for each species; add species
    cog <- calculateCenterOfGravity(p, speciesList[s])
    
    # Linear regressions between central tendency (lat and month) and environmental variables
    xaxes = c("avgsst", "avgssh", "avgsal", "avgspice")
    yaxes = c("cog.lat", "cog.month")
    xaxislabels = c("Average sea surface \n temperature [\u00b0C]", "Average sea\nsurface height [m]", "Average salinity [ppt]")
    yaxislabels = c("Center of gravity:\nlatitude [\u00b0]", "Center of gravity:\nmonth")
    colnames = c("species", "y", "x", "p", "slope", "r2")
    lms <- as.data.frame(matrix(ncol = length(colnames), nrow = length(xaxes)*length(yaxes)))
    colnames(lms) = colnames
    index = 1
    for (i in 1:length(xaxes)) {
      for (j in 1:length(yaxes)) {
        lm <- lm(data = cog, formula = as.formula(paste0(yaxes[j], "~", xaxes[i]))) %>% summary()
        lms$y[index] = yaxes[j]
        lms$x[index] = xaxes[i]
        lms$p[index] = lm$coefficients[8]
        lms$slope[index] = lm$coefficients[2]
        lms$r2[index] = lm$adj.r.squared
        lms$species[index] = speciesList[s]
        index = index+1
      }
    } 
    # Row-bind cog datasets
    cog.all <- bind_rows(cog.all, cog)
    lms.all <- bind_rows(lms.all, lms)
  }
  
  cog.all.long <- pivot_longer(cog.all, names_to = "y", cols = c("cog.month", "cog.lat"), values_to = "cog") %>% 
    pivot_longer(., names_to = "x", cols = c("avgsal", "avgssh", "avgsst", "avgspice"), values_to = "average")
  
  cog.all.long <- merge(cog.all.long, lms.all)

  cog.all <- pivot_wider(cog.all.long, names_from = y, values_from = p, names_prefix = "p_")
  
  pdf(paste0("Figures/Multispecies_COG_month_vs_latitude_.pdf"), width = 9, height = 7)
  print(ggplot(data = subset(cog.all.long, p < 0.05), mapping = aes(x = average, y = cog, col = species, fill = species)) +
    facet_grid(rows = vars(y), cols = vars(x), scales = "free", switch = "y") +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_fill_manual(values = multispeciesColors) +
    scale_color_manual(values = multispeciesColors))
  dev.off()
  
  
  # # Plot
  # plotlist <- list()
  # index = 1
  # for (j in 1:length(yaxes)) {
  #   for (i in 1:length(xaxes)) {
  #     
  #     annotations1 <- data.frame(
  #       xpos = c(-Inf), ypos =  c(-Inf), #left-bottom
  #       annotateText = paste0("p = ", round(lms$p[index], digits = 2), " , R2 = ", round(lms$r2[index], digits = 2)),
  #       hjustvar = c(-.3),   #shifts bottom left 'Text' to the right; make more negative to move it further right),  #shifts top right 'texT' to the left; make more positive to move it further left
  #       vjustvar = c(-1))    #shifts bottom left 'Text' upward; make more negative to move it further up)
  #     
  #     plotlist[[index]] <- print(ggplot(cog.all, mapping = aes_string(x = xaxes[i], y = yaxes[j])) +
  #                                  geom_point(size = 3, aes(color = species, fill = species)) +
  #                                  geom_smooth(method = "lm", aes(color = species, fill = species)) +
  #                                  theme_classic(base_size = 14) +
  #                                  coord_cartesian(expand = TRUE) + 
  #                                  scale_fill_manual(values = multispeciesColors) +
  #                                  scale_color_manual(values = multispeciesColors) +
  #                                  theme(legend.position = "none") +
  #                                  labs(x = xaxislabels[i], y = yaxislabels[j]))
  #     index = index + 1
  #   }
  # }
  # 
  # # Create legend
  # legendPlot <- ggplot(data = cog.all, aes(x = cog.lat, y = avgsst, col = species, fill = species)) +
  #   geom_point() + 
  #   geom_smooth() +
  #   scale_color_manual(values = multispeciesColors) +
  #   scale_fill_manual(values = multispeciesColors)
  # 
  # # Create some space to the left of the legend
  # legend <- get_legend(legendPlot + theme(legend.box.margin = margin(0, 0, 0, 12)))
  # 
  # Create plot
  # prow = plot_grid(plotlist = plotlist,  nrow = 2, byrow = F)
  # pdf(paste0("Figures/Multispecies_COG_month_vs_latitude_.pdf"), width = 12, height = 10)
  # print(plot_grid(prow, legend, rel_widths = c(3, 0.4), ncol = 2))
  # dev.off()
  
  # Save linear regression dataframe
  write.xlsx(lms.all, file = "Results/Multispecies_COG_lms.xlsx")
  
}
# FIG #: EMPIRICAL DATA --------------------------------------------------------
empiricalData <- function(response, species) {
  
  # Empirical catch by timeblock map
  pdf(paste0("C://KDale/Projects/Phenology/Figures/", species, "/logN1_Empirical.pdf"), width = 7, height = 4)
  print(ggplot(northAmerica) +
          facet_wrap(~timeblock, nrow = 1) +
          geom_point(data = subset(data, logN1 == 0), aes(X*1000, Y*1000), col = "gray80", size = 0.2, pch = 20, alpha = 0.5, inherit.aes = FALSE) +
          geom_point(data = subset(data, logN1 > 0), aes(X*1000, Y*1000, col = logN1), size = 1.5, alpha = 1, pch = 20, inherit.aes = FALSE) +
          geom_sf(fill = "gray20") +
          scale_color_gradient(low = "lightgoldenrod2", high = "firebrick4") +
          xlim(min(data$X)*1000-1000, max(data$X)*1000+1000) +
          ylim(min(data$Y)*1000-1000, max(data$Y)*1000+1000) +
          theme_classic(base_size = 8) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude") +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
  # Empirical catch overall map
  pdf(paste0("C://KDale/Projects/Phenology/Figures/", species, "/logN1_Empirical_overall.pdf"), width = 7, height = 4)
  print(ggplot(northAmerica) +
          geom_point(data = subset(data, logN1 == 0), aes(X*1000, Y*1000), col = "gray80", pch = 20, alpha = 0.5, inherit.aes = FALSE) +
          geom_point(data = subset(data, logN1 > 0), aes(X*1000, Y*1000, col = logN1), pch = 20, inherit.aes = FALSE) +
          geom_sf(fill = "gray20") +
          scale_color_gradient(low = "lightgoldenrod2", high = "firebrick") +
          xlim(min(data$X)*1000-1000, max(data$X)*1000+1000) +
          ylim(min(data$Y)*1000-1000, max(data$Y)*1000+1000) +
          theme_classic(base_size = 8) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude") +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
}

# FIG #: PREDICTED -------------------------------------------------------------
predictedData <- function(response, species, modelType) {
  
  # Predicted catch across months (map)
  predictedCatch.grid.month <- predictedCatch.grid %>% group_by_at(., c("gridid", "X", "Y", "month")) %>% 
    summarize(., mean_est = mean(mean_est))
  
  predictedCatch.grid.timeblock <- predictedCatch.grid %>% group_by_at(., c("gridid", "X", "Y", "timeblock")) %>% 
    summarize(., mean_est = mean(mean_est))
  
  # Predicted catch abundance for each month, across all years
  pdf(file = paste0("Figures/", species, "/Predicted_", modelType,".pdf"), width = 7, height = 4)
  print(ggplot(northAmerica) +
          facet_wrap(~month, nrow = 1) +
          geom_point(data = predictedCatch.grid.month, mapping = aes(x = (X*1000), y = (Y*1000), col = mean_est),  size = 1) +
          geom_sf(fill = "gray20") +
          scale_color_gradient(response, low = "lightgoldenrod2", high =  "firebrick") +
          xlim(min(data$X)*1000-1000, max(data$X)*1000+1000) +
          ylim(min(data$Y)*1000-1000, max(data$Y)*1000+1000) +
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          theme_classic(base_size = 8) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude"))
  dev.off()
  
  # Predicted density for each timeblock, across all years
  pdf(file = paste0("Figures/", species, "/Predicted_", modelType,".pdf"), width = 7, height = 4)
  print(ggplot(northAmerica) +
          facet_wrap(~timeblock, nrow = 1) +
          geom_point(data = predictedCatch.grid.timeblock, mapping = aes(x = (X*1000), y = (Y*1000), col = mean_est),  size = 1) +
          geom_sf(fill = "gray20") +
          scale_color_gradient(response, low = "lightgoldenrod2", high =  "firebrick") +
          xlim(min(data$X)*1000-1000, max(data$X)*1000+1000) +
          ylim(min(data$Y)*1000-1000, max(data$Y)*1000+1000) +
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          theme_classic(base_size = 8) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude"))
  dev.off()
  
  # Predicted density for each timeblock on the original data
  predictedCatch.original.timeblock <- p.original %>% 
    ungroup() %>% group_by_at(c("X", "Y", "timeblock")) %>% 
    summarize(est_retransform = mean(est_retransform)) # average across each lat/lon/timeblock combination
  
  # Prediction on original stations
  pdf(file = paste0("Figures/", species, "/Predicted_original", modelType,".pdf"), width = 7, height = 4)
  print(ggplot(northAmerica) +
          facet_wrap(~timeblock, nrow = 1) +
          geom_point(data = subset(p.original, logN1 < 0.001), aes(X*1000, Y*1000), col = "gray80", size = 0.2, pch = 20, alpha = 0.5, inherit.aes = FALSE) +
          geom_point(data = subset(p.original, logN1 > 0.001), aes(X*1000, Y*1000, col = logN1), size = 1.5, alpha = 1, pch = 20, inherit.aes = FALSE) +
          geom_sf(fill = "gray20") +
          scale_color_gradient(response, low = "lightgoldenrod2", high =  "firebrick4") +
          xlim(min(data$X)*1000-1000, max(data$X)*1000+1000) +
          ylim(min(data$Y)*1000-1000, max(data$Y)*1000+1000) +
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          theme_classic(base_size = 8) +
          theme(axis.text.x = element_text(angle = 60,hjust=1)) +
          labs(x = "Longitude", y = "Latitude"))
  dev.off()
  
  # Plot empirical vs predicted data (only on original stations)
  if (response == "logN1") {
    pdf(file = paste0("Figures/", species, "/Empirical_vs_predicted_", modelType, ".pdf"), width = 4, height = 3)
    print(ggplot(p.original, aes(x = logN1, y = est_retransform)) +
            geom_point() +
            geom_smooth(method = "lm", fill = "gray50", color = "gray50") +
            geom_abline(slope = 1, intercept = 0, lty = 2, color = "black") +
            labs(x = "Empirical abundance (log(N+1))", y = "Predicted abundance (log(N+1))"))
    dev.off()
  } else if (response == "presence") {
    pdf(file = paste0("Figures/", species, "/Empirical_vs_predicted_", modelType, ".pdf"), width = 4, height = 3)
    print(ggplot(p.original) +
            geom_boxplot(aes(x = as.factor(presence), y = est_retransform)) +
            labs(x = "Empirical presence", y = "Predicted presence"))
    dev.off()
  }
  
  
}

# PREDICTIVE PERFORMANCE -------------------------------------------------------

predictivePerformance <- function(fit, response, p.original) {
  
  if (response == "presence") {
    
    #Quantifying predictive performance for presence/absence model
    speciesData$p <- predict(fit)$est
    rocr <- ROCR::prediction(exp(speciesData$p), speciesData$presence)
    ROCR::performance(rocr, measure = "auc")@y.values[[1]] # Values near 0.5 ~ random, want close to 1
    
  } else if (response == "logN1") {
    
    # Spearman correlation coefficient - test if estimated values are close to predicted values
    print("Spearman correlation coefficient: ")
    cor.test(p.original$est_retransform, p.original$logN1) %>% print()
    
    # RAE - realized absolute error for the model - values should be near 0 ideally
    print("Realized absolute error (values should be near 0):")
    print(sum(abs(p.original$logN1 - p.original$est_retransform)) / sum(abs(p.original$logN1- mean(p.original$logN1))))
    
  }
  
  sanity(fit)
}

# LINEAR REGRESSIONS -----------------------------------------------------------
pointlinearRegressions <- function(species) {
  
  # Linear regressions for each grid point for COG and environmental variables
  for (i in 1:length(grid.df$gridid)) {
    data = subset(centralTendency.grid.year, gridid == grid.df$gridid[i])
    
    if (nrow(data) == 0) { next }
    model <- lm(data = data, formula = centralTendency ~ year) %>% summary()
    grid.df$cog_slope[i] = model$coefficients[2]
    grid.df$cog_pvalue[i] = model$coefficients[8]
  }
  # Regressions between mean estimated abundance and environmental variables, per grid cell
  predictedCatch.grid.year <- predictedCatch.grid %>%
    group_by_at(., c("gridid", "year", "latitude", "longitude","X", "Y")) %>%
    summarize(., mean_est = mean(mean_est),
              mean_sst = mean(sst_roms),
              mean_salinity = mean(salinity_roms),
              mean_ssh = mean(ssh_roms),
              var_sst = var(sst_roms),
              var_salinity = var(salinity_roms),
              var_ssh = var(ssh_roms))
  
  for (i in 1:length(grid.df$gridid)) {
    
    data <- subset(predictedCatch.grid.year, gridid == grid.df$gridid[i])
    
    model <- lm(data = data, log(mean_est) ~ (mean_sst)) %>% summary()
    grid.df$sst_pvalue[i] = model$coefficients[8]
    grid.df$sst_slope[i] = model$coefficients[2]
    
    model <- lm(data = data, log(mean_est) ~ (mean_ssh))%>% summary()
    grid.df$ssh_pvalue[i] = model$coefficients[8]
    grid.df$ssh_slope[i] = model$coefficients[2]
    
    model <- lm(data = data, log(mean_est) ~ (mean_salinity))%>% summary()
    grid.df$salinity_pvalue[i] = model$coefficients[8]
    grid.df$salinity_slope[i] = model$coefficients[2]
    
    model <- lm(data = data, log(mean_est) ~ (var_sst))%>% summary()
    grid.df$sst_var_pvalue[i] = model$coefficients[8]
    grid.df$sst_var_slope[i] = model$coefficients[2]
    
    model <- lm(data = data, log(mean_est) ~ (var_ssh))%>% summary()
    grid.df$ssh_var_pvalue[i] = model$coefficients[8]
    grid.df$ssh_var_slope[i] = model$coefficients[2]
    
    model <- lm(data = data, log(mean_est) ~ (var_salinity))%>% summary()
    grid.df$salinity_var_pvalue[i] = model$coefficients[8]
    grid.df$salinity_var_slope[i] = model$coefficients[2]
  }
  
  return(grid.df)
}
#-----------------------------------------------------------------------------#  
regionlinearRegressions <- function(species, formula) {
  
  regressionResults = data.frame(species = species, formula = formula, region = unique(data$region))
  
  # Linear regression for regions across years
  for (i in 1:nrow(regressionResults)) {
    model <- lm(data = subset(centralTendency.region.year, region == regressionResults$region[i]), centralTendency ~ year ) %>% summary()
    regressionResults$slope_year[i] = model$coefficients[2]
    regressionResults$p_year[i] = model$coefficients[8]
  }
  
  # Linear regression for regions across timeblocks
  for (i in 1:nrow(regressionResults)) {
    subset <- subset(centralTendency.region.timeblock, region == regressionResults$region[i])
    subset$timeblock_factor = as.numeric(subset$timeblock)
    if (nrow(subset) == 1) next
    model <- lm(data = subset, centralTendency ~ timeblock_factor) %>% summary()
    regressionResults$slope_timeblock[i] = model$coefficients[2]
    regressionResults$p_timeblock[i] = model$coefficients[8]
  }
  
  return(regressionResults)
}


# ENVIRO POSTERIOR ----------------------------------------------------
envPosterior <- function() {
  
  data_fit <- tidy(fit,"fixed",conf.int = TRUE) %>% subset(term != '(Intercept)')
  
  pdf(paste0("Figures/",species, "/EnvCovariate_Posteriors_", modelType,".pdf"), width = 5, height = 6)
  print(ggplot(data_fit) +
          aes(x = estimate, y = term) +
          geom_vline(xintercept = 0, lty = 2, color = "gray40") +
          geom_segment(aes(x = conf.low, xend = conf.high, yend = term), linewidth = 1, show.legend = FALSE) +
          geom_point(size = 3, pch = 21, stroke = 1, fill = "black") +
          labs(x = "Coefficient Estimate", y = "", title = bquote(~italic(.(species)))) +
          theme_classic(base_size = 12) +
          #scale_y_discrete(labels = c("Gear - Cobb MWT", "Gear - Manta", "Sea surface height [m]","Distance from shore [km]", "Salinity [psu]", "Distance from shore [km]", "Sea surface height [m]", "Sea surface temperature [\u00b0C]")) +
          theme(
            axis.line.x = element_line(color = "gray40", linewidth = 0.2),
            axis.line.y = element_line(color = "gray40", linewidth = 0.2),
            axis.title = element_text(color = "gray40"),
            panel.border = element_rect(colour = "gray10", fill=NA, size=2),
            strip.text.x =  element_text(face = "bold", hjust = 0.5),
            strip.background = element_blank(),
            strip.placement = "none",
            axis.text = element_text(size  = 12),
            strip.text = element_text(size = 12),
            panel.grid = element_blank(),
            panel.spacing = unit(1, "lines"),
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.ticks = element_line(linewidth = 1),
            plot.margin = margin(r = 10, b = 5,  t = 5),
            legend.margin = margin(t=-0,l=0,b=0,r=0.2, unit='cm'),
            legend.background = element_rect(size = 0.5)))
  dev.off()
  
}

# FIG #: ENVIRO TRACKING FIGURES -----------------------------------------------
envTracking_figures <- function(species, modelType) {
  
  plotlist <- list()
  
  variables = c("sst", "ssh", "salinity")
  variableNames = c("Sea surface\ntemperature", "Sea surface height", "Salinity")
  
  # Get scale limits from the significant data subset
  minSlope = min(bind_cols(grid.df$sst_slope, grid.df$ssh_slope, grid.df$salinity_slope))
  maxSlope = max(bind_cols(grid.df$sst_slope, grid.df$ssh_slope, grid.df$salinity_slope))
  
  for (i in 1:length(variables)) {
    
    slopeName = paste0(variables[i], "_slope")
    pName = paste0(variables[i], "_pvalue")
    
    data = subset(grid.df, eval(as.symbol(pName)) < 0.05)
    data.ns = subset(grid.df, eval(as.symbol(pName)) >= 0.05)
    plotlist[[i]] <-  
      print(ggplot(northAmerica) +
              geom_sf() +
              geom_point(data = data, mapping = aes(x = (X*1000), y = (Y*1000), col = .data[[slopeName]]), size = 2, alpha = 1) +
              scale_color_gradient2(mid = "gray70", na.value = "red", midpoint = 0, low = "deepskyblue2", high = "tomato3", limits = c(minSlope, maxSlope)) +
              geom_point(data = data.ns, mapping = aes(x = (X*1000), y = (Y*1000)), col = "gray90", size = 1, alpha = 1) +
              xlim(min(prediction_grid_roms$X)*1000-1000, max(prediction_grid_roms$X)*1000+1000) +
              ylim(min(prediction_grid_roms$Y)*1000-1000, max(prediction_grid_roms$Y)*1000+1000) +
              labs(x = "", y = "") +
              theme_classic(base_size = 11) +
              theme(legend.position = "bottom", axis.text.x = element_text(angle = 60,hjust=1)) +
              ggtitle(bquote(~italic(.(species)))) +
              ggtitle(paste0(variableNames[i], ": Mean")))
    
  }
  
  legendPlot <- ggplot() +
    geom_point(data = grid.df, mapping = aes(x = (X*1000), y = (Y*1000), col = .data[[slopeName]]), size = 2) + #
    scale_color_gradient2("Slope", na.value = "red", mid = "gray70", midpoint = 0, low = "deepskyblue2", high =  "tomato3", limits = c(minSlope, maxSlope))
  
  legend <- get_legend(legendPlot + theme(legend.box.margin = margin(0, 0, 0, 12)))
  # create some space to the left of the legend
  
  prow = plot_grid(plotlist = plotlist,  align = 'vh', nrow = 1)
  pdf(paste0("Figures/",species, "/EnvTracking.significant.changes_", modelType,".pdf"), width = 10, height = 6)
  print(plot_grid(prow, legend, rel_widths =  c(3, 0.4), nrow = 1))
  dev.off()
}

# FIG #: ENVIRO VARIANCE FIGURES -----------------------------------------------
envVariance_figures <- function(species, modelType) {
  
  plotlist <- list()
  
  variables = c("sst", "ssh", "salinity")
  variableNames = c("Sea surface\ntemperature", "Sea surface height", "Salinity")
  
  # Get scale limits from the significant data subset
  minSlope = min(bind_cols(grid.df$sst_var_slope, grid.df$ssh_var_slope, grid.df$salinity_var_slope))
  maxSlope = max(bind_cols(grid.df$sst_var_slope, grid.df$ssh_var_slope, grid.df$salinity_var_slope))
  
  for (i in 1:length(variables)) {
    
    # Get names for column headers
    slopeName = paste0(variables[i], "_var_slope")
    pName = paste0(variables[i], "_var_pvalue")
    
    # Subset significant and non-significant regressions
    data = subset(grid.df, eval(as.symbol(pName)) < 0.05)
    data.ns = subset(grid.df, eval(as.symbol(pName)) >= 0.05)
    
    # Map showing slope values at each grid point. Use .data[[<variable>]] to use string as variable
    plotlist[[i]] <-  
      print(ggplot(northAmerica) +
              geom_sf() +
              geom_point(data = data, mapping = aes(x = (X*1000), y = (Y*1000), col = .data[[slopeName]]), size = 2, alpha = 1) +
              scale_color_gradient2(mid = "gray70", na.value = "red", midpoint = 0, low = "deepskyblue2", high = "tomato3", limits = c(minSlope, maxSlope)) +
              geom_point(data = data.ns, mapping = aes(x = (X*1000), y = (Y*1000)), col = "gray90", size = 1, alpha = 1) +
              xlim(min(prediction_grid_roms$X)*1000-1000, max(prediction_grid_roms$X)*1000+1000) +
              ylim(min(prediction_grid_roms$Y)*1000-1000, max(prediction_grid_roms$Y)*1000+1000) +
              labs(x = "", y = "") +
              theme_classic(base_size = 11) +
              theme(legend.position = "bottom", axis.text.x = element_text(angle = 60,hjust=1)) +
              ggtitle(bquote(~italic(.(species)))) +
              ggtitle(paste0(variableNames[i], ": variance")))
    
  }
  
  # Create legend for the plot
  legendPlot <- ggplot() +
    geom_point(data = grid.df, mapping = aes(x = (X*1000), y = (Y*1000), col = .data[[slopeName]]), size = 2) + #
    scale_color_gradient2("Slope", na.value = "red", mid = "gray70", midpoint = 0, low = "deepskyblue2", high =  "tomato3", limits = c(minSlope, maxSlope))
  
  # create some space to the left of the legend
  legend <- get_legend(legendPlot + theme(legend.box.margin = margin(0, 0, 0, 12)))
  
  # Combine plots and save
  prow = plot_grid(plotlist = plotlist,  align = 'vh', nrow = 1)
  pdf(paste0("Figures/",species, "/EnvVariance.significant.changes_", modelType,".pdf"), width = 10, height = 6)
  print(plot_grid(prow, legend, rel_widths =  c(3, 0.4), nrow = 1))
  dev.off()
}

# FIG #: Effects ---------------------------------------------------------------
effectsFigures <- function(species, modelType) {
  
  # Spatial effect per timeblock
  if("epsilon_st" %in% colnames(p)) {
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
  }
  # Fixed effects per timeblock
  pdf(file = paste0("Figures/", species, "/Fixedeffects_map_", modelType,".pdf"), width = 7, height = 4)
  print(ggplot(northAmerica) +
          geom_sf() +
          facet_wrap(~timeblock, nrow = 1) +
          geom_point(data = predictedCatch.grid.timeblock, mapping = aes(x = (X*1000), y = (Y*1000), col = mean_fixed),  size = 2) +
          scale_color_gradient("Fixed effects", low = "gray90", high =  "firebrick") +
          xlim(min(data$X)*1000-1000, max(data$X)*1000+1000) +
          ylim(min(data$Y)*1000-1000, max(data$Y)*1000+1000) +
          theme(axis.text.x = element_text(angle = 60,hjust=1))+
          theme_classic(base_size = 8) +
          ggtitle(bquote(~italic(.(species)))) +
          labs(x = "Longitude", y = "Latitude"))
  dev.off()
  
  # Fixed effects per latitude per timeblock
  pdf(file = paste0("Figures/", species, "/Fixedeffects_", modelType,".pdf"), width = 5, height = 4)
  print(ggplot(predictedCatch.lat.timeblock, mapping = aes(x = latitude, y = mean_fixed, fill = timeblock, col = timeblock)) +
          geom_point() +
          geom_smooth() +
          labs(x = "Latitude", y = "Mean fixed effect") +
          scale_color_brewer(palette = "Spectral") +
          ggtitle(bquote(~italic(.(species)))) +
          scale_fill_brewer(palette = "Spectral"))
  dev.off()
  
  # Random effects per latitude per timeblock
  pdf(file = paste0("Figures/", species, "/Randomeffects_", modelType,".pdf"), width = 5, height = 4)
  print(ggplot(predictedCatch.lat.timeblock, mapping = aes(x = latitude, y = mean_rf, fill = timeblock, col = timeblock)) +
          geom_point() +
          geom_smooth() +
          labs(x = "Latitude", y = "Mean random effect") +
          scale_color_brewer(palette = "Spectral") +
          ggtitle(bquote(~italic(.(species)))) +
          scale_fill_brewer(palette = "Spectral"))
  dev.off()
}

# QUOTIENT CURVES --------------------------------------------------------------
quotientCurves <- function(species, modelType) {
  
  variables = c("sst_roms", "ssh_roms", "salinity_roms")
  xaxisLabels = c("Sea surface \ntemperature (\u00b0C)", "Sea surface \nheight (m)", "Salinity (ppt)")
  digits = c(0,1,0)
  plotlist <- list()
  
  for(i in 1:length(variables)) {
    
    # Create bins
    bins <- seq(floor(min(get("p")[,variables[i]])), ceiling(max(get("p")[,variables[i]])), 2)
    p$bins <- data.frame(bins = round(x = get("p")[,variables[i]], digits = digits[i]))[,1]
    
    fish.per.bin <- group_by_at(p, c("bins")) %>%
      mutate(towID = paste(latitude, longitude, year, month, sep = "_")) %>% 
      dplyr::summarize(., est_retransform = sum(est_retransform), numTows = length(unique(towID))) %>%
      mutate(., propFish = est_retransform/sum(est_retransform), propTows = numTows/sum(numTows)) %>% 
      mutate(quotient = propFish/propTows)
    
    # Plot quotient curves
    plotlist[[i]] <- print(ggplot() + 
                             geom_line(subset(fish.per.bin), mapping = aes(x = bins, y = quotient), linewidth = 1.5) +
                             geom_hline(yintercept = 1, lty = 2) +
                             scale_y_continuous(limits = c(0, 2.5)) +
                             theme_classic(base_size = 12) +
                             labs(x = xaxisLabels[i], y = "Quotient"))
    
    
  }
  # Plot quotient curves
  pdf(file = paste0("Figures/", species, "/QuotientCurves_", modelType,".pdf"), width = 9, height = 4)
  print(cowplot::plot_grid(plotlist = plotlist, align = "hv", nrow = 1))
  dev.off()
  
}
# NICHE BREADTH ----------------------------------------------------------------------------------------------

calculateNicheBreadth <- function(species, modelType) {
  
  predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", modelType, ".rdata")
  if(!file.exists(predictionObjectName)) {
    createPredictionObjects(species = species, modelType = modelType, makeNewGrid = F, makeNewPrediction = T)
  } 
  
  load(predictionObjectName)
  
  variables = c("sst_roms", "ssh_roms", "salinity_roms", "bottom_depth", "distance_from_shore_m")
  xaxisLabels = c("Sea surface \ntemperature (\u00b0C)", "Sea surface \nheight (m)", "Salinity (ppt)", "Bottom \ndepth (m)", "Distance from \nshore (m)")
  
  nicheBreadth = data.frame(species = species, variableNames = xaxisLabels, variable = variables, hurlbert = NA, smith = NA, smith.lower = NA, smith.upper = NA)
  
  for(i in 1:length(variables)) { 
    
    # Create bins
    digits = c(0,1,0,0,0) # number of digits for each variable type
    bins <- seq(floor(min(get("p")[,variables[i]])), ceiling(max(get("p")[,variables[i]])), 2)
    p$bins <- data.frame(bins = round(x = get("p")[,variables[i]], digits = digits[i]))[,1]
    
    fish.per.bin <- group_by_at(p, c("bins")) %>%
      mutate(towID = paste(latitude, longitude, year, month, sep = "_")) %>% 
      dplyr::summarize(., est_retransform = sum(est_retransform), numTows = length(unique(towID))) %>%
      mutate(., propFish = est_retransform/sum(est_retransform), propTows = numTows/sum(numTows)) %>% 
      mutate(quotient = propFish/propTows)
    
    # Levin's measure of niche breadth, standardized to 0-1
    hurlbert = 1/sum(fish.per.bin$propFish ^2 / fish.per.bin$propTows) 
    hurlbert.standardized = (hurlbert-min(fish.per.bin$propTows))/(1-min(fish.per.bin$propTows))
    nicheBreadth$hurlbert[i] = hurlbert.standardized
    # Smith's measure 
    # Varies from 0 (minimal) to 1 (maximal)
    # Trig functions are in radians
    smith = sum(sqrt(fish.per.bin$propFish * fish.per.bin$propTows))
    nicheBreadth$smith[i] = smith
    nicheBreadth$smith.lower[i] = sin(asin(smith) - (1.96/(2*sqrt(sum(fish.per.bin$est_retransform)))))
    nicheBreadth$smith.upper[i] = sin(asin(smith) + (1.96/(2*sqrt(sum(fish.per.bin$est_retransform)))))
    sin(asin(smith) - (1.96/2*sqrt(1000)))
    sin(asin(smith) + (1.96/2*sqrt(1000)))
  }
  
  # Plot niche breadth
  print(ggplot(nicheBreadth) +
          geom_point(aes(x = variableNames, y = smith), size = 4) +
          labs(x = "Environmental covariate", y = "Smith's measure") +
          theme_classic(base_size = 14) +
          geom_errorbar(aes(x = variableNames, ymin = smith.lower, ymax = smith.upper), width = 0.2, linewidth = 1.5))
  
  write.csv(nicheBreadth, file = paste0("Results/", species, "/", species, "_", modelType,  "nicheBreadth.csv"))
  
}


# NICHE BREADTH MULTISPECIES -------------------------------------------------------------------------------#
nicheBreadthMultiSpecies <- function(speciesList, modelTypes) {
  
  nicheBreadthCombined = NA
  recalculate = T
  
  for (i in 1:length(speciesList)) {
    
    # Load in niche breadth csv (or create new one)
    file = paste0("Results/", speciesList[i], "/", speciesList[i], "_", modelTypes[i],  "nicheBreadth.csv")
    if(!file.exists(file) | recalculate == T) {
      calculateNicheBreadth(speciesList[i], modelTypes[i])
    } 
    nicheBreadth = read.csv(file)
    
    # Bind niche breadth datasets together
    if (i == 1) {
      nicheBreadthCombined = nicheBreadth
    } else {
      nicheBreadthCombined = bind_rows(nicheBreadthCombined, nicheBreadth)
    }
    
  }
  
  # Plot niche breadth
  pdf("Figures/Multispecies_NicheBreadth.pdf", width = 6, height = 6)
  print(ggplot(nicheBreadthCombined) +
          geom_point(aes(x = variableNames, y = smith, color = species), size = 4) +
          geom_errorbar(aes(x = variableNames, ymin = smith.lower, ymax = smith.upper, color = species), width = 0.2, linewidth = 1.5) +
          labs(x = "Environmental covariate", y = "Smith's measure") +
          theme_classic(base_size = 14) +
          theme(legend.position = c(1,0.2), legend.justification = "right") +
          scale_color_manual("Species", values = multispeciesColors))
  dev.off()
}

# PLOT MAPS ----------------------------------------------------------------------------
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed()
}
# LEADING TRAILING ------------------------------------------------------------
leadingTrailing <- function(species) {
  
  # Get a list of unique timeblocks and set list of range percentiles
  timeblocks <- unique(p$timeblock)
  rangePercentiles = c("upper10", "middle80", "lower10")
  
  # Set up a dataframe to store info about quartiles
  # quartiles = data.frame(timeblock = unique(p$timeblock), lead.lat = NA, trail.lat = NA, max.lat = NA, min.lat = NA)
  # 
  # for (i in 1:length(timeblocks)) {
  #   
  #   p.subset <- subset(p, timeblock == timeblocks[i] & est_retransform > 0.001)
  #   # Calculate where 10% of the abundance is
  #   totalAbundance = sum(p.subset$est_retransform)
  #   tenPercentLead = totalAbundance * 0.1
  #   tenPercentTrail = totalAbundance - (totalAbundance * 0.1)
  #   p.ordered = p.subset[order(p.subset$latitude, decreasing = T),]
  #   cumulative = as.data.frame(cumsum(p.ordered$est_retransform))
  #   
  #   # Find indexes of the 10% low and high cumulative abundance
  #   lead.index <- which(cumulative$`cumsum(p.ordered$est_retransform)` >= tenPercentLead, arr.ind=TRUE)[1]
  #   trail.index <- which(cumulative$`cumsum(p.ordered$est_retransform)` >= tenPercentTrail, arr.ind=TRUE)[1]
  #   
  #   # Get latitudes for 10% low and high abundances
  #   lead.lat = p.ordered$latitude[lead.index]
  #   trail.lat = p.ordered$latitude[trail.index]
  #   
  #   quartiles$lead.lat[i] = lead.lat
  #   quartiles$trail.lat[i] = trail.lat
  #   quartiles$max.lat[i] = max(p.subset$latitude)
  #   quartiles$min.lat[i] = min(p.subset$latitude)
  #   
  # }
  
  # print(ggplot(quartiles) +
  #         geom_linerange(aes(ymin = min.lat, ymax = max.lat, x= timeblock), linewidth = 5, col = "lightblue") +
  #         geom_linerange(aes(ymin = trail.lat, ymax = lead.lat, x= timeblock), linewidth = 5) +
  #         geom_hline(yintercept =  max(p$latitude), lty = 2) +
  #         geom_hline(yintercept =  min(p$latitude), lty = 2)) 
  
  # Calculate top 10% and bottom 10% of range
  maxLat = max(p$latitude) 
  minLat = min(p$latitude)
  latRange10pct = (maxLat - minLat) * 0.1
  
  # Subset data into three sections to add range percentile
  upper10 = subset(p, latitude >= (maxLat-latRange10pct)) %>% mutate(rangePercentile = "upper10")
  lower10 = subset(p, latitude <= (minLat+latRange10pct)) %>% mutate(rangePercentile = "lower10")
  middle80 = subset(p, latitude < (maxLat-latRange10pct) & latitude > (minLat+latRange10pct)) %>% mutate(rangePercentile = "middle80")
  
  # Rebind sections
  p.dat <- rbind(upper10, lower10, middle80)
  
  # Summarize across sections and years
  p.rangePercentiles = p.dat %>% group_by_at(c("rangePercentile", "year")) %>%
    summarize(est_retransform = mean(est_retransform))
  
  # Order sections
  p.rangePercentiles$rangePercentile = factor(x = p.rangePercentiles$rangePercentile, levels = c("upper10", "middle80", "lower10"))
  
  # Set up dataframe for chi square tests
  observedExpected <- p.dat %>% group_by_at(., c("timeblock", "rangePercentile")) %>% 
    summarize(observed = sum(est_retransform)) %>% 
    ungroup %>% group_by_at(c( "timeblock")) %>% 
    mutate(totalAbundanceYear = sum(observed), expected = NA)
  
  # Calculate expected values
  for(i in 1:nrow(observedExpected)) {
    if(observedExpected$rangePercentile[i] == "upper10") {
      observedExpected$expected[i] = observedExpected$totalAbundanceYear[i] * nrow(upper10)/nrow(p)
    } else if (observedExpected$rangePercentile[i] == "lower10") {
      observedExpected$expected[i] = observedExpected$totalAbundanceYear[i] * nrow(lower10)/nrow(p)
    } else {
      observedExpected$expected[i] = observedExpected$totalAbundanceYear[i] * nrow(middle80)/nrow(p)
    }
  }
  
  # Run chi square tests
  for (i in 1:length(timeblocks)) {
    dat <- subset(observedExpected, timeblock == timeblocks[i])
    print(chisq.test(dat$expected, dat$observed))
  }
  
  for (i in 1:3) {
    dat <- subset(observedExpected, rangePercentile == rangePercentiles[i])
    print(chisq.test(dat$expected, dat$observed))
  }
  
  observedExpected.edges <- subset(observedExpected, rangePercentile != "middle80")
  chisq.test(observedExpected.edges$observed, observedExpected.edges$expected)
  
  # Chi square scatterplot
  pdf(file = paste0("Figures/", species, "/ChiSquare_", modelType,".pdf"), width = 5, height = 4)
  print(ggplot(observedExpected.edges, aes(x = expected, y = observed, color = rangePercentile)) +
          geom_abline(slope = 1, intercept = 0, linewidth = 0.75)+
          geom_point(size = 3) +
          #geom_text(data =  as.data.frame(chisquared$p.value), aes(label = chisquared$p.value)) +
          scale_color_manual("Range", values = c("upper10" = "midnightblue", "middle80" = "steelblue4", "lower10" = "skyblue1")) +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
  edges <- subset(p.rangePercentiles, rangePercentile != "middle80")
  
  # Linear regression of abundance in each area versus year
  pdf(file = paste0("Figures/", species, "/AbundanceRegressions_", modelType,".pdf"), width = 5, height = 4)
  print(ggplot(edges, mapping = aes(x = year, y = est_retransform, color = rangePercentile, fill = rangePercentile), size = 2.5) +
          geom_line(linewidth = 1) +
          #geom_smooth() +
          labs(x = "Year", y = "Estimated abundance log(N+1)") +
          scale_fill_manual(values = c("midnightblue", "steelblue4", "skyblue1")) +
          scale_color_manual("Range", values = c("upper10" = "midnightblue", "middle80" = "steelblue4", "lower10" = "skyblue1")) +
          ggtitle(bquote(~italic(.(species)))))
  dev.off()
  
}
# CUMULATIVE CATCH ------------------------------------------
cumulativeCatch <- function(species) {
  
  predictedCatch.time.cumulative <- predictedCatch.time %>% ungroup() %>% group_by_at(c("year", "month")) %>%
    summarize(mean_est_retransform = mean(mean_est_retransform)) %>% 
    mutate(cum_catch = cumsum(mean_est_retransform))
  
  cumulativeAbundance <- predictedCatch.time.cumulative %>% 
    summarize(fifteen = max(cum_catch)*0.15, fifty = max(cum_catch)*0.5, eightyfive = max(cum_catch)*0.85)
  
  timesteps = unique(cumulativeAbundance$year)
  
  for (i in 1:length(timesteps)) {
    
    # Subset timesteps
    timeSubset = subset(predictedCatch.time.cumulative, year == timesteps[i])
    
    # Find month at which 15%, 50%, and 85% are passed
    fifteen.month = max(timeSubset$month[timeSubset$cum_catch < cumulativeAbundance$fifteen[i]])
    slope = (timeSubset$cum_catch[fifteen.month]-timeSubset$cum_catch[fifteen.month-1])/(1) # differences between months are always 1
    b = timeSubset$cum_catch[fifteen.month] - slope * fifteen.month #y-mx
    cumulativeAbundance$fifteen.month[i] = (cumulativeAbundance$fifteen[i] - b) / slope #y-b/m
    
    fifty.month = max(timeSubset$month[timeSubset$cum_catch < cumulativeAbundance$fifty[i]])
    slope = (timeSubset$cum_catch[fifty.month]-timeSubset$cum_catch[fifty.month-1])/(1) # differences between months are always 1
    b = timeSubset$cum_catch[fifty.month] - slope * fifty.month #y-mx
    cumulativeAbundance$fifty.month[i] = (cumulativeAbundance$fifty[i] - b) / slope #y-b/m
    
    eightyfive.month = max(timeSubset$month[timeSubset$cum_catch < cumulativeAbundance$eightyfive[i]])
    slope = (timeSubset$cum_catch[eightyfive.month]-timeSubset$cum_catch[eightyfive.month-1])/(1) # differences between months are always 1
    b = timeSubset$cum_catch[eightyfive.month] - slope * eightyfive.month #y-mx
    cumulativeAbundance$eightyfive.month[i] = (cumulativeAbundance$eightyfive[i] - b) / slope #y-b/m
    
  }
  
  pdf(file = paste0("Figures/", species, "/", "Cumulative_Thresholds_.pdf"))
  print(ggplot(cumulativeAbundance) +
          geom_errorbarh(mapping = aes(y = year, xmin = fifteen.month, xmax = eightyfive.month))+
          geom_point(mapping = aes(y = year, x = fifty.month)) +
          labs(x = "Month", y = "Year") +
          scale_x_continuous(n.breaks = 12)+
          ggtitle(bquote(~italic(.(species)))) +
          scale_y_continuous(n.breaks = length(timesteps)))
  dev.off()
}
# CREATE PREDICTION OBJECTS---------------------------------------------------------------
createPredictionObjects <- function(species, modelType, makeNewGrid, makeNewPrediction) {
  
  # Construct filenames
  fitName = paste0("Results/", species, "/Models/", modelType, ".rdata")
  predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", modelType, ".rdata")
  gridFilename = paste0("Analysis/PredictionGrids/", species, "_grid.rdata")
  
  # Load model results
  if(file.exists(fitName)) {
    load(file = fitName)
  } else {
    print("Model object does not exist for", species, "and", modelType, "-- run modelsdmTMB.R")
    break
  }
  
  # If the prediction grid hasn't already been created, create it (or override)
  if (!file.exists(gridFilename) | makeNewGrid == T) { # Load prediction grid
    source("Analysis/Code/CreatePredictionGrid.R")
    createPredictionGrid(data = data, species = species, path = paste0("Analysis/PredictionGrids/", species, "_grid.rdata"))
  }
  
  # Then load prediction grid (and grid.df) from source
  load(gridFilename) 
  
  # Create a row for each gear type (these will ultimately be averaged together)
  gearGeneral = unique(data$gearGeneral) # Get gear categories
  prediction_grid_roms <- expand_grid(prediction_grid_roms, gearGeneral) %>% 
    subset(., !sst_roms == 0 | !salinity_roms == 0, !ssh_roms == 0)
  
  # Create factored versions of year and timeblock
  prediction_grid_roms$year_scaled = as.vector(scale(prediction_grid_roms$year, center = T, scale = T)[,1])
  prediction_grid_roms$timeblock = factor(prediction_grid_roms$timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019"))
  
  # Predict on both original data and on the prediction grid
  p.original <- predict(object = fit, newdata = data) %>%  
    mutate(est_retransform = fit$family$linkinv(est)) %>% 
    group_by_at(c("latitude", "longitude", "X", "Y", "year", "month", "timeblock", "sst_roms", "ssh_roms", "salinity_roms", "bottom_depth", "distance_from_shore_m", "spice")) %>% 
    summarize(est_retransform = sum(est_retransform), logN1 = sum(logN1))
  
  # Predict on full grid
  p.obj <- predict(object = fit, newdata = prediction_grid_roms, return_tmb_object = T)
  p <- predict(object = fit, newdata = prediction_grid_roms) %>% 
    mutate(towID = paste(latitude, longitude, year, month, sep = "_")) %>% 
    mutate(est_retransform = fit$family$linkinv(est)) %>% 
    group_by_at(c("towID", "gridid", "region", "latitude", "longitude", "X", "Y", "year", "month", "timeblock", "sst_roms", "ssh_roms", "salinity_roms", "bottom_depth", "distance_from_shore_m", "spice")) %>% 
    summarize(est_retransform = sum(est_retransform), est_rf = mean(est_rf), est_non_rf = mean(est_non_rf)) # sum for est and average for other effects across the three gear types
  
  # Get unique gridid for each grid point using lat/lon combination
  p$latlon <- paste(p$longitude, p$latitude) 
  
  # Save prediction objects
  save(p, p.obj, p.original, file = predictionObjectName)
  
}


# RUN FUNCTIONS -------------------------------------------------------------------------

northAmerica <- read_sf("C://KDale/GIS/North_South_America.shp") %>% st_union() %>% st_transform(., crs = "EPSG:5070")
programColors = c("IMECOCAL" = "firebrick3", "CalCOFI" = "coral3","RREAS" = "darkgoldenrod2", "PRS_juveniles" ="darkseagreen3","PRS_larvae" ="cornflowerblue", "NH-Line" ="deepskyblue3", "Canada" = "darkslateblue","EcoFOCI" ="darkslategray")
regionColors = c("Southern CCE" = "firebrick3", "Central CCE" = "darkgoldenrod2", "OR/WA" ="darkseagreen3","British Columbia" ="cornflowerblue", "British Columbia" = "darkslateblue","Gulf of Alaska" ="darkslategray")
multispeciesColors = c("tan3", "skyblue2","navy")
theme_set(theme_classic(base_size = 12))
timeblocks <- read_xlsx("Data/timeblocks.xlsx", sheet = 1)

tradeoffs <- read_xlsx("Results/Tradeoffs.xlsx", sheet = 1) %>% subset(., chosen_model == "x")

# Mixed taxonomic models
speciesList = tradeoffs$scientific_name
modelTypes = tradeoffs$file_name

# Single species
#speciesList = "Parophrys vetulus"
#modelTypes = c("logN1_speciesRange_allPrograms_base_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)")

makeNewPrediction = T
makeNewGrid = F

# Create prediction objects-----
  for (i in 1:length(modelTypes)) {
    
    # Set species and model
    species = speciesList[i]
    modelType = modelTypes[i]
    
    # Filenames
    fitName = paste0("Results/", species, "/Models/", modelType, ".rdata")
    predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", modelType, ".rdata")
    gridFilename = paste0("Analysis/PredictionGrids/", species, "_grid.rdata")
    
    # Create prediction objects
    if(!file.exists(predictionObjectName) | makeNewPrediction == T) {
      createPredictionObjects(
        species = species,
        modelType = modelType,
        makeNewGrid = makeNewGrid, # Change these to force a new grid or new prediction
        makeNewPrediction = makeNewPrediction
      )
    } 
    
    # Load objects
    load(fitName)
    load(predictionObjectName)
    load(gridFilename)
  }


for (j in 1:length(modelTypes)) {
  
  # Set species and model
  species = speciesList[j]
  modelType = modelTypes[j]
  
  # Filenames
  fitName = paste0("Results/", species, "/Models/", modelType, ".rdata")
  predictionObjectName = paste0("Results/", species, "/Models/prediction_objects_", modelType, ".rdata")
  gridFilename = paste0("Analysis/PredictionGrids/", species, "_grid.rdata")
  
  # Load objects
  load(fitName) # data and model
  load(predictionObjectName) # prediction objects
  load(gridFilename) # grid.df
  
  # Get response variable
  response = strsplit(modelType, "_")[[1]][1]
  
  # Predictive performance -----
  predictivePerformance(fit, response, p.original)
  
  # CT dataframes ------------
  # Summarize catch data by grid point, month, year, region, timeblock
  predictedCatch.grid <- p %>% group_by_at(., c("region", "gridid", "month", "year", "timeblock", "sst_roms", "latitude", "longitude","X", "Y", "salinity_roms", "ssh_roms")) %>%
    dplyr::summarize(., mean_est = mean(est_retransform)) %>% # summarize across grid cells
    mutate(., est_x_month = mean_est * month) %>% 
    mutate(., timeblock = factor(timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019"))) %>% 
    mutate(., region = factor(region, levels = c("Gulf of Alaska", "British Columbia", "OR/WA", "Central CCE", "Southern CCE")))
  
  predictedCatch.time <- p %>% group_by_at(., c("month", "year", "timeblock")) %>%
    dplyr::summarize(., mean_est_retransform = mean(est_retransform)) %>% # summarize across grid cells
    mutate(., est_x_month = mean_est_retransform * month) %>% 
    mutate(., timeblock = factor(timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")))
  
  if("epsilon_st" %in% colnames(p)) { # Spatial models
    # Summarize catch data by grid point, month, timeblock
    predictedCatch.grid.timeblock <- p %>% group_by_at(., c("gridid", "timeblock",  "latitude", "longitude","X", "Y")) %>%
      summarize(.,
                mean_est = mean(est_retransform),
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
                mean_est = mean(est_retransform),
                mean_rf = mean(fit$family$linkinv(est_rf)),
                mean_fixed = mean(fit$family$linkinv(est_non_rf)),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      )
  } else {
    predictedCatch.grid.timeblock <- p %>% group_by_at(., c("gridid", "timeblock", "latitude", "longitude","X", "Y")) %>%
      summarize(.,
                mean_est = mean(est_retransform),
                mean_rf = mean(fit$family$linkinv(est_rf)),
                mean_fixed = mean(fit$family$linkinv(est_non_rf)),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      )
    # Summarize catch data by latitude, month, timeblock
    predictedCatch.lat.timeblock <- p %>% group_by_at(., c("latitude", "timeblock")) %>% 
      summarize(.,
                mean_est = mean(est_retransform),
                mean_rf = mean(fit$family$linkinv(est_rf)),
                mean_fixed = mean(fit$family$linkinv(est_non_rf)),
                mean_sst = mean(sst_roms),
                mean_salinity = mean(salinity_roms),
                mean_ssh = mean(ssh_roms)
      )
  }
  # calculate central tendency for each grid cell (month # * mean abundance in month)/(sum of mean abundances across all months)
  centralTendency.grid <- group_by_at(predictedCatch.grid, c("gridid","latitude", "longitude", "X", "Y")) %>% 
    summarize(mean_est_grid = mean(mean_est), centralTendency = sum(est_x_month)/sum(mean_est))
  
  # Calculate central tendency for each grid point per timeblock
  centralTendency.grid.timeblock <- group_by_at(predictedCatch.grid, c("gridid","timeblock", "latitude", "longitude", "X", "Y")) %>% 
    summarize(mean_est_grid = mean(mean_est), centralTendency = sum(est_x_month)/sum(mean_est))
  
  # Calculate central tendency for each grid point per year
  centralTendency.grid.year <- group_by_at(predictedCatch.grid, c("gridid","year",  "latitude", "longitude", "X", "Y")) %>% 
    summarize(mean_est_grid = mean(mean_est), centralTendency = sum(est_x_month)/sum(mean_est))
  
  # Calculate central tendency for each region per timeblock
  centralTendency.region.timeblock <- group_by_at(predictedCatch.grid, c("region", "timeblock")) %>% 
    summarize(mean_est_region = mean(mean_est), centralTendency = sum(est_x_month)/sum(mean_est))
  
  # Calculate central tendency for each region per year
  centralTendency.region.year <- group_by_at(predictedCatch.grid, c("region", "year")) %>% 
    summarize(mean_est_region = mean(mean_est), centralTendency = sum(est_x_month)/sum(mean_est))
  
  # Linear regressions -----
  regressionResults <- regionlinearRegressions(species = species, formula = as.character(fit$formula))
  grid.df <- pointlinearRegressions(species = species)
  grid.df$days = grid.df$cog_slope*(365/12) # convert to days
  write.xlsx(regressionResults, file = paste0("Results/", species, "/Regressions/", modelType, "_linearRegressionResults.xlsx"))
  
  # Frequency of sampling for species----
  #frequencyOfSampling(speciesData = data)
  
  # Env tracking figures----
  #envTracking_figures(species = species, modelType = modelType)
  #envVariance_figures(species = species, modelType = modelType)
  
  # Env posterior figures----
  envPosterior()
  
  # Effects figures -----
  effectsFigures(species, modelType = modelType)
  
  # CT figures -------
  centralTendency_figures(response = response, species = species, modelType = modelType)
  
  # Empirical data ---------
  empiricalCatch.timeblock <- data %>%
    group_by_at(c("latitude", "longitude", "X", "Y", "year", "month", "timeblock")) %>%
    summarize(logN1 = sum(logN1)) %>% # sum across the three gear types
    ungroup() %>% group_by_at(c("X", "Y", "timeblock")) %>%
    summarize(logN1 = mean(logN1))
  empiricalData(response = response, species = species)
  
  # Predicted catch or presence --------
  predictedData(response = response, species = species, modelType = modelType)
  
  #Quotient curves----------
  quotientCurves(species, modelType)
  
  #Leading and trailing -------
  leadingTrailing(species)
  
  #Cumulative curves --------
  cumulativeCatch(species)
}


# Central tendency (multi-species plot) --------------
# Mixed taxonomic models
speciesList = c("Engraulis mordax", "Sardinops sagax", "Sebastes jordani", "Citharichthys sordidus", "Glyptocephalus zachirus", "Merluccius productus", "Parophyrus vetulus", "Sebastes paucispinis", "Stenobrachius leucopsarus", "Tarletonbeania crenularis", "Triphoturus mexicanus", "Vinciguerria lucetia")
modelTypes = c("logN1_speciesRange_allPrograms_both_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_geo_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_geo_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_both_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_base_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_both_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_base_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_both_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_both_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_base_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_both_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_geo_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)")


modelTypes = c("logN1_speciesRange_allPrograms_geo_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_base_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)",
               "logN1_speciesRange_allPrograms_base_spice+sst+ssh+salinity+dfs+bd+month_as.factor(gearGeneral)")

centralTendencyMultiSpecies(speciesList, modelTypes)

# Niche breadth (multi-species plot) --------------
nicheBreadthMultiSpecies(speciesList, modelTypes = modelTypes)

# Tradeoffs bar plot
tradeoffs.summary <- tradeoffs %>%  
  group_by_at(c("life_history", "model")) %>% 
  summarize(n = sum(!is.na(chosen_model)))
tradeoffs.summary$model <- factor(tradeoffs.summary$model, levels = c("base", "geo", "pheno", "both"))

pdf(file = "Figures/Tradeoff_Summary.pdf", width = 6, height = 5)
ggplot(tradeoffs.summary) +
  geom_col(aes(x = Model, y = n, fill = Life.history)) +
  labs(x = "Model", y = "N") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.05,0.8), legend.justification = "left") +
  scale_fill_manual("Life history", values = c("Groundfish" = "tan3", "Coastal pelagic" = "skyblue2", "Mesopelagic" = "navy"))
dev.off()

# Test for quadratic relationship between estimated abundance and environmental covariates
mod <- mgcv::gam(fit$family$linkinv(p$est_retransform) ~ s(p$sst_roms, k = 3), method = "REML", family = gaussian())
summary(mod)
ggplot(data = p, aes(y = est_retransform, x = sst_roms)) +
  geom_point() +
  labs(x = "Sea surface temperature (degC)", y = "Estimated abundance log(N+1)")

mod <- mgcv::gam(fit$family$linkinv(p$est_retransform) ~ s(p$salinity_roms, k = 3), method = "REML", family = gaussian())
summary(mod)
ggplot(data = p, aes(y = est_retransform, x = salinity_roms)) +
  geom_point() +
  labs(x = "Sea surface salinity (ppt)", y = "Estimated abundance log(N+1)")

mod <- mgcv::gam(fit$family$linkinv(p$est_retransform) ~ s(p$ssh_roms, k = 3), method = "REML", family = gaussian())
summary(mod)
ggplot(data = p, aes(y = est_retransform, x = ssh_roms)) +
  geom_point() +
  labs(x = "Sea surface height (m)", y = "Estimated abundance log(N+1)")


#---------------------------------------------------------------------------#

ggplot(predictedCatch.time) +
  geom_line(mapping = aes(x = month, y = cum_))

# View confidence intervals and extract parameters as a dataframe
tidy(fit, effects = "ran_pars", conf.int = TRUE)

# Plot smooother on variables along the response scale
# Doesn't work for geography model
visreg::visreg(fit, xvar = "sst_scaled", scale = "response")
visreg::visreg(fit, xvar = "ssh_scaled", scale = "response")
visreg::visreg(fit, xvar = "salinity_scaled", scale = "response")

ggeffects::ggeffect(fit,  "sst_scaled") %>% plot()

# Center of gravity --------
# works on object returned from predict()
cog = get_cog(p.obj, format = "long") # Does not work on prediction without an sdmTMB object returned
cogx = subset(cog, coord == "X") %>% dplyr::rename(., upr_x = upr, lwr_x = lwr, est_x = est, se_x = se)
cogy = subset(cog, coord == "Y") %>% dplyr::rename(., upr_y = upr, lwr_y = lwr, est_y = est, se_y = se)
cog = merge(cogx, cogy, by = "timeblock")

# Line plot for COG (month)
cog$month = as.factor(cog$month)
ggplot(cog, aes(est_x, est_y, color = month)) +
  geom_linerange(aes(xmin = lwr_x, xmax = upr_x), lwd = 1) +
  geom_linerange(aes(ymin = lwr_y, ymax = upr_y), lwd = 1) + 
  scale_color_discrete("Month", type = c("black", "#4421af", "#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#76c68f", "#a6d75b", "#c9e52f", "#d0ee11", "#d0f400", "#ffee65")) +
  theme_classic() +
  labs(x = "Latitude \n coordinates", y = "Longitude coordinates") +
  theme(axis.text.x = element_text(angle = 60,hjust=1))+
  coord_equal()

# Line plot for COG (timeblock)
cog$timeblock = as.factor(cog$timeblock)
ggplot(cog, aes(est_x, est_y, color = timeblock)) +
  geom_linerange(aes(xmin = lwr_x, xmax = upr_x), linewidth = 1) +
  geom_linerange(aes(ymin = lwr_y, ymax = upr_y), linewidth = 1) + 
  scale_color_brewer(palette = "Spectral") +
  theme_classic() +
  labs(x = "Latitude \n coordinates", y = "Longitude coordinates") +
  theme(axis.text.x = element_text(angle = 60,hjust=1))+
  coord_equal()

# SANDBOX------------------------------------------------------------------------------------------
# Plot smooother on variables in link space with randomized quantile partial residuals
# Takes a few min to run. Requires mesh.
fit <- geo

visreg(fit, xvar = "sst_scaled")
visreg(fit, xvar = "ssh_scaled")
visreg(fit, xvar = "salinity_scaled")

variables = c("sst_roms", "ssh_roms", "salinity_roms")
variables = c("latitude")

for (i in 1:length(variables)) {
  
  # Weight environmental covariates by biomass
  weighted_var <- p %>% 
    group_by(year) %>% 
    reframe("Median" = hutils::weighted_quantile(v = .data[[variables[i]]], w = est_retransform, p = 0.5),
            "1st" = hutils::weighted_quantile(v = .data[[variables[i]]], w = est_retransform, p = 0.25),
            "3rd" = hutils::weighted_quantile(v = .data[[variables[i]]], w = est_retransform, p = 0.75)) %>% 
    pivot_longer(cols = c("Median", "1st", "3rd"),
                 names_to = "series", values_to = "var") %>% 
    group_by(series) %>% 
    mutate(var_centered = var - mean(var)) %>% 
    ungroup()
  
  # Plot weighted values
  print(ggplot(weighted_var, aes(year, var, color = series, group = series, fill = series)) +
          stat_smooth(method = "gam", formula = y ~ s(x, k = 4), se = FALSE, size = 0.8) +
          geom_point(size = 1.3, alpha = 0.8) +
          scale_color_manual(values = c("1st" = "darkblue", "Median" = "blue", "3rd" = "lightblue"), name = "Quantiles") +
          scale_fill_manual(values = c("1st" = "darkblue", "Median" = "blue", "3rd" = "lightblue"), name = "Quantiles") +
          #guides(fill = "none", color = "none") +
          labs(x = "Year", y = variables[i]) +
          scale_y_reverse() + 
          theme_classic(base_size = 11) +
          theme(plot.margin = unit(c(0.9, 0, 0, 0), "cm")))
}

pdf(file = "Figures/regions.pdf", width = 4, height = 6)
ggplot() +
  geom_sf(data = a) +
  geom_sf(data = northAmerica) +
  xlim(min(prediction_grid_roms$X)*1000-1000, max(prediction_grid_roms$X)*1000+1000) +
  ylim(min(prediction_grid_roms$Y)*1000-1000, max(prediction_grid_roms$Y)*1000+1000) +
  theme_linedraw()
dev.off()

# CT over the year
ggplot(centralTendency.grid.year) +
  geom_line(mapping = aes(x = year, y = centralTendency, group = gridid, color = gridid), show.legend = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Year", y = "Central tendency (month)") +
  scale_color_gradient("Central tendency\n(month)", low = "lightgoldenrodyellow", high =  "darkred") 

# southernCCE <- list(rbind(c(-170, 20),c(-110, 20), c( -110, 34.5),c( -170, 34.5),c(-170, 20)))
# centralCCE <- list(rbind(c(-170, 34.5),c(-110,34.5),c(-110,42),c(-170,42),c(-170,34.5)))
# or.wa <- list(rbind(c(-170,42),c(-110,42),c(-110,48.3),c(-170,48.3), c(-170,42)))
# britishColumbia <- list(rbind(c(-170,48.3),c(-110,48.3),c(-110,54.4),c(-170,54.4),c(-170,48.3)))
# alaska <- list(rbind(c(-170,54.4), c(-110,54.4),c(-110,60),c(-170,60), c(-170,54.4)))
# 
# regions <- st_multipolygon(list(southernCCE, centralCCE, or.wa, britishColumbia, alaska)) %>%
#   st_combine() %>% 
#   st_set_crs(4326) %>%
#   st_transform(5070) %>% 
#   as.data.frame()


