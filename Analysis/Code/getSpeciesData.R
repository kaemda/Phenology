library(tidyr)
library(readxl)

setwd("C://KDale/Projects/Phenology/")
all_tows_roms <- read.csv("Data/AllTows_200nm.csv") %>%
  subset(., year >= 1995 & year <= 2019) 

# GET SPECIES DATA -------------------------------------------------------------
getspeciesData <- function(species, speciesRangeSubset, allgear = F) {
  
  speciesData <- read.csv(file = "Data/AllCruises_Combined_200nm.csv") %>%
    subset(scientific_name == species & year >= 1995 & year <= 2019) %>%
    mutate(., presence = 1) %>% 
    mutate(., larvae_10m2_logN1 = log(larvae_10m2+1), # log-transform each data type
           larvae_m3_logN1 = log(larvae_m3+1),
           larvae_1000m3_logN1 = log(larvae_1000m3+1),
           larvae_count_logN1 = log(larvae_count+1)) %>% 
    mutate(., abundance_logN1_scaled = coalesce( # create coalesced scaled column of log transformed individual columns
      scale(larvae_10m2_logN1, center = F)[,1],
      scale(larvae_m3_logN1, center = F)[,1],
      scale(larvae_1000m3_logN1, center = F)[,1],
      scale(larvae_count_logN1, center = F)[,1]
    )) %>% 
    mutate(., abundance_scaled = coalesce( # create coalesced scaled column of raw values
      scale(larvae_10m2, center = F)[,1],
      scale(larvae_m3, center = F)[,1],
      scale(larvae_1000m3, center = F)[,1],
      scale(larvae_count, center = F)[,1]
    )) %>% 
    mutate(., abundance = coalesce(  # create coalesced column of raw values
      larvae_10m2,
      larvae_m3,
      larvae_1000m3,
      larvae_count
    ))
  
  # Get maximum extents of positive tows
  maxLat = max(speciesData$latitude)
  minLat = min(speciesData$latitude)
  maxLon = max(speciesData$longitude)
  minLon = min(speciesData$longitude)
  
  # Merge with ROMS data
  speciesData <- merge(all_tows_roms, speciesData[c("towID", "gear",  "scientific_name", "abundance","abundance_logN1_scaled", "abundance_scaled", "presence")], all.x = TRUE)
  
  # Replace all NA abundance values with zero (these are true zeroes)
  speciesData <-  mutate(speciesData, presence = replace_na(presence, 0),
                         abundance = replace_na(abundance, 0),
                         abundance_logN1_scaled = replace_na(abundance_logN1_scaled, 0),
                         abundance_scaled = replace_na(abundance_scaled, 0),
                         scientific_name = species) %>% 
    subset(., !is.na(sst_roms) & !is.na(ssh_roms) & !is.na(salinity_roms)& !is.na(gear)) %>% # Remove any tows without ROMS data
    subset(., sst_roms != 0 & ssh_roms != 0 & salinity_roms != 0) %>%  # Remove any tows where ROMS has null values
    mutate(., sst_scaled = scale(sst_roms)[,1], ssh_scaled = scale(ssh_roms)[,1], salinity_scaled = scale(salinity_roms)[,1], spice_scaled = scale(spice)[,1]) %>%  # Center and scale enviro data
    mutate(., logN1 = log(abundance_scaled+1))
  
  if (allgear == F) {
    # subset to only programs that use rings, bongos, and mantas (i.e., larval sampling)
    speciesData = subset(speciesData, gearGeneral == "Bongo/Ring" | gearGeneral == "Manta")
  } 
  
  if(speciesRangeSubset == "speciesRange") { # Subset to within area bounded by positive tows
    speciesData <- subset(speciesData, latitude <= maxLat & latitude >= minLat & longitude >= minLon & longitude <= maxLon) 
  }
  
  # Add Conus Albers coordinates
  coords= cbind(speciesData$longitude, speciesData$latitude)
  scale = 1000
  albert_equal_area <- sf::sf_project(from = "EPSG:4326", to = 'EPSG:5070', pts = coords)/scale
  speciesData$X = albert_equal_area[,1]
  speciesData$Y = albert_equal_area[,2]
  
  # Add time block
  timeblocks <- read_xlsx(path = "Data/timeblocks.xlsx", sheet = 1)
  speciesData <- merge(speciesData, timeblocks, by = "year")
  speciesData$timeblock = factor(speciesData$timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019"))
  
  # Add latitudinal region
  speciesData$region = 0
  for (i in 1:nrow(speciesData)) {
    latitude = speciesData$latitude[i] 
    
    if(latitude < 34.5) {
      speciesData$region[i] = "Southern CCE"
    } else if (latitude < 42 ) {
      speciesData$region[i] = "Central CCE"
    } else if (latitude < 48.3) {
      speciesData$region[i] = "OR/WA"
    } else if (latitude < 54.4) {
      speciesData$region[i] = "British Columbia"
    } else {
      speciesData$region[i] = "Gulf of Alaska"
    } 
  }
  
  # Calculate top 10% and bottom 10% of range
  maxLat = max(speciesData$latitude) 
  minLat = min(speciesData$latitude)
  latRange10pct = (maxLat - minLat) * 0.1
  
  # Subset data into three sections to add range percentile
  upper10 = subset(speciesData, latitude >= (maxLat-latRange10pct)) %>% mutate(rangePercentile = "upper10")
  lower10 = subset(speciesData, latitude <= (minLat+latRange10pct)) %>% mutate(rangePercentile = "lower10")
  middle80 = subset(speciesData, latitude < (maxLat-latRange10pct) & latitude > (minLat+latRange10pct)) %>% mutate(rangePercentile = "middle80")
  
  # Rebind sections
  speciesData <- rbind(upper10, lower10, middle80)
  
  # Order sections
  speciesData$rangePercentile = factor(x = speciesData$rangePercentile, levels = c("upper10", "middle80", "lower10"))
  speciesData$region <- factor(speciesData$region, levels = c("Southern CCE", "Central CCE", "OR/WA", "British Columbia", "Gulf of Alaska"))
  
  return(speciesData)
}

