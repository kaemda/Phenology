# LOAD PACKAGES ------------
library(tidyr)
library(dplyr)
library(readxl)
library(sf)

# LOAD DATASETs ------------
all_tows_roms <- read.csv("Data/AllTows_200nm_ROMS.csv") %>%
  subset(., year >= 1995 & year <= 2019) 

positiveTows = read.csv("Data/PositiveCatches_200nm_allColumns.csv")

# GET SPECIES DATA -----------
getspeciesData <- function(species, speciesRangeSubset = "speciesRange", allgear = F) {
  
  # Extract positive tows first and center/scale larval data
  speciesData <- positiveTows %>%
    subset(scientific_name == species & year >= 1995 & year <= 2019) %>%
    mutate(., presence = 1) %>% 
    mutate(., larvae_10m2_logN1 = log(larvae_10m2+1), # Log(N+1) transform each data type
           larvae_m3_logN1 = log(larvae_m3+1)) %>% 
    mutate(., abundance_logN1_scaled = coalesce( # create coalesced scaled column of log transformed individual columns
      scale(larvae_10m2_logN1, center = F)[,1],
      scale(larvae_m3_logN1, center = F)[,1]
    )) %>% 
    mutate(., abundance = coalesce(  # create coalesced column of raw values
      larvae_10m2,
      larvae_m3
    ))
  
  # Get maximum extents of positive tows
  maxLat = max(speciesData$latitude)
  minLat = min(speciesData$latitude)
  maxLon = max(speciesData$longitude)
  minLon = min(speciesData$longitude)
  
  # Merge with tows dataset (includes zero tows and ROMS data)
  speciesData <- merge(all_tows_roms, speciesData[c("towID", "scientific_name", "gear",  "abundance","abundance_logN1_scaled")], all.x = T)
  
  # Replace all "NA" abundance values with zero (these are the true zeroes)
  speciesData <-  mutate(speciesData, 
                         abundance = replace_na(abundance, 0),
                         abundance_logN1_scaled = replace_na(abundance_logN1_scaled, 0),
                         scientific_name = species) %>% 
    subset(., !is.na(sst_roms) & !is.na(ssh_roms) & !is.na(salinity_roms)) %>% # Remove any tows without ROMS data
    subset(bottom_depth < 0) %>% 
    subset(., sst_roms > 1 & ssh_roms != 0 & salinity_roms != 0) %>%  # Remove any tows where ROMS has null values
    mutate(., sst_scaled = scale(sst_roms)[,1],
           ssh_scaled = scale(ssh_roms)[,1],
           salinity_scaled = scale(salinity_roms)[,1],
           spice_scaled = scale(spice_roms)[,1],
           bottom_depth_scaled = scale(bottom_depth)[,1])  # Center and scale enviro data

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
  
  # Add region
  speciesData <- speciesData %>%
    mutate(region = case_when(
      latitude < 34.5 ~ "Southern CCE",
      latitude < 42 ~ "Central CCE",
      latitude < 48.3 & longitude > -141 ~ "OR/WA",
      latitude < 54.4 ~ "British Columbia",
      TRUE ~ "Gulf of Alaska"
    ))
  
  return(speciesData)
}

