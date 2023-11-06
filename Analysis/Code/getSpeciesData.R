# GET SPECIES DATA -------------------------------------------------------------
getspeciesData <- function(species, speciesRangeSubset) {
  
  speciesData <- read.csv(file = "Data/AllCruises_Combined_200nm.csv") %>%
    subset(scientific_name == species & year >= 1995 & year <= 2019) %>%
    subset(gearGeneral != "MOCNESS" & gearGeneral != "Cobb MWT") %>% 
    mutate(., presence = 1) %>% 
    mutate(., abundance_scaled = coalesce( 
      scale(larvae_10m2, center = F)[,1],
      scale(larvae_m3, center = F)[,1],
      scale(larvae_1000m3, center = F)[,1],
      scale(larvae_count, center = F)[,1]
    )) %>% 
    mutate(., abundance = coalesce( 
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
  speciesData <- merge(all_tows_roms, speciesData[c("towID", "gear",  "scientific_name", "abundance", "abundance_scaled", "presence")], all.x = TRUE)
  
  # Replace all NA abundance values with zero (these are true zeroes)
  speciesData <-  mutate(speciesData, presence = replace_na(presence, 0),
                         # catch_anomaly_positive = replace_na(catch_anomaly_positive, 0),
                         abundance = replace_na(abundance, 0),
                         abundance_scaled = replace_na(abundance_scaled, 0),
                         # catch_anomaly = replace_na(catch_anomaly, 0),
                         scientific_name = species) %>% 
    subset(., !is.na(sst_roms) & !is.na(ssh_roms) & !is.na(salinity_roms)& !is.na(gear)) %>% # Remove any tows without ROMS data
    subset(., sst_roms != 0 & ssh_roms != 0 & salinity_roms != 0) %>%  # Remove any tows where ROMS has null values
    mutate(., sst_scaled = scale(sst_roms)[,1], ssh_scaled = scale(ssh_roms)[,1], salinity_scaled = scale(salinity_roms)[,1], spice_scaled = scale(spice)[,1]) %>%  # Center and scale enviro data
    mutate(., logN1 = log(abundance_scaled+1))
  

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
  
  speciesData$region <- factor(speciesData$region, levels = c("Southern CCE", "Central CCE", "OR/WA", "British Columbia", "Gulf of Alaska"))
  
  return(speciesData)
}

