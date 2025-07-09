# Load packages
library(tidyr)
library(sf)
library(marmap)
library(ggplot2)

# Scale variables based on empirical data (so that scaled values are equivalent)
scaleCovariates <- function(x, name) {
  center = attributes(scale(data[name]))$`scaled:center`
  scale = attributes(scale(data[name]))$`scaled:scale`
  (x - center) / scale # autoreturns
}

# Create prediction grid ---------------------------------------------------
createPredictionGrid <- function(species, path) {

  # load species data
  source("Analysis/Code/getSpeciesData.R")
  data =  getspeciesData(as.character(species), speciesRangeSubset = "speciesRange", allgear = F) %>% 
    subset()

  # load shapefile
  northAmerica <- read_sf("Data/Shapefiles/North_South_America/North_South_America.shp") %>% st_union() %>%
    st_transform(., crs = "EPSG:5070")

  # Get species dataset, cast to a spatial object
  data.sf <- st_as_sf(data, coords = c("longitude", "latitude")) %>%
    st_cast(., to = "MULTIPOINT") %>%
    st_combine(.) %>%
    st_set_crs(4326) %>% # WGS 84 - geographic reference system
    st_transform(5070) # Albers equal area conic (aka conic albers) - projected

  # Create study area, with an offshore buffer of 50 km.
  studyArea <-
    st_convex_hull(data.sf) %>%
    st_make_valid() %>%
    st_difference(., northAmerica) %>%
    st_buffer(., dist = -10000) # buffer distance of 10 km offshore

  # Create grid that covers just the study area. Cell size in EPSG 5070 is in meters. 50000 cell size = 50 km
  grid <- st_make_grid(studyArea, crs = st_crs(northAmerica), cellsize = 50000, square = FALSE, what = "centers") %>%
    st_make_valid() %>%
    st_intersection(., studyArea)

  # Create grid dataframe and add distance to shore
  grid.df <- fortify(grid) %>% 
    mutate(., distance_from_shore_m = as.vector(st_distance(., northAmerica)[,1])) %>%
    mutate(., distance_from_shore_scaled = scale(distance_from_shore_m)[,1])

  # Convert grid to projection that uses decimal degrees and add dd lat/long to dataframe
  grid <- st_transform(grid, 4326) %>% st_cast(., to = "MULTIPOINT")
  grid.df <- mutate(grid.df, latitude = st_coordinates(grid)[,2],
                    longitude = st_coordinates(grid)[,1],
                    gridid = seq(1:nrow(grid.df)))

  # Add Albers equal area
  coords= cbind(grid.df$longitude, grid.df$latitude)
  scale = 1000
  albert_equal_area <- sf::sf_project(from = "EPSG:4326", to = 'EPSG:5070', pts = coords)/scale
  grid.df$X = albert_equal_area[,1]
  grid.df$Y = albert_equal_area[,2]

  # Create prediction grid of all combinations of points, years, and months
  year <- seq(from = min(data$year), to = max(data$year), by = 1)
  month <- seq(from = 1, to = 12, by = 1)
  day <- 15
  prediction_grid <- expand_grid(grid.df, year, month, day) %>%
    mutate(date = strptime(paste(year,month,day, sep = "-"), format = "%Y-%m-%d")) %>% # construct date column %>% 
    mutate(day_of_year = lubridate::yday(date)) %>% # add day of year
    mutate(daylength = geosphere::daylength(lat = latitude, doy = day_of_year)) # Add daylength
    
  # Link ROMS data (takes a long time!)
  source("Analysis/Code/linkRoms.R")
  prediction_grid_roms <- linkroms(tows = prediction_grid)

  # Scale/center ROMS data based on *empirical* data values
  prediction_grid_roms <- prediction_grid_roms %>% 
    mutate(., sst_scaled = scaleCovariates(x = sst_roms, "sst_roms"),
           ssh_scaled = scaleCovariates(x = ssh_roms, "ssh_roms"),
           salinity_scaled = scaleCovariates(x = salinity_roms, "salinity_roms"),
           spice_scaled = scaleCovariates(x = spice_roms, "spice_roms")) %>% # Center and scale enviro data
    subset(., !is.na(sst_scaled) & !is.na(salinity_scaled) & !is.na(ssh_scaled)) # Subset to only locations where we have data
  
  # Get bathymetry data
  # Scale based on empirical dataset, like other covariates
  bathy <- getNOAA.bathy(lon1 = min(data$longitude)-1, lon2 = max(data$longitude)+1, lat1 = min(data$latitude)-1, lat2 = max(data$latitude)+1, resolution = 2)
  prediction_grid_roms$bottom_depth <- get.depth(bathy, x = prediction_grid_roms$longitude, y = prediction_grid_roms$latitude, locator = FALSE)$depth
  prediction_grid_roms$bottom_depth_scaled <- scaleCovariates(prediction_grid_roms$bottom_depth, "bottom_depth")
  
  save(prediction_grid_roms, grid.df, file = path)
}
