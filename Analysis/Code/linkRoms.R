# Data management
library(xlsx)
library(readxl) 
library(tidyr)
library(dplyr)
library(lubridate)

# Data access
library(ncdf4)

setwd("C://KDale/Projects/Phenology/")

# ROMS -------------------------------------------------------------------------
linkroms <- function(tows) {
  
  roms <- nc_open(filename = "Data/ROMS/nep_srf_1995-2019.nc")
  spiciness <- nc_open(filename = "Data/Environment/Spiciness/nep_avg_spice_dens_combined.nc")

  longitude=ncvar_get(roms, varid = "lon_rho")
  latitude=ncvar_get(roms, varid = "lat_rho")
  sst=ncvar_get(roms, varid = "temp")
  ssh=ncvar_get(roms, varid = "zeta")
  salinity = ncvar_get(roms, varid = "salt")
  #spice = ncvar_get(spiciness, varid = "spice")
  #spice = apply(spice, c(1,2,4), FUN=mean, na.rm=T) # average across depths
  #save(spice, file = "Data/Environment/Spiciness/spice.rdata")
  load(file = "Data/Environment/Spiciness/spice.rdata")
  
  # Get time information from the NC file
  # Create a year_month factor to link tows
  time=data.frame(date = as.Date(ncvar_get(roms, varid = "ocean_time")/86400, origin = '1900-01-01')) %>% 
    mutate(., year = year(date), month = month(date)) %>% 
    mutate(., year_month = paste(year,month, sep = "-"))
  
  time.spice = data.frame(date = as.Date(ncvar_get(spiciness, varid = "time")/86400, origin = '1900-01-01')) %>% 
    mutate(., year = year(date), month = month(date)) %>% 
    mutate(., year_month = paste(year,month, sep = "-"))
  
  nc_close(roms) # Close the netCDF file
  nc_close(spiciness)
  
  # Create a "year-month" column (ROMS output indexed by month)
  tows$year_month = paste(tows$year, tows$month, sep = "-") 
  
  # Set up progress bar
  pb <- txtProgressBar(min = 0, max = nrow(tows), char = "=", style = 3)
  tows$sst_roms = 0
  tows$ssh_roms = 0
  tows$salinity_roms = 0
  tows$spice = 0
  
  ### Match satellite data to sampling locations
  for (i in 1:nrow(tows)) {
    
    # Cancel current iteration if year is out of range
    if (tows$year[i] < min(time$year) | tows$year[i] > max(time$year)) 
      next
    
    targetLat = tows$latitude[i]
    targetLon = tows$longitude[i]
    
    # Get the best date match (ROMS output is monthly)
    dateIndex <-
      which(tows$year_month[i] == time$year_month)
    
    # Find all ROMS lat/lons within 0.5 deg of the sampling point
    nearbyLatitudeIndices <-
      which(abs(latitude -  targetLat) <= min(abs(latitude -  targetLat)+0.5))
    nearbyLongitudeIndices <-
      which(abs(longitude - targetLon) <= min(abs(longitude - targetLon)+0.5))
    
    # Find potential nearby points with valid latitude/longitudes
    matchIndices <-
      which(nearbyLatitudeIndices %in% nearbyLongitudeIndices == TRUE)
    
    # Calculate the best match (assumed to be the one with the closest latitude)
    differences = latitude[nearbyLatitudeIndices[matchIndices]] - targetLat
    bestMatch = matchIndices[which(abs(differences) == min(abs(differences)))]
    
    # Get the matching latitude
    matchLat = latitude[nearbyLatitudeIndices[bestMatch]]
    
    # Get row/column info (necessary for accessing the matrices)
    row = which(latitude == matchLat, arr.ind = TRUE)[1]
    column = which(latitude == matchLat, arr.ind = TRUE)[2]
    
    # Extract sst, ssh, salinity values
    tows$sst_roms[i] = sst[row, column, dateIndex]
    tows$ssh_roms[i] = ssh[row, column, dateIndex]
    tows$salinity_roms[i] = salinity[row, column, dateIndex]
    tows$spice[i] = spice[row, column, dateIndex]
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb) # Close progress bar
  
  # Return
  return(tows)
  
}
