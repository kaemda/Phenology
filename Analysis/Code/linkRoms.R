# Data management
library(tidyr)
library(dplyr)
library(lubridate)
library(terra)

# Data access
library(ncdf4)

# ROMS -------------------------------------------------------------------------
linkroms <- function(tows) {
  
  roms <- nc_open(filename = "Data/ROMS/nep_srf_1995-2019.nc")
  spiciness <- nc_open(filename = "Data/Environment/Spiciness/nep_avg_spice_dens_combined.nc")
  
  longitude=ncvar_get(roms, varid = "lon_rho")
  latitude=ncvar_get(roms, varid = "lat_rho")
  sst=ncvar_get(roms, varid = "temp")
  ssh=ncvar_get(roms, varid = "zeta")
  salinity = ncvar_get(roms, varid = "salt")
  
  # Access spiciness dataset and average -- long runtime (access directly if already created)
  #spice = ncvar_get(spiciness, varid = "spice")
  #spice = apply(spice, c(1,2,4), FUN=mean, na.rm=T) # average across depths (the third dimension of the array)
  #save(spice, file = "Data/Environment/Spiciness/spice.rdata")
  load(file = "Data/Environment/Spiciness/spice.rdata")
  
  # Standardize time
  time=data.frame(date = as.Date(ncvar_get(roms, varid = "ocean_time")/86400, origin = '1900-01-01')) %>% 
    mutate(., year = year(date), month = month(date)) %>% 
    mutate(., year_month = paste(year,month, sep = "-"))
  
  # Close the netCDF files
  nc_close(roms) 
  nc_close(spiciness)
  
  # Compute year_month for tows
  tows$year_month <- paste(tows$year, sprintf("%02d", tows$month), sep = "-")
  time$year_month <- paste(time$year, sprintf("%02d", time$month), sep = "-")
  tows$roms_time_index <- match(tows$year_month, time$year_month)
  dateIndices <- unique(tows$roms_time_index)
  
  # Remove any tows that fall outside of the dataset
  tows = subset(tows, !is.na(tows$roms_time_index))
  
  covariates = c("sst_roms", "ssh_roms", "salinity_roms", "spice_roms")
  
  # Create storage dataset
  tows_linked = NA
  
  # Set up progress bar
  pb <- txtProgressBar(min = 0, max = length(dateIndices), char = "=", style = 3)
  
  # Bilinear interpolation
  for (i in 1:length(dateIndices)) {
    # for (i in 1:100) {
    
    # Get year/month index (corresponding to a "slice" of the ROMS matrix)
    dateIndex <- dateIndices[i]
    
    tows_subset = subset(tows, roms_time_index == dateIndex)
    
    # Create a raster from ROMS slice at that year/month combination
    sst_slice <- sst[,,dateIndex]
    ssh_slice <- ssh[,,dateIndex]
    salinity_slice <- salinity[,,dateIndex]
    spice_slice <- spice[,,dateIndex]
    
    # Create vectors
    lon_vec <- as.vector(t(longitude))
    lat_vec <- as.vector(t(latitude))
    sst_vec <- as.vector(t(sst_slice))
    ssh_vec <- as.vector(t(ssh_slice))
    salinity_vec <- as.vector(t(salinity_slice))
    spice_vec <- as.vector(t(spice_slice))
    
    # Build spatial data frame
    grid_df <- data.frame(lon = lon_vec, lat = lat_vec, sst_roms = sst_vec, ssh_roms = ssh_vec, salinity_roms = salinity_vec, spice_roms = spice_vec)
    
    # Convert to SpatVector, then rasterize
    grid_pts <- terra::vect(grid_df, geom = c("lon", "lat"))
    crs(grid_pts) = "EPSG:4326"
    
    # Tow points vector
    tow_pts <- vect(tows_subset[, c("longitude", "latitude")], 
                    geom = c("longitude", "latitude"), 
                    crs = "EPSG:4326")
    
    # Loop through covariates, calculating interpolated values at each tow point
    
    
    for (j in 1:length(covariates)) {
      res_x <- mean(diff(longitude[,1])) # get resolution in x and y directions
      res_y <- mean(diff(latitude[1,]))
      rast <- terra::rasterize(grid_pts, terra::rast(ext=ext(grid_pts), resolution=c(res_x, res_y)), field=covariates[j])
      crs(rast) = "EPSG:4326"
      
      # Extract interpolated values 
      interp_vals <- terra::extract(rast, tow_pts, method = "bilinear")
      
      # Add results back to original tows (2nd column is the interpolated value)
      if(nrow(tows_subset) == 1) {
        tows_subset[covariates[j]] <- interp_vals[1,2]
        
      } else {
        tows_subset[covariates[j]] <- interp_vals[, 2]
      }
    }
    
    # Bind into overall
    if(i == 1) {
      tows_linked = tows_subset
    } else {
      tows_linked = bind_rows(tows_linked, tows_subset)
    }
    
    setTxtProgressBar(pb, i)
    
  }
  close(pb) # Close progress bar
  
  return(tows_linked)
  
} 
