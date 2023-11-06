# Load packages
# data management
library(xlsx)
library(readxl) 
library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)
library(RANN)

# data access
library(ncdf4)
library(rerddap)
library(rerddapXtracto)
library(marmap)

# analysis
library(qdap)

# mapping and plotting
library(sf)
library(ggplot2)
library(raster)

setwd("C://KDale/Projects/Phenology/")

speciesInfo <- read_xlsx(path = "Data/species_of_interest.xlsx", sheet = 1)
speciesNames <- as.vector(speciesInfo$scientific_name)

colors = c("IMECOCAL" = "firebrick4", "CalCOFI" = "coral2", "RREAS" = "darkgoldenrod1", "PRS_juveniles" = "darkseagreen2", "PRS_larvae" = "cornflowerblue", "NH-Line" = "lightblue1", "Canada" =  "mediumpurple3", "EcoFOCI" = "darkslategray")
regionColors = c("Southern CCE" = "firebrick3", "Central CCE" = "darkgoldenrod2", "OR/WA" ="darkseagreen3","British Columbia" ="cornflowerblue", "British Columbia" = "darkslateblue","Gulf of Alaska" ="darkslategray")

# Load North America shapefile
northAmerica <- st_read("C://KDale/GIS/North_South_America.shp") %>%
  st_union() %>% st_make_valid()
gulfOfCalifornia <- st_read("C://KDale/GIS/World_Seas_IHO_v3/World_Seas_IHO_v3.shp") %>% 
  subset(MRGID == 4314) %>% dplyr::select(., geometry)

# ERDDAP function ---------------------------------------------------------------
getERDDAP <- function(name) {
  
  data <- tabledap(x = name, url = "https://coastwatch.pfeg.noaa.gov/erddap") %>% 
    separate(., col = time, into = c("date", "time"), sep = c("T")) %>% 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date))
  
  data$cruise <- as.character(data$cruise)
  data$date <- as.Date(data$date)
  
  return(data)
}

# IMECOCAL ---------------------------------------------------------------------
getIMECOCAL <- function() {
  #species_codes <- read.csv("Species_codes.csv", as.is = TRUE)
  
  sheet_names = as.character(excel_sheets(path = "Data/OriginalDatasets/IMECOCAL_Nov 2022.xlsx"))
  nsheets = length(sheet_names)
  
  for(i in 1:nsheets) {
    
    sheet <- read_xlsx("Data/OriginalDatasets/IMECOCAL_Nov 2022.xlsx", sheet = i, skip = 1)
    
    # Go from wide to long
    sheet_long <- pivot_longer(sheet, cols = -colnames(sheet[1:9]), names_to = "scientific_name", values_to = "larvae_10m2") %>% 
      subset(., !is.na(larvae_10m2) & larvae_10m2 > 0) %>% # Remove any blank columns (artifact of Excel)
      mutate(., .before = "STATION", cruise = sheet_names[i]) # Add in cruise code
    
    # Make colummn names lowercase
    colnames(sheet_long) = tolower(colnames(sheet_long))
    
    # Use character date/time formats -- I noticed these are different
    sheet_long$hour = as.character(sheet_long$hour)
    sheet_long$date <- as.Date(sheet_long$date, tryFormats = c("%d/%m/%Y", "%m/%d/%Y"))
    
    if (i == 1) {
      all_data = sheet_long #If this is the first time around the loop, make this all_data
    } else {
      all_data = bind_rows(all_data, sheet_long) # Otherwise, add the new sheet to all_data
    }
  }
  
  # Cast the integer value in species names to a character (necessary to merge and keep non-numeric character labels)
  #species_codes$no = as.character(species_codes$no)
  
  # R creates new column names for duplicate columns (i.e., duplicate CalCOFI species codes) in the format "870...1".
  # So to return to the original CalCOFI number in order to link with species names, we need to remove the extra suffix
  #all_data$species_code = vapply(strsplit(all_data$species_code, "...", fixed = TRUE), "[", "", 1)
  
  # Link dataset to species names, convert to scientific names, add program, separate into line/station, standardize columns, extract month/day/year 
  all_data =  mutate(all_data, .before = 1, program = "IMECOCAL") %>% 
    #merge(., species_codes, by.x = "species_code", by.y = "no", all.x = TRUE) %>%
    rename(., 
           latitude = lat,
           longitude = long,
           day_night = "hour (0=night; 1= day)",
           time = hour,
           surface_temp_oC = 'sup temp',
           surface_sal_psu = 'sup sal',
           tow_depth_m = 'tow depth (m)') %>%
    separate(., col = station, into = c("line", "station"), sep = c("[.]")) %>%   # the brackets are necessary because "." is a special character 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date), across(c("line", "station"), as.character)) %>% 
    mutate(., day_night = replace(day_night, day_night == 1, "D"), gear = "CB", gearGeneral = "Bongo/Ring") %>% 
    mutate(., day_night = replace(day_night, day_night == 0, "N"))
  
  scientificNames <- vapply(strsplit(all_data$scientific_name, " \\("), "[", "", 1) 
  scientificNames <- vapply(strsplit(scientificNames, " Houttuyn"), "[", "", 1) 
  scientificNames <- vapply(strsplit(scientificNames, " Girard"), "[", "", 1) 
  scientificNames <- vapply(strsplit(scientificNames, " Jordan & Gilbert"), "[", "", 1) 
  all_data$scientific_name <- vapply(strsplit(scientificNames, ","), "[", "", 1)
  
  # Write file
  write.csv(x = all_data, file = "Data/IMECOCAL.csv", row.names = FALSE)
  
  return(all_data)
}

# CalCOFI ----------------------------------------------------------------------
getcalcofi <- function() {
  
  calcofi <- getERDDAP(name = "erdCalCOFIlrvcntpos") %>% 
    as.data.frame(.) %>%
    mutate(., .before = 1, program = "CalCOFI") %>% 
    rename(volume_sampled_m3 = volume_sampled, gear = net_type, species_code = calcofi_species_code) %>% 
    mutate(., across(c(volume_sampled_m3, latitude, longitude, larvae_10m2, larvae_1000m3), as.numeric) , across(c(cruise, line, station), as.character)) %>% 
    subset(., gear != "CV" & gear != "PV") %>% # Remove egg tow nets
    mutate(., gearGeneral = ifelse(gear == "MT", yes = "Manta", no = "Bongo/Ring"))
  
  hydros <- getERDDAP(name = "erdNOAAhydros") %>%
    subset(., standard_depth <= 210) %>%
    mutate(., across(c(cruise, line, station), as.character), across(c(standard_depth, temperature, salinity, density, oxygen, dynamic_height, percent_saturation, latitude, longitude), as.numeric)) %>%
    group_by_at(., c("cruise","line", "station", "date")) %>%
    summarise(., surface_sal_psu = mean(salinity),
              surface_temp_oC = mean(temperature),
              density_kg_m3 = mean(density),
              dissolved_oxygen_mL_L = mean(oxygen),
              dynamic_height_m = mean(dynamic_height))
  
  calcofi <- merge(
    x = calcofi,
    y = hydros,
    by.x = c("date", "line", "station", "cruise"),
    by.y = c("date", "line", "station", "cruise"),
    all.x = TRUE) 
  
  write.csv(x = calcofi, file = "Data/CalCOFI.csv", row.names = FALSE)
  
  return(calcofi)
}

# RREAS ------------------------------------------------------------------------
getrreas <- function() {
  
  # Catch data
  rreas <- getERDDAP(name = "FED_Rockfish_Catch") %>% 
    # subset(species_group == "Clupeoid" | species_group == "Cottid" | species_group == "Deep-Sea Smelt" |
    #          species_group == "Fish" | species_group == "Flatfish" | species_group == "Myctophid" | species_group == "Other Groundfish"|
    #          species_group == "Rockfish" | species_group == "Salmonid" | species_group == "Smelt") %>% 
    # mutate(., .before = 1, program = "RREAS") %>% 
    rename(., tow_number = haul_no, scientific_name = sci_name, larvae_count = catch) %>% 
    mutate(., across(c(latitude, longitude, bottom_depth), as.numeric) , across(c(station), as.character), day_night = "N") %>%   # all tows in RREAS occur at night
    mutate(., gear = "Cobb MWT", gearGeneral = "Cobb MWT") # all sampling is via a cobb midwater trawl
  
  # Environmental data
  rreas.env <- getERDDAP(name = "erdFedRockfishCtd") %>% 
    subset(.,depth <= 30) %>%
    mutate(., across(c(cruise, station), as.character),
           across(c(depth,
                    temperature,
                    salinity,
                    density,
                    transmissivity,
                    fluor_volt,
                    oxygen,
                    dyn_hgt,
                    chlorophyll,
                    latitude,
                    longitude),as.numeric)) %>%
    group_by_at(., c("cruise","station", "date")) %>%
    summarise(.,
              surface_sal_psu = mean(salinity),
              surface_temp_oC = mean(temperature),
              density_kg_m3 = mean(density),
              dissolved_oxygen_mL_L = mean(oxygen),
              dynamic_height_m = mean(dyn_hgt),
              irradiance_umol_m2_s = mean(irradiance),
              fluor_volt = mean(fluor_volt),
              oxygen_volt = mean(oxygen_volt),
              trans_percent = mean(transmissivity))
  
  rreas <- merge(rreas, rreas.env, all.x = TRUE)
  
  write.csv(x = rreas, file = "Data/RREAS.csv", row.names = FALSE)
  
  return (rreas)
}

# PRERECRUIT -------------------------------------------------------------------
getPrerecruit <- function() {
  
  # LARVAE --------------
  # bongo data
  # Pivot original dataset to long format -> reformat dates -> add program column
  prs.bongo <- read_xlsx("Data/OriginalDatasets/PreRecruit Bongo.xlsx", sheet = "Density Master") %>% 
    pivot_longer(data = ., cols = c(20:ncol(.)), values_to = "larvae_1000m3", names_to = "scientific_name") %>%
    mutate(., .before = 1, program = "PRS_larvae")
  
  # Standardize columns
  colnames(prs.bongo) = tolower(colnames(prs.bongo))
  prs.bongo <-
    rename(prs.bongo,
           station = 'station (new)',
           line = 'transect (new)',
           volume_sampled_m3 = 'volume filtered (m3)',
           tow_depth_m = 'haul depth (m)',
           bottom_depth = 'sta depth (m)',
           original_station = 'original station',
           original_transect = 'n/a transect') %>%
    mutate(., across(c("original_station", "line", "station"), as.character)) %>% 
    mutate(., .after = gear, gearGeneral = "Bongo/Ring") 

  # JUVENILES --------------
  # Midwater trawl data
  prs.mwt = read_xlsx(path = "Data/OriginalDatasets/PreRecruit MWT.xlsx", sheet = "Catch")
  prs.mwt.haul = read_xlsx(path = "Data/OriginalDatasets/PreRecruit MWT.xlsx", sheet = "Haul")
  
  colnames(prs.mwt) = tolower(colnames(prs.mwt))
  colnames(prs.mwt.haul) = tolower(colnames(prs.mwt.haul))
  
  prs.mwt <- merge(prs.mwt, prs.mwt.haul[c("distance from shore (km)", "start latitude", "start longitude", "start time", "year", "month", "day", "station (new)", "transect (new)", "total revs")]) %>% 
    rename( .,
            latitude = 'start latitude', longitude = 'start longitude',
            line = 'transect (new)',
            time = 'start time',
            scientific_name = taxon,
            larvae_count = number,
            tow_depth_m = 'start depth (m)',
            station_notes = comments,
            original_station = 'original station',
            station = 'station (new)',
            tow_number = 'new haul #',
            surface_temp_oC = 'surface temp (oc)') %>%
    mutate(., .before = 1, program = "PRS_juveniles", across(c("line", "station"), as.character), gear = "Cobb MWT", gearGeneral = "Cobb MWT")
  
  # Combine both datasets
  prerecruit <- bind_rows(prs.bongo, prs.mwt) %>% mutate(., day_night = "N")
  
  # Add date column
  prerecruit <- mutate(prerecruit, .before = year, date = as.Date(paste(prerecruit$year, prerecruit$month, prerecruit$day, sep = "-"), "%Y-%m-%d"))
  
  # Import environmental data
  prerecruit.env <- read_xlsx("Data/Environment/Prerecruit_CTD.xlsx", sheet =  "CTD Master")
  colnames(prerecruit.env) = tolower(colnames(prerecruit.env)) 
  
  # Add date column, rename columns
  prerecruit.env <- mutate(prerecruit.env, .before = year, date = as.Date(paste(prerecruit.env$year, prerecruit.env$month, prerecruit.env$day, sep = "-"), "%Y-%m-%d")) %>% 
    rename(
      line = 'transect (new)',
      station = 'station (new)',
      depth_m = 'depth (m)',
      surface_temp_oC = 'temperature (oc)',
      surface_sal_psu = 'salinity (psu)',
      density_kg_m3 = 'density (sigma-theta, kg/m3)',
      trans_percent = 'beam transmission (%)',
      fluor_volt = 'fluorescence (v)',
      fluor_mg_m3 = 'fluorescence (mg/m3)',
      dissolved_oxygen_mL_L = 'dissolved oxygen (ml/l)'
    ) %>% 
    subset(., depth_m < 100) %>% 
    group_by_at(., c("date", "line", "station")) %>% 
    summarize(surface_temp_oC = mean(surface_temp_oC))
  
  prerecruit <- merge(prerecruit, prerecruit.env, by = c("date", "line", "station"), all.x = TRUE)
  
  write.csv(prerecruit, file = "Data/Prerecruit_full.csv")
  
  return(prerecruit)
}

# NH LINE ----------------------------------------------------------------------
getNHLine <- function() {
  
  nhline <- read_xlsx("Data/OriginalDatasets/NHLine.xlsx", sheet = 1)
  
  colnames(nhline) = tolower(colnames(nhline))
  
  nhline <- rename(nhline, line = transect, volume_sampled_m3 = volume_m3_flowmeter, scientific_name = species, bottom_depth = 'station depth (m)', larvae_m3 = number_per_m3) %>% 
    mutate(., program = "NH-Line", .before = 1, larvae_count = larvae_m3 * volume_sampled_m3) %>% 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date)) %>% 
    mutate(., across(c("line"), as.character)) %>% 
    mutate(., .after = gear, gearGeneral = ifelse(gear == "MOCNESS (1.2 m2)", yes = "MOCNESS", no = ifelse(gear == "TT", yes = "Tucker trawl", no = "Bongo/Ring")))
  
  # Environmental data
  nhline.env <- read.csv(file = "Data/Environment/NHLine_CTD_lessthan100m.csv")
  colnames(nhline.env) <- tolower(colnames(nhline.env))
  nhline.env <- 
    dplyr::rename(nhline.env,
           'station code' = station.code,
           cruise = cruiseid,
           line = transect,
           date = sample.date,
           latitude = lat,
           longitude = long,
           time = time_local) %>%
    mutate(., date = as.POSIXct(date, tryFormats = "%m/%d/%y")) %>% 
    group_by_at(., c("station code")) %>% 
    summarize(., latitude = max(latitude),
              longitude = max(longitude),
              surface_sal_psu = mean(salinity),
              surface_temp_oC = mean(temperature),
              density_kg_m3 = mean(density),
              dissolved_oxygen_mL_L = mean(oxygen),
              fluor_volt = mean(fluor_volt),
              trans_percent = mean(tran_percent)) %>% 
    ungroup(.)
  
  nhline <- merge(nhline, nhline.env, by = c("station code"), all.x = TRUE)
  
  write.csv(nhline, file = "Data/NHLine.csv", row.names = FALSE)
  
  return(nhline)
  
}

# CANADA -----------------------------------------------------------------------
getCanada <- function() {
  
  canada <- read_xlsx("Data/OriginalDatasets/Canada.xlsx", sheet = "Fish")
  colnames(canada) = tolower(colnames(canada))
  canada <- rename(canada, time = stn_time, volume_sampled_m3 = 'volume filtered(m3)', latitude = lat, longitude = lon, scientific_name = name, bottom_depth = 'bottom depth(m)', day_night = twilight, gear = net_type, depth_strata = depth_strt1, larvae_m3 = 'abundance(#/m3)', phylum = 'phylum:', class = 'class:', order = 'order:', family = 'family:') %>% 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date)) %>% 
    mutate(., day_night = replace(day_night, day_night == "Daylight", "D") , day_night = replace(day_night, day_night == "Night", "N")) %>% 
    separate_wider_delim(., cols = scientific_name, names = c("genus", "species", "maturity"), delim = " ", too_many = "merge", cols_remove = F) %>% 
    mutate(., species = replace(species, species == "*sp.", NA)) %>% 
    unite(., col = scientific_name, genus, species, na.rm = TRUE, sep  = " ") %>% 
    mutate(., across(c(latitude, longitude, volume_sampled_m3, depth_strata, bottom_depth, depth_end1), as.numeric) , across(c(time), as.character)) %>% 
    mutate(., .before = larvae_m3, larvae_count = larvae_m3 * volume_sampled_m3) %>%
    mutate(., .before = 1, program = "Canada") %>% 
    mutate(., .after = gear, gearGeneral = "Bongo/Ring")
  

  write.csv(x = canada, file = "Data/Canada.csv", row.names = FALSE)
  
  return(canada)
}

# EcoFCOI ----------------------------------------------------------------------
getEcoFOCI <- function() {
  
  ecofoci <- read_xlsx("Data/OriginalDatasets/EcoFOCI.xlsx", sheet = "DuplicatesRemoved")
  
  colnames(ecofoci) <- tolower(colnames(ecofoci))
  
  ecofoci <- separate(ecofoci, col = haul_id, sep = c(" "), into = c(NA, "station", NA, "gear", "net_num"), remove = TRUE) %>% 
    rename(., tow_number = haul_name, latitude = lat, longitude = lon, date = gmt_date_time, scientific_name = species_name, region_name = geographic_area, larvae_10m2 = larvalcatchper10m2, larvae_1000m3 = larvalcatchper1000m3, larvae_count = number_measured_counted) %>% 
    mutate(., date = as.Date(date), .before = 1, program = "EcoFOCI") %>%
    mutate(., .after = date, month = month(date), day = day(date)) %>% 
    mutate(., .after = gear, gearGeneral = "Bongo/Ring")
  
  write.csv(ecofoci, "Data/EcoFOCI.csv", row.names = FALSE)
  
  return(ecofoci)
}
# GET TOWS ---------------------------------------------------------------------
getTows <- function(samples) {
  
  calcofi_tows <- getERDDAP(name = "erdCalCOFItows") %>%
    mutate(across(c(latitude, longitude), as.numeric)) %>%
    mutate(., .before = 1, program = "CalCOFI", across(c(percent_sorted, latitude, longitude), as.numeric)) %>%
    rename(gear = net_type) %>% 
    subset(., gear != "CV" & gear != "PV") %>% # Remove egg tow nets
    mutate(.,  across(c("line", "station"), as.character), program = "CalCOFI", gearGeneral = ifelse(gear == "MT", yes = "Manta", no = "Bongo/Ring"))
  
  # rreas_tows <- group_by_at(rreas, c("program", "date", "station", "gear")) %>%
  #    summarize(., latitude = max(latitude), longitude = max(longitude), total_larvae = sum(larvae_count))

  canada_tows <- read_xlsx("Data/OriginalDatasets/Canada.xlsx", sheet = "Station Headers")
  colnames(canada_tows) = tolower(colnames(canada_tows))
  canada_tows <- rename(canada_tows, latitude = lat, longitude = long) %>% 
    rename(., bottom_depth = 'max_depth1', day_night = twilight, gear = net_type, depth_strata = depth_strt1) %>% 
    separate_wider_position(., cols = key, widths = c(10, "tow_key" = 10), too_few =  "align_start" ) %>% 
    mutate(., across(c(latitude, longitude, depth_strata, bottom_depth), as.numeric)) %>% 
    mutate(., .before = 1, program = "Canada") %>% 
    mutate(., .after = gear, gearGeneral = "Bongo/Ring")
  
  canada_combined <- merge(canada_tows, canada[c("cruise", "tow_key", "latitude", "longitude", "day", "month", "year")], all = T) %>% 
    mutate(., towID = paste(
      program,
      cruise,
      date,
      station,
      gearGeneral,
      tow_key,
      sep = "_")) %>% 
    distinct(., towID, .keep_all = T)
  
  # For cruises without "line" or "tow_number" or "tow_key" columns, remove NAs in the towID
  canada_combined$towID = gsub(pattern = "_NA_", replacement = "_", x = canada_combined$towID)
  canada_combined$towID = gsub(pattern = "_NA", replacement = "", x = canada_combined$towID)
  
  all_tows <- subset(samples, program != "CalCOFI" & program != "Canada") %>% 
    bind_rows(., canada_combined, calcofi_tows) %>% 
    group_by_at(., c("towID", "program", "date", "latitude", "longitude", "year", "month", "day", "gear", "gearGeneral")) %>% 
    summarize(
      surface_temp_oC = mean(surface_temp_oC),
      surface_sal_psu = mean(surface_sal_psu),
      dissolved_oxygen_mL_L = mean(dissolved_oxygen_mL_L),
      dynamic_height_m = mean(dynamic_height_m),
      fluor_volt = mean(fluor_volt),
      density_kg_m3 = mean(density_kg_m3),
      trans_percent = mean(trans_percent)) %>% 
    subset(., !is.na(latitude)) 
  
  all_tows$date = as.POSIXct(all_tows$date, format = "%m/%d/%Y")
  all_tows$day_of_year = yday(all_tows$date)
  
  write.csv(x = all_tows, file = "Data/AllTows.csv", row.names = FALSE)
  
  return(all_tows)
}

# SELECT WITHIN EEZ --------------------------------------------------------------
selectWithinEEZ <- function(data) {
  
  # Load EEZ shapefile
  eez <- read_sf("C://KDale/GIS/NorthAmerica_EEZ.shp") %>% st_buffer(., dist = 100000)
  
  # Subset tows to within 200 km of land (EEZ)
  data_sf <- subset(data, !is.na(latitude)) %>% st_as_sf(., coords = c("longitude", "latitude"), remove = FALSE)
  data_sf$program = factor(data_sf$program, levels = c("CalCOFI","IMECOCAL","RREAS","PRS_juveniles","PRS_larvae", "NH-Line","Canada", "EcoFOCI"))
  st_crs(data_sf) <- st_crs(eez)
  eez_data <- st_intersection(data_sf, eez)
  eez_data <- subset(eez_data, latitude > 23)
  
  return(eez_data)
  
}

# MAKE MAP ----------------------------------------------------------------------
createMap <- function(data) { 
  
  data = data %>% subset(., year >= 1995 & year <= 2019) %>% 
    subset(., !is.na(latitude) & !is.na(longitude)) %>% 
    subset(., program != "RREAS" & program != "PRS_juveniles") %>% 
    st_as_sf(., coords = c("longitude", "latitude"))
  
  st_crs(data) = st_crs(northAmerica)
  
  data$program = factor(data$program, levels = c("CalCOFI",  "Canada","PRS_larvae", "PRS_juveniles", "NH-Line", "RREAS",  "IMECOCAL", "EcoFOCI"))
  nhline <- subset(data, program == "NH-Line")
  canada <- subset(data, program == "Canada")
  calcofi <- subset(data, program == "CalCOFI")
  
  jpeg("Figures/All_tows.jpg",units = "in", width = 8, height = 8, res = 500)
  print(ggplot() +
          geom_sf(data = calcofi, mapping = aes(color = program), alpha = 0.3) +
          geom_sf(data = subset(data, program != "CalCOFI" & program != "Canada" & program != "NH-Line"), mapping = aes(color = program), alpha = 0.3) +
          geom_sf(data = nhline, mapping = aes(color = program), alpha = 0.3) +
          geom_sf(data = canada, mapping = aes(color = program), alpha = 0.5) +
          geom_sf(data = northAmerica) +
          coord_sf(expand = FALSE) +
          scale_color_manual("Program", values = colors) +
          coord_sf(xlim = c(-165, -106), ylim = c(19,61)) +
          theme_bw(base_size = 14))
  dev.off()
  
  jpeg("Figures/Stations_All_tows_Legend.jpg", width = 8, height = 8, res = 400, units = "in")
  print(ggplot() +
          geom_sf(data = data, mapping = aes(color = program), alpha = 1) +
          geom_sf(data = nhline, mapping = aes(color = program), alpha = 1) +
          geom_sf(data = northAmerica) +
          coord_sf(expand = FALSE) +
          scale_color_manual("Program", values = colors) +
          coord_sf(xlim = c(-170, -105), ylim = c(18,61)) +
          theme_bw(base_size = 14))
  dev.off()
  
  jpeg("Figures/Stations_Canada.jpg", width = 8, height = 8, res = 400, units = "in")
  print(ggplot() +
          geom_sf(data = subset(data, program == "Canada"), mapping = aes(color = program), alpha = 1) +
          geom_sf(data = northAmerica) +
          coord_sf(expand = FALSE) +
          scale_color_manual("Program", values = colors) +
          coord_sf(xlim = c(-170, -105), ylim = c(18,61)) +
          theme_bw(base_size = 14))
  dev.off()
  
}

# CREATE SPECIES TABLE ---------------------------------------------------------
createSpeciesTable <- function(data) {
  
  data = data %>% as.data.frame()
  
  speciesOfInterest <- filter(data, scientific_name %in% speciesNames) %>% subset(., larvae_count > 0)
  
  speciesAbundances <-
    group_by_at(speciesOfInterest, c("program","scientific_name")) %>%
    summarize(., n = sum(larvae_count)) %>% 
    pivot_wider(., id_cols = "scientific_name", names_from = "program", values_from = "n", names_prefix = "n_")
  speciesAbundances$nMin = min(speciesAbundances[,which(colnames(speciesAbundances)=="n_IMECOCAL"):which(colnames(speciesAbundances)=="n_EcoFOCI")])
  speciesAbundances$nMax = max(speciesAbundances[,which(colnames(speciesAbundances)=="n_IMECOCAL"):which(colnames(speciesAbundances)=="n_EcoFOCI")])
  
  years <- group_by_at(data, c("program", "year")) %>%
    summarise(n = n()) %>%
    group_by(., program) %>%
    summarise(program_years = n())
  
  speciesFrequency <-
    group_by_at(speciesOfInterest, c("program", "year", "scientific_name")) %>%
    summarise(n = n()) %>%
    group_by_at(., c("program","scientific_name")) %>%
    summarise(., seen_years = n()) %>% 
    merge(., years, by = "program") %>% 
    mutate(., freq_across_years = seen_years / program_years) %>% 
    pivot_wider(., id_cols = "scientific_name", names_from = "program", values_from = "freq_across_years", names_prefix = "freq_")
  speciesFrequency$freqMin = min(speciesFrequency[,which(colnames(speciesFrequency)=="freq_IMECOCAL"):which(colnames(speciesFrequency)=="freq_EcoFOCI")])
  speciesFrequency$freqMax = max(speciesFrequency[,which(colnames(speciesFrequency)=="freq_IMECOCAL"):which(colnames(speciesFrequency)=="freq_EcoFOCI")])
  
  speciesTable <- merge(speciesInfo,speciesAbundances,  by = "scientific_name") %>% 
    merge(., speciesFrequency, by = "scientific_name") %>%
    replace(is.na(.), 0)
  
  write.xlsx(speciesTable, file = "Data/Species_table.xlsx")
  
  return(speciesTable)
}

# RUN FUNCTIONS ----------------------------------------------------------------
## Standardize datasets  -----
imecocal <- getIMECOCAL()
calcofi <- getcalcofi()
rreas <- getrreas()
prerecruit <- getPrerecruit()
nhline <- getNHLine()
canada <- getCanada()
ecofoci <- getEcoFOCI()

## Load datasets --------
# As an alternative to the above, import already-cleaned CSV versions (use mutate/across to assign correct data types)
imecocal <- read.csv("Data/IMECOCAL.csv") %>% mutate(., across(date, as.Date)) %>% mutate(., across(c(cruise, line, station), as.character))
calcofi <- read.csv("Data/CalCOFI.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across(c("cruise", "line", "station"), as.character))
rreas <- read.csv("Data/RREAS.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across(c("cruise", "station"), as.character))
prerecruit <- read.csv("Data/Prerecruit_full.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across(c("line", "station"), as.character))
nhline <- read.csv("Data/NHLine.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across("station", as.character))
canada <- read.csv("Data/canada.csv") %>% mutate(., across("date", as.Date))
ecofoci <- read.csv("Data/ecofoci.csv") %>%  mutate(., across("gear", as.character))
ecofoci$date = as.Date(ecofoci$date, format = c("%m/%d/%Y"))

## Combine datasets, select relevant columns, add unique tow ID -----
all_data <- bind_rows(list(imecocal, calcofi, rreas, prerecruit, nhline, canada, ecofoci)) %>%
  dplyr::select(., program, cruise, line, station, date, year, month, day, time,
         day_night, gear, gearGeneral, tow_number, tow_key, latitude, longitude, scientific_name, 
         volume_sampled_m3, larvae_count, larvae_10m2, larvae_m3, larvae_1000m3, maturity,
         tow_depth_m, bottom_depth,
         surface_temp_oC, surface_sal_psu, dissolved_oxygen_mL_L, dynamic_height_m, fluor_volt, density_kg_m3, trans_percent) %>% 
  mutate(., .after = program, towID = paste(
    program,
    cruise,
    date,
    line,
    station,
    gearGeneral,
    tow_number,
    tow_key,
    sep = "_"))

# For cruises without "line" or "tow_number" or "tow_key" columns, remove NAs in the towID
all_data$towID = gsub(pattern = "_NA_", replacement = "_", x = all_data$towID)
all_data$towID = gsub(pattern = "_NA", replacement = "", x = all_data$towID)

## Select positive samples within 300nm, remove Gulf of CA -----
all_data_eez <- selectWithinEEZ(all_data) %>%
  st_difference(., gulfOfCalifornia) %>% 
  st_drop_geometry(.) %>% 
  subset(., larvae_count > 0 | larvae_10m2 > 0 | larvae_m3 > 0 | larvae_1000m3 > 0) %>% 
  mutate(., .after = day, day_of_year = yday(date)) # Add a day of the year column

## Get tows -----
all_tows <- getTows(all_data)

## Link ROMS to tows, subset data ------
source("Analysis/Code/linkRoms.R")
all_tows_roms <- selectWithinEEZ(all_tows) %>%
  st_difference(., gulfOfCalifornia) %>% 
  mutate(., distance_from_shore_m = as.vector(st_distance(., northAmerica)[,1])) %>% 
  mutate(., distance_from_shore_scaled = scale(distance_from_shore_m)[,1]) %>% 
  st_drop_geometry(.) %>% 
  linkroms(.) %>% 
  mutate(daylength = geosphere::daylength(lat = latitude, doy = day_of_year)) # daylength

# Get bathymetry data
bathy <- getNOAA.bathy(lon1 = -180, lon2 = -100, lat1 = 15, lat2 = 70, resolution = 4)
all_tows_roms$bottom_depth <- get.depth(bathy, x = all_tows_roms$longitude, y = all_tows_roms$latitude, locator = FALSE)$depth
all_tows_roms$bottom_depth_scaled <- scale(all_tows_roms$bottom_depth)[,1]

## Write data -----
write.csv(all_data, file = "Data/AllCruises_Combined.csv")
write.csv(all_data_eez, file = "Data/AllCruises_Combined_200nm.csv", row.names = F)
write.csv(all_tows, file = "Data/AllTows.csv")
write.csv(all_tows_roms, file = "Data/AllTows_200nm.csv", row.names=F) 

## Read data -----
all_data <- read.csv(file = "Data/AllCruises_Combined.csv")
all_data_eez <- read.csv(file = "Data/AllCruises_Combined_200nm.csv")
all_tows <- read.csv(file = "Data/AllTows.csv")
all_tows_roms <- read.csv(file = "Data/AllTows_200nm.csv")

## Species table ----
speciesTable = createSpeciesTable(all_data_eez)

## Map of all tows ----
createMap(data = all_tows_roms)

# SPATIAL OPERATIONS ----------------------------------------------------------

# Create EEZ buffer for North America
# eez <- read_sf("C://KDale/GIS/World_EEZ_v11_20191118/World_EEZ_v11_20191118/eez_v11.shp") %>% 
#   st_make_valid(.) %>%  
#   st_crop(., xmin = -170, xmax = -109,
#           ymin = 20, ymax = 62,
#           expand=FALSE) %>% 
#   subset(., TERRITORY1 != "Hawaii" & TERRITORY1 != "Johnston Atoll") %>% 
#   st_union(.) %>% 
#   write_sf(., "C://KDale/GIS/NorthAmerica_EEZ.shp")

# northAmerica <- read_sf("C://KDale/GIS/North_South_America.shp") %>% 
#   st_crop(., xmin = -100, xmax = -180,
#           ymin = 15, ymax = 65,
#           expand=TRUE) %>% 
#   st_union(northAmerica)
# write_sf(northAmerica, "C://KDale/GIS/NorthAmerica.shp")

# SANDBOX --------------------------

