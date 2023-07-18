# Load libraries
library(xlsx)
library(readxl) 
library(tidyr)
library(dplyr)
library(lubridate)

# Set working directory

# ---------------------------------------------------------------------
getIMECOCAL <- function() {

  sheet_names = as.character(excel_sheets(path = "OriginalDatasets/IMECOCAL_Nov 2022.xlsx"))
  nsheets = length(sheet_names)
  
  for(i in 1:nsheets) {
    
    sheet <- read_xlsx("OriginalDatasets/IMECOCAL_Nov 2022.xlsx", sheet = i, skip = 1)
    colnames(sheet) = tolower(colnames(sheet))
    
    # Go from wide to long
    sheet_long <- pivot_longer(sheet, cols = -colnames(sheet[1:9]), names_to = "scientific_name", values_to = "larvae_10m2") %>% 
      subset(., !is.na(larvae_10m2) & larvae_10m2 > 0) %>% # Remove any blank columns (artifact of Excel)
      mutate(., .before = "station", cruise = as.character(sheet_names[i])) # Add in cruise code, taken from spreadsheet name
    
    # Use character date/time formats -- I noticed these are different
    sheet_long$hour = as.character(sheet_long$hour)
    sheet_long$date <- as.Date(sheet_long$date, tryFormats = c("%d/%m/%Y", "%m/%d/%Y"))
    
    if (i == 1) {
      all_data = sheet_long #If this is the first time around the loop, make this all_data
    } else {
      all_data = bind_rows(all_data, sheet_long) # Otherwise, add the new sheet to all_data
    }
  }
  
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
  
  # Write file
  write.csv(x = all_data, file = "IMECOCAL.csv", row.names = FALSE)
  
  return(all_data)
}
#-------------------------------------------------------------------------------
imecocal <- getIMECOCAL()
