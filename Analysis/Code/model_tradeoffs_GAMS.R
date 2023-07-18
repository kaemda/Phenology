library(mgcv)
library(xlsx)
library(ggplot2)
library(visreg)
library(tidyr)
library(sf)

setwd("C://KDale/Projects/Phenology/")

# Load data
all_tows_roms <- read.csv("Data/AllTows_200nm.csv") %>%
  subset(., year >= 1995 & year <= 2019) %>% 
  subset(gearGeneral == "Bongo/Ring" | gearGeneral == "Cobb MWT" | gearGeneral == "Manta")

northAmerica <- read_sf("C://KDale/GIS/North_South_America.shp")
northAmerica <- sf::st_transform(northAmerica, crs = 5070)
# Select species
species = "Sardinops sagax"
source("Analysis/Code/getSpeciesData.R")
speciesRange = "speciesRange"
data <- getspeciesData(species = species, speciesRangeSubset = speciesRange)
data$year_scaled = scale(data$year)

dependentVar = "logN1"

# Run base model
base.gam <- gam(formula = logN1 ~ s(sst_scaled) + s(ssh_scaled) + s(salinity_scaled) + s(distance_from_shore_scaled) + s(latitude,longitude) + as.factor(gearGeneral) + s(month, bs = "cc", k = 12),
            family = tw(link = "log"),
            data = data)

# Run geography model
geo.gam <- gam(formula = logN1 ~ s(sst_scaled) + s(ssh_scaled) + s(salinity_scaled) + s(distance_from_shore_scaled) + s(latitude,longitude, by=as.factor(year)) + s(month, bs = "cc", k = 12) + as.factor(gearGeneral),
            family = tw(link = "log"),
            data = data)

# Run seasonality model
pheno.gam <- gam(formula = logN1 ~ s(sst_scaled) + s(ssh_scaled) + s(salinity_scaled) + s(distance_from_shore_scaled) + s(latitude,longitude) + as.factor(gearGeneral) + s(month, by=factor(year), bs = "cc", k = 12),
            family = tw(link = "log"),
            data = data)

# Compare
base.gam %>% summary()
geo.gam %>% summary()
pheno.gam %>% summary()

base.gam$aic
geo.gam$aic
pheno.gam$aic

save(data, base.gam, geo.gam, pheno.gam, file = paste0("Results/", species, "/Models/", dependentVar,"_", speciesRange, "GAM_tradeoff_models.rdata"))

