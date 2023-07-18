# Data management
library(xlsx)
library(readxl)
library(tidyverse)
library(dplyr)

setwd("C://KDale/Projects/Phenology/")
programColors = c("IMECOCAL" = "firebrick3", "CalCOFI" = "coral3","RREAS" = "darkgoldenrod2", "PRS_juveniles" ="darkseagreen3","PRS_larvae" ="cornflowerblue", "NH-Line" ="deepskyblue3", "Canada" = "darkslateblue","EcoFOCI" ="darkslategray")

# Load model and data -------
load(file = paste0("Results/Sardinops sagax/Models/presence_speciesRange_allPrograms_ar1_s(sst_scaled) + s(ssh_scaled) +s(salinity_scaled) +s(distance_from_shore_scaled)__as.factor(gearGeneral).rdata"))

# Bin temperatures
tempBins <- seq(floor(min(data$sst_roms)), ceiling(max(data$sst_roms)), 2)

data$tempBins <- round(data$sst_roms)

# Summarize number of fish in each bin
fish.per.temp.bin <- group_by(data, gearGeneral, tempBins) %>%
  summarize(., presence_sum = sum(presence), numTows = length(unique(towID))) %>%
  subset(., presence_sum > 0) %>%
  mutate(., percFish = presence_sum/sum(presence_sum), percTows = numTows/sum(numTows)) %>% 
  mutate(quotient = percFish/percTows)

fish.per.temp.bin <- group_by(data, gearGeneral, tempBins) %>%
  summarize(., catch_anomaly_positive_sum = sum(catch_anomaly_positive), numTows = length(unique(towID))) %>%
  subset(., catch_anomaly_positive_sum > 0) %>%
  mutate(., percFish = catch_anomaly_positive_sum/sum(catch_anomaly_positive_sum), percTows = numTows/sum(numTows)) %>% 
  mutate(quotient = percFish/percTows)

# Plot from empirical catch data
gearColors = c("IMECOCAL" = "firebrick3", "CalCOFI" = "coral3","RREAS" = "darkgoldenrod2", "PRS_juveniles" ="darkseagreen3","PRS_larvae" ="cornflowerblue", "NH-Line" ="deepskyblue3", "Canada" = "darkslateblue","EcoFOCI" ="darkslategray")

ggplot() + 
  geom_line(fish.per.temp.bin, mapping = aes(x = tempBins, y = quotient, color = gearGeneral), linewidth = 1.5) +
  scale_color_manual("Gear", values = c("goldenrod", "dodgerblue", "slateblue4")) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic(base_size = 12) +
  labs(x = "Temperature (degC)", y = "Quotient ")
