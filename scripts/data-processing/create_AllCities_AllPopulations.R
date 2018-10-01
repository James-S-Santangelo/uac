
# Load required packages
library(dplyr)

# Load in data with presence/absence of HCN for every plant from every population
datAllPlants <- read.csv("data-clean/AllCities_AllPlants.csv")

# Summarise HCN, Ac and Li frequency from each population from each city.
datPops <- datAllPlants %>%
  group_by(City, Transect) %>%
  mutate(min_dist = min(Distance),
         max_dist = max(Distance)) %>%
  ungroup() %>%
  group_by(City, Transect, Population) %>%
  mutate(std_distance = Distance / max_dist) %>%
  ungroup() %>%
  group_by(City, Transect, Population, Distance, std_distance) %>%
  summarize(n_HCN = sum(!is.na(HCN_Result)), sumC = sum(HCN_Result, na.rm = T), FreqC = (sumC/n_HCN),
            n_Ac = sum(!is.na(Locus.Ac)), sumAc = sum(Locus.Ac, na.rm = T), FreqAc = (sumAc/n_Ac), 
            n_Li = sum(!is.na(Locus.Li)), sumLi = sum(Locus.Li, na.rm = T), FreqLi = (sumLi/n_Li)) %>%
  mutate(Distance_squared = Distance^2, 
         std_distance_squared = std_distance^2) %>%
  na_if("NaN")

# Write Population dataset to disk
write.csv(datPops, "data-clean/AllCities_AllPopulations.csv", row.names = FALSE)
