# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script to summarize individual-plant level phenotype data and generate
# population-means

# Load in data with presence/absence of HCN for every plant from every population
datAllPlants <- read_csv("data-clean/AllCities_AllPlants.csv",
                         col_types = "ccncddnnnd")

# Summarise HCN, Ac and Li frequency from each population from each city.
datPops <- datAllPlants %>%
  group_by(City, Transect) %>%
  
  # Min and Max distances
  mutate(min_dist = min(Distance),
         max_dist = max(Distance)) %>%
  ungroup() %>%
  
  # Calculate standardized distances
  group_by(City, Transect, Population) %>%
  mutate(std_distance = (Distance - min_dist)/ (max_dist - min_dist)) %>%
  ungroup() %>%
  
  # Calculate phenotype/allele frequencies
  group_by(City, Transect, Population, Distance, std_distance) %>%
  summarise(n_HCN = sum(!is.na(HCN_Result)), 
            sumHCN = sum(HCN_Result, na.rm = T), 
            freqHCN = (sumHCN/n_HCN),
            
            n_Ac = sum(!is.na(Locus.Ac)), 
            sum_Ac = sum(Locus.Ac, na.rm = T), 
            freqAc_marker = (sum_Ac/n_Ac), 
            sum_acac = n_Ac - sum_Ac,
            freq_acac = (sum_acac / n_Ac),
            acHWE = sqrt(freq_acac),
            AcHWE = 1 - acHWE,
            
            n_Li = sum(!is.na(Locus.Li)), 
            sum_Li = sum(Locus.Li, na.rm = T), 
            freqLi_marker = (sum_Li / n_Li),
            sum_lili = n_Li - sum_Li,
            freq_lili = (sum_lili / n_Li),
            liHWE = sqrt(freq_lili),
            LiHWE = 1 - liHWE) %>%
  
  # Add squared distance terms for quadratic regressions
  mutate(Distance_squared = Distance^2, 
         std_distance_squared = std_distance^2) %>%
  na_if("NaN") %>% 
  ungroup() %>% 
  mutate(Population = as.character(Population))

# Write Population dataset to disk
write.csv(datPops, "data-clean/AllCities_AllPopulations.csv", row.names = FALSE)

