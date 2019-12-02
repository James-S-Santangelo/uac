# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script to perform power analyses to examine the minimum number of plants per
# population, and number of populations, to sample to identify urban-rural
# clines in HCN without loss of power.

set.seed(42)

# Load in data for Ken's cities
# Coltypes required to convert population to character
# due to presence of letters for some populations after first 1000 rows.
KT_plants <- read_csv("data-clean/AllCities_AllPlants.csv",
                      col_types = "ccncddnnnd") %>% 
  filter(City %in% c("Toronto", "NewYork", "Boston", "Montreal"))

# Load in population level data
KT_pops <- read_csv("data-clean/AllCities_AllPopulations.csv") %>%
  filter(City %in% c("Toronto", "NewYork", "Boston", "Montreal"))

#### ANALYSIS OF BOSTON ####

# Model for Boston
KT_pops_Bos <- KT_pops %>% filter(City == "Boston")
ModBos <- lm(freqHCN ~ Distance, data = KT_pops_Bos)
SumBos <- summary(ModBos)
obs_BosSlope <- SumBos$coefficients[2]

# Get plant-level data for Boston, which shows the weakest cline overall
KT_plants_Bos <- KT_plants %>% filter(City == "Boston")

### SAMPLING PLANTS ###

nreps <- 10

# Slopes for Boston
S1B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 20))
S2B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 19))
S3B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 18))
S4B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 17))
S5B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 16))
S6B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 15))
S7B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 14))
S8B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 13))
S9B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 12))
S10B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 11))
S11B <- replicate(nreps, FunSlope(KT_plants_Bos, size = 10))

#Bind results into table#
SlopeBosPlants <- cbind(S1B, S2B, S3B, S4B, S5B, S6B, S7B, S8B, S9B, S10B, S11B)
SlopeBosPlants <- as.data.frame(SlopeBosPlants)


