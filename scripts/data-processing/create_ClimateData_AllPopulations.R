# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script to generate combined dataset with all climatic variables

#### LOAD CLIMATE DATASETS ####

# Aridity index
annual_AI <-
  read_csv("data-raw/enviro-data/extractAnnualAridityIndex.csv") %>%
  mutate(AI_actual = AridityInd * 0.0001) # Multiply AI values according to CGIAR README file.

# Monthly precpipitation
monthly_precip <-
  read_csv("data-raw/enviro-data/extractedPrecipValues.csv") %>%
  rename(
    sept_precip = September,
    aug_precip = August_Pre,
    july_precip = July_Preci,
    june_precip = June_Preci,
    may_precip = May_Precip
  ) %>%
  mutate(meanMonthlyPrecip = rowMeans(
    select(
      .,
      sept_precip,
      aug_precip,
      july_precip,
      june_precip,
      may_precip
    ),
    na.rm = TRUE
  ))

# Monthly PET
monthly_pet <-
  read_csv("data-raw/enviro-data/extractMonthlyPET.csv") %>%
  rename(
    sept_pet = Septempter,
    aug_pet = August_PET,
    july_pet = July_PET,
    june_pet = June_PET,
    may_pet = May_PET
  ) %>%
  mutate(meanMonthlyPET = rowMeans(
    select(.,
           sept_pet,
           aug_pet,
           july_pet,
           june_pet,
           may_pet),
    na.rm = TRUE
  ))

# Annual PET
annual_pet <- read_csv("data-raw/enviro-data/extractAnnualPET.csv")

# Bioclimatic variables
bioclim <- read_csv("data-raw/enviro-data/extractBioClim.csv")

#### MERGE AND WRITE DATASETS ####

# Merge datasets
allClimateData <- left_join(annual_AI, annual_pet, by = c("City", "Population", "Transect", "Lat.pop", "Long.pop")) %>%
  left_join(., monthly_pet, by = c("City", "Population", "Transect", "Lat.pop", "Long.pop")) %>%
  left_join(., monthly_precip, by = c("City", "Population", "Transect", "Lat.pop", "Long.pop")) %>%
  left_join(., bioclim, by = c("City", "Population", "Transect", "Lat.pop", "Long.pop")) %>%
  mutate(meanSMD = meanMonthlyPrecip - meanMonthlyPET)

# Write climate data to disk
write_csv(allClimateData, "data-clean/ClimateData_AllPopulations.csv")

