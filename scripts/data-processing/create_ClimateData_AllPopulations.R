# Load required packages
library(tidyverse)

#### LOAD CLIMATE DATASETS ####

# Aridity index
annual_AI <-
  read_csv("GIS/extract-climate-data/AI_annual/extractAnnualAridityIndex.csv") %>%
  mutate(AI_actual = AridityInd * 0.0001) # Multiply AI values according to CGIAR README file.

# Monthly precpipitation
monthly_precip <-
  read_csv("GIS/extract-climate-data/WorldClim_Precip/extractedPrecipValues.csv") %>%
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
  read_csv("GIS/extract-climate-data/Global_PET/extractMonthlyPET.csv") %>%
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
annual_pet <- read_csv("GIS/extract-climate-data/Global_PET/extractAnnualPET.csv")

# Bioclimatic variables
bioclim <- read_csv("GIS/extract-climate-data/BioClim/extractBioClim.csv")

#### MERGE AND WRITE DATASETS ####

# Merge datasets
allClimateData <- left_join(annual_AI, annual_pet) %>%
  left_join(., monthly_pet) %>%
  left_join(., monthly_precip) %>%
  left_join(., bioclim) %>%
  mutate(meanSMD = meanMonthlyPrecip - meanMonthlyPET)

# Write climate data to disk
write_csv(allClimateData, "data-clean/ClimateData_AllPopulations.csv")

