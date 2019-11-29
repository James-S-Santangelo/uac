# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Generates tables and figures for maanuscript

# Load required packages
library(tidyverse)

################
#### TABLES ####
################
## TABLE 1 

# Linear cline models only. Part of table 1
linearClines <- read_csv("analysis/clineModelOutput.csv") %>% 
  select(City, betaLinOnly, betaLin_AcHWE, betaLin_LiHWE)

# Extract number of populations and plants per city. Merge with linear clines model output above
table1 <- read.csv("data-clean/AllCities_AllPlants.csv") %>%
  group_by(City) %>%
  summarise(numPops = n_distinct(Population),
            numPlants = n()) %>%
  merge(., linearClines, by = "City")

# Write table 1 to disk
write_csv(table1, "analysis/tables/main-text/Table1_cityClineSummary.csv")

## TABLE S1

weather_data <- read_csv("data-clean/DailyNormals_AllCities_Filtered.csv")

# Summarize number of observations for each city
weather_summ <- weather_data %>% 
  group_by(City, STATION_NAME) %>%
  summarise(min_year = min(Year),
            max_year = max(Year),
            count = n()) %>%
  mutate(STATION_NAME = str_to_title(STATION_NAME))

# Write summarized data to disk
write_csv(weather_summ, "analysis/tables/supplemental/TableS1_DailyNormals_Summary.csv")

## TABLE S2

# Written in "CORRELATIONS AMONG CLIMATIC VARIABLES" section of masterAnalysisScript.R
# It can be loaded below
tableS2_corr_mat <- read_csv("analysis/tables/supplemental/TableS2_weatherCorrMat.csv")

## TABLE S3

# Source functions
source("scripts/functions.R")

# Load population-mean dataset
datPops <- read_csv("data-clean/AllCities_AllPopulations.csv")

# Get the best fit cline model for each city and locus
# Possible values for the response variables
possibleResponses <- c("freqHCN", "AcHWE", "LiHWE")

tbl_out <- c() # vector to hold resulting dataframes

# Loop over response variables
for(res in possibleResponses){

  tbl_out[[res]] <- datPops %>% 
    split(.$City) %>%  # Separate list for each city
    map_dfr(., ~tidy(getBestFitClineModelOrder(res, .)), .id = "City") %>% 
    mutate(response = res) %>% 
    spread(key = response, value = x)
  
}

# Left join all dataframes with model orders
allModelOutputs <- reduce(tbl_out, left_join, by = "City")

# Write model order table to disk
write_csv(allModelOutputs, "analysis/tables/supplemental/TableS3_clineOrderData.csv")

## TABLE S4

# Written in "ANALYSIS PREDICTING MEAN HCN FREQUENCIES" section of masterAnalysisScript.R
# It can be loaded below
tableS4_HCN_dredge_models <- read_csv("analysis/tables/supplemental/TableS4_HCN_dredge_output.csv")

## TABLE S5

# Written in "ANALYSIS PREDICTING MEAN HCN FREQUENCIES" section of masterAnalysisScript.R
# It can be loaded below
tableS5_fullCoeffs <- read_csv("analysis/tables/supplemental/TableS5_freqHCN_FullModelAvg.csv")


