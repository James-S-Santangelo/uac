# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Main script to reproduce manuscript results

###############
#### SETUP ####
###############

## LOAD PACKAGES ##

library(tidyverse)
library(broom)
library(Hmisc)
library(factoextra)
library(vegan)
library(MuMIn)
library(car)
library(patchwork)
source('scripts/functions.R')

##################################
#### STEP 1: PROCESS RAW DATA ####
##################################

# 1.1: Merge Thompson (2016) data with mine

source("scripts/data-processing/merge-clean_JSS-KAT_AllPlants.R")

# 1.2: Create population-mean dataset

source("scripts/data-processing/create_AllCities_AllPopulations.R")

# 1.3: Combine climatic variables

source("scripts/data-processing/create_ClimateData_AllPopulations.R")

# 1.4: Clean and filter Daily Normals data

source("scripts/data-processing/process_DailyNormals.R")

# 1.5: Create clines and environment summary by city

## 1.5.1 Create summary (e.g., coefficients and p-values of linear
## and logistic clines by city

source("scripts/createClineModelOutput.R")

## 1.5.2: Create city-level summary of environmental data and clines

source("scripts/data-processing/create_citySummaryData.R")

## 1.6: Clean raw-haplotype data

source("scripts/data-processing/create-haplotypeData.R")

## 1.7: Process latitude and longitude coordinates for creating maps
## (Resulting file not used for analyses, only for creating supplemental maps)

source("scripts/data-processing/create_latLong_forMaps_CSVs.R")







