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

## CREATE DIRECTORIES ##

paths <- c("analysis/figures/individualCline_biplots/Ac/",
           "analysis/figures/individualCline_biplots/HCN/",
           "analysis/figures/individualCline_biplots/Li/",
           "analysis/figures/main-text/",
           "analysis/figures/supplemental/",
           "analysis/figures/supplemental/",
           "data-clean/",
           "analysis/tables/main-text/",
           "analysis/tables/supplemental/",
           "analysis/inividual-cline-models/AcHWE/",
           "analysis/inividual-cline-models/LiHWE/",
           "analysis/inividual-cline-models/freqHCN/")

purrr::walk(paths, dir.create, recursive = T, showWarnings = T)

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

####################################
#### STEP 2: RUN POWER ANALYSES ####
####################################

# Will run power analyses. Output only used for plotting
source("scripts/powerAnalyses.R")

###################################
#### STEP 3: RUN MAIN ANALYSES ####
###################################

# Run all main analyses in the manuscript. 
# Will create dataframes and objects used for
# creating most figures and tables.
source("scripts/analysisScript.R")

############################################
### STEP 4: GENERATE TABLES AND FIGURES ####
############################################

source("scripts/figures_tables.R")






