# Script to generate a summary of clines models and environmental variables for each city

# Load required packages
library(tidyverse)

#### LOAD IN DATASETS ####

# Load dataset with slopes of clines
slopes <- read_csv("analysis/clineModelOutput.csv") %>%
  select(-pvalCyanSlopeLin, -pvalSlopeQuad, -pvalAcSlopeLin,
         -pvalAcSlopeQuad, -pvalLoSlopeLin, -pvalLiSlopeQuad)



# Load climate data and take mean by city 
climate_data <- read_csv("data-clean/ClimateData_AllPopulations.csv") %>%
  group_by(City) %>%
  summarise(annualAI = mean(AI_actual, na.rm = TRUE),
            monthlyPET = mean(meanMonthlyPET, na.rm = TRUE),
            annualPET = mean(annualPET, na.rm = TRUE),
            monthlyPrecip = mean(meanMonthlyPrecip, na.rm = TRUE),
            mwtBio = mean(MWT, na.rm = TRUE),
            mstBio = mean(MST, na.rm = TRUE),
            smd = mean(meanSMD, na.rm = TRUE))

# Load city centre lats and longs
lat_long <- read_csv("data-raw/Lat-longs_City-centres.csv") %>%
  select(-Comments)

# Load in all population data and get mean HCN and alleles Freqs by city
datPops <- read_csv("data-clean/AllCities_AllPopulations.csv") 
gene_freqs <- datPops %>%
  group_by(City) %>%
  summarise(freqHCN = mean(freqHCN, na.rm = TRUE),
            AcHWE = mean(AcHWE, na.rm = TRUE), 
            LiHWE = mean(LiHWE, na.rm = TRUE))

# Add slopes from linear model of HCN against standardized distance
# Same as 'cyanSlopeLin' in slopes data for most cities. 
# Different if best fit model was quadratic. Required to test
# for predictors of clines strength
slopes <- datPops %>%
  group_by(City) %>%
  do(mod = lm(freqHCN ~ std_distance, data = .)) %>%
  broom::tidy(., mod) %>%
  filter(term == "std_distance") %>%
  select(estimate, p.value) %>%
  rename(cyanSlopeForAnalysis = estimate, cyanPvalSlopeAnalysis = p.value) %>%
  merge(., slopes, by = "City")

# Load in daily normals (i.e. weather data)
daily_normals <- read_csv("data-clean/DailyNormals_AllCities_Filtered.csv") %>%
  mutate(negNoSnow = ifelse(TMIN < 0 & SNWDcm == 0, 1, 0)) 

# Days below 0 with no snow summary dataset
daysNegNoSnow <- daily_normals %>%
  group_by(City, Year) %>%
  summarise(daysNegNoSnow = sum(negNoSnow)) %>%
  ungroup() %>%
  group_by(City) %>%
  summarise(daysNegNoSnow = mean(daysNegNoSnow))

# Final daily normal dataset
daily_normals_merged <- daily_normals %>%
  group_by(City) %>%
  summarise(snow_depth = mean(SNWDcm, na.rm = TRUE),
            snowfall = mean(SNOWcm, na.rm = TRUE),
            # mstWea = max(TMAX, na.rm = TRUE), 
            mwtWea = min(TMIN, na.rm = TRUE)) %>%
  left_join(., daysNegNoSnow)

 #### MERGE DATASETS TO FINAL CITY SUMMARY DATASET ####

citySummaryData <- left_join(slopes, gene_freqs) %>%
  left_join(., lat_long) %>%
  left_join(., climate_data) %>%
  left_join(., daily_normals_merged)
  
# Write dataset 
write_csv(citySummaryData, "data-clean/citySummaryData.csv")


