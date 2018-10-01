
# Here is will filter the daily temerature and snow normals for use in predicting the strength of clines
# The daily normals were downloaded on September 6th, 2016 from the National Centers for Environmental
# Information, which are part of NOAA (https://www.ncdc.noaa.gov/cdo-web/datatools/selectlocation).

# Load in required packages
library(tidyverse)

# Load in daily summaries
weather_data <- read_csv("data-raw/DailyNormals_AllCities_Unfiltered.csv", col_names = TRUE)
glimpse(weather_data)

# Filter data
weather_data_filtered <- weather_data %>%
  
  # Only include January and February from 1980 to 2015
  filter(Year >= 1980 & Year <= 2015, 
         Month == 1 | Month == 2) %>%
  group_by(City, Year) %>%
  
  # Remove observations if both months aren't represented
  filter(n_distinct(Month) == 2) %>%
  ungroup() %>%
  group_by(City, Year, Month) %>%
  
  # Remove observations if there isn't at least 10 observations in a month
  filter(n() > 10) %>%
  ungroup() %>%
  
  # Remove rows if any values are -9999 or NA
  filter(!rowSums(. == -9999)) %>% 
  na.omit() 

# Write filtered data to disk
write_csv(weather_data_filtered, "data-clean/DailyNormals_AllCities_Filtered.csv")
