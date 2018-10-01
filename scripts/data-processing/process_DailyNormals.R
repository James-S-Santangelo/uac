
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
  filter(Year >= 1980 & Year <= 2015,  # Only include January and February from 1980 to 2015
         Month == 1 | Month == 2) %>%
  group_by(City, Year) %>%
  filter(n_distinct(Month) == 2) %>% # Remove observations if both months aren't represented
  ungroup() %>%
  group_by(City, Year, Month) %>%
  filter(n() > 10) %>%  # Remove observations if there isn't at least 10 observations in a month
  ungroup() %>%
  filter(!rowSums(. == -9999)) %>% # Remove rows if any values are -9999 (i.e. missing data)
  na.omit() # Omit NAs (if any)

# Write filtered data to disk
write_csv(weather_data_filtered, "data-clean/DailyNormals_AllCities_Filtered.csv")
