
# Load required packages
library(dplyr)

# Load in Master data with presence/absence of HCN for 12 cities sampled by JSS
datJSS <- read.csv("data-raw/AllPlants_JSS-Cities.csv", na.strings = "NA") %>%
  select(-Comments.Ac.Li, -Comments.HCN, -Pop.density, -Check.HCN, -Distance) %>% # Select required columns
  mutate(Transect = "NA") # Add transect column as NA for later merging with Ken's data

# Load in KAT's Population data to extract Lat and Long
datKATPops <- read.csv("data-raw/AllPops_KAT-Cities.csv", na.strings = "NA") %>%
  select(City, Pop, Transect, Latitude, Longitude) %>% # Select required columns
  rename(Population = Pop, Lat.pop = Latitude, Long.pop = Longitude) %>% # Rename column so they're consistent with my data
  mutate(Transect = ifelse(City == "B" | City == "M" | City == "Y", "NA", as.character(Transect)))

# Load in KAT's Plant data, which will be merged with JSS's data
datKATPlants <- read.csv("data-raw/AllPlants_KAT-Cities.csv", na.strings = "NA") %>%
  select(CITY, Pop_Num, Plant_Num, Transect, HCN_Result, FINAL_AC, FINAL_Li) %>% # Select required columns
  rename(City = CITY, Population = Pop_Num, Plant = Plant_Num, Locus.Ac = FINAL_AC, Locus.Li = FINAL_Li) %>% # Rename columns so they're consistent with my data
  mutate(Dmg.1 = "NA", Dmg.2 = "NA", Dmg.Avg = "NA") %>%
  mutate(Transect = ifelse(City == "B" | City == "M" | City == "Y", "NA", as.character(Transect))) %>%
  mutate(Locus.Ac = ifelse(City != "Toronto" & Transect != "A", "NA", Locus.Ac), 
         Locus.Li = ifelse(City != "Toronto" & Transect != "A", "NA", Locus.Li))

# Add Latitude and Longitude for KAT's Pop data to KAT's AllPlant dataset.
# Done because JSS's AllPlant data contains Lat and Long coordinates
datKATPlants <- merge(datKATPlants, datKATPops, 
                      by = c("City", "Population", "Transect"), 
                      all.x = TRUE)

# Add Damage columns to KAT's AllPlants data. Required because JSS measured damaged for
# some plants and these should be included in the final merged dataset. Since KAT did not measure damage, 
# these will be input as "NA"
datKATPlants <- datKATPlants %>%
  mutate(City = ifelse(City == "T", "Toronto", 
                       ifelse(City == "M", "Montreal", 
                              ifelse(City == "B", "Boston", "NewYork")))) 

# Merge AllPlant datasets for JSS and KAT. 
datAllPlants <- rbind(datJSS, datKATPlants) 

# Load dataset with latitudes and longitudes for each population and for each city centre
datLatLong <- read.csv("data-raw/Lat-Longs_AllCities_AllPops.csv")

# Source R script with Haversine formula for distance calculation
source(file = "scripts/haversine.R")

# Add distance to Lat Long dataset
datLatLong = datLatLong %>%
  mutate(Distance = haversine(Long.pop, Lat.pop, Long.City, Lat.City)) %>%
  select(City, Population, Transect, Distance) 

# Merge Distance with AllPlants dataset
datAllPlants <- merge(datAllPlants, datLatLong, 
        by = c("City", "Population", "Transect"), 
        all.x = TRUE)

# Write merged AllPlant dataset to disk.
write.csv(datAllPlants, "data-clean/AllCities_AllPlants.csv", na = "NA", row.names = FALSE)
