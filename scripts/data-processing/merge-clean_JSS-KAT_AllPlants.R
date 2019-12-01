# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script to merge the data collected by Thompson et al (2016) with the additional 12 cities
# collected by James Santangelo. Merges individual plant-level phenotype data. 

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
  select(CITY, Pop_Num, Plant_Num, Transect, HCN_Result, LINAMARIN_RESULT, LINAMARASE_RESULT) %>% # Select required columns
  rename(City = CITY, Population = Pop_Num, Plant = Plant_Num, Locus.Ac = LINAMARIN_RESULT, Locus.Li = LINAMARASE_RESULT) %>% # Rename columns so they're consistent with my data
  mutate(Dmg.1 = "NA", Dmg.2 = "NA", Dmg.Avg = "NA") %>%
  mutate(Transect = ifelse(City == "B" | City == "M" | City == "Y", "NA", as.character(Transect))) %>%
  mutate(Locus.Ac = ifelse(City != "T", "NA", Locus.Ac),
         Locus.Li = ifelse(City != "T", "NA", Locus.Li)) %>%
  group_by(City, Population) %>%
  mutate(Locus.Ac = ifelse(City == "T" & any(!is.na(Locus.Ac)) & HCN_Result == 1, 1, Locus.Ac),
         Locus.Li = ifelse(City == "T" & any(!is.na(Locus.Li)) & HCN_Result == 1, 1, Locus.Li))

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

# Replace Locus.Ac and Locus.Li values for New York with data for the individual genes
# that was collected by Jibran Syed after the publication of Thompson (2016). 
datNY_genes <- read.csv("data-raw/AllPlants_NY_IndGenes_Jibran.csv") %>%
  select(City, Population, Plant, Transect, HCN_Result, Locus.Li, Locus.Ac) %>% # Select required columns
  mutate(Transect = "NA") %>%
  mutate(City = "NewYork")
datNY_genes <- datAllPlants %>%
  filter(City == "NewYork") %>%
  select(-Locus.Li, -Locus.Ac) %>%
  merge(., datNY_genes,
        by = c("City", "Population", "Plant", "HCN_Result", "Transect"))
datAllPlants <- datAllPlants %>%
  filter(City != "NewYork") %>%
  rbind(., datNY_genes) 

# Load dataset with latitudes and longitudes for each population. Merge with
# dataframe containing lat/longs for city centres.
datLatLong <- read.csv("data-raw/Lat-Longs_AllCities_AllPops.csv")
datCityCentres <- read.csv("data-raw/Lat-longs_City-centres.csv") %>%
  select(City, Latitude, Longitude) %>%
  rename(Lat.City = Latitude, Long.City = Longitude)
datLatLong <- merge(datLatLong, datCityCentres, by = "City", all.x = TRUE)

# Source R script with Haversine formula for distance calculation
source(file = "scripts/haversine.R")

# Add distance to Lat Long dataset
datLatLong <- datLatLong %>%
  mutate(Distance = haversine(Long.pop, Lat.pop, Long.City, Lat.City)) %>%
  select(City, Population, Transect, Distance) 

# Merge Distance with AllPlants dataset
datAllPlants <- merge(datAllPlants, datLatLong, 
        by = c("City", "Population", "Transect"), 
        all.x = TRUE) %>% 
  select(-contains("Dmg")) %>% 
  mutate(Population = as.character(case_when(City =="Boston" & Population == "4A" ~ "44",
                                             City =="Boston" & Population == "4B" ~ "45",
                                             City =="NewYork" & Population == "11-I" ~ "11",
                                             City =="NewYork" & Population == "11-II" ~ "12",
                                             TRUE ~ Population)))

# Write merged AllPlant dataset to disk.
write.csv(datAllPlants, "data-clean/AllCities_AllPlants.csv", na = "NA", row.names = FALSE)
