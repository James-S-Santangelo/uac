# Load in data containing all plants from all populations
datAllCities_AllPlants <- read.csv("data-raw/20170415_UAC_AllCities_AllPlants.csv")

# Load dataset with population level data for all cities
datAllCities_AllPops <- read.csv("data-raw/20170415_UAC_AllCities_AllPopulations.csv")
datAllCities_AllPops <- datAllCities_AllPops[-c(91, 290),] # Duplicated Rows.
# datAllCities_AllPops[which(duplicated(datAllCities_AllPops[,c("City", "Population")])),]

# Load packages
library(dplyr)

# Columns required for merging
to_merge <- datAllCities_AllPops %>%
  select(City, Population, Lat.pop, Long.pop)

# Temporary dataset used for manipulating
tmp <- datAllCities_AllPlants

# Merge dataframes by Lat and Long
tmp <- merge(tmp, to_merge, 
      by.x = c("City", "Population"), 
      by.y = c("City", "Population"), all = TRUE)

# Remove unnecessary columns and rename
tmp <- tmp %>%
  select(-Comments.Ac.Li, -Comments.HCN, -Check.HCN, -Lat.pop.x, -Long.pop.x) %>%
  rename(Lat.pop = Lat.pop.y, Long.pop = Long.pop.y)

# Write clean dataset
write.csv(tmp, "data-clean/AllCities_AllPlants.csv", na = "NA", row.names = FALSE)
