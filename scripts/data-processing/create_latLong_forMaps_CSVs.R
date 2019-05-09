# Load required packages
library(tidyverse)

datPops <- read_csv("data-clean/AllCities_AllPopulations.csv")
latLongs <- read_csv("data-raw/Lat-Longs_AllCities_AllPops.csv") %>%
  select(-Transect)
city_centres <- read_csv("data-raw/Lat-longs_City-centres.csv") %>%
  select(-Comments)

latLongs_forMap <- function(df, city_centres, latLongs){
  
  city_name <- df$City[1]

  centre <- city_centres %>% filter(City == city_name)
  df_out <- datPops %>%
    filter(City == city_name) %>%
    select(City, Transect, Population, freqHCN) %>%
    left_join(latLongs %>% filter(City == city_name)) %>%
    rbind(c(centre$City, "NA","City Centre", "NA", centre$Latitude, centre$Longitude))
  
  return(df_out)
}

# Create list with city dataframes as elements
city_df_list <- datPops %>% split(.$City)

# Generate single dataframe with populations and City centres
allPops_withCentres <- purrr::map_dfr(city_df_list, latLongs_forMap, 
               city_centres = city_centres, 
               latLongs = latLongs)
write_csv(allPops_withCentres, "GIS/city_maps/latLongs_CSVs/allPops_withCentres.csv")
